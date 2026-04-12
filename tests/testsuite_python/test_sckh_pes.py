"""SCKH_PES testsuite: spectrum from PES surfaces and dynamics sampling.

Runs classical trajectory dynamics on PES surfaces, computes the SCKH
spectrum via FFT, and compares against the Fortran reference output.

This test mirrors tests/testsuite/SCKH_PES/ and also serves as an example
of the full dynamics-to-spectrum workflow.

Example usage (as a script)::

    from pathlib import Path
    from dynamics_1d import (
        load_full_config,
        DynamicsRunner,
        SpectrumCalculator,
    )

    config = load_full_config("dynamics.yaml")
    config.spectrum.compatibility_mode = "fortran"
    config.dynamics.sampling.compatibility_mode = "fortran"

    # Run dynamics
    runner = DynamicsRunner(config)
    result_dyn = runner.run()

    # Compute spectrum
    calculator = SpectrumCalculator(config)
    calculator.load_surfaces()
    result_sp = calculator.compute_spectrum(result_dyn.trajectories)

    # result_sp.omega -> frequency grid (eV)
    # result_sp.sigma_tot -> total cross-section
"""

from pathlib import Path

import numpy as np
import pytest

from dynamics_1d import (
    SpectrumCalculator,
    DynamicsRunner,
)
from dynamics_1d.config import DynamicsConfig, GridConfig, TimeConfig, SamplingConfig
from dynamics_1d.spectrum_config import InterpolationConfig, SpectrumConfig, FullConfig

from .conftest import TESTSUITE_FORTRAN_DIR, assert_spectrum_matches_ref


SCKH_PES_DIR = TESTSUITE_FORTRAN_DIR / "SCKH_PES"
INPUT_DIR = SCKH_PES_DIR / "input"
REF_DIR = SCKH_PES_DIR / "ref"


def _make_config() -> FullConfig:
    """Build configuration matching SCKH_PES/input/input.txt.

    Parameters are taken directly from the Fortran input file.
    """
    dynamics = DynamicsConfig(
        mu=1.0078825,
        grid=GridConfig(start=0.5, dx=0.025, npoints=77),
        time=TimeConfig(dt=0.1, nsteps=512),
        sampling=SamplingConfig(
            mode=1,
            npoints_x=10,
            npoints_mom=10,
            compatibility_mode="fortran",
        ),
        pes_initial=INPUT_DIR / "test_energy_tot.dat",
        pes_dynamics=INPUT_DIR / "test_energy_tot_exc.dat",
        units="angstrom",
        outfile="testout",
    )

    # Final state PES files (order from pes_file_list_f.txt)
    pes_final_list = [
        INPUT_DIR / f"pesfile_{i}.dat" for i in [9, 8, 7, 6, 5, 4]
    ]

    # Dipole files (order from dipole_file_list_f.txt)
    dipole_final_list = [
        INPUT_DIR / f"dipolefile_{i}.dat" for i in [9, 8, 7, 6, 5, 4]
    ]

    interpolation = InterpolationConfig(
        dipole_mode="DIPOLE",
        pes_final_list=pes_final_list,
        dipole_final_list=dipole_final_list,
        pes_lp_corr=INPUT_DIR / "test_lp_lp_energy.dat",
        compatibility_mode="fortran",
    )

    spectrum = SpectrumConfig(gamma_fwhm=0.18)

    return FullConfig(dynamics1d=dynamics, interpolation1d=interpolation, spectrum=spectrum)


@pytest.fixture
def config():
    """Configuration for SCKH_PES test."""
    return _make_config()


@pytest.fixture
def spectrum_result(config):
    """Run dynamics and compute spectrum (cached for test class)."""
    from dynamics_1d.spectrum import SurfaceInterpolator, SCKHSpectrumCalculator

    runner = DynamicsRunner(config)
    result_dyn = runner.run()

    interp = SurfaceInterpolator(config)
    interp.load_surfaces()
    sckh_trajs = interp.trajectories_to_sckh(result_dyn.trajectories)
    E_mean = interp.compute_mean_transition_energy(result_dyn.trajectories)
    D_ni = interp.compute_D_ni()

    calc = SCKHSpectrumCalculator(config)
    return calc.compute_spectrum(sckh_trajs, E_mean=E_mean, D_ni=D_ni)


class TestSCKHPESDynamics:
    """Test dynamics setup and trajectory generation."""

    def test_number_of_trajectories(self, config):
        """Should generate npoints_x * npoints_mom trajectories."""
        runner = DynamicsRunner(config)
        result_dyn = runner.run()

        expected = (
            config.dynamics1d.sampling.npoints_x
            * config.dynamics1d.sampling.npoints_mom
        )
        assert len(result_dyn.trajectories) == expected

    def test_trajectory_length(self, config):
        """Each trajectory should have nsteps time points."""
        runner = DynamicsRunner(config)
        result_dyn = runner.run()

        for traj in result_dyn.trajectories:
            assert len(traj.time) == config.dynamics1d.time.nsteps
            assert len(traj.x) == config.dynamics1d.time.nsteps


class TestSCKHPESSpectrum:
    """Test spectrum computation against Fortran reference."""

    def test_spectrum_matches_fortran_reference(self, spectrum_result):
        """Computed spectrum should match Fortran reference output.

        This is the main validation test, using the same comparison
        metric as the Fortran testsuite (my_diff_sum_ref).
        """
        assert_spectrum_matches_ref(
            spectrum_result.omega,
            spectrum_result.sigma_tot,
            REF_DIR / "testout_sigma.dat",
        )

    def test_frequency_grid_matches_reference(self, spectrum_result):
        """Frequency grid should match the Fortran reference."""
        ref = np.loadtxt(REF_DIR / "testout_sigma.dat")
        np.testing.assert_allclose(
            spectrum_result.omega, ref[:, 0], rtol=1e-4
        )

    def test_spectrum_is_normalized(self, spectrum_result):
        """Spectrum integral should be close to 1."""
        integral = np.trapezoid(
            spectrum_result.sigma_tot, spectrum_result.omega
        )
        assert 0.99 < integral < 1.01, f"Integral={integral:.4f}"

    def test_e_mean_in_expected_range(self, spectrum_result):
        """E_mean should be around 519 eV for this system."""
        assert 518 < spectrum_result.E_mean < 521

    def test_peak_position(self, spectrum_result):
        """Peak should be in the expected XES energy range."""
        peak_idx = np.argmax(spectrum_result.sigma_tot)
        peak_energy = spectrum_result.omega[peak_idx]
        assert 520 < peak_energy < 530, f"Peak at {peak_energy:.2f} eV"


class TestStandardVsFortranMode:
    """Compare standard and fortran compatibility modes.

    Both modes should produce physically reasonable spectra that are
    similar but not identical (different numerical algorithms).
    """

    def test_both_modes_produce_valid_spectra(self):
        """Both modes should give normalized spectra with correct peak."""
        from dynamics_1d.spectrum import SurfaceInterpolator, SCKHSpectrumCalculator

        for mode in ("standard", "fortran"):
            config = _make_config()
            config.interpolation1d.compatibility_mode = mode
            config.dynamics1d.sampling.compatibility_mode = mode

            runner = DynamicsRunner(config)
            result_dyn = runner.run()

            interp = SurfaceInterpolator(config)
            interp.load_surfaces()
            sckh_trajs = interp.trajectories_to_sckh(result_dyn.trajectories)
            E_mean = interp.compute_mean_transition_energy(result_dyn.trajectories)
            D_ni = interp.compute_D_ni()

            calc = SCKHSpectrumCalculator(config)
            result = calc.compute_spectrum(sckh_trajs, E_mean=E_mean, D_ni=D_ni)

            # Normalized
            integral = np.trapezoid(result.sigma_tot, result.omega)
            assert 0.99 < integral < 1.01, f"{mode}: integral={integral:.4f}"

            # Peak in range
            peak = result.omega[np.argmax(result.sigma_tot)]
            assert 520 < peak < 530, f"{mode}: peak at {peak:.2f} eV"
