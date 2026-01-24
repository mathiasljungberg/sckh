"""Integration tests comparing Python implementation to Fortran reference.

These tests use the exact input files from the SCKH_PES test case and compare
the Python output against the Fortran reference output.
"""

import numpy as np
import pytest
from pathlib import Path

from python_scripts.dynamics_1d import (
    load_full_config,
    SpectrumCalculator,
    trajectory,
)


# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"
REFERENCE_DIR = FIXTURES_DIR / "reference"


def _fix_paths(config):
    """Convert relative paths in config to absolute paths based on fixtures dir."""
    # Fix dynamics paths
    config.dynamics.pes_initial = FIXTURES_DIR / config.dynamics.pes_initial.name
    config.dynamics.pes_dynamics = FIXTURES_DIR / config.dynamics.pes_dynamics.name

    # Fix spectrum paths
    config.spectrum.pes_final_list = [
        FIXTURES_DIR / p.name for p in config.spectrum.pes_final_list
    ]
    config.spectrum.dipole_final_list = [
        FIXTURES_DIR / p.name for p in config.spectrum.dipole_final_list
    ]
    if config.spectrum.pes_lp_corr:
        config.spectrum.pes_lp_corr = FIXTURES_DIR / config.spectrum.pes_lp_corr.name

    return config


class TestFortranCompatibility:
    """Integration tests comparing Python (fortran mode) to Fortran reference."""

    @pytest.fixture
    def config_fortran_mode(self):
        """Load config and set to fortran compatibility mode."""
        config = load_full_config(FIXTURES_DIR / "dynamics.yaml")
        config = _fix_paths(config)
        config.spectrum.compatibility_mode = "fortran"
        config.dynamics.sampling.compatibility_mode = "fortran"
        return config

    @pytest.fixture
    def fortran_reference(self):
        """Load Fortran reference spectrum."""
        data = np.loadtxt(REFERENCE_DIR / "fortran_sigma.dat")
        return {"omega": data[:, 0], "sigma": data[:, 1]}

    def test_spectrum_matches_fortran_reference(self, config_fortran_mode, fortran_reference):
        """Test that Python (fortran mode) matches Fortran output."""
        # Run dynamics
        runner = trajectory.DynamicsRunner(config_fortran_mode)
        result_dyn = runner.run(verbose=False)

        # Compute spectrum
        calculator = SpectrumCalculator(config_fortran_mode)
        calculator.load_surfaces()
        result_sp = calculator.compute_spectrum(result_dyn.trajectories)

        # Compare frequency grids
        assert len(result_sp.omega) == len(fortran_reference["omega"])
        np.testing.assert_allclose(
            result_sp.omega,
            fortran_reference["omega"],
            rtol=1e-4,
            err_msg="Frequency grids do not match",
        )

        # Compare peak positions
        py_peak_idx = np.argmax(result_sp.sigma_tot)
        f_peak_idx = np.argmax(fortran_reference["sigma"])

        py_peak_omega = result_sp.omega[py_peak_idx]
        f_peak_omega = fortran_reference["omega"][f_peak_idx]

        assert abs(py_peak_omega - f_peak_omega) < 0.01, (
            f"Peak position mismatch: Python={py_peak_omega:.4f} eV, "
            f"Fortran={f_peak_omega:.4f} eV"
        )

        # Compare peak intensities (within 1%)
        py_peak_sigma = result_sp.sigma_tot[py_peak_idx]
        f_peak_sigma = fortran_reference["sigma"][f_peak_idx]

        rel_diff = abs(py_peak_sigma - f_peak_sigma) / f_peak_sigma
        assert rel_diff < 0.01, (
            f"Peak intensity mismatch: Python={py_peak_sigma:.6f}, "
            f"Fortran={f_peak_sigma:.6f}, diff={rel_diff*100:.2f}%"
        )

        # Compare integrated intensity
        py_integral = np.trapezoid(result_sp.sigma_tot, result_sp.omega)
        f_integral = np.trapezoid(fortran_reference["sigma"], fortran_reference["omega"])

        rel_diff_integral = abs(py_integral - f_integral) / f_integral
        assert rel_diff_integral < 0.001, (
            f"Integrated intensity mismatch: Python={py_integral:.6f}, "
            f"Fortran={f_integral:.6f}, diff={rel_diff_integral*100:.3f}%"
        )

        # Compare RMS difference in significant region
        # Interpolate to common grid for detailed comparison
        omega_common = np.linspace(
            max(result_sp.omega[0], fortran_reference["omega"][0]),
            min(result_sp.omega[-1], fortran_reference["omega"][-1]),
            500,
        )
        py_interp = np.interp(omega_common, result_sp.omega, result_sp.sigma_tot)
        f_interp = np.interp(omega_common, fortran_reference["omega"], fortran_reference["sigma"])

        # Only compare where signal is significant (> 1% of max)
        mask = f_interp > f_interp.max() * 0.01
        rms_diff = np.sqrt(np.mean((py_interp[mask] - f_interp[mask]) ** 2))
        rel_rms = rms_diff / np.max(f_interp)

        assert rel_rms < 0.005, (
            f"RMS difference too large: {rel_rms*100:.2f}% of max intensity"
        )

    def test_e_mean_matches_expected(self, config_fortran_mode):
        """Test that E_mean is computed correctly."""
        runner = trajectory.DynamicsRunner(config_fortran_mode)
        result_dyn = runner.run(verbose=False)

        calculator = SpectrumCalculator(config_fortran_mode)
        calculator.load_surfaces()
        E_mean = calculator.compute_mean_transition_energy(result_dyn.trajectories)

        # Expected E_mean from Fortran (approximately 519.3 eV based on previous runs)
        assert 518 < E_mean < 521, f"E_mean={E_mean:.2f} eV is outside expected range"

    def test_number_of_trajectories(self, config_fortran_mode):
        """Test that correct number of trajectories is generated."""
        runner = trajectory.DynamicsRunner(config_fortran_mode)
        result_dyn = runner.run(verbose=False)

        expected_n_traj = (
            config_fortran_mode.dynamics.sampling.npoints_x
            * config_fortran_mode.dynamics.sampling.npoints_mom
        )
        assert len(result_dyn.trajectories) == expected_n_traj

    def test_trajectory_length(self, config_fortran_mode):
        """Test that trajectories have correct number of time steps."""
        runner = trajectory.DynamicsRunner(config_fortran_mode)
        result_dyn = runner.run(verbose=False)

        expected_nsteps = config_fortran_mode.dynamics.time.nsteps
        for traj in result_dyn.trajectories:
            assert len(traj.time) == expected_nsteps
            assert len(traj.x) == expected_nsteps


class TestStandardVsFortranMode:
    """Tests comparing standard mode to fortran mode."""

    @pytest.fixture
    def configs(self):
        """Load config for both modes."""
        config_standard = load_full_config(FIXTURES_DIR / "dynamics.yaml")
        config_standard = _fix_paths(config_standard)
        config_standard.spectrum.compatibility_mode = "standard"
        config_standard.dynamics.sampling.compatibility_mode = "standard"

        config_fortran = load_full_config(FIXTURES_DIR / "dynamics.yaml")
        config_fortran = _fix_paths(config_fortran)
        config_fortran.spectrum.compatibility_mode = "fortran"
        config_fortran.dynamics.sampling.compatibility_mode = "fortran"

        return {"standard": config_standard, "fortran": config_fortran}

    def test_both_modes_produce_valid_spectra(self, configs):
        """Both modes should produce physically reasonable spectra."""
        results = {}

        for mode, config in configs.items():
            runner = trajectory.DynamicsRunner(config)
            result_dyn = runner.run(verbose=False)

            calculator = SpectrumCalculator(config)
            calculator.load_surfaces()
            result_sp = calculator.compute_spectrum(result_dyn.trajectories)
            results[mode] = result_sp

        for mode, result in results.items():
            # Check spectrum is normalized (integral ≈ 1)
            integral = np.trapezoid(result.sigma_tot, result.omega)
            assert 0.99 < integral < 1.01, f"{mode} mode: integral={integral:.4f}"

            # Check E_mean is reasonable
            assert 500 < result.E_mean < 550, f"{mode} mode: E_mean={result.E_mean:.2f}"

            # Check peak is in expected range
            peak_idx = np.argmax(result.sigma_tot)
            peak_omega = result.omega[peak_idx]
            assert 520 < peak_omega < 530, f"{mode} mode: peak at {peak_omega:.2f} eV"

    def test_modes_give_similar_results(self, configs):
        """Standard and fortran modes should give similar (but not identical) results."""
        results = {}

        for mode, config in configs.items():
            runner = trajectory.DynamicsRunner(config)
            result_dyn = runner.run(verbose=False)

            calculator = SpectrumCalculator(config)
            calculator.load_surfaces()
            result_sp = calculator.compute_spectrum(result_dyn.trajectories)
            results[mode] = result_sp

        # Peak positions should be within 1 eV
        peak_standard = results["standard"].omega[np.argmax(results["standard"].sigma_tot)]
        peak_fortran = results["fortran"].omega[np.argmax(results["fortran"].sigma_tot)]

        assert abs(peak_standard - peak_fortran) < 1.0, (
            f"Peak positions differ by {abs(peak_standard - peak_fortran):.2f} eV"
        )

        # Peak intensities should be within 10%
        intensity_standard = np.max(results["standard"].sigma_tot)
        intensity_fortran = np.max(results["fortran"].sigma_tot)

        rel_diff = abs(intensity_standard - intensity_fortran) / intensity_fortran
        assert rel_diff < 0.10, f"Peak intensities differ by {rel_diff*100:.1f}%"
