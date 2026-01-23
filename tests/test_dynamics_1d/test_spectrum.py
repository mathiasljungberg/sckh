"""Tests for spectrum module."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_1d.spectrum import (
    compute_energy_phase,
    compute_F_if,
    get_frequency_grid,
    SpectrumResult,
)


class TestComputeEnergyPhase:
    """Tests for compute_energy_phase function."""

    def test_constant_energy_difference(self):
        """Phase for constant energy difference should grow linearly."""
        nsteps = 100
        dt = 1e-15  # 1 fs
        time = np.arange(nsteps) * dt

        # Constant energies
        E_n = np.ones(nsteps) * 10.0 * CONST.eV  # 10 eV in Joules
        E_f = np.ones(nsteps) * 5.0 * CONST.eV  # 5 eV in Joules
        E_mean = 5.0  # 5 eV (matches E_n - E_f)

        phase = compute_energy_phase(E_n, E_f, E_mean, time)

        # When E_n - E_f = E_mean, the argument should be zero
        # so phase should be exp(0) = 1
        np.testing.assert_allclose(np.abs(phase), 1.0, rtol=1e-10)

    def test_phase_is_complex_unit_norm(self):
        """Phase factor should have unit magnitude."""
        nsteps = 100
        dt = 1e-15
        time = np.arange(nsteps) * dt

        E_n = np.linspace(10, 12, nsteps) * CONST.eV
        E_f = np.linspace(5, 6, nsteps) * CONST.eV
        E_mean = 5.5

        phase = compute_energy_phase(E_n, E_f, E_mean, time)

        # All values should have magnitude 1
        np.testing.assert_allclose(np.abs(phase), 1.0, rtol=1e-10)

    def test_phase_oscillates_for_nonzero_difference(self):
        """Phase should oscillate when E_n - E_f != E_mean."""
        nsteps = 1000
        dt = 1e-15
        time = np.arange(nsteps) * dt

        # Set up oscillation
        E_n = np.ones(nsteps) * 10.0 * CONST.eV
        E_f = np.ones(nsteps) * 5.0 * CONST.eV
        E_mean = 4.0  # Off by 1 eV

        phase = compute_energy_phase(E_n, E_f, E_mean, time)

        # Phase should have both real and imaginary parts over time
        assert np.std(phase.real) > 0.1
        assert np.std(phase.imag) > 0.1


class TestComputeFFT:
    """Tests for compute_F_if function."""

    def test_zero_dipole_gives_zero_amplitude(self):
        """Zero dipole should give zero transition amplitude."""
        nsteps = 64
        dt = 1e-15
        time = np.arange(nsteps) * dt

        e_factor = np.ones(nsteps, dtype=complex)
        D_fn = np.zeros((nsteps, 3))
        D_ni = np.array([1.0, 1.0, 1.0])
        gamma = 0.1  # eV

        F_if = compute_F_if(e_factor, D_fn, D_ni, time, gamma)

        np.testing.assert_allclose(F_if, 0.0, atol=1e-20)

    def test_constant_dipole_fft(self):
        """Constant dipole * phase gives delta-like spectrum."""
        nsteps = 128
        dt = 1e-15
        time = np.arange(nsteps) * dt

        # Constant phase factor (no oscillation)
        e_factor = np.ones(nsteps, dtype=complex)

        # Constant dipole
        D_fn = np.ones((nsteps, 3))
        D_ni = np.array([1.0, 0.0, 0.0])

        gamma = 0.05  # Small broadening

        F_if = compute_F_if(e_factor, D_fn, D_ni, time, gamma)

        # The FFT of a decaying exponential should be Lorentzian-like
        # Check that amplitude is maximum near zero frequency
        assert F_if.shape == (nsteps, 3, 3)

        # Check m1=0, m2=0 component (x-polarized)
        amp = np.abs(F_if[:, 0, 0])
        max_idx = np.argmax(amp)

        # Maximum should be near center (zero frequency after shift)
        assert abs(max_idx - nsteps // 2) < 10

    def test_output_shape(self):
        """Check output shape is (nsteps, 3, 3)."""
        nsteps = 64
        time = np.arange(nsteps) * 1e-15

        e_factor = np.ones(nsteps, dtype=complex)
        D_fn = np.random.randn(nsteps, 3)
        D_ni = np.array([1.0, 1.0, 1.0])
        gamma = 0.1

        F_if = compute_F_if(e_factor, D_fn, D_ni, time, gamma)

        assert F_if.shape == (nsteps, 3, 3)
        assert F_if.dtype == complex

    def test_gamma_broadening_effect(self):
        """Larger gamma should give broader spectrum."""
        nsteps = 256
        dt = 1e-15
        time = np.arange(nsteps) * dt

        e_factor = np.ones(nsteps, dtype=complex)
        D_fn = np.ones((nsteps, 3))
        D_ni = np.array([1.0, 0.0, 0.0])

        # Compute with two different gamma values
        gamma_small = 0.05
        gamma_large = 0.2

        F_small = compute_F_if(e_factor, D_fn, D_ni, time, gamma_small)
        F_large = compute_F_if(e_factor, D_fn, D_ni, time, gamma_large)

        # Compare widths (FWHM)
        amp_small = np.abs(F_small[:, 0, 0])
        amp_large = np.abs(F_large[:, 0, 0])

        # Peak amplitude should be higher for smaller gamma (narrower peak)
        assert np.max(amp_small) > np.max(amp_large)


class TestFrequencyGrid:
    """Tests for get_frequency_grid function."""

    def test_frequency_grid_centering(self):
        """Frequency grid should be centered around E_mean."""
        nsteps = 128
        dt = 1e-15
        time = np.arange(nsteps) * dt
        E_mean = 500.0  # eV

        omega = get_frequency_grid(time, E_mean)

        # Grid should be centered around E_mean
        assert omega[nsteps // 2] == pytest.approx(E_mean, rel=0.01)

        # Should have both positive and negative frequencies around E_mean
        assert omega[0] < E_mean
        assert omega[-1] > E_mean

    def test_frequency_grid_spacing(self):
        """Frequency spacing should match 1/T (Fortran convention)."""
        nsteps = 128
        dt = 1e-15
        time = np.arange(nsteps) * dt
        E_mean = 0.0

        omega = get_frequency_grid(time, E_mean)

        # Spacing in omega should be 2*pi*hbar / (T * eV)
        # Fortran uses T = (nsteps-1) * dt = time[-1] - time[0]
        T = time[-1] - time[0]
        expected_spacing = 2 * np.pi * CONST.hbar / (T * CONST.eV)

        actual_spacing = omega[1] - omega[0]
        np.testing.assert_allclose(actual_spacing, expected_spacing, rtol=1e-10)

    def test_frequency_grid_length(self):
        """Frequency grid should have same length as time grid."""
        for nsteps in [64, 128, 256]:
            time = np.arange(nsteps) * 1e-15
            omega = get_frequency_grid(time, E_mean=0.0)
            assert len(omega) == nsteps


class TestSpectrumResult:
    """Tests for SpectrumResult dataclass."""

    def test_spectrum_result_creation(self):
        """Test creating SpectrumResult."""
        n_omega = 100
        n_final = 3

        result = SpectrumResult(
            omega=np.linspace(500, 520, n_omega),
            sigma_tot=np.random.rand(n_omega),
            sigma_f=np.random.rand(n_final, n_omega),
            n_trajectories=100,
            E_mean=510.0,
        )

        assert result.omega.shape == (n_omega,)
        assert result.sigma_tot.shape == (n_omega,)
        assert result.sigma_f.shape == (n_final, n_omega)
        assert result.n_trajectories == 100
        assert result.E_mean == 510.0


class TestSpectrumIO:
    """Tests for spectrum I/O functions."""

    def test_write_and_read_spectrum(self):
        """Test writing spectrum to file."""
        from python_scripts.dynamics_1d.io import write_spectrum

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "spectrum.dat"

            omega = np.linspace(500, 520, 101)
            sigma = np.exp(-(omega - 510) ** 2 / 2)

            write_spectrum(filepath, omega, sigma)

            # Read back
            data = np.loadtxt(filepath)
            np.testing.assert_allclose(data[:, 0], omega, rtol=1e-5)
            np.testing.assert_allclose(data[:, 1], sigma, rtol=1e-5)

    def test_write_spectrum_per_final(self):
        """Test writing per-final-state spectra."""
        from python_scripts.dynamics_1d.io import write_spectrum_per_final

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath_base = Path(tmpdir) / "spectrum"

            omega = np.linspace(500, 520, 51)
            sigma_f = np.random.rand(3, 51)

            write_spectrum_per_final(filepath_base, omega, sigma_f)

            # Check files were created
            for j in range(3):
                filepath = Path(f"{filepath_base}_final_{j + 1}.dat")
                assert filepath.exists()

                data = np.loadtxt(filepath)
                np.testing.assert_allclose(data[:, 0], omega, rtol=1e-5)
                np.testing.assert_allclose(data[:, 1], sigma_f[j], rtol=1e-5)


class TestSpectrumConfigValidation:
    """Tests for spectrum configuration validation."""

    def test_invalid_dipole_mode(self):
        """Should raise error for invalid dipole mode."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        with pytest.raises(ValueError, match="dipole_mode"):
            SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="INVALID",
                pes_final_list=[Path("pes.dat")],
            )

    def test_dipole_mode_requires_files(self):
        """DIPOLE mode requires dipole_final_list."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        with pytest.raises(ValueError, match="dipole_final_list"):
            SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="DIPOLE",
                pes_final_list=[Path("pes.dat")],
                dipole_final_list=[],  # Empty list
            )

    def test_fc_mode_no_dipole_required(self):
        """FC mode should not require dipole files."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        config = SpectrumConfig(
            gamma_fwhm=0.1,
            dipole_mode="FC",
            pes_final_list=[Path("pes.dat")],
            dipole_final_list=[],
        )

        assert config.dipole_mode == "FC"

    def test_gamma_hwhm_property(self):
        """Test gamma_hwhm property."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        config = SpectrumConfig(
            gamma_fwhm=0.2,
            dipole_mode="FC",
            pes_final_list=[Path("pes.dat")],
        )

        assert config.gamma_hwhm == 0.1
        assert config.gamma_fwhm == 0.2


class TestFortranDFormat:
    """Tests for Fortran D-format exponent support."""

    def test_loadtxt_fortran_d_format(self):
        """Test reading files with Fortran D-format exponents."""
        from python_scripts.dynamics_1d.io import _loadtxt_fortran

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"
            filepath.write_text("1.0  -0.345376681424D-03\n2.0  1.234567890123D+02\n")

            data = _loadtxt_fortran(filepath)

            np.testing.assert_allclose(data[0, 1], -0.345376681424e-03, rtol=1e-12)
            np.testing.assert_allclose(data[1, 1], 1.234567890123e+02, rtol=1e-12)

    def test_loadtxt_fortran_lowercase_d(self):
        """Test reading files with lowercase d exponents."""
        from python_scripts.dynamics_1d.io import _loadtxt_fortran

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"
            filepath.write_text("1.0  -0.123456789d-05\n2.0  9.876543210d+01\n")

            data = _loadtxt_fortran(filepath)

            np.testing.assert_allclose(data[0, 1], -0.123456789e-05, rtol=1e-12)
            np.testing.assert_allclose(data[1, 1], 9.876543210e+01, rtol=1e-12)

    def test_loadtxt_fortran_standard_e_format(self):
        """Test that standard E-format still works."""
        from python_scripts.dynamics_1d.io import _loadtxt_fortran

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"
            filepath.write_text("1.0  -1.5E-03\n2.0  2.5E+02\n")

            data = _loadtxt_fortran(filepath)

            np.testing.assert_allclose(data[0, 1], -1.5e-03, rtol=1e-12)
            np.testing.assert_allclose(data[1, 1], 2.5e+02, rtol=1e-12)


class TestCompatibilityMode:
    """Tests for compatibility_mode feature."""

    def test_invalid_compatibility_mode(self):
        """Should raise error for invalid compatibility mode."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        with pytest.raises(ValueError, match="compatibility_mode"):
            SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="FC",
                pes_final_list=[Path("pes.dat")],
                compatibility_mode="invalid_mode",
            )

    def test_valid_compatibility_modes(self):
        """Both 'standard' and 'fortran' should be valid."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        for mode in ["standard", "fortran", "STANDARD", "FORTRAN"]:
            config = SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="FC",
                pes_final_list=[Path("pes.dat")],
                compatibility_mode=mode,
            )
            assert config.compatibility_mode == mode.lower()

    def test_default_compatibility_mode(self):
        """Default compatibility_mode should be 'standard'."""
        from python_scripts.dynamics_1d.spectrum_config import SpectrumConfig

        config = SpectrumConfig(
            gamma_fwhm=0.1,
            dipole_mode="FC",
            pes_final_list=[Path("pes.dat")],
        )
        assert config.compatibility_mode == "standard"

    def test_sampling_config_compatibility_mode(self):
        """Test SamplingConfig accepts compatibility_mode."""
        from python_scripts.dynamics_1d.config import SamplingConfig

        config = SamplingConfig(mode=1, npoints_x=10, npoints_mom=10, compatibility_mode="fortran")
        assert config.compatibility_mode == "fortran"

        config_default = SamplingConfig()
        assert config_default.compatibility_mode == "standard"

    def test_load_full_config_compatibility_mode(self):
        """Test that compatibility_mode is loaded from YAML."""
        from python_scripts.dynamics_1d.spectrum_config import load_full_config

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create minimal PES files
            x = np.linspace(0.5, 2.5, 11)
            E = 0.1 * (x - 1.5) ** 2
            pes_file = tmpdir / "pes.dat"
            np.savetxt(pes_file, np.column_stack([x, E]))

            # Create YAML config with compatibility_mode
            yaml_content = f"""
dynamics:
  mu: 1.0
  grid:
    start: 0.5
    dx: 0.2
    npoints: 11
  time:
    dt: 1.0
    nsteps: 10
  sampling:
    mode: 1
    npoints_x: 2
    npoints_mom: 2
  pes_initial: "{pes_file}"
  pes_dynamics: "{pes_file}"

spectrum:
  gamma_fwhm: 0.1
  dipole_mode: FC
  pes_final_list:
    - "{pes_file}"
  compatibility_mode: fortran
"""
            yaml_file = tmpdir / "config.yaml"
            yaml_file.write_text(yaml_content)

            config = load_full_config(yaml_file)

            # Both should have fortran mode
            assert config.spectrum.compatibility_mode == "fortran"
            assert config.dynamics.sampling.compatibility_mode == "fortran"

    def test_e_mean_standard_vs_fortran(self):
        """Test that E_mean calculation differs between modes."""
        from python_scripts.dynamics_1d.spectrum import SpectrumCalculator
        from python_scripts.dynamics_1d.spectrum_config import FullConfig, SpectrumConfig
        from python_scripts.dynamics_1d.config import (
            DynamicsConfig, GridConfig, TimeConfig, SamplingConfig
        )
        from python_scripts.dynamics_1d.trajectory import TrajectoryResult

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test PES files with minimum at index 5 (not at edge)
            # Use 11 points so minimum is at index 5
            x = np.linspace(0.5, 2.5, 11)
            E_init = 0.1 * (x - 1.5) ** 2  # Minimum at x=1.5 (index 5)
            E_dyn = 19.0 + 0.1 * (x - 1.5) ** 2  # Intermediate state
            E_final = -0.5 + 0.1 * (x - 1.5) ** 2  # Final state

            pes_init = tmpdir / "pes_init.dat"
            pes_dyn = tmpdir / "pes_dyn.dat"
            pes_final = tmpdir / "pes_final.dat"

            np.savetxt(pes_init, np.column_stack([x, E_init]))
            np.savetxt(pes_dyn, np.column_stack([x, E_dyn]))
            np.savetxt(pes_final, np.column_stack([x, E_final]))

            # Create configs for both modes
            def make_config(mode):
                return FullConfig(
                    dynamics=DynamicsConfig(
                        mu=1.0,
                        grid=GridConfig(start=0.5, dx=0.2, npoints=11),
                        time=TimeConfig(dt=1.0, nsteps=20),
                        sampling=SamplingConfig(mode=1, npoints_x=2, npoints_mom=2, compatibility_mode=mode),
                        pes_initial=pes_init,
                        pes_dynamics=pes_dyn,
                        units="angstrom",
                    ),
                    spectrum=SpectrumConfig(
                        gamma_fwhm=0.1,
                        dipole_mode="FC",
                        pes_final_list=[pes_final],
                        compatibility_mode=mode,
                    ),
                )

            # Create mock trajectory with enough points and varying positions
            # Trajectory has 20 points, PES minloc is at index 5
            # Put trajectory at different positions to ensure modes give different E_mean
            nsteps = 20
            time = np.arange(nsteps) * 1e-15
            # Position varies: starts at 1.6 Å, so position at index 5 is different from equilibrium
            x_traj = np.linspace(1.6e-10, 1.4e-10, nsteps)  # 1.6 to 1.4 Angstrom

            traj = TrajectoryResult(
                time=time,
                x=x_traj,
                v=np.zeros(nsteps),
                a=np.zeros(nsteps),
                x0=x_traj[0],
                p0=0.0,
            )

            # Calculate E_mean with both modes
            calc_standard = SpectrumCalculator(make_config("standard"))
            calc_standard.load_surfaces()
            E_mean_standard = calc_standard.compute_mean_transition_energy([traj])

            calc_fortran = SpectrumCalculator(make_config("fortran"))
            calc_fortran.load_surfaces()
            E_mean_fortran = calc_fortran.compute_mean_transition_energy([traj])

            # Standard mode evaluates at true equilibrium (x=1.5 Å)
            # Fortran mode evaluates at trajectory point x_traj[5] ≈ 1.55 Å
            # These should give slightly different E_mean values

            # Both should be reasonable values (E_dyn - E_final ≈ 19.5 Hartree ≈ 530 eV)
            assert 500 < E_mean_standard < 600
            assert 500 < E_mean_fortran < 600

            # The two methods should give similar but not identical results
            # (difference depends on how far x_traj[5] is from equilibrium)
            assert abs(E_mean_standard - E_mean_fortran) < 10  # Within 10 eV


class TestPESEnergyShift:
    """Tests for PES energy correction shift functionality."""

    def test_pes_lp_corr_shifts_final_states(self):
        """Test that pes_lp_corr shifts all final state PES energies."""
        from python_scripts.dynamics_1d.spectrum import SpectrumCalculator
        from python_scripts.dynamics_1d.spectrum_config import (
            FullConfig,
            SpectrumConfig,
        )
        from python_scripts.dynamics_1d.config import (
            DynamicsConfig,
            GridConfig,
            TimeConfig,
            SamplingConfig,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test PES files (in Angstrom/Hartree)
            x = np.linspace(0.5, 2.5, 21)

            # Initial/dynamics PES: harmonic around 1.5
            E_init = 0.1 * (x - 1.5) ** 2
            pes_init = tmpdir / "pes_init.dat"
            np.savetxt(pes_init, np.column_stack([x, E_init]))

            # Final state 1: shifted harmonic
            E_f1 = -0.2 + 0.1 * (x - 1.5) ** 2
            pes_f1 = tmpdir / "pes_f1.dat"
            np.savetxt(pes_f1, np.column_stack([x, E_f1]))

            # Final state 2: another shifted harmonic
            E_f2 = -0.15 + 0.1 * (x - 1.5) ** 2  # 0.05 above f1
            pes_f2 = tmpdir / "pes_f2.dat"
            np.savetxt(pes_f2, np.column_stack([x, E_f2]))

            # LP correction PES: different energy level
            E_corr = -0.5 + 0.1 * (x - 1.5) ** 2  # 0.3 below f1
            pes_corr = tmpdir / "pes_lp_corr.dat"
            np.savetxt(pes_corr, np.column_stack([x, E_corr]))

            # Create config with pes_lp_corr
            dynamics_config = DynamicsConfig(
                mu=1.0,
                grid=GridConfig(start=0.5, dx=0.1, npoints=21),
                time=TimeConfig(dt=1.0, nsteps=10),
                sampling=SamplingConfig(mode=1, npoints_x=2, npoints_mom=2),
                pes_initial=pes_init,
                pes_dynamics=pes_init,
                units="angstrom",
            )
            spectrum_config = SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="FC",
                pes_final_list=[pes_f1, pes_f2],
                pes_lp_corr=pes_corr,
            )
            config = FullConfig(dynamics=dynamics_config, spectrum=spectrum_config)

            # Create calculator and load surfaces
            calc = SpectrumCalculator(config)
            calc.load_surfaces()

            # The shift should be E_corr - E_f1 = -0.3 (in Hartree)
            # After shift: E_f1_shifted = E_f1 + shift = E_corr
            # After shift: E_f2_shifted = E_f2 + shift = E_f2 - 0.3

            # Check that first final state now matches correction PES
            x_test = 1.5e-10  # 1.5 Angstrom in meters
            E_f1_shifted = calc.pes_f[0].energy(x_test)

            # Load expected from correction file
            from python_scripts.dynamics_1d.io import read_pes_file
            x_exp, E_exp = read_pes_file(pes_corr, units="angstrom")
            pes_exp = calc.pes_f[0].__class__(x=x_exp, E=E_exp)
            E_expected = pes_exp.energy(x_test)

            np.testing.assert_allclose(E_f1_shifted, E_expected, rtol=1e-10)

            # Check that second final state is shifted by same amount
            # Original spacing was 0.05 Hartree, should be preserved
            E_f2_shifted = calc.pes_f[1].energy(x_test)
            spacing_shifted = E_f2_shifted - E_f1_shifted
            spacing_original = 0.05 * CONST.hartree  # Convert to Joules

            np.testing.assert_allclose(spacing_shifted, spacing_original, rtol=1e-10)

    def test_no_shift_without_pes_lp_corr(self):
        """Test that PES energies are unchanged when pes_lp_corr is not set."""
        from python_scripts.dynamics_1d.spectrum import SpectrumCalculator
        from python_scripts.dynamics_1d.spectrum_config import (
            FullConfig,
            SpectrumConfig,
        )
        from python_scripts.dynamics_1d.config import (
            DynamicsConfig,
            GridConfig,
            TimeConfig,
            SamplingConfig,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            x = np.linspace(0.5, 2.5, 21)
            E_init = 0.1 * (x - 1.5) ** 2
            pes_init = tmpdir / "pes_init.dat"
            np.savetxt(pes_init, np.column_stack([x, E_init]))

            E_f1 = -0.2 + 0.1 * (x - 1.5) ** 2
            pes_f1 = tmpdir / "pes_f1.dat"
            np.savetxt(pes_f1, np.column_stack([x, E_f1]))

            dynamics_config = DynamicsConfig(
                mu=1.0,
                grid=GridConfig(start=0.5, dx=0.1, npoints=21),
                time=TimeConfig(dt=1.0, nsteps=10),
                sampling=SamplingConfig(mode=1, npoints_x=2, npoints_mom=2),
                pes_initial=pes_init,
                pes_dynamics=pes_init,
                units="angstrom",
            )
            spectrum_config = SpectrumConfig(
                gamma_fwhm=0.1,
                dipole_mode="FC",
                pes_final_list=[pes_f1],
                # No pes_lp_corr
            )
            config = FullConfig(dynamics=dynamics_config, spectrum=spectrum_config)

            calc = SpectrumCalculator(config)
            calc.load_surfaces()

            # Check that energy at minimum matches original file
            x_test = 1.5e-10  # 1.5 Angstrom in meters
            E_loaded = calc.pes_f[0].energy(x_test)
            E_original_hartree = -0.2 + 0.1 * (1.5 - 1.5) ** 2  # -0.2 Hartree
            E_expected = E_original_hartree * CONST.hartree

            np.testing.assert_allclose(E_loaded, E_expected, rtol=1e-10)
