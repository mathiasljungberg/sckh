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
        """Frequency spacing should match 1/T."""
        nsteps = 128
        dt = 1e-15
        time = np.arange(nsteps) * dt
        E_mean = 0.0

        omega = get_frequency_grid(time, E_mean)

        # Spacing in omega should be 2*pi*hbar / (T * eV)
        T = time[-1] + dt
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
