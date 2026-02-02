"""Tests for KH scattering amplitude computation."""

import numpy as np
import pytest
from pathlib import Path

from python_scripts.kh_1d.amplitude import (
    compute_amplitude_F_res,
    compute_amplitude_F_nonres,
    compute_amplitude_F,
    compute_cross_section_from_F,
    compute_XES_nonres,
    compute_XES_per_final_state,
)


class TestAmplitudeBasic:
    """Basic tests for amplitude computation."""

    def test_single_intermediate_state_lorentzian(self):
        """Single intermediate state should give Lorentzian lineshape."""
        # Single intermediate state
        E_i = 0.0
        E_n = np.array([520.0])  # One intermediate state at 520 eV
        E_f = 0.5  # Final state at 0.5 eV
        gamma = 0.2  # HWHM

        # Unit dipoles in z-direction
        D_ni = np.array([[0.0, 0.0, 1.0]])  # (1, 3)
        D_fn = np.array([[0.0, 0.0, 1.0]])  # (1, 3)

        # Peak should be at omega = E_n - E_f = 519.5 eV
        omega_peak = E_n[0] - E_f

        # Compute amplitude at peak
        F_peak = compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega_peak, gamma)

        # At resonance, denominator is purely imaginary: i*gamma
        # F = D_fn * D_ni / (i*gamma) = 1 / (i*0.2) = -5i
        expected = 1.0 / (1j * gamma)
        assert abs(F_peak[2, 2] - expected) < 1e-10

    def test_amplitude_shape(self):
        """Amplitude should be 3x3 complex array."""
        E_n = np.array([520.0, 521.0, 522.0])
        D_ni = np.random.randn(3, 3)
        D_fn = np.random.randn(3, 3)

        F = compute_amplitude_F_res(0.0, E_n, 0.0, D_ni, D_fn, 520.0, 0.2)

        assert F.shape == (3, 3)
        assert F.dtype == np.complex128

    def test_nonres_no_pole(self):
        """Non-resonant amplitude should have no imaginary pole."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = 0.5
        D_ni = np.array([[0.0, 0.0, 1.0]])
        D_fn = np.array([[0.0, 0.0, 1.0]])

        # Non-resonant term has real denominator: omega + (E_n - E_i)
        # At omega = 0: denom = 520
        F_nonres = compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, 0.0, 0.2)

        # Should be real (no imaginary part in denominator)
        expected = -1.0 / (0.0 + (520.0 - 0.0))
        assert abs(F_nonres[2, 2] - expected) < 1e-10
        assert abs(F_nonres[2, 2].imag) < 1e-15


class TestCrossSection:
    """Tests for cross-section computation."""

    def test_cross_section_positive(self):
        """Cross-section should always be positive."""
        F = np.random.randn(3, 3) + 1j * np.random.randn(3, 3)
        sigma = compute_cross_section_from_F(F)
        assert sigma > 0

    def test_cross_section_from_real_F(self):
        """For real F, cross-section is sum of squares."""
        F = np.array([[1.0, 2.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 4.0]])
        sigma = compute_cross_section_from_F(F)
        expected = 1 + 4 + 4 + 9 + 16  # Sum of |F_ij|^2
        assert abs(sigma - expected) < 1e-10


class TestXESSpectrum:
    """Tests for full XES spectrum computation."""

    def test_spectrum_peak_position(self):
        """Spectrum peak should be at expected energy."""
        E_i = 0.0
        E_n = np.array([520.0])  # Single intermediate
        E_f = np.array([0.0])  # Single final vibrational state
        gamma = 0.1

        D_ni = np.array([[0.0, 0.0, 1.0]])
        D_fn = np.array([[[0.0, 0.0, 1.0]]])  # (1, 1, 3)

        # Peak at omega = E_n - E_f = 520 eV
        omega_grid = np.linspace(515, 525, 1001)

        sigma = compute_XES_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma)

        # Find peak
        peak_idx = np.argmax(sigma)
        peak_omega = omega_grid[peak_idx]

        # Peak should be near 520 eV
        assert abs(peak_omega - 520.0) < 0.1

    def test_spectrum_fwhm(self):
        """Spectrum FWHM should match input broadening."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = np.array([0.0])
        gamma = 0.2  # HWHM, so FWHM = 0.4

        D_ni = np.array([[0.0, 0.0, 1.0]])
        D_fn = np.array([[[0.0, 0.0, 1.0]]])

        omega_grid = np.linspace(515, 525, 10001)

        sigma = compute_XES_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma)

        # Find FWHM
        max_val = np.max(sigma)
        half_max = max_val / 2

        # Find where spectrum crosses half-maximum
        above_half = sigma > half_max
        left_idx = np.argmax(above_half)
        right_idx = len(above_half) - np.argmax(above_half[::-1]) - 1

        fwhm = omega_grid[right_idx] - omega_grid[left_idx]

        # FWHM should be approximately 2*gamma for Lorentzian
        expected_fwhm = 2 * gamma
        assert abs(fwhm - expected_fwhm) < 0.02  # Within 2% of grid spacing

    def test_multiple_intermediate_states(self):
        """Multiple intermediate states should produce multiple peaks."""
        E_i = 0.0
        E_n = np.array([518.0, 520.0, 522.0])  # Three intermediate states
        E_f = np.array([0.0])
        gamma = 0.15

        # Different transition strengths
        D_ni = np.array([[0, 0, 1.0], [0, 0, 0.8], [0, 0, 0.5]])
        D_fn = np.array([[[0, 0, 1.0], [0, 0, 0.8], [0, 0, 0.5]]])

        omega_grid = np.linspace(515, 525, 2001)

        sigma = compute_XES_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma)

        # Should have 3 peaks
        # Find local maxima (simple approach)
        peaks = []
        for i in range(1, len(sigma) - 1):
            if sigma[i] > sigma[i - 1] and sigma[i] > sigma[i + 1]:
                peaks.append(omega_grid[i])

        # Should find peaks near 518, 520, 522
        assert len(peaks) >= 3
        peak_positions = sorted(peaks[-3:])  # Take largest 3

        # Check peak positions
        expected_peaks = [518.0, 520.0, 522.0]
        for p_found, p_expected in zip(peak_positions, expected_peaks):
            assert abs(p_found - p_expected) < 0.3

    def test_per_final_state_spectrum(self):
        """Test spectrum resolved by final vibrational state."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = np.array([0.0, 0.2, 0.4])  # Three final vibrational states
        gamma = 0.1

        D_ni = np.array([[0, 0, 1.0]])
        D_fn = np.array([
            [[0, 0, 1.0]],
            [[0, 0, 0.7]],
            [[0, 0, 0.4]],
        ])

        omega_grid = np.linspace(515, 525, 501)

        sigma_f = compute_XES_per_final_state(E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma)

        # Shape should be (3, 501)
        assert sigma_f.shape == (3, 501)

        # Each spectrum should have a peak at different positions
        peaks = []
        for j in range(3):
            peak_idx = np.argmax(sigma_f[j])
            peaks.append(omega_grid[peak_idx])

        # Peaks should be at omega = E_n - E_f[j] = 520, 519.8, 519.6
        assert abs(peaks[0] - 520.0) < 0.1
        assert abs(peaks[1] - 519.8) < 0.1
        assert abs(peaks[2] - 519.6) < 0.1


class TestAmplitudeSumRules:
    """Tests for sum rules and physical constraints."""

    def test_resonant_dominates_near_pole(self):
        """Near resonance, resonant term should dominate non-resonant."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = 0.5
        gamma = 0.1

        D_ni = np.array([[0, 0, 1.0]])
        D_fn = np.array([[0, 0, 1.0]])

        omega_peak = E_n[0] - E_f

        F_res = compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega_peak, gamma)
        F_nonres = compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_peak, gamma)

        # At resonance, |F_res| >> |F_nonres|
        assert np.abs(F_res[2, 2]) > 10 * np.abs(F_nonres[2, 2])

    def test_combined_amplitude(self):
        """Combined amplitude should sum resonant and non-resonant."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = 0.5
        gamma = 0.1
        omega = 519.0

        D_ni = np.array([[0.1, 0.2, 0.3]])
        D_fn = np.array([[0.4, 0.5, 0.6]])

        F_res = compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega, gamma)
        F_nonres = compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega, gamma)
        F_total = compute_amplitude_F(
            E_i, E_n, E_f, D_ni, D_fn, omega, gamma,
            include_resonant=True, include_nonresonant=True
        )

        np.testing.assert_allclose(F_total, F_res + F_nonres, rtol=1e-10)

    def test_resonant_only_flag(self):
        """include_resonant=False should exclude resonant term."""
        E_i = 0.0
        E_n = np.array([520.0])
        E_f = 0.5
        gamma = 0.1
        omega = 519.0

        D_ni = np.array([[0.1, 0.2, 0.3]])
        D_fn = np.array([[0.4, 0.5, 0.6]])

        F_nonres = compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega, gamma)
        F_total = compute_amplitude_F(
            E_i, E_n, E_f, D_ni, D_fn, omega, gamma,
            include_resonant=False, include_nonresonant=True
        )

        np.testing.assert_allclose(F_total, F_nonres, rtol=1e-10)
