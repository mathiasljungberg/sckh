"""Tests for vibrational module (Fourier grid solver)."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.vibrational import (
    build_kinetic_matrix_fast,
    build_kinetic_matrix_vectorized,
    solve_vibrational,
    harmonic_eigenvalues,
    harmonic_omega_from_k,
)
from python_scripts.dynamics_1d.constants import CONST


class TestBuildKineticMatrix:
    """Tests for kinetic energy matrix construction."""

    def test_matrix_is_symmetric(self, position_grid):
        """Kinetic energy matrix should be symmetric."""
        nstates = len(position_grid)
        dx = position_grid[1] - position_grid[0]
        mass = 1.0 * CONST.u

        H_kin = build_kinetic_matrix_vectorized(nstates, dx, mass)

        np.testing.assert_allclose(H_kin, H_kin.T, rtol=1e-12)

    def test_fast_and_vectorized_match(self, position_grid):
        """Fast and vectorized implementations should give same result."""
        nstates = len(position_grid)
        dx = position_grid[1] - position_grid[0]
        mass = 1.0 * CONST.u

        H_fast = build_kinetic_matrix_fast(nstates, dx, mass)
        H_vec = build_kinetic_matrix_vectorized(nstates, dx, mass)

        np.testing.assert_allclose(H_fast, H_vec, rtol=1e-10)

    def test_kinetic_matrix_toeplitz_structure(self, position_grid):
        """Kinetic matrix should have Toeplitz structure (depends only on |i-j|)."""
        nstates = len(position_grid)
        dx = position_grid[1] - position_grid[0]
        mass = 1.0 * CONST.u

        H_kin = build_kinetic_matrix_vectorized(nstates, dx, mass)

        # Check that H[i,j] = H[i+k, j+k] for all valid k
        for offset in range(1, min(5, nstates - 1)):
            for i in range(nstates - offset):
                np.testing.assert_allclose(
                    H_kin[i, i + offset],
                    H_kin[0, offset],
                    rtol=1e-12,
                )


class TestSolveVibrational:
    """Tests for vibrational Schrödinger equation solver."""

    def test_harmonic_oscillator_eigenvalues(self):
        """Eigenvalues should match E_n = ℏω(n + 1/2) for harmonic oscillator."""
        # Use a fine grid for accuracy
        nstates = 201
        x_range = 2.0e-10  # 2 Angstrom total range
        x = np.linspace(-x_range / 2, x_range / 2, nstates)

        mass = 1.0 * CONST.u
        k = 500.0  # N/m
        omega = np.sqrt(k / mass)

        # Harmonic potential
        V = 0.5 * k * x**2

        eigenvalues, eigenvectors = solve_vibrational(x, V, mass, n_states=5)

        # Expected eigenvalues
        E_expected = harmonic_eigenvalues(omega, 5)

        # Compare first few eigenvalues (numerical accuracy decreases for higher states)
        np.testing.assert_allclose(eigenvalues[:3], E_expected[:3], rtol=1e-3)

    def test_harmonic_oscillator_spacing(self):
        """Energy spacing should be constant ℏω for harmonic oscillator."""
        nstates = 201
        x = np.linspace(-1.5e-10, 1.5e-10, nstates)

        mass = 1.0 * CONST.u
        k = 400.0
        omega = np.sqrt(k / mass)

        V = 0.5 * k * x**2
        eigenvalues, _ = solve_vibrational(x, V, mass, n_states=6)

        # Energy spacing
        dE = np.diff(eigenvalues)
        dE_expected = CONST.hbar * omega

        # First few spacings should be approximately constant
        np.testing.assert_allclose(dE[:3], dE_expected, rtol=1e-2)

    def test_eigenvectors_normalized(self):
        """Eigenvectors should be normalized."""
        nstates = 101
        x = np.linspace(-1e-10, 1e-10, nstates)
        dx = x[1] - x[0]

        mass = 1.0 * CONST.u
        V = 0.5 * 500 * x**2

        _, eigenvectors = solve_vibrational(x, V, mass, n_states=3)

        for i in range(3):
            norm = np.sum(eigenvectors[:, i] ** 2) * dx
            np.testing.assert_allclose(norm, 1.0, rtol=1e-6)

    def test_eigenvectors_orthogonal(self):
        """Eigenvectors should be orthogonal."""
        nstates = 101
        x = np.linspace(-1e-10, 1e-10, nstates)
        dx = x[1] - x[0]

        mass = 1.0 * CONST.u
        V = 0.5 * 500 * x**2

        _, eigenvectors = solve_vibrational(x, V, mass, n_states=5)

        # Check orthogonality
        for i in range(5):
            for j in range(i + 1, 5):
                overlap = np.sum(eigenvectors[:, i] * eigenvectors[:, j]) * dx
                np.testing.assert_allclose(overlap, 0.0, atol=1e-10)

    def test_ground_state_no_nodes(self):
        """Ground state wavefunction should have no nodes."""
        nstates = 101
        x = np.linspace(-1e-10, 1e-10, nstates)

        mass = 1.0 * CONST.u
        V = 0.5 * 500 * x**2

        _, eigenvectors = solve_vibrational(x, V, mass, n_states=1)
        psi_0 = eigenvectors[:, 0]

        # Ground state should be all positive or all negative (no sign changes)
        # Count sign changes, excluding near-zero values at boundaries
        interior = psi_0[10:-10]
        sign_changes = np.sum(np.diff(np.sign(interior)) != 0)
        assert sign_changes == 0


class TestHarmonicHelperFunctions:
    """Tests for harmonic oscillator helper functions."""

    def test_harmonic_eigenvalues(self):
        """Test analytical harmonic eigenvalue formula."""
        omega = 1e14  # rad/s
        E = harmonic_eigenvalues(omega, 4)

        E_expected = CONST.hbar * omega * np.array([0.5, 1.5, 2.5, 3.5])
        np.testing.assert_allclose(E, E_expected, rtol=1e-12)

    def test_harmonic_omega_from_k(self):
        """Test omega = sqrt(k/m) calculation."""
        k = 500.0  # N/m
        mass = 1.0 * CONST.u

        omega = harmonic_omega_from_k(k, mass)
        omega_expected = np.sqrt(k / mass)

        np.testing.assert_allclose(omega, omega_expected, rtol=1e-12)
