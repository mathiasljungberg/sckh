"""Tests for dipole matrix element computation."""

import numpy as np
import pytest
from pathlib import Path

from python_scripts.kh_1d.dipole_matrix import (
    compute_fc_overlap,
    compute_dipole_overlap,
    compute_transition_dipoles_fc,
    compute_transition_dipoles_full,
    compute_dipole_matrix_elements,
)
from python_scripts.dynamics_1d.vibrational import solve_vibrational
from python_scripts.dynamics_1d.pes import create_harmonic_pes
from python_scripts.dynamics_1d.constants import CONST


class TestFCOverlap:
    """Tests for Franck-Condon overlap computation."""

    def test_self_overlap_is_unity(self, harmonic_vibrational_states):
        """Wavefunction should be normalized: <psi|psi> = 1."""
        eigenvalues, eigenvectors, x = harmonic_vibrational_states
        dx = x[1] - x[0]

        # Test first few states
        for i in range(min(5, eigenvectors.shape[1])):
            psi = eigenvectors[:, i]
            overlap = compute_fc_overlap(psi, psi, dx)
            assert abs(overlap - 1.0) < 1e-10, f"State {i} not normalized: {overlap}"

    def test_orthogonality(self, harmonic_vibrational_states):
        """Different states should be orthogonal: <psi_i|psi_j> = 0."""
        eigenvalues, eigenvectors, x = harmonic_vibrational_states
        dx = x[1] - x[0]

        # Test orthogonality between first few states
        n_test = min(5, eigenvectors.shape[1])
        for i in range(n_test):
            for j in range(i + 1, n_test):
                psi_i = eigenvectors[:, i]
                psi_j = eigenvectors[:, j]
                overlap = compute_fc_overlap(psi_i, psi_j, dx)
                assert abs(overlap) < 1e-8, f"States {i},{j} not orthogonal: {overlap}"

    def test_overlap_different_potentials(self, harmonic_grid, hydrogen_mass):
        """Overlaps between states of shifted potentials should be non-trivial."""
        x, dx = harmonic_grid
        k = 500.0  # Force constant

        # Two harmonic oscillators with different equilibrium positions
        pes1 = create_harmonic_pes(x, x0=0.96e-10, k=k)
        pes2 = create_harmonic_pes(x, x0=1.00e-10, k=k)  # Shifted by 0.04 Angstrom

        V1 = pes1.energy(x)
        V2 = pes2.energy(x)

        eig1, vec1 = solve_vibrational(x, V1, hydrogen_mass, n_states=10)
        eig2, vec2 = solve_vibrational(x, V2, hydrogen_mass, n_states=10)

        # Ground state overlap should be less than 1 (shifted oscillators)
        overlap_00 = compute_fc_overlap(vec1[:, 0], vec2[:, 0], dx)
        assert overlap_00 < 1.0
        assert overlap_00 > 0.5  # But still significant


class TestDipoleOverlap:
    """Tests for dipole matrix element computation."""

    def test_constant_dipole_gives_fc(self, harmonic_vibrational_states):
        """With constant dipole, dipole matrix element equals FC overlap * dipole."""
        eigenvalues, eigenvectors, x = harmonic_vibrational_states
        dx = x[1] - x[0]

        d_const = 2.5  # Constant dipole value
        dipole = np.zeros((len(x), 3))
        dipole[:, 2] = d_const  # z-component only

        psi_0 = eigenvectors[:, 0]
        psi_1 = eigenvectors[:, 1]

        fc_01 = compute_fc_overlap(psi_0, psi_1, dx)
        d_01 = compute_dipole_overlap(psi_0, psi_1, dipole, dx)

        # z-component should equal FC * d_const
        assert abs(d_01[2] - fc_01 * d_const) < 1e-10
        # x and y components should be zero
        assert abs(d_01[0]) < 1e-12
        assert abs(d_01[1]) < 1e-12

    def test_symmetric_matrix_element(self, harmonic_vibrational_states):
        """Dipole matrix element should be symmetric for real wavefunctions."""
        eigenvalues, eigenvectors, x = harmonic_vibrational_states
        dx = x[1] - x[0]

        # Linear dipole: d(x) = x
        dipole = np.zeros((len(x), 3))
        dipole[:, 2] = x * 1e10  # Scale for reasonable values

        psi_0 = eigenvectors[:, 0]
        psi_1 = eigenvectors[:, 1]

        d_01 = compute_dipole_overlap(psi_0, psi_1, dipole, dx)
        d_10 = compute_dipole_overlap(psi_1, psi_0, dipole, dx)

        np.testing.assert_allclose(d_01, d_10, rtol=1e-10)


class TestTransitionDipoles:
    """Tests for batch computation of transition dipoles."""

    def test_fc_matrix_shape(self, harmonic_grid, hydrogen_mass):
        """FC overlap matrix should have correct shape."""
        x, dx = harmonic_grid
        k = 500.0

        pes1 = create_harmonic_pes(x, x0=0.96e-10, k=k)
        pes2 = create_harmonic_pes(x, x0=1.00e-10, k=k)

        _, vec1 = solve_vibrational(x, pes1.energy(x), hydrogen_mass, n_states=5)
        _, vec2 = solve_vibrational(x, pes2.energy(x), hydrogen_mass, n_states=10)

        fc_matrix = compute_transition_dipoles_fc(vec1, vec2, dx)

        assert fc_matrix.shape == (10, 5)

    def test_full_dipole_matrix_shape(self, harmonic_grid, hydrogen_mass):
        """Full dipole matrix should have correct shape."""
        x, dx = harmonic_grid
        k = 500.0

        pes = create_harmonic_pes(x, x0=0.96e-10, k=k)
        _, vec = solve_vibrational(x, pes.energy(x), hydrogen_mass, n_states=5)

        dipole = np.zeros((len(x), 3))
        dipole[:, 2] = 1.0

        D = compute_transition_dipoles_full(vec, vec, dipole, dx)

        assert D.shape == (5, 5, 3)

    def test_fc_matrix_orthonormality(self, harmonic_vibrational_states):
        """FC matrix for same potential should give identity."""
        eigenvalues, eigenvectors, x = harmonic_vibrational_states
        dx = x[1] - x[0]

        n_states = min(5, eigenvectors.shape[1])
        vec = eigenvectors[:, :n_states]

        fc_matrix = compute_transition_dipoles_fc(vec, vec, dx)

        # Should be close to identity
        np.testing.assert_allclose(fc_matrix, np.eye(n_states), atol=1e-8)


class TestDipoleMatrixElements:
    """Tests for the combined dipole matrix element computation."""

    def test_fc_mode_shapes(self, harmonic_grid, hydrogen_mass):
        """Test shapes in FC mode."""
        x, dx = harmonic_grid
        k = 500.0

        pes_i = create_harmonic_pes(x, x0=0.96e-10, k=k)
        pes_n = create_harmonic_pes(x, x0=1.00e-10, k=k)
        pes_f = create_harmonic_pes(x, x0=0.98e-10, k=k)

        _, vec_i = solve_vibrational(x, pes_i.energy(x), hydrogen_mass, n_states=1)
        _, vec_n = solve_vibrational(x, pes_n.energy(x), hydrogen_mass, n_states=10)
        _, vec_f = solve_vibrational(x, pes_f.energy(x), hydrogen_mass, n_states=8)

        dipole = np.zeros((len(x), 3))
        dipole[:, 2] = 1.0

        D_ni, D_fn = compute_dipole_matrix_elements(
            vib_i=vec_i,
            vib_n=vec_n,
            vib_f=vec_f,
            dipole_on_grid=dipole,
            mode="FC",
            dx=dx,
        )

        # D_ni: (n_states_n, n_states_i, 3) = (10, 1, 3)
        assert D_ni.shape == (10, 1, 3)

        # D_fn: (n_states_f, n_states_n, 3) = (8, 10, 3)
        assert D_fn.shape == (8, 10, 3)

    def test_dipole_mode_consistency(self, harmonic_grid, hydrogen_mass):
        """DIPOLE mode should give same D_ni as FC mode for constant dipole."""
        x, dx = harmonic_grid
        k = 500.0

        pes_i = create_harmonic_pes(x, x0=0.96e-10, k=k)
        pes_n = create_harmonic_pes(x, x0=1.00e-10, k=k)
        pes_f = create_harmonic_pes(x, x0=0.98e-10, k=k)

        _, vec_i = solve_vibrational(x, pes_i.energy(x), hydrogen_mass, n_states=1)
        _, vec_n = solve_vibrational(x, pes_n.energy(x), hydrogen_mass, n_states=5)
        _, vec_f = solve_vibrational(x, pes_f.energy(x), hydrogen_mass, n_states=5)

        # Constant dipole
        dipole = np.ones((len(x), 3)) * np.array([0.0, 0.0, 1.0])

        D_ni_fc, D_fn_fc = compute_dipole_matrix_elements(
            vib_i=vec_i,
            vib_n=vec_n,
            vib_f=vec_f,
            dipole_on_grid=dipole,
            mode="FC",
            dx=dx,
        )

        D_ni_dip, D_fn_dip = compute_dipole_matrix_elements(
            vib_i=vec_i,
            vib_n=vec_n,
            vib_f=vec_f,
            dipole_on_grid=dipole,
            mode="DIPOLE",
            dx=dx,
        )

        # D_ni should be the same (both use FC for absorption)
        np.testing.assert_allclose(D_ni_fc, D_ni_dip, rtol=1e-10)

        # D_fn should also be the same for constant dipole
        np.testing.assert_allclose(D_fn_fc, D_fn_dip, rtol=1e-10)
