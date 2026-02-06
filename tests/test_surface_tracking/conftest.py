"""Pytest fixtures for surface tracking tests."""

import numpy as np
import pytest


@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return np.random.default_rng(42)


@pytest.fixture
def smooth_2state_3x3():
    """3x3 grid with 2 smoothly varying states, no crossings.

    State 0: low energy, dipole along x.
    State 1: high energy, dipole along y.
    """
    n_x1, n_x2, n_states = 3, 3, 2
    E = np.zeros((n_x1, n_x2, n_states))
    D = np.zeros((n_x1, n_x2, n_states, 3))

    for i in range(n_x1):
        for j in range(n_x2):
            E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
            E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
            D[i, j, 0] = [1.0 + 0.01 * i, 0.05 * j, 0.0]
            D[i, j, 1] = [0.05 * i, 1.0 + 0.01 * j, 0.0]

    return E, D


@pytest.fixture
def crossing_2state_3x3():
    """3x3 grid with 2 states that cross at one point.

    Two surfaces with distinct dipole character (x vs y) that approach
    each other in energy and cross at (2,2). Energy ordering at (2,2)
    swaps the state labels relative to the character-based ordering.
    """
    n_x1, n_x2, n_states = 3, 3, 2
    E = np.zeros((n_x1, n_x2, n_states))
    D = np.zeros((n_x1, n_x2, n_states, 3))

    # "True" surfaces (by character):
    # Surface A: dipole along x, energy rises across grid
    # Surface B: dipole along y, energy falls across grid
    # They cross at (2,2).
    E_A = np.zeros((n_x1, n_x2))
    E_B = np.zeros((n_x1, n_x2))
    D_A = np.zeros((n_x1, n_x2, 3))
    D_B = np.zeros((n_x1, n_x2, 3))

    for i in range(n_x1):
        for j in range(n_x2):
            # Surfaces approach each other: A rises, B falls
            E_A[i, j] = 2.0 + 0.15 * i + 0.15 * j
            E_B[i, j] = 2.5 - 0.15 * i - 0.15 * j
            # Strong, distinct dipole characters
            D_A[i, j] = [1.0, 0.0, 0.0]
            D_B[i, j] = [0.0, 1.0, 0.0]

    # At (2,2): E_A = 2+0.6+0.6 = 3.2, E_B = 4-0.6-0.6 = 2.8
    # So E_A > E_B at (2,2): crossing! Energy ordering puts B first.

    # Build energy-ordered arrays
    for i in range(n_x1):
        for j in range(n_x2):
            if E_A[i, j] <= E_B[i, j]:
                E[i, j, 0] = E_A[i, j]
                E[i, j, 1] = E_B[i, j]
                D[i, j, 0] = D_A[i, j]
                D[i, j, 1] = D_B[i, j]
            else:
                E[i, j, 0] = E_B[i, j]
                E[i, j, 1] = E_A[i, j]
                D[i, j, 0] = D_B[i, j]
                D[i, j, 1] = D_A[i, j]

    return E, D, E_A, E_B, D_A, D_B


@pytest.fixture
def random_permuted_5x5(rng):
    """5x5 grid with 3 smooth states, randomly permuted at each point.

    Returns (E_permuted, D_permuted, E_true, D_true, perms)
    where perms[i,j] is the permutation applied.
    """
    n_x1, n_x2, n_states = 5, 5, 3
    E_true = np.zeros((n_x1, n_x2, n_states))
    D_true = np.zeros((n_x1, n_x2, n_states, 3))

    for i in range(n_x1):
        for j in range(n_x2):
            E_true[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
            E_true[i, j, 1] = 2.0 + 0.01 * i - 0.01 * j
            E_true[i, j, 2] = 4.0 - 0.02 * i + 0.01 * j
            D_true[i, j, 0] = [1.0, 0.1 * j, 0.0]
            D_true[i, j, 1] = [0.0, 1.0, 0.1 * i]
            D_true[i, j, 2] = [0.1 * i, 0.0, 1.0]

    E_perm = np.zeros_like(E_true)
    D_perm = np.zeros_like(D_true)
    perms = np.zeros((n_x1, n_x2, n_states), dtype=int)

    for i in range(n_x1):
        for j in range(n_x2):
            p = rng.permutation(n_states)
            perms[i, j] = p
            E_perm[i, j] = E_true[i, j, p]
            D_perm[i, j] = D_true[i, j, p]

    return E_perm, D_perm, E_true, D_true, perms


@pytest.fixture
def sign_flipped_3x3():
    """3x3 grid with consistent states but random sign flips on dipoles.

    Returns (E, D_flipped, D_true, flips) where flips[i,j,k] = +/-1.
    """
    rng = np.random.default_rng(123)
    n_x1, n_x2, n_states = 3, 3, 2
    E = np.zeros((n_x1, n_x2, n_states))
    D_true = np.zeros((n_x1, n_x2, n_states, 3))

    for i in range(n_x1):
        for j in range(n_x2):
            E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
            E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
            D_true[i, j, 0] = [1.0, 0.1, 0.0]
            D_true[i, j, 1] = [0.0, 0.1, 1.0]

    flips = rng.choice([-1, 1], size=(n_x1, n_x2, n_states))
    D_flipped = D_true * flips[:, :, :, None]

    return E, D_flipped, D_true, flips
