"""Dipole matrix element computation for KH calculations.

Computes transition dipole matrix elements between vibrational states,
supporting Franck-Condon, full dipole, and dipole-at-x0 approximations.
"""

from typing import Tuple

import numpy as np


def compute_fc_overlap(
    psi_a: np.ndarray,
    psi_b: np.ndarray,
    dx: float,
) -> float:
    """Compute Franck-Condon overlap integral.

    <psi_a|psi_b> = integral(psi_a(x) * psi_b(x) dx)

    Args:
        psi_a: First wavefunction on grid
        psi_b: Second wavefunction on grid
        dx: Grid spacing in meters

    Returns:
        Overlap integral (dimensionless)
    """
    return np.sum(psi_a * psi_b) * dx


def compute_dipole_overlap(
    psi_a: np.ndarray,
    psi_b: np.ndarray,
    dipole: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Compute dipole matrix element.

    <psi_a|d|psi_b> = integral(psi_a(x) * d(x) * psi_b(x) dx)

    Args:
        psi_a: First wavefunction on grid
        psi_b: Second wavefunction on grid
        dipole: Dipole moment on grid, shape (npoints, 3)
        dx: Grid spacing in meters

    Returns:
        Dipole matrix element (3,) array in atomic units
    """
    # psi_a * psi_b -> (npoints,)
    # dipole -> (npoints, 3)
    # integrand -> (npoints, 3)
    integrand = (psi_a * psi_b)[:, np.newaxis] * dipole
    return np.sum(integrand, axis=0) * dx


def compute_transition_dipoles_fc_loop(
    vib_i: np.ndarray,
    vib_n: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Loop-based FC overlaps (reference implementation).

    Args:
        vib_i: Initial state eigenvectors, shape (npoints, n_states_i)
        vib_n: Intermediate state eigenvectors, shape (npoints, n_states_n)
        dx: Grid spacing in meters

    Returns:
        FC overlaps, shape (n_states_n, n_states_i)
    """
    n_states_n = vib_n.shape[1]
    n_states_i = vib_i.shape[1]

    overlaps = np.zeros((n_states_n, n_states_i))

    for i in range(n_states_i):
        for n in range(n_states_n):
            overlaps[n, i] = compute_fc_overlap(vib_n[:, n], vib_i[:, i], dx)

    return overlaps


def compute_transition_dipoles_fc(
    vib_i: np.ndarray,
    vib_n: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Compute FC overlaps between initial and intermediate vibrational states.

    Vectorized version using matrix multiply.

    Args:
        vib_i: Initial state eigenvectors, shape (npoints, n_states_i)
        vib_n: Intermediate state eigenvectors, shape (npoints, n_states_n)
        dx: Grid spacing in meters

    Returns:
        FC overlaps, shape (n_states_n, n_states_i)
    """
    # overlaps[n, i] = sum_x vib_n[x, n] * vib_i[x, i] * dx
    return (vib_n.T @ vib_i) * dx


def compute_transition_dipoles_full_loop(
    vib_a: np.ndarray,
    vib_b: np.ndarray,
    dipole: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Loop-based full dipole matrix elements (reference implementation).

    Args:
        vib_a: First set of eigenvectors, shape (npoints, n_states_a)
        vib_b: Second set of eigenvectors, shape (npoints, n_states_b)
        dipole: Dipole moment on grid, shape (npoints, 3)
        dx: Grid spacing in meters

    Returns:
        Dipole matrix elements, shape (n_states_a, n_states_b, 3)
    """
    n_states_a = vib_a.shape[1]
    n_states_b = vib_b.shape[1]

    D = np.zeros((n_states_a, n_states_b, 3))

    for a in range(n_states_a):
        for b in range(n_states_b):
            D[a, b, :] = compute_dipole_overlap(vib_a[:, a], vib_b[:, b], dipole, dx)

    return D


def compute_transition_dipoles_full(
    vib_a: np.ndarray,
    vib_b: np.ndarray,
    dipole: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Compute full dipole matrix elements between vibrational states.

    Vectorized version using einsum: D[a,b,m] = sum_x vib_a[x,a] * vib_b[x,b] * dipole[x,m] * dx

    Args:
        vib_a: First set of eigenvectors, shape (npoints, n_states_a)
        vib_b: Second set of eigenvectors, shape (npoints, n_states_b)
        dipole: Dipole moment on grid, shape (npoints, 3)
        dx: Grid spacing in meters

    Returns:
        Dipole matrix elements, shape (n_states_a, n_states_b, 3)
    """
    return np.einsum('xa,xb,xm->abm', vib_a, vib_b, dipole) * dx


def compute_dipole_matrix_elements(
    vib_i: np.ndarray,
    vib_n: np.ndarray,
    vib_f: np.ndarray,
    dipole_on_grid: np.ndarray,
    mode: str,
    dx: float,
    x: np.ndarray = None,
    x0: float = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute dipole matrix elements for KH calculation.

    Computes:
    - D_ni: transition dipole from initial to intermediate states
    - D_fn: transition dipole from intermediate to final states

    Args:
        vib_i: Initial state eigenvectors, shape (npoints, n_states_i)
        vib_n: Intermediate state eigenvectors, shape (npoints, n_states_n)
        vib_f: Final state eigenvectors, shape (npoints, n_states_f)
        dipole_on_grid: Dipole moment on grid, shape (npoints, 3)
        mode: Calculation mode - "FC", "DIPOLE", or "DIPOLE_X0"
        dx: Grid spacing in meters
        x: Grid positions in meters (needed for DIPOLE_X0 mode)
        x0: Reference position for DIPOLE_X0 mode (meters)

    Returns:
        D_ni: shape (n_states_n, n_states_i, 3) - initial to intermediate
        D_fn: shape (n_states_f, n_states_n, 3) - intermediate to final
    """
    n_states_i = vib_i.shape[1]
    n_states_n = vib_n.shape[1]
    n_states_f = vib_f.shape[1]

    if mode == "FC":
        # Franck-Condon approximation: use overlaps with unit dipole in z-direction
        # This is the simplest approximation where dipole is assumed constant

        # Initial to intermediate: D_ni = <psi_n|psi_i> * (0, 0, 1)
        fc_ni = compute_transition_dipoles_fc(vib_i, vib_n, dx)
        D_ni = np.zeros((n_states_n, n_states_i, 3))
        D_ni[:, :, 2] = fc_ni  # Put overlap in z-component

        # Intermediate to final: D_fn = <psi_f|psi_n> * d(x)
        # For FC mode, we still use the dipole for the emission step
        D_fn = compute_transition_dipoles_full(vib_f, vib_n, dipole_on_grid, dx)

    elif mode == "DIPOLE":
        # Full dipole: use actual dipole surface for both transitions
        # For absorption (i->n), we might need a different dipole surface
        # Here we assume the same dipole surface for simplicity
        # (In Fortran code, absorption uses FC and emission uses full dipole)

        # For KH non-resonant, typically:
        # - D_ni uses FC (absorption)
        # - D_fn uses full dipole (emission)
        fc_ni = compute_transition_dipoles_fc(vib_i, vib_n, dx)
        D_ni = np.zeros((n_states_n, n_states_i, 3))
        D_ni[:, :, 2] = fc_ni

        D_fn = compute_transition_dipoles_full(vib_f, vib_n, dipole_on_grid, dx)

    elif mode == "DIPOLE_X0":
        # Dipole at reference position times FC overlap
        if x is None or x0 is None:
            raise ValueError("DIPOLE_X0 mode requires x and x0 arguments")

        # Find dipole at x0 by interpolation
        from scipy.interpolate import interp1d

        d_interp = interp1d(x, dipole_on_grid, axis=0, kind="cubic")
        d_x0 = d_interp(x0)  # Shape: (3,)

        # D_ni = <psi_n|psi_i> * d(x0)
        fc_ni = compute_transition_dipoles_fc(vib_i, vib_n, dx)
        D_ni = np.zeros((n_states_n, n_states_i, 3))
        for m in range(3):
            D_ni[:, :, m] = fc_ni * d_x0[m]

        # D_fn = <psi_f|psi_n> * d(x0)
        fc_fn = compute_transition_dipoles_fc(vib_n, vib_f, dx)
        D_fn = np.zeros((n_states_f, n_states_n, 3))
        for m in range(3):
            D_fn[:, :, m] = fc_fn * d_x0[m]

    else:
        raise ValueError(f"Unknown dipole mode: {mode}")

    return D_ni, D_fn
