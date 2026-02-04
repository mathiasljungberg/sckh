"""Kramers-Heisenberg scattering amplitude computation.

Implements the KH scattering amplitude F for resonant and non-resonant
X-ray emission spectroscopy (XES).

Based on src/m_KH_utils.F90:1108-1185.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Loop-based reference implementations
# ---------------------------------------------------------------------------

def compute_amplitude_F_res_loop(
    E_i: float,
    E_n: np.ndarray,
    E_f: float,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega: float,
    gamma: float,
) -> np.ndarray:
    """Loop-based resonant KH amplitude (reference implementation)."""
    n_states_n = len(E_n)
    F = np.zeros((3, 3), dtype=complex)

    for n in range(n_states_n):
        denom = omega - (E_n[n] - E_f) + 1j * gamma
        for m1 in range(3):
            for m2 in range(3):
                F[m1, m2] += D_fn[n, m1] * D_ni[n, m2] / denom

    return F


def compute_amplitude_F_nonres_loop(
    E_i: float,
    E_n: np.ndarray,
    E_f: float,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega: float,
    gamma: float,
) -> np.ndarray:
    """Loop-based non-resonant KH amplitude (reference implementation)."""
    n_states_n = len(E_n)
    F = np.zeros((3, 3), dtype=complex)

    for n in range(n_states_n):
        denom = omega + (E_n[n] - E_i)
        if abs(denom) < 1e-12:
            continue
        for m1 in range(3):
            for m2 in range(3):
                F[m1, m2] -= D_fn[n, m1] * D_ni[n, m2] / denom

    return F


def compute_XES_per_final_state_loop(
    E_i: float,
    E_n: np.ndarray,
    E_f: np.ndarray,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega_grid: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Loop-based per-final-state XES spectrum (reference implementation)."""
    n_omega = len(omega_grid)
    n_states_f = len(E_f)
    n_states_n = len(E_n)

    sigma_f = np.zeros((n_states_f, n_omega))

    for f in range(n_states_f):
        for i_omega in range(n_omega):
            omega = omega_grid[i_omega]

            F = np.zeros((3, 3), dtype=complex)
            for n in range(n_states_n):
                denom = omega - (E_n[n] - E_f[f]) + 1j * gamma
                for m1 in range(3):
                    for m2 in range(3):
                        F[m1, m2] += D_fn[f, n, m1] * D_ni[n, m2] / denom

            sigma_f[f, i_omega] = np.sum(np.abs(F) ** 2)

    return sigma_f


# ---------------------------------------------------------------------------
# Vectorized implementations
# ---------------------------------------------------------------------------

def compute_amplitude_F_res(
    E_i: float,
    E_n: np.ndarray,
    E_f: float,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega: float,
    gamma: float,
) -> np.ndarray:
    """Compute resonant KH scattering amplitude.

    F_{m1,m2}(omega) = sum_n [D_fn^{m1} * D_ni^{m2}] / [omega - (E_n - E_f) + i*gamma]

    Args:
        E_i: Initial state energy (eV) - not used in resonant term
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energy (eV)
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_n, 3)
        omega: Output photon energy (eV)
        gamma: Broadening half-width (eV)

    Returns:
        F: Scattering amplitude (3, 3) complex array
    """
    # denom shape: (n_states_n,)
    denom = omega - (E_n - E_f) + 1j * gamma
    weights = 1.0 / denom  # (n_states_n,)
    # F[m1, m2] = sum_n weights[n] * D_fn[n, m1] * D_ni[n, m2]
    F = np.einsum('n,nm,np->mp', weights, D_fn, D_ni)
    return F


def compute_amplitude_F_nonres(
    E_i: float,
    E_n: np.ndarray,
    E_f: float,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega: float,
    gamma: float,
) -> np.ndarray:
    """Compute non-resonant KH scattering amplitude.

    F_{m1,m2}(omega) = -sum_n [D_fn^{m1} * D_ni^{m2}] / [omega + (E_n - E_i)]

    Args:
        E_i: Initial state energy (eV)
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energy (eV) - not used in non-resonant term
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_n, 3)
        omega: Output photon energy (eV)
        gamma: Broadening half-width (eV) - not used in standard non-resonant

    Returns:
        F: Scattering amplitude (3, 3) complex array
    """
    denom = omega + (E_n - E_i)
    # Mask out near-zero denominators
    mask = np.abs(denom) >= 1e-12
    weights = np.zeros_like(denom, dtype=complex)
    weights[mask] = -1.0 / denom[mask]
    F = np.einsum('n,nm,np->mp', weights, D_fn, D_ni)
    return F


def compute_amplitude_F(
    E_i: float,
    E_n: np.ndarray,
    E_f: float,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega: float,
    gamma: float,
    include_resonant: bool = True,
    include_nonresonant: bool = True,
) -> np.ndarray:
    """Compute total KH scattering amplitude.

    Combines resonant and non-resonant contributions.

    Args:
        E_i: Initial state energy (eV)
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energy (eV)
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_n, 3)
        omega: Output photon energy (eV)
        gamma: Broadening half-width (eV)
        include_resonant: Include resonant contribution
        include_nonresonant: Include non-resonant contribution

    Returns:
        F: Scattering amplitude (3, 3) complex array
    """
    F = np.zeros((3, 3), dtype=complex)

    if include_resonant:
        F += compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega, gamma)

    if include_nonresonant:
        F += compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega, gamma)

    return F


def compute_cross_section_from_F(F: np.ndarray) -> float:
    """Compute cross-section from scattering amplitude.

    sigma = sum_{m1,m2} |F_{m1,m2}|^2

    Args:
        F: Scattering amplitude (3, 3) complex array

    Returns:
        Cross-section (arbitrary units)
    """
    return np.sum(np.abs(F) ** 2)


def compute_XES_nonres(
    E_i: float,
    E_n: np.ndarray,
    E_f: np.ndarray,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega_grid: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Compute non-resonant XES spectrum for a single final electronic state.

    sigma(omega) = sum_f |F_f(omega)|^2

    Args:
        E_i: Initial state energy (eV) - single vibrational ground state
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energies (eV), shape (n_states_f,)
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_f, n_states_n, 3)
        omega_grid: Output frequency grid (eV), shape (n_omega,)
        gamma: Broadening half-width (eV)

    Returns:
        sigma: Cross-section array, shape (n_omega,)
    """
    # Use per-final-state and sum
    sigma_f = compute_XES_per_final_state(
        E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma
    )
    return np.sum(sigma_f, axis=0)


def compute_XES_nonres_vectorized(
    E_i: float,
    E_n: np.ndarray,
    E_f: np.ndarray,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega_grid: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Vectorized non-resonant XES spectrum computation.

    Alias for compute_XES_nonres (now always vectorized).

    Args:
        E_i: Initial state energy (eV)
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energies (eV), shape (n_states_f,)
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_f, n_states_n, 3)
        omega_grid: Output frequency grid (eV), shape (n_omega,)
        gamma: Broadening half-width (eV)

    Returns:
        sigma: Cross-section array, shape (n_omega,)
    """
    return compute_XES_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_grid, gamma)


def compute_XES_per_final_state(
    E_i: float,
    E_n: np.ndarray,
    E_f: np.ndarray,
    D_ni: np.ndarray,
    D_fn: np.ndarray,
    omega_grid: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Compute XES spectrum resolved by final vibrational state.

    Vectorized over the omega grid for each final state.

    Args:
        E_i: Initial state energy (eV)
        E_n: Intermediate state energies (eV), shape (n_states_n,)
        E_f: Final state energies (eV), shape (n_states_f,)
        D_ni: Transition dipoles initial->intermediate, shape (n_states_n, 3)
        D_fn: Transition dipoles intermediate->final, shape (n_states_f, n_states_n, 3)
        omega_grid: Output frequency grid (eV), shape (n_omega,)
        gamma: Broadening half-width (eV)

    Returns:
        sigma_f: Cross-section per final state, shape (n_states_f, n_omega)
    """
    n_omega = len(omega_grid)
    n_states_f = len(E_f)

    sigma_f = np.zeros((n_states_f, n_omega))

    for f in range(n_states_f):
        # denom shape: (n_omega, n_states_n)
        denom = omega_grid[:, None] - (E_n[None, :] - E_f[f]) + 1j * gamma
        # weights shape: (n_omega, n_states_n)
        weights = 1.0 / denom
        # F[omega, m1, m2] = sum_n weights[omega, n] * D_fn[f, n, m1] * D_ni[n, m2]
        F = np.einsum('on,nm,np->omp', weights, D_fn[f], D_ni)
        # sigma_f[f, omega] = sum_{m1,m2} |F|^2
        sigma_f[f, :] = np.sum(np.abs(F) ** 2, axis=(1, 2))

    return sigma_f
