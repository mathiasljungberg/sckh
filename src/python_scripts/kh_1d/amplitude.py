"""Kramers-Heisenberg scattering amplitude computation.

Implements the KH scattering amplitude F for resonant and non-resonant
X-ray emission spectroscopy (XES).

Based on src/m_KH_utils.F90:1108-1185.
"""

import numpy as np


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

    This is the resonant contribution where the emitted photon energy
    matches the intermediate-to-final transition energy.

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
    n_states_n = len(E_n)
    F = np.zeros((3, 3), dtype=complex)

    for n in range(n_states_n):
        # Denominator: omega - (E_n - E_f) + i*gamma
        denom = omega - (E_n[n] - E_f) + 1j * gamma

        # F(m1, m2) += D_fn(m1) * D_ni(m2) / denom
        for m1 in range(3):
            for m2 in range(3):
                F[m1, m2] += D_fn[n, m1] * D_ni[n, m2] / denom

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

    This is the non-resonant contribution (no imaginary part in denominator).

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
    n_states_n = len(E_n)
    F = np.zeros((3, 3), dtype=complex)

    for n in range(n_states_n):
        # Denominator: omega + (E_n - E_i)
        denom = omega + (E_n[n] - E_i)

        # Avoid division by zero
        if abs(denom) < 1e-12:
            continue

        # F(m1, m2) -= D_fn(m1) * D_ni(m2) / denom
        for m1 in range(3):
            for m2 in range(3):
                F[m1, m2] -= D_fn[n, m1] * D_ni[n, m2] / denom

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

    For unpolarized light, we sum over all polarization components.
    Alternatively, can use trace of |F|^2 matrix.

    Args:
        F: Scattering amplitude (3, 3) complex array

    Returns:
        Cross-section (arbitrary units)
    """
    # Sum over all components: |F|^2
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

    where the sum is over final vibrational states.

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
    n_omega = len(omega_grid)
    n_states_f = len(E_f)
    sigma = np.zeros(n_omega)

    for i_omega, omega in enumerate(omega_grid):
        for f in range(n_states_f):
            # Get dipoles for this final vibrational state
            D_fn_f = D_fn[f, :, :]  # Shape: (n_states_n, 3)

            # Compute amplitude for this final state
            F = compute_amplitude_F_res(
                E_i, E_n, E_f[f], D_ni, D_fn_f, omega, gamma
            )

            # Add to cross-section
            sigma[i_omega] += compute_cross_section_from_F(F)

    return sigma


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

    More efficient implementation using numpy broadcasting.

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
    n_omega = len(omega_grid)
    n_states_n = len(E_n)
    n_states_f = len(E_f)

    sigma = np.zeros(n_omega)

    # For each final vibrational state
    for f in range(n_states_f):
        # For each frequency point
        for i_omega in range(n_omega):
            omega = omega_grid[i_omega]

            # Compute scattering amplitude for this (omega, f) pair
            F = np.zeros((3, 3), dtype=complex)

            for n in range(n_states_n):
                # Resonant denominator
                denom = omega - (E_n[n] - E_f[f]) + 1j * gamma

                # Outer product of dipole vectors
                for m1 in range(3):
                    for m2 in range(3):
                        F[m1, m2] += D_fn[f, n, m1] * D_ni[n, m2] / denom

            # Add |F|^2 to cross-section
            sigma[i_omega] += np.sum(np.abs(F) ** 2)

    return sigma


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
