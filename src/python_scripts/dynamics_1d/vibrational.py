"""Fourier grid method for solving 1D vibrational Schrödinger equation.

Based on the Fourier grid method from src/m_fourier_grid.F90.
Solves: -ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ
"""

from typing import Tuple

import numpy as np

from .constants import CONST


def build_kinetic_matrix_fast(
    nstates: int,
    dx: float,
    mass: float,
) -> np.ndarray:
    """Build kinetic energy matrix using Fourier grid method (fast version).

    The kinetic energy matrix elements H_kin(i,j) depend only on |i-j|,
    exploiting translational symmetry. This allows O(N²) construction
    instead of O(N³).

    Uses sum over Fourier components k = 1 to (N-1)/2:
        prefac(k) = -2 * (k * 2π/(dx*N))² / N
        T_ij = -ℏ²/(2m) * Σ_k prefac(k) * cos(2π(i-j)k/N)

    Based on src/m_fourier_grid.F90:391-429

    Args:
        nstates: Number of grid points
        dx: Grid spacing (m)
        mass: Particle mass (kg)

    Returns:
        H_kin: Kinetic energy matrix (nstates x nstates) in Joules
    """
    recip_vec = 2.0 * CONST.pi / (dx * nstates)
    ncoeff = (nstates - 1) // 2  # Number of Fourier coefficients (excluding zero)

    # Precompute phase_and_prefac for each possible |i-j| value
    # Index 0 corresponds to |i-j| = 0, index 1 to |i-j| = 1, etc.
    phase_and_prefac = np.zeros(nstates)

    for i in range(nstates):
        for k in range(1, ncoeff + 1):
            prefac = -2.0 * (k * recip_vec) ** 2 / nstates
            phase_and_prefac[i] += prefac * np.cos(2.0 * CONST.pi * i * k / nstates)

    # Build symmetric matrix using precomputed values
    H_kin = np.zeros((nstates, nstates))
    for i in range(nstates):
        H_kin[i, i] = phase_and_prefac[0]
        for j in range(i + 1, nstates):
            H_kin[i, j] = phase_and_prefac[j - i]
            H_kin[j, i] = H_kin[i, j]

    # Apply kinetic energy prefactor: -ℏ²/(2m)
    prefac_kin = -CONST.hbar**2 / (2.0 * mass)
    H_kin *= prefac_kin

    return H_kin


def build_kinetic_matrix_vectorized(
    nstates: int,
    dx: float,
    mass: float,
) -> np.ndarray:
    """Build kinetic energy matrix using vectorized operations.

    More efficient numpy implementation of the Fourier grid kinetic energy.

    Args:
        nstates: Number of grid points
        dx: Grid spacing (m)
        mass: Particle mass (kg)

    Returns:
        H_kin: Kinetic energy matrix (nstates x nstates) in Joules
    """
    recip_vec = 2.0 * CONST.pi / (dx * nstates)
    ncoeff = (nstates - 1) // 2

    # Create arrays for vectorized computation
    k_vals = np.arange(1, ncoeff + 1)
    i_vals = np.arange(nstates)

    # Compute phase_and_prefac for each |i-j| value using broadcasting
    # prefac[k] = -2 * (k * recip_vec)^2 / nstates
    prefac_k = -2.0 * (k_vals * recip_vec) ** 2 / nstates

    # phase_and_prefac[i] = sum_k prefac[k] * cos(2*pi*i*k/nstates)
    # Shape: (nstates, ncoeff) -> sum over k -> (nstates,)
    angles = 2.0 * CONST.pi * np.outer(i_vals, k_vals) / nstates
    phase_and_prefac = np.sum(prefac_k * np.cos(angles), axis=1)

    # Build Toeplitz matrix (symmetric, depends only on |i-j|)
    from scipy.linalg import toeplitz

    H_kin = toeplitz(phase_and_prefac)

    # Apply kinetic energy prefactor: -ℏ²/(2m)
    prefac_kin = -CONST.hbar**2 / (2.0 * mass)
    H_kin *= prefac_kin

    return H_kin


def solve_vibrational(
    x: np.ndarray,
    V: np.ndarray,
    mass: float,
    n_states: int = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve 1D vibrational Schrödinger equation using Fourier grid method.

    Solves: -ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ

    Args:
        x: Position grid (SI: meters), must be uniformly spaced
        V: Potential energy on grid (SI: Joules)
        mass: Particle mass (SI: kg)
        n_states: Number of eigenstates to return (default: all)

    Returns:
        eigenvalues: Energy eigenvalues (SI: Joules), sorted ascending
        eigenvectors: Wavefunctions on grid (columns), normalized
    """
    nstates = len(x)
    dx = x[1] - x[0]

    # Build Hamiltonian: H = T + V
    # Kinetic energy matrix (off-diagonal from Fourier expansion)
    H_kin = build_kinetic_matrix_vectorized(nstates, dx, mass)

    # Potential energy matrix (diagonal)
    H_pot = np.diag(V)

    # Total Hamiltonian
    H = H_kin + H_pot

    # Diagonalize (symmetric matrix)
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    # Normalize wavefunctions (they should already be normalized by eigh,
    # but we ensure proper normalization for the grid)
    for i in range(nstates):
        norm = np.sqrt(np.sum(eigenvectors[:, i] ** 2) * dx)
        eigenvectors[:, i] /= norm

    # Return requested number of states
    if n_states is not None:
        eigenvalues = eigenvalues[:n_states]
        eigenvectors = eigenvectors[:, :n_states]

    return eigenvalues, eigenvectors


def harmonic_eigenvalues(omega: float, n_max: int) -> np.ndarray:
    """Compute analytical eigenvalues for harmonic oscillator.

    E_n = ℏω(n + 1/2)

    Args:
        omega: Angular frequency (rad/s)
        n_max: Maximum quantum number (0 to n_max-1)

    Returns:
        Eigenvalues in Joules
    """
    n = np.arange(n_max)
    return CONST.hbar * omega * (n + 0.5)


def harmonic_omega_from_k(k: float, mass: float) -> float:
    """Compute angular frequency from force constant.

    ω = √(k/m)

    Args:
        k: Force constant (N/m)
        mass: Particle mass (kg)

    Returns:
        Angular frequency (rad/s)
    """
    return np.sqrt(k / mass)
