"""Product approximation solver for 2D vibrational ground state.

Uses the product wavefunction approximation:
    psi(x1, x2) ~ psi1(x1) * psi2(x2)

This reduces the 2D problem to two independent 1D problems, which can be
solved using the Fourier grid method from dynamics_1d.
"""

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from python_scripts.dynamics_1d.vibrational import solve_vibrational


@dataclass
class ProductGroundState:
    """Result from product approximation ground state calculation.

    Attributes:
        E1: Ground state energy for x1 coordinate (J)
        E2: Ground state energy for x2 coordinate (J)
        E_total: Total ground state energy E1 + E2 (J)
        psi1: Ground state wavefunction for x1
        psi2: Ground state wavefunction for x2
        x1_grid: Grid points for x1 (m)
        x2_grid: Grid points for x2 (m)
    """

    E1: float
    E2: float
    E_total: float
    psi1: np.ndarray
    psi2: np.ndarray
    x1_grid: np.ndarray
    x2_grid: np.ndarray

    def psi_squared_1(self) -> np.ndarray:
        """Return |psi1|^2 probability density."""
        return np.abs(self.psi1) ** 2

    def psi_squared_2(self) -> np.ndarray:
        """Return |psi2|^2 probability density."""
        return np.abs(self.psi2) ** 2

    def psi_2d(self) -> np.ndarray:
        """Return 2D product wavefunction psi(x1, x2).

        Returns:
            2D array of shape (len(x1_grid), len(x2_grid))
        """
        return np.outer(self.psi1, self.psi2)

    def psi_squared_2d(self) -> np.ndarray:
        """Return 2D probability density |psi(x1, x2)|^2.

        Returns:
            2D array of shape (len(x1_grid), len(x2_grid))
        """
        return np.outer(self.psi_squared_1(), self.psi_squared_2())


def solve_product_ground_state(
    x1_grid: np.ndarray,
    V1: np.ndarray,
    mass1: float,
    x2_grid: np.ndarray,
    V2: np.ndarray,
    mass2: float,
) -> ProductGroundState:
    """Solve 2D vibrational problem using product approximation.

    psi(x1, x2) ~ psi1(x1) * psi2(x2)

    Solves two independent 1D Schrodinger equations and returns the
    product state.

    Args:
        x1_grid: Position grid for x1 (SI: meters)
        V1: Potential energy for x1 on grid (SI: Joules)
        mass1: Particle mass for x1 (SI: kg)
        x2_grid: Position grid for x2 (SI: meters)
        V2: Potential energy for x2 on grid (SI: Joules)
        mass2: Particle mass for x2 (SI: kg)

    Returns:
        ProductGroundState with wavefunctions and energies
    """
    # Solve 1D problem for x1 coordinate
    E1_vals, psi1_vecs = solve_vibrational(x1_grid, V1, mass1, n_states=1)
    E1 = E1_vals[0]
    psi1 = psi1_vecs[:, 0]

    # Solve 1D problem for x2 coordinate
    E2_vals, psi2_vecs = solve_vibrational(x2_grid, V2, mass2, n_states=1)
    E2 = E2_vals[0]
    psi2 = psi2_vecs[:, 0]

    return ProductGroundState(
        E1=E1,
        E2=E2,
        E_total=E1 + E2,
        psi1=psi1,
        psi2=psi2,
        x1_grid=x1_grid,
        x2_grid=x2_grid,
    )


def solve_product_states(
    x1_grid: np.ndarray,
    V1: np.ndarray,
    mass1: float,
    x2_grid: np.ndarray,
    V2: np.ndarray,
    mass2: float,
    n1_states: int = 1,
    n2_states: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve 2D vibrational problem for multiple product states.

    Args:
        x1_grid: Position grid for x1 (SI: meters)
        V1: Potential energy for x1 on grid (SI: Joules)
        mass1: Particle mass for x1 (SI: kg)
        x2_grid: Position grid for x2 (SI: meters)
        V2: Potential energy for x2 on grid (SI: Joules)
        mass2: Particle mass for x2 (SI: kg)
        n1_states: Number of states to compute for x1
        n2_states: Number of states to compute for x2

    Returns:
        E1: Eigenvalues for x1 (n1_states,)
        psi1: Eigenvectors for x1 (len(x1_grid), n1_states)
        E2: Eigenvalues for x2 (n2_states,)
        psi2: Eigenvectors for x2 (len(x2_grid), n2_states)
    """
    # Solve 1D problems
    E1, psi1 = solve_vibrational(x1_grid, V1, mass1, n_states=n1_states)
    E2, psi2 = solve_vibrational(x2_grid, V2, mass2, n_states=n2_states)

    return E1, psi1, E2, psi2
