"""Potential Energy Surface with spline interpolation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Union

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class PES1D:
    """1D Potential Energy Surface with cubic spline interpolation.

    Uses scipy.interpolate.CubicSpline with natural boundary conditions
    (second derivative = 0 at endpoints), matching Fortran's natural splines.

    Attributes:
        x: Grid points (SI: meters)
        E: Energies at grid points (SI: Joules)
    """

    x: np.ndarray
    E: np.ndarray
    _spline: CubicSpline = field(init=False, repr=False)

    def __post_init__(self):
        """Initialize spline interpolation."""
        # Natural splines: bc_type='natural' gives zero second derivative at boundaries
        self._spline = CubicSpline(self.x, self.E, bc_type="natural")

    def energy(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate energy at position(s) x.

        Args:
            x: Position(s) in meters

        Returns:
            Energy in Joules
        """
        return self._spline(x)

    def force(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate force at position(s) x.

        Force = -dV/dx, computed analytically from spline derivative.

        Args:
            x: Position(s) in meters

        Returns:
            Force in Newtons (kg*m/s^2)
        """
        return -self._spline(x, 1)  # First derivative

    def energy_and_force(
        self, x: Union[float, np.ndarray]
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        """Evaluate both energy and force efficiently.

        Args:
            x: Position(s) in meters

        Returns:
            Tuple of (energy, force)
        """
        return self.energy(x), self.force(x)

    def second_derivative(
        self, x: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Evaluate second derivative of potential at position(s) x.

        d²V/dx², useful for computing harmonic frequencies.

        Args:
            x: Position(s) in meters

        Returns:
            Second derivative in J/m^2
        """
        return self._spline(x, 2)

    @property
    def x_min(self) -> float:
        """Minimum x value of the grid."""
        return float(self.x[0])

    @property
    def x_max(self) -> float:
        """Maximum x value of the grid."""
        return float(self.x[-1])

    @property
    def dx(self) -> float:
        """Grid spacing (assumes uniform grid)."""
        return float(self.x[1] - self.x[0])

    @property
    def npoints(self) -> int:
        """Number of grid points."""
        return len(self.x)

    def find_minimum(self) -> Tuple[float, float]:
        """Find the minimum of the PES.

        Returns:
            Tuple of (x_min, E_min) at the minimum
        """
        # Start from grid minimum
        idx_min = np.argmin(self.E)
        x_guess = self.x[idx_min]

        # Refine using Newton-Raphson on the derivative
        from scipy.optimize import minimize_scalar

        result = minimize_scalar(
            self.energy, bounds=(self.x_min, self.x_max), method="bounded"
        )
        return result.x, result.fun


def create_pes_from_file(
    filepath: Path,
    units: str = "angstrom",
    energy_column: int = 1,
) -> PES1D:
    """Factory function to create PES from file.

    Args:
        filepath: Path to PES file
        units: "angstrom" or "bohr" for x coordinate
        energy_column: Which column to use for energy (1-indexed)

    Returns:
        PES1D object with spline interpolation
    """
    from .io import read_pes_file

    x, E = read_pes_file(filepath, units, energy_column)
    return PES1D(x=x, E=E)


def create_harmonic_pes(
    x: np.ndarray,
    x0: float,
    k: float,
    E0: float = 0.0,
) -> PES1D:
    """Create a harmonic PES for testing.

    V(x) = E0 + 0.5 * k * (x - x0)^2

    Args:
        x: Grid points (SI: meters)
        x0: Equilibrium position (SI: meters)
        k: Force constant (SI: N/m = kg/s^2)
        E0: Energy offset (SI: Joules)

    Returns:
        PES1D object
    """
    E = E0 + 0.5 * k * (x - x0) ** 2
    return PES1D(x=x, E=E)
