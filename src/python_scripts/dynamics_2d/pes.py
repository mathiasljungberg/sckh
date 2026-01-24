"""2D Potential Energy Surface with spline interpolation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Union

import numpy as np
from scipy.interpolate import RectBivariateSpline


@dataclass
class PES2D:
    """2D Potential Energy Surface with bicubic spline interpolation.

    Uses scipy.interpolate.RectBivariateSpline for smooth interpolation
    of the 2D potential and its derivatives.

    Attributes:
        x1: Grid points for first coordinate (SI: meters)
        x2: Grid points for second coordinate (SI: meters)
        E: Energies on 2D grid (SI: Joules), shape (len(x1), len(x2))
    """

    x1: np.ndarray
    x2: np.ndarray
    E: np.ndarray
    _spline: RectBivariateSpline = field(init=False, repr=False)

    def __post_init__(self):
        """Initialize spline interpolation."""
        # RectBivariateSpline requires sorted, monotonic increasing grids
        # kx=3, ky=3 for cubic spline in both directions
        self._spline = RectBivariateSpline(
            self.x1, self.x2, self.E, kx=3, ky=3, s=0
        )

    def energy(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Evaluate energy at position(s).

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Energy in Joules
        """
        result = self._spline(x1, x2, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    def force_x1(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Force in x1 direction: F1 = -dV/dx1.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Force in Newtons (kg*m/s^2)
        """
        result = -self._spline(x1, x2, dx=1, dy=0, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    def force_x2(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Force in x2 direction: F2 = -dV/dx2.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Force in Newtons (kg*m/s^2)
        """
        result = -self._spline(x1, x2, dx=0, dy=1, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    def forces(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        """Both force components: (F1, F2) = (-dV/dx1, -dV/dx2).

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Tuple of (F1, F2) forces in Newtons
        """
        return self.force_x1(x1, x2), self.force_x2(x1, x2)

    def energy_and_forces(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Tuple[
        Union[float, np.ndarray],
        Union[float, np.ndarray],
        Union[float, np.ndarray],
    ]:
        """Evaluate energy and both forces efficiently.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Tuple of (energy, F1, F2)
        """
        return self.energy(x1, x2), self.force_x1(x1, x2), self.force_x2(x1, x2)

    def second_derivative_x1(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Second derivative d^2V/dx1^2.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Second derivative in J/m^2
        """
        result = self._spline(x1, x2, dx=2, dy=0, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    def second_derivative_x2(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Second derivative d^2V/dx2^2.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Second derivative in J/m^2
        """
        result = self._spline(x1, x2, dx=0, dy=2, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    def mixed_derivative(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Mixed derivative d^2V/dx1dx2.

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Mixed derivative in J/m^2
        """
        result = self._spline(x1, x2, dx=1, dy=1, grid=False)
        if np.isscalar(x1) and np.isscalar(x2):
            return float(result)
        return result

    @property
    def x1_min(self) -> float:
        """Minimum x1 value of the grid."""
        return float(self.x1[0])

    @property
    def x1_max(self) -> float:
        """Maximum x1 value of the grid."""
        return float(self.x1[-1])

    @property
    def x2_min(self) -> float:
        """Minimum x2 value of the grid."""
        return float(self.x2[0])

    @property
    def x2_max(self) -> float:
        """Maximum x2 value of the grid."""
        return float(self.x2[-1])

    @property
    def dx1(self) -> float:
        """Grid spacing in x1 (assumes uniform grid)."""
        return float(self.x1[1] - self.x1[0])

    @property
    def dx2(self) -> float:
        """Grid spacing in x2 (assumes uniform grid)."""
        return float(self.x2[1] - self.x2[0])

    @property
    def npoints_x1(self) -> int:
        """Number of grid points in x1."""
        return len(self.x1)

    @property
    def npoints_x2(self) -> int:
        """Number of grid points in x2."""
        return len(self.x2)

    def find_minimum(self) -> Tuple[float, float, float]:
        """Find the minimum of the PES.

        Returns:
            Tuple of (x1_min, x2_min, E_min) at the minimum
        """
        from scipy.optimize import minimize

        # Start from grid minimum
        idx = np.unravel_index(np.argmin(self.E), self.E.shape)
        x0 = [self.x1[idx[0]], self.x2[idx[1]]]

        # Refine using optimization
        result = minimize(
            lambda x: self.energy(x[0], x[1]),
            x0,
            bounds=[
                (self.x1_min, self.x1_max),
                (self.x2_min, self.x2_max),
            ],
            method="L-BFGS-B",
        )
        return result.x[0], result.x[1], result.fun


def create_pes_from_file_2d(
    filepath: Path,
    units: str = "angstrom",
) -> PES2D:
    """Factory function to create PES2D from file.

    Args:
        filepath: Path to 2D PES file
        units: "angstrom" or "bohr" for coordinates

    Returns:
        PES2D object with spline interpolation
    """
    from .io import read_pes_file_2d

    x1, x2, E = read_pes_file_2d(filepath, units)
    return PES2D(x1=x1, x2=x2, E=E)


def create_harmonic_pes_2d(
    x1: np.ndarray,
    x2: np.ndarray,
    x1_0: float,
    x2_0: float,
    k1: float,
    k2: float,
    E0: float = 0.0,
) -> PES2D:
    """Create a separable 2D harmonic PES for testing.

    V(x1, x2) = E0 + 0.5 * k1 * (x1 - x1_0)^2 + 0.5 * k2 * (x2 - x2_0)^2

    Args:
        x1: Grid points for first coordinate (SI: meters)
        x2: Grid points for second coordinate (SI: meters)
        x1_0: Equilibrium position for x1 (SI: meters)
        x2_0: Equilibrium position for x2 (SI: meters)
        k1: Force constant for x1 (SI: N/m = kg/s^2)
        k2: Force constant for x2 (SI: N/m = kg/s^2)
        E0: Energy offset (SI: Joules)

    Returns:
        PES2D object
    """
    # Create meshgrid for 2D potential
    X1, X2 = np.meshgrid(x1, x2, indexing="ij")
    E = E0 + 0.5 * k1 * (X1 - x1_0) ** 2 + 0.5 * k2 * (X2 - x2_0) ** 2
    return PES2D(x1=x1, x2=x2, E=E)
