"""2D Dipole moment surface with spline interpolation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

import numpy as np
from scipy.interpolate import RectBivariateSpline


@dataclass
class Dipole2D:
    """2D dipole moment surface with bicubic spline interpolation.

    Uses scipy.interpolate.RectBivariateSpline for smooth interpolation
    of the 2D dipole components and their derivatives.

    Attributes:
        x1: Grid points for first coordinate (SI: meters)
        x2: Grid points for second coordinate (SI: meters)
        d: Dipole components at grid points (n_x1, n_x2, 3) for x, y, z components
           in atomic units (e * bohr)
    """

    x1: np.ndarray
    x2: np.ndarray
    d: np.ndarray  # Shape: (n_x1, n_x2, 3)
    _splines: List[RectBivariateSpline] = field(init=False, repr=False)

    def __post_init__(self):
        """Initialize spline interpolation for each component."""
        if self.d.ndim != 3 or self.d.shape[2] != 3:
            raise ValueError(
                f"Dipole array must have shape (n_x1, n_x2, 3), got {self.d.shape}"
            )
        if len(self.x1) != self.d.shape[0]:
            raise ValueError(
                f"x1 and d must have matching first dimension: "
                f"{len(self.x1)} vs {self.d.shape[0]}"
            )
        if len(self.x2) != self.d.shape[1]:
            raise ValueError(
                f"x2 and d must have matching second dimension: "
                f"{len(self.x2)} vs {self.d.shape[1]}"
            )

        # Create a spline for each component (x, y, z)
        # kx=3, ky=3 for cubic spline in both directions
        self._splines = [
            RectBivariateSpline(self.x1, self.x2, self.d[:, :, i], kx=3, ky=3, s=0)
            for i in range(3)
        ]

    def dipole(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> np.ndarray:
        """Interpolate dipole moment at position(s) (x1, x2).

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Dipole moment(s) in atomic units. Shape (3,) for scalar inputs,
            or (n, 3) for array inputs.
        """
        scalar_input = np.isscalar(x1) and np.isscalar(x2)
        x1 = np.atleast_1d(x1)
        x2 = np.atleast_1d(x2)

        result = np.column_stack([
            spline(x1, x2, grid=False) for spline in self._splines
        ])

        if scalar_input:
            return result[0]
        return result

    def dipole_derivative_x1(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> np.ndarray:
        """Evaluate dipole derivative with respect to x1 at position(s).

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Dipole derivative(s). Shape (3,) for scalar inputs,
            or (n, 3) for array inputs.
        """
        scalar_input = np.isscalar(x1) and np.isscalar(x2)
        x1 = np.atleast_1d(x1)
        x2 = np.atleast_1d(x2)

        result = np.column_stack([
            spline(x1, x2, dx=1, dy=0, grid=False) for spline in self._splines
        ])

        if scalar_input:
            return result[0]
        return result

    def dipole_derivative_x2(
        self, x1: Union[float, np.ndarray], x2: Union[float, np.ndarray]
    ) -> np.ndarray:
        """Evaluate dipole derivative with respect to x2 at position(s).

        Args:
            x1: Position(s) along first coordinate (meters)
            x2: Position(s) along second coordinate (meters)

        Returns:
            Dipole derivative(s). Shape (3,) for scalar inputs,
            or (n, 3) for array inputs.
        """
        scalar_input = np.isscalar(x1) and np.isscalar(x2)
        x1 = np.atleast_1d(x1)
        x2 = np.atleast_1d(x2)

        result = np.column_stack([
            spline(x1, x2, dx=0, dy=1, grid=False) for spline in self._splines
        ])

        if scalar_input:
            return result[0]
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
    def npoints_x1(self) -> int:
        """Number of grid points in x1."""
        return len(self.x1)

    @property
    def npoints_x2(self) -> int:
        """Number of grid points in x2."""
        return len(self.x2)


def create_dipole_from_file_2d(
    filepath: Path,
    position_units: str = "angstrom",
    index_order: str = "C",
    dipole_components: int = 3,
) -> Dipole2D:
    """Factory function to create Dipole2D from file.

    Args:
        filepath: Path to 2D dipole file
            Format for 3 components: x1 x2 d_x d_y d_z
            Format for 1 component: x1 x2 |d|^2 (dipole squared)
        position_units: "angstrom" or "bohr" for coordinates
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
        dipole_components: Number of dipole components in file (1 or 3).
            If 1, the file contains |d|^2. The square root is taken and
            placed in the z-component (index 2), with x and y set to 0.

    Returns:
        Dipole2D object with spline interpolation
    """
    from .io import read_dipole_file_2d

    x1, x2, d = read_dipole_file_2d(
        filepath, position_units, index_order, dipole_components
    )
    return Dipole2D(x1=x1, x2=x2, d=d)


def create_constant_dipole_2d(
    x1: np.ndarray,
    x2: np.ndarray,
    d_value: np.ndarray,
) -> Dipole2D:
    """Create a constant dipole moment surface (Franck-Condon approximation).

    Args:
        x1: Grid points for first coordinate (SI: meters)
        x2: Grid points for second coordinate (SI: meters)
        d_value: Constant dipole value (3,) in atomic units

    Returns:
        Dipole2D object with constant dipole
    """
    d_value = np.atleast_1d(d_value)
    if len(d_value) != 3:
        raise ValueError(f"d_value must have 3 components, got {len(d_value)}")

    # Create 2D array with constant dipole at all grid points
    d = np.zeros((len(x1), len(x2), 3))
    d[:, :, :] = d_value  # Broadcast to all grid points

    return Dipole2D(x1=x1, x2=x2, d=d)
