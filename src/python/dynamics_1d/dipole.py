"""Dipole moment surface with spline interpolation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class Dipole1D:
    """1D dipole moment surface with cubic spline interpolation.

    Uses scipy.interpolate.CubicSpline with natural boundary conditions
    (second derivative = 0 at endpoints), matching Fortran's natural splines.

    Attributes:
        x: Grid points (SI: meters)
        d: Dipole components at grid points (n_points, 3) for x, y, z components
           in atomic units (e * bohr)
    """

    x: np.ndarray
    d: np.ndarray  # Shape: (n_points, 3)
    _splines: List[CubicSpline] = field(init=False, repr=False)

    def __post_init__(self):
        """Initialize spline interpolation for each component."""
        if self.d.ndim != 2 or self.d.shape[1] != 3:
            raise ValueError(
                f"Dipole array must have shape (n_points, 3), got {self.d.shape}"
            )
        if len(self.x) != self.d.shape[0]:
            raise ValueError(
                f"x and d must have same number of points: {len(self.x)} vs {self.d.shape[0]}"
            )

        # Create a spline for each component (x, y, z)
        self._splines = [
            CubicSpline(self.x, self.d[:, i], bc_type="natural") for i in range(3)
        ]

    def dipole(self, x: Union[float, np.ndarray]) -> np.ndarray:
        """Interpolate dipole moment at position(s) x.

        Args:
            x: Position(s) in meters

        Returns:
            Dipole moment(s) in atomic units. Shape (3,) for scalar x,
            or (len(x), 3) for array x.
        """
        x = np.atleast_1d(x)
        result = np.column_stack([spline(x) for spline in self._splines])

        if result.shape[0] == 1:
            return result[0]
        return result

    def dipole_derivative(self, x: Union[float, np.ndarray]) -> np.ndarray:
        """Evaluate dipole derivative at position(s) x.

        Args:
            x: Position(s) in meters

        Returns:
            Dipole derivative(s). Shape (3,) for scalar x, or (len(x), 3) for array x.
        """
        x = np.atleast_1d(x)
        result = np.column_stack([spline(x, 1) for spline in self._splines])

        if result.shape[0] == 1:
            return result[0]
        return result

    @property
    def x_min(self) -> float:
        """Minimum x value of the grid."""
        return float(self.x[0])

    @property
    def x_max(self) -> float:
        """Maximum x value of the grid."""
        return float(self.x[-1])

    @property
    def npoints(self) -> int:
        """Number of grid points."""
        return len(self.x)


def create_dipole_from_file(filepath: Path, units: str = "angstrom") -> Dipole1D:
    """Factory function to create Dipole1D from file.

    Args:
        filepath: Path to dipole file (format: x d_x d_y d_z)
        units: "angstrom" or "bohr" for x coordinate

    Returns:
        Dipole1D object with spline interpolation
    """
    from .io import read_dipole_file

    x, d = read_dipole_file(filepath, units)
    return Dipole1D(x=x, d=d)


def create_constant_dipole(x: np.ndarray, d_value: np.ndarray) -> Dipole1D:
    """Create a constant dipole moment surface (Franck-Condon approximation).

    Args:
        x: Grid points (SI: meters)
        d_value: Constant dipole value (3,) in atomic units

    Returns:
        Dipole1D object with constant dipole
    """
    d_value = np.atleast_1d(d_value)
    if len(d_value) != 3:
        raise ValueError(f"d_value must have 3 components, got {len(d_value)}")

    d = np.tile(d_value, (len(x), 1))
    return Dipole1D(x=x, d=d)
