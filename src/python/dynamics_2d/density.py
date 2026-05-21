"""Trajectory density on the dynamics 2D grid.

Provides elementary functions to place broadened trajectory points onto
the same 2D grid used for PES interpolation, plus composite functions that
operate on an EnsembleResult2D for a chosen set of time steps.

All grid arrays here are in SI units (meters), matching EnsembleResult2D.
The composite functions accept FWHM in user position units
(Angstrom or Bohr) and convert internally.
"""

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from dynamics_1d.constants import CONST

from .trajectory import EnsembleResult2D, TrajectoryResult2D


_FWHM_TO_SIGMA_FACTOR = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


@dataclass
class TrajectoryDensityResult2D:
    """Gridded trajectory density at a set of time steps."""

    density: np.ndarray  # (nsteps, nx1, nx2); sum * dx1 * dx2 == 1 per slice
    time: np.ndarray  # (nsteps,), SI seconds
    step_indices: np.ndarray  # (nsteps,), int
    x1_grid: np.ndarray  # SI meters
    x2_grid: np.ndarray  # SI meters
    fwhm_SI: float  # FWHM used, in meters


def fwhm_to_sigma(fwhm: float) -> float:
    """Convert FWHM of a Gaussian to its standard deviation sigma."""
    return fwhm * _FWHM_TO_SIGMA_FACTOR


def _position_unit_to_meters(position_units: str) -> float:
    """Return the conversion factor from `position_units` to meters."""
    u = position_units.lower()
    if u == "angstrom":
        return 1e-10
    if u == "bohr":
        return CONST.bohr
    raise ValueError(f"Unknown position units: {position_units}")


def gaussian_2d_on_grid(
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
    x1_0: float,
    x2_0: float,
    sigma: float,
) -> np.ndarray:
    """Normalized isotropic 2D Gaussian evaluated on a rectangular grid.

    Args:
        x1_grid: 1D x1 grid (SI: meters)
        x2_grid: 1D x2 grid (SI: meters)
        x1_0, x2_0: Gaussian center (SI: meters)
        sigma: Standard deviation (SI: meters)

    Returns:
        Array of shape (nx1, nx2) with the continuous normalization
        1 / (2 pi sigma^2). For a grid that contains the Gaussian bulk,
        out.sum() * dx1 * dx2 == 1 to good approximation.
    """
    inv_two_sigma2 = 0.5 / (sigma * sigma)
    gx = np.exp(-((x1_grid - x1_0) ** 2) * inv_two_sigma2)
    gy = np.exp(-((x2_grid - x2_0) ** 2) * inv_two_sigma2)
    prefactor = 1.0 / (2.0 * np.pi * sigma * sigma)
    return prefactor * np.outer(gx, gy)


def positions_to_grid_density(
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
    x1_points: np.ndarray,
    x2_points: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Sum isotropic 2D Gaussians at given points and normalize to 1.

    Args:
        x1_grid: 1D x1 grid (SI: meters)
        x2_grid: 1D x2 grid (SI: meters)
        x1_points: 1D array of x1 positions (SI: meters)
        x2_points: 1D array of x2 positions (SI: meters); same length as x1_points
        sigma: Gaussian standard deviation (SI: meters)

    Returns:
        Array of shape (nx1, nx2). Normalized so that
        out.sum() * dx1 * dx2 == 1.
    """
    x1_points = np.asarray(x1_points)
    x2_points = np.asarray(x2_points)
    if x1_points.shape != x2_points.shape:
        raise ValueError(
            f"x1_points and x2_points must have the same shape, "
            f"got {x1_points.shape} vs {x2_points.shape}"
        )

    inv_two_sigma2 = 0.5 / (sigma * sigma)
    # gx[i, k] = exp(-0.5 ((x1_grid[k] - x1_points[i])/sigma)^2)
    gx = np.exp(-((x1_grid[None, :] - x1_points[:, None]) ** 2) * inv_two_sigma2)
    gy = np.exp(-((x2_grid[None, :] - x2_points[:, None]) ** 2) * inv_two_sigma2)
    density = gx.T @ gy  # (nx1, nx2)

    dx1 = x1_grid[1] - x1_grid[0]
    dx2 = x2_grid[1] - x2_grid[0]
    total = density.sum() * dx1 * dx2
    if total <= 0.0:
        raise ValueError(
            "Total density is non-positive; trajectory points may lie far "
            "outside the grid relative to sigma."
        )
    return density / total


def trajectories_density_at_step(
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
    trajectories: Sequence[TrajectoryResult2D],
    step: int,
    fwhm: float,
    position_units: str = "angstrom",
) -> np.ndarray:
    """Gridded broadened density of a list of trajectories at one time step.

    Args:
        x1_grid: 1D x1 grid (SI: meters).
        x2_grid: 1D x2 grid (SI: meters).
        trajectories: List of TrajectoryResult2D (positions in SI meters).
        step: Time-step index into each trajectory.
        fwhm: Gaussian FWHM in `position_units`.
        position_units: "angstrom" or "bohr".

    Returns:
        Density on (x1_grid, x2_grid), shape (nx1, nx2), normalized so
        sum * dx1 * dx2 == 1.
    """
    fwhm_SI = fwhm * _position_unit_to_meters(position_units)
    sigma = fwhm_to_sigma(fwhm_SI)
    x1_points = np.array([t.x1[step] for t in trajectories])
    x2_points = np.array([t.x2[step] for t in trajectories])
    return positions_to_grid_density(x1_grid, x2_grid, x1_points, x2_points, sigma)


def trajectories_density_timeseries(
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
    trajectories: Sequence[TrajectoryResult2D],
    fwhm: float,
    step_indices: Optional[Sequence[int]] = None,
    position_units: str = "angstrom",
) -> TrajectoryDensityResult2D:
    """Gridded broadened density at a chosen set of time steps.

    Same semantics as :func:`ensemble_density_timeseries` but takes the
    trajectory list and grid arrays directly; suitable for trajectories
    loaded from disk via :func:`dynamics_2d.io.read_trajectories_2d`.

    Args:
        x1_grid: 1D x1 grid (SI: meters).
        x2_grid: 1D x2 grid (SI: meters).
        trajectories: Non-empty list of TrajectoryResult2D.
        fwhm: Gaussian FWHM in `position_units`.
        step_indices: Indices into the trajectory time axis. If None,
            uses every step (0 .. n_steps - 1).
        position_units: "angstrom" or "bohr".

    Returns:
        TrajectoryDensityResult2D with one density slice per selected step.
    """
    if not trajectories:
        raise ValueError("trajectories list is empty")

    time_full = trajectories[0].time
    if step_indices is None:
        step_indices_arr = np.arange(len(time_full), dtype=int)
    else:
        step_indices_arr = np.asarray(step_indices, dtype=int)

    fwhm_SI = fwhm * _position_unit_to_meters(position_units)
    sigma = fwhm_to_sigma(fwhm_SI)

    nx1 = len(x1_grid)
    nx2 = len(x2_grid)
    density = np.empty((len(step_indices_arr), nx1, nx2), dtype=float)
    for out_i, step in enumerate(step_indices_arr):
        step_i = int(step)
        x1_points = np.array([t.x1[step_i] for t in trajectories])
        x2_points = np.array([t.x2[step_i] for t in trajectories])
        density[out_i] = positions_to_grid_density(
            x1_grid, x2_grid, x1_points, x2_points, sigma
        )

    time = time_full[step_indices_arr]

    return TrajectoryDensityResult2D(
        density=density,
        time=time,
        step_indices=step_indices_arr,
        x1_grid=x1_grid,
        x2_grid=x2_grid,
        fwhm_SI=fwhm_SI,
    )


def ensemble_density_at_step(
    ensemble: EnsembleResult2D,
    step: int,
    fwhm: float,
    position_units: str = "angstrom",
) -> np.ndarray:
    """Gridded broadened density of an ensemble at a single time step."""
    return trajectories_density_at_step(
        ensemble.x1_grid,
        ensemble.x2_grid,
        ensemble.trajectories,
        step,
        fwhm,
        position_units=position_units,
    )


def ensemble_density_timeseries(
    ensemble: EnsembleResult2D,
    fwhm: float,
    step_indices: Optional[Sequence[int]] = None,
    position_units: str = "angstrom",
) -> TrajectoryDensityResult2D:
    """Gridded broadened density at a chosen set of time steps."""
    return trajectories_density_timeseries(
        ensemble.x1_grid,
        ensemble.x2_grid,
        ensemble.trajectories,
        fwhm,
        step_indices=step_indices,
        position_units=position_units,
    )
