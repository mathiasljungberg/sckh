"""
dynamics_2d: Classical dynamics on 2D potential energy surfaces.

This package provides tools for running classical molecular dynamics
simulations on 2D potential energy surfaces with two degrees of freedom
(potentially different masses). Initial conditions are sampled from
quantum mechanical ground state distributions using a product wavefunction
approximation.

Spectrum computation itself lives in the ``sckh`` package.
"""

from .config import (
    GridConfig2D,
    TimeConfig,
    SamplingConfig2D,
    DynamicsConfig2D,
    InterpolationConfig2D,
    SpectrumConfig2D,  # backwards compat alias
    load_config,
    save_config,
)
from .pes import PES2D, create_pes_from_file_2d, create_harmonic_pes_2d
from .dipole import Dipole2D, create_dipole_from_file_2d, create_constant_dipole_2d
from .vibrational import solve_product_ground_state
from .sampling import create_initial_conditions_2d
from .integrators import velocity_verlet_step_2d, run_trajectory_2d
from .trajectory import (
    TrajectoryResult2D,
    EnsembleResult2D,
    DynamicsRunner2D,
)
from .interpolation import SurfaceInterpolator2D, SpectrumCalculator2D
from .io import read_pes_file_2d, read_dipole_file_2d

# Re-export constants from dynamics_1d
from dynamics_1d.constants import CONST

# Re-export shared types from sckh package
from sckh import SCKHTrajectory, SpectrumResult, SCKHSpectrumCalculator, FullConfig

__all__ = [
    # Constants
    "CONST",
    # Config
    "GridConfig2D",
    "TimeConfig",
    "SamplingConfig2D",
    "DynamicsConfig2D",
    "InterpolationConfig2D",
    "SpectrumConfig2D",
    "FullConfig",
    "load_config",
    "save_config",
    # PES
    "PES2D",
    "create_pes_from_file_2d",
    "create_harmonic_pes_2d",
    # Dipole
    "Dipole2D",
    "create_dipole_from_file_2d",
    "create_constant_dipole_2d",
    # Vibrational
    "solve_product_ground_state",
    # Sampling
    "create_initial_conditions_2d",
    # Integrators
    "velocity_verlet_step_2d",
    "run_trajectory_2d",
    # Trajectory
    "TrajectoryResult2D",
    "EnsembleResult2D",
    "DynamicsRunner2D",
    # Surface interpolation
    "SurfaceInterpolator2D",
    "SpectrumCalculator2D",
    # Spectrum (re-exported from sckh)
    "SCKHSpectrumCalculator",
    "SCKHTrajectory",
    "SpectrumResult",
    # I/O
    "read_pes_file_2d",
    "read_dipole_file_2d",
]
