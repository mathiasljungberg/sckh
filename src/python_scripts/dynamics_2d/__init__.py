"""
dynamics_2d: Classical dynamics on 2D potential energy surfaces.

This package provides tools for running classical molecular dynamics
simulations on 2D potential energy surfaces with two degrees of freedom
(potentially different masses). Initial conditions are sampled from
quantum mechanical ground state distributions using a product wavefunction
approximation.
"""

from .config import (
    GridConfig2D,
    TimeConfig,
    SamplingConfig2D,
    DynamicsConfig2D,
    FullConfig2D,
    load_config,
    save_config,
)
from .pes import PES2D, create_pes_from_file_2d, create_harmonic_pes_2d
from .vibrational import solve_product_ground_state
from .sampling import create_initial_conditions_2d
from .integrators import velocity_verlet_step_2d, run_trajectory_2d
from .trajectory import (
    TrajectoryResult2D,
    EnsembleResult2D,
    DynamicsRunner2D,
)
from .io import read_pes_file_2d

# Re-export constants from dynamics_1d
from python_scripts.dynamics_1d.constants import CONST

__all__ = [
    # Constants
    "CONST",
    # Config
    "GridConfig2D",
    "TimeConfig",
    "SamplingConfig2D",
    "DynamicsConfig2D",
    "FullConfig2D",
    "load_config",
    "save_config",
    # PES
    "PES2D",
    "create_pes_from_file_2d",
    "create_harmonic_pes_2d",
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
    # I/O
    "read_pes_file_2d",
]
