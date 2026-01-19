"""
dynamics_1d: Classical dynamics on 1D potential energy surfaces.

This package provides tools for running classical molecular dynamics
simulations on 1D potential energy surfaces, with initial conditions
sampled from quantum mechanical ground state distributions.
"""

from .constants import CONST
from .config import DynamicsConfig, load_config
from .pes import PES1D
from .integrators import velocity_verlet_step, run_trajectory
from .vibrational import solve_vibrational
from .sampling import create_initial_conditions
from .trajectory import DynamicsRunner, TrajectoryResult, EnsembleResult

__all__ = [
    "CONST",
    "DynamicsConfig",
    "load_config",
    "PES1D",
    "velocity_verlet_step",
    "run_trajectory",
    "solve_vibrational",
    "create_initial_conditions",
    "DynamicsRunner",
    "TrajectoryResult",
    "EnsembleResult",
]
