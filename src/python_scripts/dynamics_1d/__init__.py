"""
dynamics_1d: Classical dynamics on 1D potential energy surfaces.

This package provides tools for running classical molecular dynamics
simulations on 1D potential energy surfaces, with initial conditions
sampled from quantum mechanical ground state distributions. It also
includes spectrum calculation via FFT of trajectory data.
"""

from .constants import CONST
from .config import DynamicsConfig, load_config
from .pes import PES1D
from .integrators import velocity_verlet_step, run_trajectory
from .vibrational import solve_vibrational
from .sampling import create_initial_conditions
from .trajectory import DynamicsRunner, TrajectoryResult, EnsembleResult
from .dipole import Dipole1D, create_dipole_from_file, create_constant_dipole
from .spectrum_config import SpectrumConfig, FullConfig, load_full_config
from .spectrum import (
    SpectrumResult,
    SpectrumCalculator,
    compute_energy_phase,
    compute_F_if,
    get_frequency_grid,
)

__all__ = [
    # Constants
    "CONST",
    # Config
    "DynamicsConfig",
    "load_config",
    "SpectrumConfig",
    "FullConfig",
    "load_full_config",
    # PES
    "PES1D",
    # Dipole
    "Dipole1D",
    "create_dipole_from_file",
    "create_constant_dipole",
    # Integrators
    "velocity_verlet_step",
    "run_trajectory",
    # Vibrational
    "solve_vibrational",
    # Sampling
    "create_initial_conditions",
    # Trajectory
    "DynamicsRunner",
    "TrajectoryResult",
    "EnsembleResult",
    # Spectrum
    "SpectrumResult",
    "SpectrumCalculator",
    "compute_energy_phase",
    "compute_F_if",
    "get_frequency_grid",
]
