"""
dynamics_1d: Classical dynamics on 1D potential energy surfaces.

This package provides tools for running classical molecular dynamics
simulations on 1D potential energy surfaces, with initial conditions
sampled from quantum mechanical ground state distributions. It also
provides surface interpolation (SurfaceInterpolator) to produce
SCKHTrajectory objects from dynamics results.

Spectrum computation itself lives in the ``sckh`` package.
"""

from .constants import CONST
from .config import DynamicsConfig, load_config
from .pes import PES1D
from .integrators import velocity_verlet_step, run_trajectory
from .vibrational import solve_vibrational
from .sampling import create_initial_conditions
from .trajectory import DynamicsRunner, TrajectoryResult, EnsembleResult
from .dipole import Dipole1D, create_dipole_from_file, create_constant_dipole
from .spectrum_config import InterpolationConfig
from .sckh_trajectory import sckh_trajectory_from_dynamics_1d
from .interpolation import SurfaceInterpolator, SpectrumCalculator
from .workflow import SCKHInputs1D, compute_sckh_trajectories_1d

# Lazy re-exports from the sckh package (backwards compatibility).
# Defined via module __getattr__ so that ``import sckh`` does not trigger
# a circular import when sckh.spectrum loads dynamics_1d.constants.
_SCKH_REEXPORTS = {
    "SCKHTrajectory": ("sckh.trajectory", "SCKHTrajectory"),
    "SpectrumConfig": ("sckh.config", "SpectrumConfig"),
    "FullConfig": ("sckh.config", "FullConfig"),
    "load_full_config": ("sckh.config", "load_full_config"),
    "save_full_config": ("sckh.config", "save_full_config"),
    "SpectrumResult": ("sckh.spectrum", "SpectrumResult"),
    "SCKHSpectrumCalculator": ("sckh.spectrum", "SCKHSpectrumCalculator"),
    "compute_energy_phase": ("sckh.spectrum", "compute_energy_phase"),
    "compute_F_if": ("sckh.spectrum", "compute_F_if"),
    "get_frequency_grid": ("sckh.spectrum", "get_frequency_grid"),
    "compute_spectrum_from_sckh": ("sckh.spectrum", "compute_spectrum_from_sckh"),
    "read_sckh_trajectory": ("sckh.io", "read_sckh_trajectory"),
    "write_sckh_trajectory": ("sckh.io", "write_sckh_trajectory"),
    "read_sckh_trajectory_list": ("sckh.io", "read_sckh_trajectory_list"),
}


def __getattr__(name):
    if name in _SCKH_REEXPORTS:
        import importlib
        module_name, attr_name = _SCKH_REEXPORTS[name]
        value = getattr(importlib.import_module(module_name), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'dynamics_1d' has no attribute {name!r}")


__all__ = [
    # Constants
    "CONST",
    # Config
    "DynamicsConfig",
    "load_config",
    "InterpolationConfig",
    "SpectrumConfig",
    "FullConfig",
    "load_full_config",
    "save_full_config",
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
    # SCKH Trajectory (re-exported from sckh)
    "SCKHTrajectory",
    "sckh_trajectory_from_dynamics_1d",
    "read_sckh_trajectory",
    "write_sckh_trajectory",
    "read_sckh_trajectory_list",
    # Surface interpolation
    "SurfaceInterpolator",
    "SpectrumCalculator",
    # Workflow
    "SCKHInputs1D",
    "compute_sckh_trajectories_1d",
    # Spectrum (re-exported from sckh)
    "SpectrumResult",
    "SCKHSpectrumCalculator",
    "compute_energy_phase",
    "compute_F_if",
    "get_frequency_grid",
    "compute_spectrum_from_sckh",
]
