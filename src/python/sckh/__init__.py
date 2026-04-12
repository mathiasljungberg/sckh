"""
sckh: Semi-Classical Kramers-Heisenberg spectrum computation.

This package provides the shared infrastructure for SCKH spectrum
calculation: the SCKHTrajectory data type, FFT-based spectrum
computation, unified configuration, and SCKH-specific I/O.

Typical workflow::

    from sckh import SCKHSpectrumCalculator, SpectrumConfig, read_sckh_trajectory

    trajs = [read_sckh_trajectory(p) for p in paths]
    calc = SCKHSpectrumCalculator(spectrum_config=SpectrumConfig(gamma_fwhm=0.18))
    result = calc.compute_spectrum(trajs)
"""

from .trajectory import SCKHTrajectory
from .config import (
    SpectrumConfig,
    SCKHRunConfig,
    FullConfig,
    load_full_config,
    save_full_config,
)
from .spectrum import (
    SpectrumResult,
    SCKHSpectrumCalculator,
    compute_energy_phase,
    compute_F_if,
    get_frequency_grid,
    compute_spectrum_from_sckh,
)
from .io import (
    read_sckh_trajectory,
    read_sckh_trajectories,
    write_sckh_trajectory,
    read_sckh_trajectory_list,
    write_spectrum,
    write_spectrum_per_final,
    write_spectrum_result,
)

__all__ = [
    # Trajectory
    "SCKHTrajectory",
    # Config
    "SpectrumConfig",
    "SCKHRunConfig",
    "FullConfig",
    "load_full_config",
    "save_full_config",
    # Spectrum computation
    "SpectrumResult",
    "SCKHSpectrumCalculator",
    "compute_energy_phase",
    "compute_F_if",
    "get_frequency_grid",
    "compute_spectrum_from_sckh",
    # I/O
    "read_sckh_trajectory",
    "read_sckh_trajectories",
    "write_sckh_trajectory",
    "read_sckh_trajectory_list",
    "write_spectrum",
    "write_spectrum_per_final",
    "write_spectrum_result",
]
