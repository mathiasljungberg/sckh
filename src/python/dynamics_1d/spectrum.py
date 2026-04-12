"""Backwards-compatibility re-exports.

Surface interpolation has moved to dynamics_1d.interpolation.
Spectrum computation has moved to sckh.spectrum.
"""

# Surface interpolation
from .interpolation import SurfaceInterpolator, SpectrumCalculator

# Spectrum computation (from sckh package)
from sckh.spectrum import (
    SpectrumResult,
    SCKHSpectrumCalculator,
    compute_energy_phase,
    compute_F_if,
    get_frequency_grid,
    compute_spectrum_from_sckh,
)

# Re-export SCKHTrajectory for callers that import it from here
from sckh.trajectory import SCKHTrajectory

# Re-export configs for callers that import them from here
from sckh.config import FullConfig, SpectrumConfig
from .spectrum_config import InterpolationConfig

__all__ = [
    "SurfaceInterpolator",
    "SpectrumCalculator",
    "SpectrumResult",
    "SCKHSpectrumCalculator",
    "compute_energy_phase",
    "compute_F_if",
    "get_frequency_grid",
    "compute_spectrum_from_sckh",
    "SCKHTrajectory",
    "FullConfig",
    "SpectrumConfig",
    "InterpolationConfig",
]
