"""Backwards-compatibility re-exports.

Surface interpolation has moved to dynamics_2d.interpolation.
Spectrum computation has moved to sckh.spectrum.
"""

from .interpolation import SurfaceInterpolator2D, SpectrumCalculator2D

__all__ = [
    "SurfaceInterpolator2D",
    "SpectrumCalculator2D",
]
