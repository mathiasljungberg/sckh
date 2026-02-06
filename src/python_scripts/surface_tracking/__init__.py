"""Surface tracking for 2D potential energy surfaces.

Track states by character continuity across a 2D grid using energies and
transition dipole vectors, producing smooth, consistently labeled surfaces.
"""

from .config import TrackingConfig
from .tracking import TrackingResult, track_surfaces

__all__ = ["track_surfaces", "TrackingConfig", "TrackingResult"]
