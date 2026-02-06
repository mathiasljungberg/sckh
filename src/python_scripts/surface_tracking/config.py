"""Configuration dataclass for surface tracking."""

from dataclasses import dataclass

import numpy as np


@dataclass
class TrackingConfig:
    """Configuration for surface tracking."""

    w_E: float = 1.0
    w_D: float = 1.0
    w_f: float = 0.0
    sigma_E: float = 0.1
    sigma_f: float = 0.01
    E_gate: float = np.inf
    D_ref: float = 0.1
    eps: float = 1e-12
    ref_point: tuple[int, int] | None = None
