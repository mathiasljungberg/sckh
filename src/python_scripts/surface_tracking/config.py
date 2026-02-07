"""Configuration dataclass for surface tracking."""

from dataclasses import dataclass

import numpy as np


@dataclass
class TrackingConfig:
    """Configuration for surface tracking."""

    # --- Original parameters ---
    w_E: float = 1.0
    w_D: float = 1.0
    w_f: float = 0.0
    sigma_E: float = 0.1
    sigma_f: float = 0.01
    E_gate: float = np.inf
    D_ref: float = 0.1
    eps: float = 1e-12
    ref_point: tuple[int, int] | str | None = None

    # --- Phase 1 toggles ---
    # A: Confidence weighting - downweight dipole cost for weak dipoles
    use_confidence: bool = False
    # B: Priority-queue BFS - visit low-cost edges first
    use_priority_queue: bool = False
    # C: Multi-neighbor consensus - average cost from all visited neighbors
    use_multi_neighbor: bool = False

    # --- Phase 2 toggles ---
    # D: Post-assignment repair pass
    n_repair_iter: int = 0
    # E: Adaptive sigma_E from local energy gaps
    adaptive_sigma_E: bool = False

    # --- Phase 3 toggles ---
    # G: Dipole-norm similarity cost
    w_norm: float = 0.0
    sigma_norm: float = 0.1
    # H: Subspace tracking for degenerate clusters
    use_subspace: bool = False
    degen_threshold: float = 0.01
    # I: Energy-derivative continuity cost
    w_dE: float = 0.0
    sigma_dE: float = 0.1
    # J: Dipole-derivative continuity cost
    w_dD: float = 0.0
    sigma_dD: float = 0.1
