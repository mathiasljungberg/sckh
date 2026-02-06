"""Public API: track_surfaces(), TrackingConfig, TrackingResult."""

from dataclasses import dataclass

import numpy as np

from .bfs import bfs_assign
from .config import TrackingConfig
from .cost import precompute_features
from .diagnostics import DiagnosticReport, compute_diagnostics
from .phase import fix_phases_bfs


@dataclass
class TrackingResult:
    """Result of surface tracking."""

    E: np.ndarray
    D: np.ndarray
    perm: np.ndarray
    cost_at_assignment: np.ndarray
    bfs_parent: np.ndarray
    diagnostics: DiagnosticReport


def _validate_inputs(E: np.ndarray, D: np.ndarray) -> None:
    """Validate input array shapes."""
    if E.ndim != 3:
        raise ValueError(
            f"E must be 3D (n_x1, n_x2, n_states), got shape {E.shape}"
        )
    if D.ndim != 4:
        raise ValueError(
            f"D must be 4D (n_x1, n_x2, n_states, 3), got shape {D.shape}"
        )
    n_x1, n_x2, n_states = E.shape
    if D.shape != (n_x1, n_x2, n_states, 3):
        raise ValueError(
            f"D shape {D.shape} incompatible with E shape {E.shape}; "
            f"expected ({n_x1}, {n_x2}, {n_states}, 3)"
        )


def track_surfaces(
    E: np.ndarray,
    D: np.ndarray,
    config: TrackingConfig | None = None,
) -> TrackingResult:
    """Track electronic states across a 2D grid by character continuity.

    Parameters
    ----------
    E : ndarray, shape (n_x1, n_x2, n_states)
        Energies at each grid point, ordered by energy.
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
        Transition dipole vectors at each grid point.
    config : TrackingConfig, optional
        Tracking parameters. Uses defaults if None.

    Returns
    -------
    TrackingResult
        Tracked energies, phase-regularized dipoles, permutation map,
        assignment costs, BFS parent map, and diagnostic report.
    """
    if config is None:
        config = TrackingConfig()

    _validate_inputs(E, D)

    features = precompute_features(E, D, config)
    E_tracked, D_tracked, perm, bfs_order, bfs_parent, cost_at_assignment = (
        bfs_assign(E, D, features, config)
    )
    fix_phases_bfs(D_tracked, bfs_order, bfs_parent)
    diagnostics = compute_diagnostics(E_tracked, D_tracked, cost_at_assignment)

    return TrackingResult(
        E=E_tracked,
        D=D_tracked,
        perm=perm,
        cost_at_assignment=cost_at_assignment,
        bfs_parent=bfs_parent,
        diagnostics=diagnostics,
    )
