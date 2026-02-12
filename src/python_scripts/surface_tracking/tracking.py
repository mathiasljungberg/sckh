"""Public API: track_surfaces(), TrackingConfig, TrackingResult."""

from dataclasses import dataclass

import numpy as np

from .bfs import bfs_assign, repair_assignments, repair_reassignment
from .config import TrackingConfig
from .cost import precompute_features
from .diagnostics import DiagnosticReport, compute_diagnostics
from .phase import fix_phases_bfs, relax_phases
from .sweeps import sweep_cols, sweep_rows


@dataclass
class TrackingResult:
    """Result of surface tracking."""

    E: np.ndarray
    D: np.ndarray
    perm: np.ndarray
    cost_at_assignment: np.ndarray
    bfs_parent: np.ndarray
    diagnostics: DiagnosticReport


def _select_ref_point(E: np.ndarray, D: np.ndarray) -> tuple[int, int]:
    """Select the best reference point for BFS (feature F).

    Scores each grid point by:
    1. Mean dipole norm (higher = better, more reliable dipoles)
    2. Energy smoothness (lower second derivative = better)

    Returns the (i, j) with the best combined score.
    """
    n_x1, n_x2, n_states = E.shape

    # Mean dipole norm at each point
    dipole_norms = np.linalg.norm(D, axis=-1)  # (n_x1, n_x2, n_states)
    mean_norm = np.mean(dipole_norms, axis=-1)  # (n_x1, n_x2)

    # Energy roughness: mean absolute second derivative across states
    roughness = np.zeros((n_x1, n_x2))
    count = 0
    if n_x1 > 2:
        d2E_x1 = np.abs(E[2:, :, :] - 2 * E[1:-1, :, :] + E[:-2, :, :])
        # Pad edges with the nearest interior value
        r_x1 = np.mean(d2E_x1, axis=-1)
        padded = np.pad(r_x1, ((1, 1), (0, 0)), mode="edge")
        roughness += padded
        count += 1
    if n_x2 > 2:
        d2E_x2 = np.abs(E[:, 2:, :] - 2 * E[:, 1:-1, :] + E[:, :-2, :])
        r_x2 = np.mean(d2E_x2, axis=-1)
        padded = np.pad(r_x2, ((0, 0), (1, 1)), mode="edge")
        roughness += padded
        count += 1
    if count > 0:
        roughness /= count

    # Normalize both to [0, 1] range, then combine
    # Higher norm is better, lower roughness is better
    norm_range = mean_norm.max() - mean_norm.min()
    if norm_range > 0:
        norm_score = (mean_norm - mean_norm.min()) / norm_range
    else:
        norm_score = np.ones_like(mean_norm)

    rough_range = roughness.max() - roughness.min()
    if rough_range > 0:
        rough_score = 1.0 - (roughness - roughness.min()) / rough_range
    else:
        rough_score = np.ones_like(roughness)

    combined = norm_score + rough_score
    best = np.unravel_index(np.argmax(combined), combined.shape)
    return int(best[0]), int(best[1])


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

    # Feature F: smart reference point selection
    if config.ref_point == "auto":
        ref_i, ref_j = _select_ref_point(E, D)
        # Create a modified config with the resolved ref_point
        from dataclasses import replace

        config = replace(config, ref_point=(ref_i, ref_j))

    features = precompute_features(E, D, config)
    E_tracked, D_tracked, perm, bfs_order, bfs_parent, cost_at_assignment = (
        bfs_assign(E, D, features, config)
    )

    # Feature D: post-assignment swap repair (fast)
    if config.n_repair_iter > 0:
        repair_assignments(E_tracked, D_tracked, perm, features, config)

    # Feature M: full Hungarian re-assignment repair
    if config.n_reassign_iter > 0:
        repair_reassignment(
            E, D, E_tracked, D_tracked, perm, features, config
        )

    # Feature N: 1D row/column sweep consistency
    for _ in range(config.n_sweep_iter):
        row_changes = sweep_rows(
            E, D, E_tracked, D_tracked, perm, features, config
        )
        col_changes = sweep_cols(
            E, D, E_tracked, D_tracked, perm, features, config
        )
        if row_changes + col_changes == 0:
            break

    fix_phases_bfs(D_tracked, bfs_order, bfs_parent)

    # Feature L: iterative neighbor-consensus phase relaxation
    if config.n_phase_iter > 0:
        relax_phases(D_tracked, config.n_phase_iter)

    diagnostics = compute_diagnostics(E_tracked, D_tracked, cost_at_assignment)

    return TrackingResult(
        E=E_tracked,
        D=D_tracked,
        perm=perm,
        cost_at_assignment=cost_at_assignment,
        bfs_parent=bfs_parent,
        diagnostics=diagnostics,
    )
