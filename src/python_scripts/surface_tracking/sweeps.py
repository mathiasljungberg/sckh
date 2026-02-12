"""1D row/column sweep consistency passes (feature N)."""

import numpy as np
from scipy.optimize import linear_sum_assignment

from .bfs import _neighbor_cost, grid_neighbors
from .config import TrackingConfig
from .cost import Features, build_cost_matrix


def _sweep_1d(
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    indices: list[tuple[int, int]],
    config: TrackingConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Track states along a 1D path of grid points.

    Parameters
    ----------
    E, D : raw energies and dipoles (full grid).
    E_tracked, D_tracked : current tracked arrays (read-only reference).
    indices : ordered list of (i, j) grid coordinates forming the 1D path.
    config : TrackingConfig

    Returns
    -------
    perms : ndarray, shape (len(indices), n_states)
        Permutation at each point along the path.
    D_sweep : ndarray, shape (len(indices), n_states, 3)
        Phase-aligned dipoles along the path.
    """
    n_states = E.shape[2]
    n_pts = len(indices)
    perms = np.empty((n_pts, n_states), dtype=int)
    E_sweep = np.empty((n_pts, n_states))
    D_sweep = np.empty((n_pts, n_states, 3))

    # Initialize first point from current tracked assignment
    i0, j0 = indices[0]
    perms[0] = np.arange(n_states)  # identity relative to tracked
    E_sweep[0] = E_tracked[i0, j0]
    D_sweep[0] = D_tracked[i0, j0]

    for p in range(1, n_pts):
        ip, jp = indices[p]
        # Build cost from previous point in sweep to raw data at current point
        C = build_cost_matrix(
            E_sweep[p - 1], D_sweep[p - 1],
            E[ip, jp], D[ip, jp],
            config,
        )
        _, col_ind = linear_sum_assignment(C)

        perms[p] = col_ind
        E_sweep[p] = E[ip, jp, col_ind]
        D_sweep[p] = D[ip, jp, col_ind]

        # Inline phase alignment with previous point
        for k in range(n_states):
            if np.dot(D_sweep[p, k], D_sweep[p - 1, k]) < 0:
                D_sweep[p, k] *= -1

    return perms, D_sweep


def sweep_rows(
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    perm: np.ndarray,
    features: Features,
    config: TrackingConfig,
) -> int:
    """Sweep along rows (varying x1, fixed x2) and reconcile.

    For each row, does forward and backward 1D sweeps. At each point,
    compares the sweep result with the current assignment and takes the
    one with lower total neighbor cost.

    Modifies E_tracked, D_tracked, perm in-place. Returns total changes.
    """
    n_x1, n_x2, n_states = E.shape
    changes = 0

    for j in range(n_x2):
        row_indices = [(i, j) for i in range(n_x1)]

        # Forward sweep (i=0 → n_x1-1)
        fwd_perms, fwd_D = _sweep_1d(
            E, D, E_tracked, D_tracked, row_indices, config
        )
        # Backward sweep (i=n_x1-1 → 0)
        bwd_perms, bwd_D = _sweep_1d(
            E, D, E_tracked, D_tracked, row_indices[::-1], config
        )

        # Reconcile: at each point, pick best among current, forward, backward
        for p, (i, jj) in enumerate(row_indices):
            neighbors = grid_neighbors(i, jj, n_x1, n_x2)
            if not neighbors:
                continue

            current_cost = _neighbor_cost(
                E_tracked[i, jj], D_tracked[i, jj],
                E_tracked, D_tracked, neighbors, config,
            )

            # Forward candidate
            fwd_perm_global = fwd_perms[p]
            fwd_E_cand = E[i, jj, fwd_perm_global]
            fwd_D_cand = fwd_D[p].copy()
            # Phase-align with neighbor majority
            for k in range(n_states):
                agreement = sum(
                    np.dot(fwd_D_cand[k], D_tracked[ni, nj, k])
                    for ni, nj in neighbors
                )
                if agreement < 0:
                    fwd_D_cand[k] *= -1
            fwd_cost = _neighbor_cost(
                fwd_E_cand, fwd_D_cand,
                E_tracked, D_tracked, neighbors, config,
            )

            # Backward candidate
            bwd_idx = n_x1 - 1 - p
            bwd_perm_global = bwd_perms[bwd_idx]
            bwd_E_cand = E[i, jj, bwd_perm_global]
            bwd_D_cand = bwd_D[bwd_idx].copy()
            for k in range(n_states):
                agreement = sum(
                    np.dot(bwd_D_cand[k], D_tracked[ni, nj, k])
                    for ni, nj in neighbors
                )
                if agreement < 0:
                    bwd_D_cand[k] *= -1
            bwd_cost = _neighbor_cost(
                bwd_E_cand, bwd_D_cand,
                E_tracked, D_tracked, neighbors, config,
            )

            # Pick the best
            best_cost = current_cost
            best_E, best_D, best_perm = None, None, None

            if fwd_cost < best_cost:
                best_cost = fwd_cost
                best_E, best_D, best_perm = fwd_E_cand, fwd_D_cand, fwd_perm_global

            if bwd_cost < best_cost:
                best_cost = bwd_cost
                best_E, best_D, best_perm = bwd_E_cand, bwd_D_cand, bwd_perm_global

            if best_perm is not None:
                E_tracked[i, jj] = best_E
                D_tracked[i, jj] = best_D
                perm[i, jj] = best_perm
                changes += 1

    return changes


def sweep_cols(
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    perm: np.ndarray,
    features: Features,
    config: TrackingConfig,
) -> int:
    """Sweep along columns (varying x2, fixed x1) and reconcile.

    Same logic as sweep_rows but along the other axis.

    Modifies E_tracked, D_tracked, perm in-place. Returns total changes.
    """
    n_x1, n_x2, n_states = E.shape
    changes = 0

    for i in range(n_x1):
        col_indices = [(i, j) for j in range(n_x2)]

        # Forward sweep (j=0 → n_x2-1)
        fwd_perms, fwd_D = _sweep_1d(
            E, D, E_tracked, D_tracked, col_indices, config
        )
        # Backward sweep (j=n_x2-1 → 0)
        bwd_perms, bwd_D = _sweep_1d(
            E, D, E_tracked, D_tracked, col_indices[::-1], config
        )

        # Reconcile
        for p, (ii, j) in enumerate(col_indices):
            neighbors = grid_neighbors(ii, j, n_x1, n_x2)
            if not neighbors:
                continue

            current_cost = _neighbor_cost(
                E_tracked[ii, j], D_tracked[ii, j],
                E_tracked, D_tracked, neighbors, config,
            )

            # Forward candidate
            fwd_perm_global = fwd_perms[p]
            fwd_E_cand = E[ii, j, fwd_perm_global]
            fwd_D_cand = fwd_D[p].copy()
            for k in range(n_states):
                agreement = sum(
                    np.dot(fwd_D_cand[k], D_tracked[ni, nj, k])
                    for ni, nj in neighbors
                )
                if agreement < 0:
                    fwd_D_cand[k] *= -1
            fwd_cost = _neighbor_cost(
                fwd_E_cand, fwd_D_cand,
                E_tracked, D_tracked, neighbors, config,
            )

            # Backward candidate
            bwd_idx = n_x2 - 1 - p
            bwd_perm_global = bwd_perms[bwd_idx]
            bwd_E_cand = E[ii, j, bwd_perm_global]
            bwd_D_cand = bwd_D[bwd_idx].copy()
            for k in range(n_states):
                agreement = sum(
                    np.dot(bwd_D_cand[k], D_tracked[ni, nj, k])
                    for ni, nj in neighbors
                )
                if agreement < 0:
                    bwd_D_cand[k] *= -1
            bwd_cost = _neighbor_cost(
                bwd_E_cand, bwd_D_cand,
                E_tracked, D_tracked, neighbors, config,
            )

            best_cost = current_cost
            best_E, best_D, best_perm = None, None, None

            if fwd_cost < best_cost:
                best_cost = fwd_cost
                best_E, best_D, best_perm = fwd_E_cand, fwd_D_cand, fwd_perm_global

            if bwd_cost < best_cost:
                best_cost = bwd_cost
                best_E, best_D, best_perm = bwd_E_cand, bwd_D_cand, bwd_perm_global

            if best_perm is not None:
                E_tracked[ii, j] = best_E
                D_tracked[ii, j] = best_D
                perm[ii, j] = best_perm
                changes += 1

    return changes
