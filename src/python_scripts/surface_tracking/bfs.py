"""BFS traversal and Hungarian assignment for surface tracking."""

from collections import deque

import numpy as np
from scipy.optimize import linear_sum_assignment

from .cost import Features, build_cost_matrix
from .config import TrackingConfig


def grid_neighbors(
    i: int, j: int, n_x1: int, n_x2: int
) -> list[tuple[int, int]]:
    """Return 4-connected neighbors of (i, j) within grid bounds.

    Parameters
    ----------
    i, j : int
        Grid indices.
    n_x1, n_x2 : int
        Grid dimensions.

    Returns
    -------
    list of (int, int)
    """
    neighbors = []
    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        ni, nj = i + di, j + dj
        if 0 <= ni < n_x1 and 0 <= nj < n_x2:
            neighbors.append((ni, nj))
    return neighbors


def bfs_assign(
    E: np.ndarray,
    D: np.ndarray,
    features: Features,
    config: TrackingConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list, np.ndarray, np.ndarray]:
    """Assign states across grid using BFS and Hungarian algorithm.

    Parameters
    ----------
    E : ndarray, shape (n_x1, n_x2, n_states)
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
    features : Features
    config : TrackingConfig

    Returns
    -------
    E_tracked : ndarray, shape (n_x1, n_x2, n_states)
    D_tracked : ndarray, shape (n_x1, n_x2, n_states, 3)
    perm : ndarray, shape (n_x1, n_x2, n_states), dtype int
    bfs_order : list of (int, int)
    bfs_parent : ndarray, shape (n_x1, n_x2, 2), dtype int
    cost_at_assignment : ndarray, shape (n_x1, n_x2)
    """
    n_x1, n_x2, n_states = E.shape

    E_tracked = np.empty_like(E)
    D_tracked = np.empty_like(D)
    perm = np.empty((n_x1, n_x2, n_states), dtype=int)
    cost_at_assignment = np.zeros((n_x1, n_x2))
    bfs_parent = np.full((n_x1, n_x2, 2), -1, dtype=int)

    visited = np.zeros((n_x1, n_x2), dtype=bool)

    # Determine reference point
    if config.ref_point is not None:
        ri, rj = config.ref_point
    else:
        ri, rj = n_x1 // 2, n_x2 // 2

    # Initialize reference point with identity permutation
    E_tracked[ri, rj] = E[ri, rj]
    D_tracked[ri, rj] = D[ri, rj]
    perm[ri, rj] = np.arange(n_states)
    visited[ri, rj] = True

    bfs_order = [(ri, rj)]
    queue = deque([(ri, rj)])

    while queue:
        ai, aj = queue.popleft()

        for bi, bj in grid_neighbors(ai, aj, n_x1, n_x2):
            if visited[bi, bj]:
                continue
            visited[bi, bj] = True

            # Build cost: tracked states at A vs raw states at B
            f_A = features.osc_strengths[ai, aj] if config.w_f > 0 else None
            f_B = features.osc_strengths[bi, bj] if config.w_f > 0 else None

            C = build_cost_matrix(
                E_tracked[ai, aj],
                D_tracked[ai, aj],
                E[bi, bj],
                D[bi, bj],
                config,
                f_A=f_A,
                f_B=f_B,
            )

            row_ind, col_ind = linear_sum_assignment(C)
            cost_at_assignment[bi, bj] = C[row_ind, col_ind].sum()

            # col_ind[k] = raw state index at B that matches tracked state k at A
            E_tracked[bi, bj] = E[bi, bj, col_ind]
            D_tracked[bi, bj] = D[bi, bj, col_ind]
            perm[bi, bj] = col_ind

            bfs_parent[bi, bj] = [ai, aj]
            bfs_order.append((bi, bj))
            queue.append((bi, bj))

    return E_tracked, D_tracked, perm, bfs_order, bfs_parent, cost_at_assignment
