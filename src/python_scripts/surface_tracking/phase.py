"""Phase/sign fixing via parallel transport."""

import numpy as np

from .bfs import grid_neighbors


def fix_phases_bfs(
    D: np.ndarray,
    bfs_order: list[tuple[int, int]],
    bfs_parent: np.ndarray,
    min_norm: float = 1e-10,
) -> np.ndarray:
    """Fix dipole sign ambiguity using parallel transport along BFS tree.

    For each state at each grid point, ensures the dipole vector has
    consistent sign relative to its BFS parent by flipping if the dot
    product with the parent's dipole is negative.

    Parameters
    ----------
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
        Tracked dipole vectors. Modified in-place.
    bfs_order : list of (int, int)
        Grid points in BFS visitation order.
    bfs_parent : ndarray, shape (n_x1, n_x2, 2)
        Parent indices for each grid point.
    min_norm : float
        Skip phase fixing if either dipole norm is below this threshold.

    Returns
    -------
    D : ndarray
        The same array, modified in-place.
    """
    n_states = D.shape[2]

    # Skip the reference point (first in BFS order)
    for bi, bj in bfs_order[1:]:
        ai, aj = bfs_parent[bi, bj]

        for k in range(n_states):
            d_A = D[ai, aj, k]
            d_B = D[bi, bj, k]

            if np.linalg.norm(d_A) < min_norm or np.linalg.norm(d_B) < min_norm:
                continue

            if np.dot(d_A, d_B) < 0:
                D[bi, bj, k] *= -1

    return D


def relax_phases(
    D: np.ndarray,
    n_iter: int,
    min_norm: float = 1e-10,
) -> int:
    """Iterative neighbor-consensus phase relaxation (feature L).

    For each grid point and state, flips the dipole sign if the total
    dot-product agreement with all grid neighbors is negative.  Repeats
    until convergence or *n_iter* iterations.

    Parameters
    ----------
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
        Tracked dipole vectors. Modified in-place.
    n_iter : int
        Maximum number of relaxation iterations.
    min_norm : float
        Skip dipoles with norm below this threshold.

    Returns
    -------
    total_flips : int
        Total number of sign flips applied across all iterations.
    """
    n_x1, n_x2, n_states = D.shape[:3]
    total_flips = 0

    for _ in range(n_iter):
        flips = 0
        for i in range(n_x1):
            for j in range(n_x2):
                for k in range(n_states):
                    d = D[i, j, k]
                    if np.linalg.norm(d) < min_norm:
                        continue
                    agreement = 0.0
                    for ni, nj in grid_neighbors(i, j, n_x1, n_x2):
                        agreement += np.dot(d, D[ni, nj, k])
                    if agreement < 0:
                        D[i, j, k] *= -1
                        flips += 1
        total_flips += flips
        if flips == 0:
            break

    return total_flips
