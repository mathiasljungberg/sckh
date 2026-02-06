"""Phase/sign fixing via parallel transport."""

import numpy as np


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
