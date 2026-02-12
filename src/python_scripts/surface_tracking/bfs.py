"""BFS traversal and Hungarian assignment for surface tracking."""

import heapq
from collections import deque

import numpy as np
from scipy.optimize import linear_sum_assignment

from .config import TrackingConfig
from .cost import Features, build_cost_matrix, build_cost_matrix_subspace


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


def _build_cost_for_pair(
    ai: int,
    aj: int,
    bi: int,
    bj: int,
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    features: Features,
    config: TrackingConfig,
    energy_gradients: np.ndarray | None = None,
    dipole_gradients: np.ndarray | None = None,
) -> np.ndarray:
    """Build cost matrix for assigning raw states at B from tracked states at A."""
    f_A = features.osc_strengths[ai, aj] if config.w_f > 0 else None
    f_B = features.osc_strengths[bi, bj] if config.w_f > 0 else None
    conf_A = features.confidence[ai, aj] if config.use_confidence else None
    conf_B = features.confidence[bi, bj] if config.use_confidence else None

    sigma_E_local = None
    if config.adaptive_sigma_E and features.sigma_E_local is not None:
        sigma_E_local = 0.5 * (
            features.sigma_E_local[ai, aj] + features.sigma_E_local[bi, bj]
        )

    # Feature I: energy gradient at A
    grad_A = None
    if config.w_dE > 0 and energy_gradients is not None:
        grad_A = energy_gradients[ai, aj]

    # Feature J: dipole gradient at A
    grad_D_A = None
    if config.w_dD > 0 and dipole_gradients is not None:
        grad_D_A = dipole_gradients[ai, aj]

    cost_fn = build_cost_matrix_subspace if config.use_subspace else build_cost_matrix
    return cost_fn(
        E_tracked[ai, aj],
        D_tracked[ai, aj],
        E[bi, bj],
        D[bi, bj],
        config,
        f_A=f_A,
        f_B=f_B,
        conf_A=conf_A,
        conf_B=conf_B,
        sigma_E_override=sigma_E_local,
        grad_A=grad_A,
        grad_D_A=grad_D_A,
    )


def _assign_point(
    bi: int,
    bj: int,
    parent_ai: int,
    parent_aj: int,
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    features: Features,
    config: TrackingConfig,
    visited: np.ndarray,
    energy_gradients: np.ndarray | None = None,
    dipole_gradients: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    """Compute the optimal permutation for point B.

    If use_multi_neighbor is True, averages cost matrices from all
    already-visited neighbors. Otherwise uses only the BFS parent.

    Returns (col_ind, total_cost).
    """
    n_x1, n_x2 = E.shape[0], E.shape[1]

    if config.use_multi_neighbor:
        # Collect cost matrices from all visited neighbors
        neighbors = grid_neighbors(bi, bj, n_x1, n_x2)
        visited_neighbors = [(ni, nj) for ni, nj in neighbors if visited[ni, nj]]
        if not visited_neighbors:
            visited_neighbors = [(parent_ai, parent_aj)]

        C_sum = None
        for ni, nj in visited_neighbors:
            C_nb = _build_cost_for_pair(
                ni, nj, bi, bj, E, D, E_tracked, D_tracked, features, config,
                energy_gradients=energy_gradients,
                dipole_gradients=dipole_gradients,
            )
            if C_sum is None:
                C_sum = C_nb.copy()
            else:
                C_sum += C_nb
        C = C_sum / len(visited_neighbors)
    else:
        C = _build_cost_for_pair(
            parent_ai, parent_aj, bi, bj,
            E, D, E_tracked, D_tracked, features, config,
            energy_gradients=energy_gradients,
            dipole_gradients=dipole_gradients,
        )

    row_ind, col_ind = linear_sum_assignment(C)
    total_cost = C[row_ind, col_ind].sum()
    return col_ind, total_cost


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

    # Feature I: energy gradients (estimated from parent->child differences)
    energy_gradients = None
    if config.w_dE > 0:
        energy_gradients = np.zeros((n_x1, n_x2, n_states))

    # Feature J: dipole gradients (estimated from parent->child differences)
    dipole_gradients = None
    if config.w_dD > 0:
        dipole_gradients = np.zeros((n_x1, n_x2, n_states, 3))

    # Determine reference point (feature F handled in tracking.py)
    if config.ref_point is not None and config.ref_point != "auto":
        ri, rj = config.ref_point
    else:
        ri, rj = n_x1 // 2, n_x2 // 2

    # Initialize reference point with identity permutation
    E_tracked[ri, rj] = E[ri, rj]
    D_tracked[ri, rj] = D[ri, rj]
    perm[ri, rj] = np.arange(n_states)
    visited[ri, rj] = True

    bfs_order = [(ri, rj)]

    def _finalize_point(bi, bj, ai, aj, col_ind, total_cost):
        """Apply assignment and update gradients for a newly visited point."""
        E_tracked[bi, bj] = E[bi, bj, col_ind]
        D_tracked[bi, bj] = D[bi, bj, col_ind]
        perm[bi, bj] = col_ind
        cost_at_assignment[bi, bj] = total_cost
        bfs_parent[bi, bj] = [ai, aj]
        bfs_order.append((bi, bj))

        # Feature K: inline phase alignment with parent
        # Align dipole signs before computing gradients so that
        # D_tracked is phase-consistent throughout BFS.
        for k in range(n_states):
            if np.dot(D_tracked[bi, bj, k], D_tracked[ai, aj, k]) < 0:
                D_tracked[bi, bj, k] *= -1

        # Feature I: update energy gradient estimate at B
        if energy_gradients is not None:
            energy_gradients[bi, bj] = E_tracked[bi, bj] - E_tracked[ai, aj]

        # Feature J: update dipole gradient estimate at B
        if dipole_gradients is not None:
            dipole_gradients[bi, bj] = D_tracked[bi, bj] - D_tracked[ai, aj]

    if config.use_priority_queue:
        # Feature B: priority queue keyed by assignment cost
        counter = 0
        heap: list[tuple[float, int, int, int, int, int]] = []
        for ni, nj in grid_neighbors(ri, rj, n_x1, n_x2):
            _, cost_tent = _assign_point(
                ni, nj, ri, rj, E, D, E_tracked, D_tracked,
                features, config, visited, energy_gradients, dipole_gradients,
            )
            heapq.heappush(heap, (cost_tent, counter, ni, nj, ri, rj))
            counter += 1

        while heap:
            _, _, bi, bj, ai, aj = heapq.heappop(heap)
            if visited[bi, bj]:
                continue
            visited[bi, bj] = True

            col_ind, total_cost = _assign_point(
                bi, bj, ai, aj, E, D, E_tracked, D_tracked,
                features, config, visited, energy_gradients, dipole_gradients,
            )
            _finalize_point(bi, bj, ai, aj, col_ind, total_cost)

            for ni, nj in grid_neighbors(bi, bj, n_x1, n_x2):
                if not visited[ni, nj]:
                    _, cost_tent = _assign_point(
                        ni, nj, bi, bj, E, D, E_tracked, D_tracked,
                        features, config, visited, energy_gradients, dipole_gradients,
                    )
                    heapq.heappush(heap, (cost_tent, counter, ni, nj, bi, bj))
                    counter += 1
    else:
        # Standard FIFO BFS
        queue = deque([(ri, rj)])

        while queue:
            ai, aj = queue.popleft()

            for bi, bj in grid_neighbors(ai, aj, n_x1, n_x2):
                if visited[bi, bj]:
                    continue
                visited[bi, bj] = True

                col_ind, total_cost = _assign_point(
                    bi, bj, ai, aj, E, D, E_tracked, D_tracked,
                    features, config, visited, energy_gradients, dipole_gradients,
                )
                _finalize_point(bi, bj, ai, aj, col_ind, total_cost)
                queue.append((bi, bj))

    return E_tracked, D_tracked, perm, bfs_order, bfs_parent, cost_at_assignment


def repair_assignments(
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    perm: np.ndarray,
    features: Features,
    config: TrackingConfig,
) -> int:
    """Post-assignment repair pass (feature D).

    Iterates over all grid points and checks if swapping any pair of
    state labels would reduce the total cost with all neighbors.
    Repeats for config.n_repair_iter iterations or until convergence.

    Modifies E_tracked, D_tracked, perm in-place.

    Parameters
    ----------
    E_tracked : ndarray, shape (n_x1, n_x2, n_states)
    D_tracked : ndarray, shape (n_x1, n_x2, n_states, 3)
    perm : ndarray, shape (n_x1, n_x2, n_states)
    features : Features
    config : TrackingConfig

    Returns
    -------
    total_swaps : int
        Total number of swaps applied across all iterations.
    """
    n_x1, n_x2, n_states = E_tracked.shape
    total_swaps = 0

    # Precompute state pair indices
    s1_idx, s2_idx = np.triu_indices(n_states, k=1)

    for _ in range(config.n_repair_iter):
        swaps_this_iter = 0

        for i in range(n_x1):
            for j in range(n_x2):
                neighbors = grid_neighbors(i, j, n_x1, n_x2)
                if not neighbors:
                    continue

                improvement = _swap_improvements_vectorized(
                    i, j, E_tracked, D_tracked, neighbors,
                    config, s1_idx, s2_idx,
                )

                best_pair = np.argmax(improvement)
                if improvement[best_pair] > 0:
                    s1, s2 = s1_idx[best_pair], s2_idx[best_pair]
                    # Apply swap
                    E_tracked[i, j, [s1, s2]] = E_tracked[i, j, [s2, s1]]
                    D_tracked[i, j, [s1, s2]] = D_tracked[i, j, [s2, s1]]
                    perm[i, j, [s1, s2]] = perm[i, j, [s2, s1]]
                    swaps_this_iter += 1

        total_swaps += swaps_this_iter
        if swaps_this_iter == 0:
            break

    return total_swaps


def _swap_improvements_vectorized(
    i: int,
    j: int,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    neighbors: list[tuple[int, int]],
    config: TrackingConfig,
    s1_idx: np.ndarray,
    s2_idx: np.ndarray,
) -> np.ndarray:
    """Compute cost improvement for all state-pair swaps at point (i,j).

    Vectorized over state pairs and neighbors.

    Returns array of shape (n_pairs,) where positive = swap reduces cost.
    """
    E_ij = E_tracked[i, j]  # (n_states,)
    D_ij = D_tracked[i, j]  # (n_states, 3)

    # Stack neighbor data: (n_neighbors, n_states) and (n_neighbors, n_states, 3)
    nb_indices = np.array(neighbors)
    E_nb = E_tracked[nb_indices[:, 0], nb_indices[:, 1]]  # (n_nb, n_states)
    D_nb = D_tracked[nb_indices[:, 0], nb_indices[:, 1]]  # (n_nb, n_states, 3)

    # Per-state energy cost between (i,j) and each neighbor: (n_nb, n_states)
    dE_per_state = np.abs(E_ij[None, :] - E_nb) / config.sigma_E

    # Per-state dipole cost between (i,j) and each neighbor: (n_nb, n_states)
    norm_ij = np.linalg.norm(D_ij, axis=-1)  # (n_states,)
    norm_nb = np.linalg.norm(D_nb, axis=-1)  # (n_nb, n_states)
    # dot product: sum over xyz dim
    dot_vals = np.abs(np.sum(D_ij[None, :, :] * D_nb, axis=-1))  # (n_nb, n_states)
    denom = norm_ij[None, :] * norm_nb + config.eps
    dD_per_state = 1.0 - dot_vals / denom  # (n_nb, n_states)

    # Total per-state cost: (n_nb, n_states)
    cost_per_state = config.w_E * dE_per_state + config.w_D * dD_per_state

    # Sum over neighbors: (n_states,)
    cost_by_state = cost_per_state.sum(axis=0)

    # For each swap (s1, s2): cost_before = cost[s1] + cost[s2]
    cost_before = cost_by_state[s1_idx] + cost_by_state[s2_idx]  # (n_pairs,)

    # Cost after swap: s1 slot gets s2 data compared to s1 at neighbors, and vice versa
    # "s1 slot with s2 data vs neighbor s1" = cross-cost(s2_data, s1_neighbor)
    # We need cross-costs: cost of state a at (i,j) vs state b at neighbor

    # Energy cross-costs for swap pairs
    # s1_idx data swapped to s2_idx position, compared to s1_idx at neighbors
    dE_s2_vs_s1 = np.abs(E_ij[s2_idx][None, :] - E_nb[:, s1_idx]) / config.sigma_E  # (n_nb, n_pairs)
    dE_s1_vs_s2 = np.abs(E_ij[s1_idx][None, :] - E_nb[:, s2_idx]) / config.sigma_E  # (n_nb, n_pairs)

    # Dipole cross-costs
    # D_ij[s2] vs D_nb[:, s1]
    dot_s2_vs_s1 = np.abs(np.sum(D_ij[s2_idx][None, :, :] * D_nb[:, s1_idx, :], axis=-1))
    denom_s2_vs_s1 = norm_ij[s2_idx][None, :] * norm_nb[:, s1_idx] + config.eps
    dD_s2_vs_s1 = 1.0 - dot_s2_vs_s1 / denom_s2_vs_s1

    dot_s1_vs_s2 = np.abs(np.sum(D_ij[s1_idx][None, :, :] * D_nb[:, s2_idx, :], axis=-1))
    denom_s1_vs_s2 = norm_ij[s1_idx][None, :] * norm_nb[:, s2_idx] + config.eps
    dD_s1_vs_s2 = 1.0 - dot_s1_vs_s2 / denom_s1_vs_s2

    cost_after_per_nb = (
        config.w_E * (dE_s2_vs_s1 + dE_s1_vs_s2)
        + config.w_D * (dD_s2_vs_s1 + dD_s1_vs_s2)
    )  # (n_nb, n_pairs)
    cost_after = cost_after_per_nb.sum(axis=0)  # (n_pairs,)

    return cost_before - cost_after


def repair_reassignment(
    E: np.ndarray,
    D: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    perm: np.ndarray,
    features: Features,
    config: TrackingConfig,
) -> int:
    """Full Hungarian re-assignment repair (feature M).

    For each grid point, builds an averaged cost matrix from ALL neighbors
    against the raw data at that point, runs Hungarian, and accepts the
    new permutation if it reduces the total neighbor cost.

    Modifies E_tracked, D_tracked, perm in-place.

    Parameters
    ----------
    E : ndarray, shape (n_x1, n_x2, n_states)
        Raw energies.
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
        Raw dipoles.
    E_tracked, D_tracked : ndarray
        Tracked arrays (modified in-place).
    perm : ndarray
        Permutation map (modified in-place).
    features : Features
    config : TrackingConfig

    Returns
    -------
    total_changes : int
        Total number of re-assignments across all iterations.
    """
    n_x1, n_x2, n_states = E.shape
    total_changes = 0

    for _ in range(config.n_reassign_iter):
        changes = 0

        for i in range(n_x1):
            for j in range(n_x2):
                neighbors = grid_neighbors(i, j, n_x1, n_x2)
                if not neighbors:
                    continue

                # Build averaged cost matrix: neighbors → raw data at (i,j)
                C_avg = np.zeros((n_states, n_states))
                for ni, nj in neighbors:
                    C_nb = build_cost_matrix(
                        E_tracked[ni, nj], D_tracked[ni, nj],
                        E[i, j], D[i, j], config,
                    )
                    C_avg += C_nb
                C_avg /= len(neighbors)

                _, new_perm = linear_sum_assignment(C_avg)

                if np.array_equal(new_perm, perm[i, j]):
                    continue

                # Evaluate: is the new assignment better?
                # Compute total per-state neighbor cost for current and new
                new_E = E[i, j, new_perm]
                new_D = D[i, j, new_perm].copy()

                # Phase-align new dipoles with neighbor majority
                for k in range(n_states):
                    agreement = sum(
                        np.dot(new_D[k], D_tracked[ni, nj, k])
                        for ni, nj in neighbors
                    )
                    if agreement < 0:
                        new_D[k] *= -1

                current_cost = _neighbor_cost(
                    E_tracked[i, j], D_tracked[i, j],
                    E_tracked, D_tracked, neighbors, config,
                )
                new_cost = _neighbor_cost(
                    new_E, new_D,
                    E_tracked, D_tracked, neighbors, config,
                )

                if new_cost < current_cost:
                    E_tracked[i, j] = new_E
                    D_tracked[i, j] = new_D
                    perm[i, j] = new_perm
                    changes += 1

        total_changes += changes
        if changes == 0:
            break

    return total_changes


def _neighbor_cost(
    E_point: np.ndarray,
    D_point: np.ndarray,
    E_tracked: np.ndarray,
    D_tracked: np.ndarray,
    neighbors: list[tuple[int, int]],
    config: TrackingConfig,
) -> float:
    """Compute total per-state cost of a point against its neighbors."""
    total = 0.0
    for ni, nj in neighbors:
        dE = np.abs(E_point - E_tracked[ni, nj]) / config.sigma_E

        norm_p = np.linalg.norm(D_point, axis=-1)
        norm_n = np.linalg.norm(D_tracked[ni, nj], axis=-1)
        dot_val = np.abs(np.sum(D_point * D_tracked[ni, nj], axis=-1))
        denom = norm_p * norm_n + config.eps
        dD = 1.0 - dot_val / denom

        total += np.sum(config.w_E * dE + config.w_D * dD)
    return total
