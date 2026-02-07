"""Feature precomputation and cost matrix construction."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import subspace_angles

from .config import TrackingConfig


@dataclass
class Features:
    """Precomputed features for cost matrix construction."""

    dipole_norms: np.ndarray  # (n_x1, n_x2, n_states)
    osc_strengths: np.ndarray  # (n_x1, n_x2, n_states)
    confidence: np.ndarray  # (n_x1, n_x2, n_states)
    sigma_E_local: np.ndarray | None = None  # (n_x1, n_x2) - feature E


def precompute_features(
    E: np.ndarray, D: np.ndarray, config: TrackingConfig
) -> Features:
    """Precompute dipole norms, oscillator strengths, and confidence weights.

    Parameters
    ----------
    E : ndarray, shape (n_x1, n_x2, n_states)
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
    config : TrackingConfig

    Returns
    -------
    Features
    """
    dipole_norms = np.linalg.norm(D, axis=-1)  # (n_x1, n_x2, n_states)

    # Oscillator strengths: f = dE * |D|^2, where dE is relative to lowest state
    dE = E - E[:, :, 0:1]  # energy relative to lowest state
    osc_strengths = dE * dipole_norms**2

    confidence = np.minimum(1.0, dipole_norms / config.D_ref)

    # Feature E: adaptive sigma_E
    sigma_E_local = None
    if config.adaptive_sigma_E:
        # Median nearest-neighbor energy gap at each grid point
        sorted_E = np.sort(E, axis=-1)
        gaps = np.diff(sorted_E, axis=-1)  # (n_x1, n_x2, n_states-1)
        sigma_E_local = np.median(gaps, axis=-1)  # (n_x1, n_x2)
        # Clip to avoid extremes
        sigma_E_local = np.clip(
            sigma_E_local, config.sigma_E / 10.0, config.sigma_E * 10.0
        )

    return Features(
        dipole_norms=dipole_norms,
        osc_strengths=osc_strengths,
        confidence=confidence,
        sigma_E_local=sigma_E_local,
    )


def build_cost_matrix(
    E_A: np.ndarray,
    D_A: np.ndarray,
    E_B: np.ndarray,
    D_B: np.ndarray,
    config: TrackingConfig,
    f_A: np.ndarray | None = None,
    f_B: np.ndarray | None = None,
    conf_A: np.ndarray | None = None,
    conf_B: np.ndarray | None = None,
    sigma_E_override: float | None = None,
    grad_A: np.ndarray | None = None,
    grad_D_A: np.ndarray | None = None,
) -> np.ndarray:
    """Build cost matrix for state matching between two grid points.

    Parameters
    ----------
    E_A, E_B : ndarray, shape (n_states,)
        Energies at points A and B.
    D_A, D_B : ndarray, shape (n_states, 3)
        Dipole vectors at points A and B.
    config : TrackingConfig
    f_A, f_B : ndarray, shape (n_states,), optional
        Oscillator strengths at points A and B.
    conf_A, conf_B : ndarray, shape (n_states,), optional
        Confidence weights for dipole cost (feature A).
    sigma_E_override : float, optional
        Local sigma_E override (feature E: adaptive sigma_E).
    grad_A : ndarray, shape (n_states,), optional
        Energy gradient estimates at A (feature I).
    grad_D_A : ndarray, shape (n_states, 3), optional
        Dipole gradient estimates at A (feature J).

    Returns
    -------
    C : ndarray, shape (n_states, n_states)
        Cost matrix where C[k, m] is the cost of assigning state k at A
        to state m at B.
    """
    n_states = len(E_A)

    sigma_E = sigma_E_override if sigma_E_override is not None else config.sigma_E

    # Energy cost: |E_A[k] - E_B[m]| / sigma_E
    dE = np.abs(E_A[:, None] - E_B[None, :])
    c_E = dE / sigma_E

    # Dipole similarity (phase-invariant): 1 - |D_A . D_B| / (|D_A| * |D_B| + eps)
    norm_A = np.linalg.norm(D_A, axis=-1)  # (n_states,)
    norm_B = np.linalg.norm(D_B, axis=-1)  # (n_states,)
    dot_AB = np.abs(D_A @ D_B.T)  # (n_states, n_states)
    denom = norm_A[:, None] * norm_B[None, :] + config.eps
    c_D = 1.0 - dot_AB / denom

    # Feature A: confidence weighting on dipole cost
    if config.use_confidence and conf_A is not None and conf_B is not None:
        conf_weight = np.minimum(conf_A[:, None], conf_B[None, :])
        c_D = c_D * conf_weight

    C = config.w_E * c_E + config.w_D * c_D

    # Feature G: dipole-norm similarity cost
    if config.w_norm > 0:
        c_norm = np.abs(norm_A[:, None] - norm_B[None, :]) / config.sigma_norm
        C += config.w_norm * c_norm

    # Optional oscillator strength cost
    if config.w_f > 0 and f_A is not None and f_B is not None:
        c_f = np.abs(f_A[:, None] - f_B[None, :]) / config.sigma_f
        C += config.w_f * c_f

    # Feature I: energy-derivative continuity cost
    if config.w_dE > 0 and grad_A is not None:
        # For each candidate assignment (k -> m), the predicted energy at B
        # is E_A[k] + grad_A[k]. Penalize deviation from prediction.
        E_predicted = E_A[:, None] + grad_A[:, None]  # (n_states_A, 1)
        c_dE = np.abs(E_predicted - E_B[None, :]) / config.sigma_dE
        C += config.w_dE * c_dE

    # Feature J: dipole-derivative continuity cost
    if config.w_dD > 0 and grad_D_A is not None:
        # Predicted dipole at B from trend at A
        D_predicted = D_A + grad_D_A  # (n_states, 3)
        # Phase-invariant distance: min of ‖pred - D_B‖ and ‖pred + D_B‖
        diff_pos = D_predicted[:, None, :] - D_B[None, :, :]  # (n_s, n_s, 3)
        diff_neg = D_predicted[:, None, :] + D_B[None, :, :]
        dist = np.minimum(
            np.linalg.norm(diff_pos, axis=-1),
            np.linalg.norm(diff_neg, axis=-1),
        )
        C += config.w_dD * dist / config.sigma_dD

    # Energy gate: block assignments with large energy difference
    if np.isfinite(config.E_gate):
        C[dE > config.E_gate] = 1e12

    return C


def _find_degenerate_clusters(
    E: np.ndarray, threshold: float
) -> list[list[int]]:
    """Find clusters of near-degenerate states.

    Parameters
    ----------
    E : ndarray, shape (n_states,)
        Sorted energies at a single grid point.
    threshold : float
        Maximum energy gap to consider states as degenerate.

    Returns
    -------
    clusters : list of list of int
        Each inner list contains indices of states in a degenerate cluster.
        Singletons are included.
    """
    n_states = len(E)
    sorted_idx = np.argsort(E)
    sorted_E = E[sorted_idx]

    clusters = []
    current = [sorted_idx[0]]
    for i in range(1, n_states):
        if sorted_E[i] - sorted_E[i - 1] < threshold:
            current.append(sorted_idx[i])
        else:
            clusters.append(current)
            current = [sorted_idx[i]]
    clusters.append(current)
    return clusters


def build_cost_matrix_subspace(
    E_A: np.ndarray,
    D_A: np.ndarray,
    E_B: np.ndarray,
    D_B: np.ndarray,
    config: TrackingConfig,
    **kwargs,
) -> np.ndarray:
    """Build cost matrix with subspace tracking for degenerate clusters (feature H).

    For near-degenerate states, replaces individual dipole costs with
    subspace angle cost. Non-degenerate states use standard cost matrix.

    Parameters
    ----------
    E_A, E_B : ndarray, shape (n_states,)
    D_A, D_B : ndarray, shape (n_states, 3)
    config : TrackingConfig
    **kwargs : passed to build_cost_matrix

    Returns
    -------
    C : ndarray, shape (n_states, n_states)
    """
    # Start with standard cost matrix
    C = build_cost_matrix(E_A, D_A, E_B, D_B, config, **kwargs)

    # Find degenerate clusters at point B (the raw point being assigned)
    clusters_B = _find_degenerate_clusters(E_B, config.degen_threshold)

    for cluster_B in clusters_B:
        if len(cluster_B) < 2:
            continue

        # Find the corresponding cluster at A: states in A whose energies
        # are within degen_threshold of any state in cluster_B
        cluster_B_energies = E_B[cluster_B]
        e_min = cluster_B_energies.min() - config.degen_threshold
        e_max = cluster_B_energies.max() + config.degen_threshold
        cluster_A = [k for k in range(len(E_A)) if e_min <= E_A[k] <= e_max]

        if len(cluster_A) < 2 or len(cluster_A) != len(cluster_B):
            continue

        # Compute subspace angle between the cluster dipole subspaces
        D_sub_A = D_A[cluster_A]  # (cluster_size, 3)
        D_sub_B = D_B[cluster_B]  # (cluster_size, 3)

        # Check norms - skip if too small
        if (np.linalg.norm(D_sub_A) < config.eps or
                np.linalg.norm(D_sub_B) < config.eps):
            continue

        try:
            angles = subspace_angles(D_sub_A.T, D_sub_B.T)
            subspace_cost = np.mean(angles) / (np.pi / 2)  # normalize to [0, 1]
        except (np.linalg.LinAlgError, ValueError):
            continue

        # Replace individual dipole costs within the cluster with
        # uniform subspace cost + small energy-ordering tiebreaker
        for ka_idx, ka in enumerate(cluster_A):
            for kb_idx, kb in enumerate(cluster_B):
                # Subspace cost (same for all pairs in cluster)
                c_sub = config.w_D * subspace_cost
                # Small tiebreaker: prefer preserving energy ordering
                tiebreaker = 0.01 * abs(ka_idx - kb_idx)
                # Replace dipole cost component, keep energy cost
                dE = abs(E_A[ka] - E_B[kb])
                sigma_E = kwargs.get("sigma_E_override", None) or config.sigma_E
                C[ka, kb] = config.w_E * dE / sigma_E + c_sub + tiebreaker

    return C
