"""Feature precomputation and cost matrix construction."""

from dataclasses import dataclass

import numpy as np

from .config import TrackingConfig


@dataclass
class Features:
    """Precomputed features for cost matrix construction."""

    dipole_norms: np.ndarray  # (n_x1, n_x2, n_states)
    osc_strengths: np.ndarray  # (n_x1, n_x2, n_states)
    confidence: np.ndarray  # (n_x1, n_x2, n_states)


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

    return Features(
        dipole_norms=dipole_norms,
        osc_strengths=osc_strengths,
        confidence=confidence,
    )


def build_cost_matrix(
    E_A: np.ndarray,
    D_A: np.ndarray,
    E_B: np.ndarray,
    D_B: np.ndarray,
    config: TrackingConfig,
    f_A: np.ndarray | None = None,
    f_B: np.ndarray | None = None,
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

    Returns
    -------
    C : ndarray, shape (n_states, n_states)
        Cost matrix where C[k, m] is the cost of assigning state k at A
        to state m at B.
    """
    n_states = len(E_A)
    C = np.zeros((n_states, n_states))

    # Energy cost: |E_A[k] - E_B[m]| / sigma_E
    dE = np.abs(E_A[:, None] - E_B[None, :])
    c_E = dE / config.sigma_E

    # Dipole similarity (phase-invariant): 1 - |D_A . D_B| / (|D_A| * |D_B| + eps)
    norm_A = np.linalg.norm(D_A, axis=-1)  # (n_states,)
    norm_B = np.linalg.norm(D_B, axis=-1)  # (n_states,)
    dot_AB = np.abs(D_A @ D_B.T)  # (n_states, n_states)
    denom = norm_A[:, None] * norm_B[None, :] + config.eps
    c_D = 1.0 - dot_AB / denom

    C = config.w_E * c_E + config.w_D * c_D

    # Optional oscillator strength cost
    if config.w_f > 0 and f_A is not None and f_B is not None:
        c_f = np.abs(f_A[:, None] - f_B[None, :]) / config.sigma_f
        C += config.w_f * c_f

    # Energy gate: block assignments with large energy difference
    if np.isfinite(config.E_gate):
        C[dE > config.E_gate] = 1e12

    return C
