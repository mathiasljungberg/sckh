"""Mismatch scores and smoothness metrics."""

from dataclasses import dataclass

import numpy as np


@dataclass
class DiagnosticReport:
    """Diagnostic metrics for tracked surfaces."""

    neighbor_energy_mismatch: float
    neighbor_dipole_mismatch: float
    per_state_energy_roughness: np.ndarray  # (n_states,)
    max_assignment_cost: float


def compute_diagnostics(
    E: np.ndarray,
    D: np.ndarray,
    cost_at_assignment: np.ndarray,
) -> DiagnosticReport:
    """Compute smoothness and mismatch diagnostics for tracked surfaces.

    Parameters
    ----------
    E : ndarray, shape (n_x1, n_x2, n_states)
        Tracked energies.
    D : ndarray, shape (n_x1, n_x2, n_states, 3)
        Tracked dipole vectors.
    cost_at_assignment : ndarray, shape (n_x1, n_x2)
        Total Hungarian assignment cost at each grid point.

    Returns
    -------
    DiagnosticReport
    """
    n_x1, n_x2, n_states = E.shape

    # Neighbor energy mismatch: avg |E[A,k] - E[B,k]| over neighbor pairs
    # Consider both x1 and x2 directions
    energy_diffs = []
    dipole_mismatches = []

    # x1 direction
    if n_x1 > 1:
        dE_x1 = np.abs(E[1:, :, :] - E[:-1, :, :])
        energy_diffs.append(dE_x1)

        # Dipole mismatch: 1 - |cos(angle)| between neighbors
        D_a = D[:-1, :, :, :]
        D_b = D[1:, :, :, :]
        norm_a = np.linalg.norm(D_a, axis=-1, keepdims=True)
        norm_b = np.linalg.norm(D_b, axis=-1, keepdims=True)
        eps = 1e-12
        cos_angle = np.sum(D_a * D_b, axis=-1) / (
            norm_a.squeeze(-1) * norm_b.squeeze(-1) + eps
        )
        dipole_mismatches.append(1.0 - np.abs(cos_angle))

    # x2 direction
    if n_x2 > 1:
        dE_x2 = np.abs(E[:, 1:, :] - E[:, :-1, :])
        energy_diffs.append(dE_x2)

        D_a = D[:, :-1, :, :]
        D_b = D[:, 1:, :, :]
        norm_a = np.linalg.norm(D_a, axis=-1, keepdims=True)
        norm_b = np.linalg.norm(D_b, axis=-1, keepdims=True)
        cos_angle = np.sum(D_a * D_b, axis=-1) / (
            norm_a.squeeze(-1) * norm_b.squeeze(-1) + eps
        )
        dipole_mismatches.append(1.0 - np.abs(cos_angle))

    if energy_diffs:
        neighbor_energy_mismatch = float(
            np.mean(np.concatenate([d.ravel() for d in energy_diffs]))
        )
    else:
        neighbor_energy_mismatch = 0.0

    if dipole_mismatches:
        neighbor_dipole_mismatch = float(
            np.mean(np.concatenate([d.ravel() for d in dipole_mismatches]))
        )
    else:
        neighbor_dipole_mismatch = 0.0

    # Per-state energy roughness: mean absolute second difference
    roughness = np.zeros(n_states)
    count = 0
    if n_x1 > 2:
        d2E_x1 = E[2:, :, :] - 2 * E[1:-1, :, :] + E[:-2, :, :]
        roughness += np.mean(np.abs(d2E_x1), axis=(0, 1))
        count += 1
    if n_x2 > 2:
        d2E_x2 = E[:, 2:, :] - 2 * E[:, 1:-1, :] + E[:, :-2, :]
        roughness += np.mean(np.abs(d2E_x2), axis=(0, 1))
        count += 1
    if count > 0:
        roughness /= count

    max_assignment_cost = float(np.max(cost_at_assignment))

    return DiagnosticReport(
        neighbor_energy_mismatch=neighbor_energy_mismatch,
        neighbor_dipole_mismatch=neighbor_dipole_mismatch,
        per_state_energy_roughness=roughness,
        max_assignment_cost=max_assignment_cost,
    )
