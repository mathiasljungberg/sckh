"""Tests for cost module."""

import numpy as np
import pytest

from python_scripts.surface_tracking.config import TrackingConfig
from python_scripts.surface_tracking.cost import (
    _find_degenerate_clusters,
    build_cost_matrix,
    build_cost_matrix_subspace,
    precompute_features,
)


class TestPrecomputeFeatures:
    def test_dipole_norms(self):
        E = np.array([[[1.0, 2.0]]])
        D = np.array([[[[3.0, 4.0, 0.0], [0.0, 0.0, 5.0]]]])
        features = precompute_features(E, D, TrackingConfig())
        np.testing.assert_allclose(features.dipole_norms[0, 0], [5.0, 5.0])

    def test_confidence_clipped(self):
        E = np.array([[[1.0, 2.0]]])
        D = np.array([[[[10.0, 0.0, 0.0], [0.01, 0.0, 0.0]]]])
        config = TrackingConfig(D_ref=0.1)
        features = precompute_features(E, D, config)
        assert features.confidence[0, 0, 0] == 1.0
        assert features.confidence[0, 0, 1] == pytest.approx(0.1)

    def test_osc_strengths(self):
        E = np.array([[[0.0, 2.0]]])
        D = np.array([[[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]]])
        features = precompute_features(E, D, TrackingConfig())
        # f = dE * |D|^2: state 0: 0*1=0, state 1: 2*1=2
        np.testing.assert_allclose(features.osc_strengths[0, 0], [0.0, 2.0])


class TestBuildCostMatrix:
    def test_identical_states_zero_diagonal(self):
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig()
        C = build_cost_matrix(E, D, E, D, config)
        # Diagonal should be zero (identical states)
        np.testing.assert_allclose(np.diag(C), 0.0, atol=1e-10)

    def test_swapped_states_offdiag_lower(self):
        E_A = np.array([1.0, 3.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        # B has states swapped
        E_B = np.array([3.0, 1.0])
        D_B = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        config = TrackingConfig()
        C = build_cost_matrix(E_A, D_A, E_B, D_B, config)
        # Off-diagonal should be lower than diagonal for swapped states
        assert C[0, 1] < C[0, 0]
        assert C[1, 0] < C[1, 1]

    def test_energy_gate_blocks(self):
        E_A = np.array([1.0, 10.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        E_B = np.array([1.1, 10.1])
        D_B = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(E_gate=0.5)
        C = build_cost_matrix(E_A, D_A, E_B, D_B, config)
        # (0,0) and (1,1): dE=0.1, within gate
        assert C[0, 0] < 1e10
        assert C[1, 1] < 1e10
        # (0,1): dE=9.1, blocked
        assert C[0, 1] >= 1e12
        assert C[1, 0] >= 1e12

    def test_sign_flip_invariance(self):
        E = np.array([1.0, 3.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        D_B_pos = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        D_B_neg = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])
        config = TrackingConfig()
        C_pos = build_cost_matrix(E, D_A, E, D_B_pos, config)
        C_neg = build_cost_matrix(E, D_A, E, D_B_neg, config)
        np.testing.assert_allclose(C_pos, C_neg, atol=1e-12)

    def test_zero_dipole_falls_back_to_energy(self):
        E_A = np.array([1.0, 3.0])
        D_A = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        E_B = np.array([1.0, 3.0])
        D_B = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        config = TrackingConfig(w_E=1.0, w_D=1.0)
        C = build_cost_matrix(E_A, D_A, E_B, D_B, config)
        # With zero dipoles, dipole cost is ~1 everywhere
        # Energy cost distinguishes: diagonal should be lower
        assert C[0, 0] < C[0, 1]
        assert C[1, 1] < C[1, 0]

    def test_oscillator_strength_term(self):
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        f_A = np.array([0.0, 2.0])
        f_B = np.array([0.0, 2.0])
        config = TrackingConfig(w_f=1.0, sigma_f=0.01)
        C = build_cost_matrix(E, D, E, D, config, f_A=f_A, f_B=f_B)
        # Diagonal should still be zero (identical)
        np.testing.assert_allclose(np.diag(C), 0.0, atol=1e-10)

    def test_confidence_weighting_zeros_out_weak_dipole(self):
        """Feature A: confidence weighting should zero out dipole cost for weak dipoles."""
        E = np.array([1.0, 3.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.001, 0.0, 0.0]])  # state 1: near-zero
        D_B = np.array([[1.0, 0.0, 0.0], [0.0, 0.001, 0.0]])  # state 1: near-zero, different dir
        conf_A = np.array([1.0, 0.01])  # low confidence for state 1
        conf_B = np.array([1.0, 0.01])

        # Without confidence: dipole cost for (1,1) is large because dirs differ
        config_off = TrackingConfig(w_E=0.0, w_D=1.0, use_confidence=False)
        C_off = build_cost_matrix(E, D_A, E, D_B, config_off)

        # With confidence: dipole cost for (1,1) should be much smaller
        config_on = TrackingConfig(w_E=0.0, w_D=1.0, use_confidence=True)
        C_on = build_cost_matrix(E, D_A, E, D_B, config_on, conf_A=conf_A, conf_B=conf_B)

        # Row 1 (weak dipole state) should have lower costs with confidence on
        assert C_on[1, 1] < C_off[1, 1]

    def test_confidence_weighting_no_effect_when_off(self):
        """Feature A: with use_confidence=False, conf_A/conf_B are ignored."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        conf = np.array([0.5, 0.5])
        config = TrackingConfig(use_confidence=False)
        C_with = build_cost_matrix(E, D, E, D, config, conf_A=conf, conf_B=conf)
        C_without = build_cost_matrix(E, D, E, D, config)
        np.testing.assert_allclose(C_with, C_without)

    def test_dipole_norm_cost(self):
        """Feature G: dipole-norm similarity cost discriminates by magnitude."""
        E = np.array([1.0, 1.0])  # same energy
        D_A = np.array([[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]])  # same dir, diff norm
        D_B = np.array([[0.9, 0.0, 0.0], [0.55, 0.0, 0.0]])
        config = TrackingConfig(w_E=0.0, w_D=0.0, w_norm=1.0, sigma_norm=0.1)
        C = build_cost_matrix(E, D_A, E, D_B, config)
        # (0,0) match: |1.0-0.9|/0.1 = 1.0; (0,1) match: |1.0-0.55|/0.1 = 4.5
        assert C[0, 0] < C[0, 1]
        # (1,1) match: |0.5-0.55|/0.1 = 0.5; (1,0) match: |0.5-0.9|/0.1 = 4.0
        assert C[1, 1] < C[1, 0]

    def test_dipole_norm_cost_off_by_default(self):
        """Feature G: w_norm=0 means no norm cost."""
        E = np.array([1.0, 1.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        D_B = np.array([[0.9, 0.0, 0.0], [0.55, 0.0, 0.0]])
        config = TrackingConfig(w_norm=0.0)
        C1 = build_cost_matrix(E, D_A, E, D_B, config)
        config2 = TrackingConfig(w_norm=0.0, sigma_norm=999.0)
        C2 = build_cost_matrix(E, D_A, E, D_B, config2)
        np.testing.assert_allclose(C1, C2)

    def test_sigma_E_override(self):
        """Feature E: sigma_E_override changes energy normalization."""
        E = np.array([1.0, 2.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(sigma_E=0.1)
        C_default = build_cost_matrix(E, D, E, D, config)
        C_override = build_cost_matrix(E, D, E, D, config, sigma_E_override=1.0)
        # With larger sigma_E, energy cost should be smaller
        # Off-diagonal: |1-2|/0.1=10 vs |1-2|/1.0=1
        assert C_override[0, 1] < C_default[0, 1]


class TestAdaptiveSigmaE:
    """Feature E: adaptive sigma_E tests."""

    def test_sigma_E_local_computed(self):
        E = np.array([[[1.0, 2.0, 5.0]]])  # gaps: 1.0, 3.0; median=2.0
        D = np.zeros((1, 1, 3, 3))
        D[0, 0] = np.eye(3)
        config = TrackingConfig(adaptive_sigma_E=True, sigma_E=1.0)
        features = precompute_features(E, D, config)
        assert features.sigma_E_local is not None
        # median gap = 2.0, clipped to [0.1, 10.0]
        np.testing.assert_allclose(features.sigma_E_local[0, 0], 2.0)

    def test_sigma_E_local_clipped(self):
        # Very small gaps: sigma_E=1.0, gap=0.001 -> clipped to sigma_E/10=0.1
        E = np.array([[[1.0, 1.001]]])
        D = np.zeros((1, 1, 2, 3))
        D[0, 0, 0] = [1, 0, 0]
        D[0, 0, 1] = [0, 1, 0]
        config = TrackingConfig(adaptive_sigma_E=True, sigma_E=1.0)
        features = precompute_features(E, D, config)
        np.testing.assert_allclose(features.sigma_E_local[0, 0], 0.1)

    def test_sigma_E_local_none_when_off(self):
        E = np.array([[[1.0, 2.0]]])
        D = np.zeros((1, 1, 2, 3))
        config = TrackingConfig(adaptive_sigma_E=False)
        features = precompute_features(E, D, config)
        assert features.sigma_E_local is None

    def test_adaptive_sigma_E_through_pipeline(self, smooth_2state_3x3):
        from python_scripts.surface_tracking import track_surfaces

        E, D = smooth_2state_3x3
        config = TrackingConfig(adaptive_sigma_E=True)
        result = track_surfaces(E, D, config)
        # Should produce smooth surfaces like default
        np.testing.assert_allclose(result.E, E)
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(result.perm[i, j], [0, 1])


class TestSubspaceTracking:
    """Feature H: subspace tracking for degenerate clusters."""

    def test_find_degenerate_clusters_single(self):
        E = np.array([1.0, 3.0, 5.0])
        clusters = _find_degenerate_clusters(E, threshold=0.1)
        assert len(clusters) == 3
        for c in clusters:
            assert len(c) == 1

    def test_find_degenerate_clusters_pair(self):
        E = np.array([1.0, 1.05, 3.0])
        clusters = _find_degenerate_clusters(E, threshold=0.1)
        assert len(clusters) == 2
        # First cluster should have 2 states
        sizes = sorted([len(c) for c in clusters])
        assert sizes == [1, 2]

    def test_find_degenerate_clusters_all_degenerate(self):
        E = np.array([1.0, 1.01, 1.02])
        clusters = _find_degenerate_clusters(E, threshold=0.1)
        assert len(clusters) == 1
        assert len(clusters[0]) == 3

    def test_subspace_cost_matrix_non_degenerate(self):
        """Non-degenerate states should use standard cost matrix."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(use_subspace=True, degen_threshold=0.01)
        C_sub = build_cost_matrix_subspace(E, D, E, D, config)
        C_std = build_cost_matrix(E, D, E, D, config)
        np.testing.assert_allclose(C_sub, C_std)

    def test_subspace_cost_matrix_degenerate_pair(self):
        """Degenerate states should use subspace cost."""
        E_A = np.array([2.0, 2.005, 5.0])
        E_B = np.array([2.001, 2.004, 5.1])
        # Rotated dipoles within the degenerate subspace
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        D_B = np.array([[0.7, 0.7, 0.0], [-0.7, 0.7, 0.0], [0.0, 0.0, 1.0]])
        config = TrackingConfig(
            use_subspace=True, degen_threshold=0.1, sigma_E=1.0
        )
        C_sub = build_cost_matrix_subspace(E_A, D_A, E_B, D_B, config)
        C_std = build_cost_matrix(E_A, D_A, E_B, D_B, config)
        # Subspace cost should differ from standard for the degenerate cluster
        # but be the same for the non-degenerate state
        assert C_sub[2, 2] == pytest.approx(C_std[2, 2], abs=1e-10)
        # The degenerate block should have been modified
        assert not np.allclose(C_sub[:2, :2], C_std[:2, :2])

    def test_subspace_off_by_default(self):
        """use_subspace=False means no subspace cost."""
        E = np.array([2.0, 2.005])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(use_subspace=False, degen_threshold=0.1)
        C = build_cost_matrix(E, D, E, D, config)
        # Standard cost matrix, no subspace modification
        np.testing.assert_allclose(np.diag(C), 0.0, atol=1e-10)


class TestEnergyDerivativeCost:
    """Feature I: energy-derivative continuity cost."""

    def test_gradient_cost_prefers_predicted_energy(self):
        """Gradient cost should prefer the state matching the energy trend."""
        E_A = np.array([1.0, 3.0])
        E_B = np.array([1.1, 2.9])  # B energies
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        grad_A = np.array([0.1, -0.1])  # state 0 rising, state 1 falling
        config = TrackingConfig(w_E=0.0, w_D=0.0, w_dE=1.0, sigma_dE=0.1)
        C = build_cost_matrix(E_A, D, E_B, D, config, grad_A=grad_A)
        # State 0: predicted E=1.1, E_B[0]=1.1 (perfect match), E_B[1]=2.9 (bad)
        assert C[0, 0] < C[0, 1]
        # State 1: predicted E=2.9, E_B[1]=2.9 (perfect match), E_B[0]=1.1 (bad)
        assert C[1, 1] < C[1, 0]

    def test_gradient_cost_off_by_default(self):
        """w_dE=0 means no gradient cost."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        grad = np.array([0.1, -0.1])
        config = TrackingConfig(w_dE=0.0)
        C1 = build_cost_matrix(E, D, E, D, config)
        C2 = build_cost_matrix(E, D, E, D, config, grad_A=grad)
        np.testing.assert_allclose(C1, C2)

    def test_gradient_cost_no_effect_without_grad(self):
        """Without grad_A, w_dE has no effect."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(w_dE=1.0)
        C = build_cost_matrix(E, D, E, D, config, grad_A=None)
        # Should be same as w_dE=0
        config_off = TrackingConfig(w_dE=0.0)
        C_off = build_cost_matrix(E, D, E, D, config_off)
        np.testing.assert_allclose(C, C_off)


class TestDipoleDerivativeCost:
    """Feature J: dipole-derivative continuity cost."""

    def test_dipole_gradient_prefers_predicted_direction(self):
        """Dipole gradient cost should prefer the state matching the dipole trend."""
        E = np.array([1.0, 3.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        # State 0 dipole trending toward [1.2, 0, 0], state 1 toward [0, 1.2, 0]
        grad_D_A = np.array([[0.2, 0.0, 0.0], [0.0, 0.2, 0.0]])
        # D_B[0] matches state 0's predicted dipole, D_B[1] matches state 1's
        D_B = np.array([[1.15, 0.05, 0.0], [0.05, 1.15, 0.0]])
        config = TrackingConfig(w_E=0.0, w_D=0.0, w_dD=1.0, sigma_dD=0.1)
        C = build_cost_matrix(E, D_A, E, D_B, config, grad_D_A=grad_D_A)
        # (0,0) should be low cost: predicted [1.2,0,0] close to D_B[0]=[1.15,0.05,0]
        # (0,1) should be high cost: predicted [1.2,0,0] far from D_B[1]=[0.05,1.15,0]
        assert C[0, 0] < C[0, 1]
        assert C[1, 1] < C[1, 0]

    def test_dipole_gradient_sign_invariant(self):
        """Dipole gradient cost should be phase-invariant (sign of D_B)."""
        E = np.array([1.0, 3.0])
        D_A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        grad_D_A = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
        D_B_pos = np.array([[1.1, 0.0, 0.0], [0.0, 1.1, 0.0]])
        D_B_neg = -D_B_pos  # sign-flipped
        config = TrackingConfig(w_E=0.0, w_D=0.0, w_dD=1.0, sigma_dD=0.1)
        C_pos = build_cost_matrix(E, D_A, E, D_B_pos, config, grad_D_A=grad_D_A)
        C_neg = build_cost_matrix(E, D_A, E, D_B_neg, config, grad_D_A=grad_D_A)
        np.testing.assert_allclose(C_pos, C_neg, atol=1e-12)

    def test_dipole_gradient_off_by_default(self):
        """w_dD=0 means no dipole gradient cost."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        grad_D = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
        config = TrackingConfig(w_dD=0.0)
        C1 = build_cost_matrix(E, D, E, D, config)
        C2 = build_cost_matrix(E, D, E, D, config, grad_D_A=grad_D)
        np.testing.assert_allclose(C1, C2)

    def test_dipole_gradient_no_effect_without_grad(self):
        """Without grad_D_A, w_dD has no effect."""
        E = np.array([1.0, 3.0])
        D = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        config = TrackingConfig(w_dD=1.0)
        C = build_cost_matrix(E, D, E, D, config, grad_D_A=None)
        config_off = TrackingConfig(w_dD=0.0)
        C_off = build_cost_matrix(E, D, E, D, config_off)
        np.testing.assert_allclose(C, C_off)
