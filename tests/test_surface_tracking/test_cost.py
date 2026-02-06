"""Tests for cost module."""

import numpy as np
import pytest

from python_scripts.surface_tracking.cost import (
    build_cost_matrix,
    precompute_features,
)
from python_scripts.surface_tracking.config import TrackingConfig


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
