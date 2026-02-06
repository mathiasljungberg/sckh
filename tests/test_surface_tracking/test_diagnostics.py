"""Tests for diagnostics module."""

import numpy as np

from python_scripts.surface_tracking.diagnostics import compute_diagnostics


class TestComputeDiagnostics:
    def test_smooth_surfaces_low_mismatch(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        cost = np.zeros((3, 3))
        report = compute_diagnostics(E, D, cost)
        # Smooth surfaces should have low mismatch
        assert report.neighbor_energy_mismatch < 0.05
        assert report.neighbor_dipole_mismatch < 0.05
        assert report.max_assignment_cost == 0.0

    def test_rough_surfaces_high_mismatch(self):
        rng = np.random.default_rng(42)
        E = rng.standard_normal((5, 5, 2))
        D = rng.standard_normal((5, 5, 2, 3))
        cost = np.zeros((5, 5))
        report = compute_diagnostics(E, D, cost)
        # Random surfaces should have high mismatch
        assert report.neighbor_energy_mismatch > 0.1
        assert report.neighbor_dipole_mismatch > 0.1

    def test_output_shapes(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        cost = np.zeros((3, 3))
        report = compute_diagnostics(E, D, cost)
        assert report.per_state_energy_roughness.shape == (2,)
        assert isinstance(report.neighbor_energy_mismatch, float)
        assert isinstance(report.neighbor_dipole_mismatch, float)
        assert isinstance(report.max_assignment_cost, float)

    def test_1x1_grid(self):
        E = np.array([[[1.0, 2.0]]])
        D = np.array([[[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]])
        cost = np.array([[0.0]])
        report = compute_diagnostics(E, D, cost)
        assert report.neighbor_energy_mismatch == 0.0
        assert report.neighbor_dipole_mismatch == 0.0

    def test_max_assignment_cost(self):
        E = np.zeros((3, 3, 2))
        D = np.zeros((3, 3, 2, 3))
        D[:, :, 0, 0] = 1.0
        D[:, :, 1, 1] = 1.0
        cost = np.array([[0.1, 0.5, 0.2], [0.3, 0.0, 0.4], [0.1, 0.2, 0.9]])
        report = compute_diagnostics(E, D, cost)
        assert report.max_assignment_cost == 0.9
