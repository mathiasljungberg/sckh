"""Tests for BFS module."""

import numpy as np

from python_scripts.surface_tracking.bfs import bfs_assign, grid_neighbors
from python_scripts.surface_tracking.cost import precompute_features
from python_scripts.surface_tracking.tracking import TrackingConfig


class TestGridNeighbors:
    def test_center(self):
        neighbors = grid_neighbors(1, 1, 3, 3)
        assert set(neighbors) == {(0, 1), (2, 1), (1, 0), (1, 2)}

    def test_corner(self):
        neighbors = grid_neighbors(0, 0, 3, 3)
        assert set(neighbors) == {(1, 0), (0, 1)}

    def test_edge(self):
        neighbors = grid_neighbors(0, 1, 3, 3)
        assert set(neighbors) == {(0, 0), (0, 2), (1, 1)}

    def test_bottom_right_corner(self):
        neighbors = grid_neighbors(2, 2, 3, 3)
        assert set(neighbors) == {(1, 2), (2, 1)}

    def test_1x1_grid(self):
        neighbors = grid_neighbors(0, 0, 1, 1)
        assert neighbors == []


class TestBfsAssign:
    def test_no_crossing_identity_perm(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1))
        features = precompute_features(E, D, config)
        E_t, D_t, perm, bfs_order, bfs_parent, cost = bfs_assign(
            E, D, features, config
        )
        # With no crossing, permutation should be identity everywhere
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(perm[i, j], [0, 1])

    def test_crossing_detected(self, crossing_2state_3x3):
        E, D, E_A, E_B, D_A, D_B = crossing_2state_3x3
        config = TrackingConfig(ref_point=(0, 0), sigma_E=1.0)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, bfs_order, bfs_parent, cost = bfs_assign(
            E, D, features, config
        )
        # At crossed points, states should be reordered to follow character.
        # (0,2) is the first crossed point reachable from (0,0).
        assert np.array_equal(perm[0, 2], [1, 0])
        # (2,2) is also crossed
        assert np.array_equal(perm[2, 2], [1, 0])

    def test_all_points_visited(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1))
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, _, _ = bfs_assign(E, D, features, config)
        assert len(bfs_order) == 9  # 3x3 grid

    def test_ref_point_default_center(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig()  # ref_point=None -> center
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, _, _ = bfs_assign(E, D, features, config)
        assert bfs_order[0] == (1, 1)

    def test_random_perms_recovered(self, random_permuted_5x5):
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(ref_point=(2, 2), sigma_E=0.1)
        features = precompute_features(E_perm, D_perm, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E_perm, D_perm, features, config)

        # Tracked energies should be smooth like E_true
        # Check that at each point, the set of energies matches
        for i in range(5):
            for j in range(5):
                np.testing.assert_allclose(
                    sorted(E_t[i, j]), sorted(E_true[i, j]), atol=1e-10
                )

        # Check tracked surfaces are smooth: neighbor differences should be small
        max_diff = 0
        for i in range(4):
            for j in range(5):
                diff = np.max(np.abs(E_t[i + 1, j] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        for i in range(5):
            for j in range(4):
                diff = np.max(np.abs(E_t[i, j + 1] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        # Smooth surfaces should have small neighbor differences
        assert max_diff < 0.1

    def test_bfs_parent_structure(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1))
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, bfs_parent, _ = bfs_assign(
            E, D, features, config
        )
        # Reference point has parent (-1, -1)
        np.testing.assert_array_equal(bfs_parent[1, 1], [-1, -1])
        # All other points should have valid parents
        for i in range(3):
            for j in range(3):
                if (i, j) != (1, 1):
                    assert bfs_parent[i, j, 0] >= 0
                    assert bfs_parent[i, j, 1] >= 0
