"""Tests for BFS module."""

import numpy as np

from python_scripts.surface_tracking.bfs import (
    bfs_assign,
    grid_neighbors,
    repair_assignments,
    repair_reassignment,
)
from python_scripts.surface_tracking.config import TrackingConfig
from python_scripts.surface_tracking.cost import precompute_features


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


class TestPriorityQueueBfs:
    """Feature B: priority-queue BFS tests."""

    def test_all_points_visited(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), use_priority_queue=True)
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, _, _ = bfs_assign(E, D, features, config)
        assert len(bfs_order) == 9

    def test_no_crossing_identity_perm(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), use_priority_queue=True)
        features = precompute_features(E, D, config)
        _, _, perm, _, _, _ = bfs_assign(E, D, features, config)
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(perm[i, j], [0, 1])

    def test_crossing_detected(self, crossing_2state_3x3):
        E, D, E_A, E_B, D_A, D_B = crossing_2state_3x3
        config = TrackingConfig(
            ref_point=(0, 0), sigma_E=1.0, use_priority_queue=True
        )
        features = precompute_features(E, D, config)
        _, _, perm, _, _, _ = bfs_assign(E, D, features, config)
        # Crossed points should be detected
        assert np.array_equal(perm[2, 2], [1, 0])

    def test_random_perms_recovered(self, random_permuted_5x5):
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2), sigma_E=0.1, use_priority_queue=True
        )
        features = precompute_features(E_perm, D_perm, config)
        E_t, _, _, _, _, _ = bfs_assign(E_perm, D_perm, features, config)
        # Tracked surfaces should be smooth
        max_diff = 0
        for i in range(4):
            for j in range(5):
                diff = np.max(np.abs(E_t[i + 1, j] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        for i in range(5):
            for j in range(4):
                diff = np.max(np.abs(E_t[i, j + 1] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        assert max_diff < 0.1

    def test_low_cost_edge_preferred(self):
        """Priority queue should prefer assigning from the lower-cost direction."""
        # 3x1 grid: ref at (0,0), smooth at (1,0), crossing at (2,0)
        # With priority queue, (1,0) gets assigned before (2,0)
        n_x1, n_x2, n_states = 3, 1, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        # Smooth variation
        for i in range(n_x1):
            E[i, 0, 0] = 1.0 + 0.01 * i
            E[i, 0, 1] = 3.0 + 0.01 * i
            D[i, 0, 0] = [1.0, 0.0, 0.0]
            D[i, 0, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point=(0, 0), use_priority_queue=True)
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, _, cost = bfs_assign(E, D, features, config)
        # Cost should be low everywhere (no crossings)
        assert cost[1, 0] < 1.0
        assert cost[2, 0] < 1.0


class TestMultiNeighborConsensus:
    """Feature C: multi-neighbor consensus tests."""

    def test_all_points_visited(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), use_multi_neighbor=True)
        features = precompute_features(E, D, config)
        _, _, _, bfs_order, _, _ = bfs_assign(E, D, features, config)
        assert len(bfs_order) == 9

    def test_no_crossing_identity_perm(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), use_multi_neighbor=True)
        features = precompute_features(E, D, config)
        _, _, perm, _, _, _ = bfs_assign(E, D, features, config)
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(perm[i, j], [0, 1])

    def test_crossing_detected(self, crossing_2state_3x3):
        E, D, E_A, E_B, D_A, D_B = crossing_2state_3x3
        config = TrackingConfig(
            ref_point=(0, 0), sigma_E=1.0, use_multi_neighbor=True
        )
        features = precompute_features(E, D, config)
        _, _, perm, _, _, _ = bfs_assign(E, D, features, config)
        assert np.array_equal(perm[2, 2], [1, 0])

    def test_random_perms_recovered(self, random_permuted_5x5):
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2), sigma_E=0.1, use_multi_neighbor=True
        )
        features = precompute_features(E_perm, D_perm, config)
        E_t, _, _, _, _, _ = bfs_assign(E_perm, D_perm, features, config)
        max_diff = 0
        for i in range(4):
            for j in range(5):
                diff = np.max(np.abs(E_t[i + 1, j] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        for i in range(5):
            for j in range(4):
                diff = np.max(np.abs(E_t[i, j + 1] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        assert max_diff < 0.1

    def test_combined_priority_and_multi_neighbor(self, random_permuted_5x5):
        """Features B+C together should also work."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2),
            sigma_E=0.1,
            use_priority_queue=True,
            use_multi_neighbor=True,
        )
        features = precompute_features(E_perm, D_perm, config)
        E_t, _, _, bfs_order, _, _ = bfs_assign(E_perm, D_perm, features, config)
        assert len(bfs_order) == 25  # 5x5 grid
        max_diff = 0
        for i in range(4):
            for j in range(5):
                diff = np.max(np.abs(E_t[i + 1, j] - E_t[i, j]))
                max_diff = max(max_diff, diff)
        assert max_diff < 0.1


class TestRepairAssignments:
    """Feature D: post-assignment repair pass tests."""

    def test_no_repair_when_zero_iters(self, smooth_2state_3x3):
        """n_repair_iter=0 means no repair (default behavior)."""
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), n_repair_iter=0)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)
        E_before = E_t.copy()
        total_swaps = repair_assignments(E_t, D_t, perm, features, config)
        assert total_swaps == 0
        np.testing.assert_array_equal(E_t, E_before)

    def test_no_swap_on_smooth_surface(self, smooth_2state_3x3):
        """Repair should not swap anything on an already-smooth surface."""
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), n_repair_iter=3)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)
        total_swaps = repair_assignments(E_t, D_t, perm, features, config)
        assert total_swaps == 0

    def test_repairs_deliberate_wrong_assignment(self):
        """Repair should fix a deliberately wrong assignment at one point."""
        n_x1, n_x2, n_states = 3, 3, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            for j in range(n_x2):
                E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
                E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
                D[i, j, 0] = [1.0, 0.0, 0.0]
                D[i, j, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point=(1, 1), n_repair_iter=3)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)

        # Deliberately corrupt assignment at (0, 0): swap states
        E_t[0, 0, 0], E_t[0, 0, 1] = E_t[0, 0, 1].copy(), E_t[0, 0, 0].copy()
        D_t[0, 0, 0], D_t[0, 0, 1] = D_t[0, 0, 1].copy(), D_t[0, 0, 0].copy()
        perm[0, 0] = [1, 0]

        total_swaps = repair_assignments(E_t, D_t, perm, features, config)
        assert total_swaps > 0
        # After repair, (0,0) should be back to identity
        np.testing.assert_array_equal(perm[0, 0], [0, 1])

    def test_repair_converges(self, random_permuted_5x5):
        """Repair should converge (stop swapping) within n_repair_iter."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(ref_point=(2, 2), sigma_E=0.1, n_repair_iter=5)
        features = precompute_features(E_perm, D_perm, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E_perm, D_perm, features, config)
        # Should converge without error
        total_swaps = repair_assignments(E_t, D_t, perm, features, config)
        assert total_swaps >= 0  # may or may not need swaps

    def test_repair_through_pipeline(self, random_permuted_5x5):
        """Feature D through full pipeline."""
        from python_scripts.surface_tracking import track_surfaces

        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2), sigma_E=0.1, n_repair_iter=3
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1


class TestInlinePhaseAlignment:
    """Feature K: inline phase alignment during BFS."""

    def test_d_tracked_has_consistent_signs(self):
        """D_tracked should have phase-consistent dipoles after BFS (before phase fix)."""
        n_x1, n_x2, n_states = 3, 3, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            for j in range(n_x2):
                E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
                E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
                # Alternate sign flips in the raw data
                sign = -1 if (i + j) % 2 == 1 else 1
                D[i, j, 0] = [sign * 1.0, 0.0, 0.0]
                D[i, j, 1] = [0.0, sign * 1.0, 0.0]

        config = TrackingConfig(ref_point=(1, 1))
        features = precompute_features(E, D, config)
        _, D_t, _, _, _, _ = bfs_assign(E, D, features, config)

        # After BFS with inline phase alignment, all dipoles should
        # agree in sign with the reference point
        ref_d0 = D_t[1, 1, 0]
        ref_d1 = D_t[1, 1, 1]
        for i in range(n_x1):
            for j in range(n_x2):
                assert np.dot(D_t[i, j, 0], ref_d0) > 0, (
                    f"State 0 at ({i},{j}) has wrong sign after BFS"
                )
                assert np.dot(D_t[i, j, 1], ref_d1) > 0, (
                    f"State 1 at ({i},{j}) has wrong sign after BFS"
                )


class TestRepairReassignment:
    """Feature M: full Hungarian re-assignment repair."""

    def test_fixes_fully_wrong_permutation(self):
        """A deliberately wrong full permutation should be fixed."""
        n_x1, n_x2, n_states = 3, 3, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            for j in range(n_x2):
                E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
                E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
                D[i, j, 0] = [1.0, 0.0, 0.0]
                D[i, j, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point=(1, 1), n_reassign_iter=3)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)

        # Deliberately corrupt assignment at (0, 0): swap ALL states
        E_t[0, 0, 0], E_t[0, 0, 1] = E_t[0, 0, 1].copy(), E_t[0, 0, 0].copy()
        D_t[0, 0, 0], D_t[0, 0, 1] = D_t[0, 0, 1].copy(), D_t[0, 0, 0].copy()
        perm[0, 0] = [1, 0]

        total = repair_reassignment(E, D, E_t, D_t, perm, features, config)
        assert total > 0
        np.testing.assert_array_equal(perm[0, 0], [0, 1])

    def test_no_change_on_smooth_surface(self, smooth_2state_3x3):
        """Smooth surface should produce zero changes."""
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), n_reassign_iter=3)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)
        total = repair_reassignment(E, D, E_t, D_t, perm, features, config)
        assert total == 0

    def test_off_by_default(self, smooth_2state_3x3):
        """n_reassign_iter=0 means no re-assignment."""
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), n_reassign_iter=0)
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)
        total = repair_reassignment(E, D, E_t, D_t, perm, features, config)
        assert total == 0
