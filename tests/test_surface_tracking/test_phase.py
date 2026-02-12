"""Tests for phase module."""

import numpy as np

from python_scripts.surface_tracking.phase import fix_phases_bfs, relax_phases


class TestFixPhasesBfs:
    def test_no_flip_needed(self):
        """Consistent signs should remain unchanged."""
        D = np.zeros((3, 3, 2, 3))
        D[:, :, 0] = [1.0, 0.0, 0.0]
        D[:, :, 1] = [0.0, 1.0, 0.0]
        D_orig = D.copy()

        bfs_order = [(1, 1), (0, 1), (2, 1), (1, 0), (1, 2),
                     (0, 0), (0, 2), (2, 0), (2, 2)]
        bfs_parent = np.full((3, 3, 2), -1, dtype=int)
        bfs_parent[0, 1] = [1, 1]
        bfs_parent[2, 1] = [1, 1]
        bfs_parent[1, 0] = [1, 1]
        bfs_parent[1, 2] = [1, 1]
        bfs_parent[0, 0] = [0, 1]
        bfs_parent[0, 2] = [0, 1]
        bfs_parent[2, 0] = [2, 1]
        bfs_parent[2, 2] = [2, 1]

        fix_phases_bfs(D, bfs_order, bfs_parent)
        np.testing.assert_array_equal(D, D_orig)

    def test_single_flip_corrected(self):
        """A flipped dipole should be corrected."""
        D = np.zeros((2, 1, 1, 3))
        D[0, 0, 0] = [1.0, 0.0, 0.0]
        D[1, 0, 0] = [-1.0, 0.0, 0.0]  # flipped

        bfs_order = [(0, 0), (1, 0)]
        bfs_parent = np.full((2, 1, 2), -1, dtype=int)
        bfs_parent[1, 0] = [0, 0]

        fix_phases_bfs(D, bfs_order, bfs_parent)
        np.testing.assert_allclose(D[1, 0, 0], [1.0, 0.0, 0.0])

    def test_near_zero_norm_skipped(self):
        """Near-zero dipoles should not be flipped."""
        D = np.zeros((2, 1, 1, 3))
        D[0, 0, 0] = [1.0, 0.0, 0.0]
        D[1, 0, 0] = [-1e-15, 0.0, 0.0]  # tiny, should be skipped
        D_orig_1 = D[1, 0, 0].copy()

        bfs_order = [(0, 0), (1, 0)]
        bfs_parent = np.full((2, 1, 2), -1, dtype=int)
        bfs_parent[1, 0] = [0, 0]

        fix_phases_bfs(D, bfs_order, bfs_parent)
        np.testing.assert_array_equal(D[1, 0, 0], D_orig_1)

    def test_magnitude_preserved(self):
        """Phase fix should preserve dipole magnitude."""
        D = np.zeros((2, 1, 1, 3))
        D[0, 0, 0] = [2.5, 1.0, 0.5]
        D[1, 0, 0] = [-2.5, -1.0, -0.5]  # flipped

        bfs_order = [(0, 0), (1, 0)]
        bfs_parent = np.full((2, 1, 2), -1, dtype=int)
        bfs_parent[1, 0] = [0, 0]

        norm_before = np.linalg.norm(D[1, 0, 0])
        fix_phases_bfs(D, bfs_order, bfs_parent)
        norm_after = np.linalg.norm(D[1, 0, 0])
        assert abs(norm_before - norm_after) < 1e-12

    def test_multiple_states_independent(self):
        """Each state's phase should be fixed independently."""
        D = np.zeros((2, 1, 2, 3))
        D[0, 0, 0] = [1.0, 0.0, 0.0]
        D[0, 0, 1] = [0.0, 1.0, 0.0]
        D[1, 0, 0] = [-1.0, 0.0, 0.0]  # flipped
        D[1, 0, 1] = [0.0, 1.0, 0.0]   # not flipped

        bfs_order = [(0, 0), (1, 0)]
        bfs_parent = np.full((2, 1, 2), -1, dtype=int)
        bfs_parent[1, 0] = [0, 0]

        fix_phases_bfs(D, bfs_order, bfs_parent)
        np.testing.assert_allclose(D[1, 0, 0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(D[1, 0, 1], [0.0, 1.0, 0.0])


class TestRelaxPhases:
    """Feature L: iterative neighbor-consensus phase relaxation."""

    def test_relax_fixes_isolated_flip(self):
        """A single flipped point surrounded by consistent neighbors should be fixed."""
        D = np.zeros((3, 3, 1, 3))
        D[:, :, 0] = [1.0, 0.0, 0.0]
        # Flip the center point
        D[1, 1, 0] = [-1.0, 0.0, 0.0]

        total_flips = relax_phases(D, n_iter=5)
        assert total_flips > 0
        # Center should now agree with neighbors
        np.testing.assert_allclose(D[1, 1, 0], [1.0, 0.0, 0.0])

    def test_relax_converges_immediately_if_consistent(self):
        """Already-consistent grid should produce zero flips."""
        D = np.zeros((3, 3, 2, 3))
        D[:, :, 0] = [1.0, 0.0, 0.0]
        D[:, :, 1] = [0.0, 1.0, 0.0]

        total_flips = relax_phases(D, n_iter=5)
        assert total_flips == 0

    def test_relax_off_when_zero_iters(self):
        """n_iter=0 means no relaxation, even if flips needed."""
        D = np.zeros((3, 3, 1, 3))
        D[:, :, 0] = [1.0, 0.0, 0.0]
        D[1, 1, 0] = [-1.0, 0.0, 0.0]  # flipped
        D_orig = D.copy()

        total_flips = relax_phases(D, n_iter=0)
        assert total_flips == 0
        np.testing.assert_array_equal(D, D_orig)

    def test_relax_skips_weak_dipoles(self):
        """Near-zero dipoles should not be flipped."""
        D = np.zeros((3, 1, 1, 3))
        D[0, 0, 0] = [1.0, 0.0, 0.0]
        D[1, 0, 0] = [-1e-15, 0.0, 0.0]  # tiny
        D[2, 0, 0] = [1.0, 0.0, 0.0]
        D_mid_orig = D[1, 0, 0].copy()

        relax_phases(D, n_iter=5)
        np.testing.assert_array_equal(D[1, 0, 0], D_mid_orig)

    def test_relax_multiple_states_independent(self):
        """Each state should be relaxed independently."""
        D = np.zeros((3, 3, 2, 3))
        D[:, :, 0] = [1.0, 0.0, 0.0]
        D[:, :, 1] = [0.0, 1.0, 0.0]
        # Flip center for state 0 only
        D[1, 1, 0] = [-1.0, 0.0, 0.0]

        relax_phases(D, n_iter=5)
        # State 0 center should be fixed
        np.testing.assert_allclose(D[1, 1, 0], [1.0, 0.0, 0.0])
        # State 1 should be unchanged
        np.testing.assert_allclose(D[1, 1, 1], [0.0, 1.0, 0.0])
