"""Tests for 1D sweep module."""

import numpy as np

from python_scripts.surface_tracking import TrackingConfig, track_surfaces
from python_scripts.surface_tracking.bfs import bfs_assign
from python_scripts.surface_tracking.config import TrackingConfig
from python_scripts.surface_tracking.cost import precompute_features
from python_scripts.surface_tracking.sweeps import sweep_cols, sweep_rows


class TestSweepRows:
    def test_fixes_isolated_hole_in_row(self):
        """Row sweep should fix a single wrong point in a row."""
        n_x1, n_x2, n_states = 5, 1, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            E[i, 0, 0] = 1.0 + 0.01 * i
            E[i, 0, 1] = 3.0 + 0.01 * i
            D[i, 0, 0] = [1.0, 0.0, 0.0]
            D[i, 0, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point=(2, 0))
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)

        # Corrupt one point in the middle
        E_t[2, 0, 0], E_t[2, 0, 1] = E_t[2, 0, 1].copy(), E_t[2, 0, 0].copy()
        D_t[2, 0, 0], D_t[2, 0, 1] = D_t[2, 0, 1].copy(), D_t[2, 0, 0].copy()
        perm[2, 0] = [1, 0]

        changes = sweep_rows(E, D, E_t, D_t, perm, features, config)
        assert changes > 0
        np.testing.assert_array_equal(perm[2, 0], [0, 1])


class TestSweepCols:
    def test_fixes_isolated_hole_in_col(self):
        """Column sweep should fix a single wrong point in a column."""
        n_x1, n_x2, n_states = 1, 5, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for j in range(n_x2):
            E[0, j, 0] = 1.0 + 0.01 * j
            E[0, j, 1] = 3.0 + 0.01 * j
            D[0, j, 0] = [1.0, 0.0, 0.0]
            D[0, j, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point=(0, 2))
        features = precompute_features(E, D, config)
        E_t, D_t, perm, _, _, _ = bfs_assign(E, D, features, config)

        # Corrupt one point in the middle
        E_t[0, 2, 0], E_t[0, 2, 1] = E_t[0, 2, 1].copy(), E_t[0, 2, 0].copy()
        D_t[0, 2, 0], D_t[0, 2, 1] = D_t[0, 2, 1].copy(), D_t[0, 2, 0].copy()
        perm[0, 2] = [1, 0]

        changes = sweep_cols(E, D, E_t, D_t, perm, features, config)
        assert changes > 0
        np.testing.assert_array_equal(perm[0, 2], [0, 1])


class TestSweepsOffByDefault:
    def test_n_sweep_iter_zero(self, smooth_2state_3x3):
        """n_sweep_iter=0 means no sweeps in the pipeline."""
        E, D = smooth_2state_3x3
        config = TrackingConfig(ref_point=(1, 1), n_sweep_iter=0)
        result = track_surfaces(E, D, config)
        # Should produce same result as without sweeps
        np.testing.assert_allclose(result.E, E)


class TestSweepPipeline:
    def test_sweep_through_pipeline(self, random_permuted_5x5):
        """Feature N through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2),
            sigma_E=0.1,
            n_sweep_iter=2,
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1
