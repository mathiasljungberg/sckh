"""Tests for the full tracking pipeline."""

import numpy as np
import pytest

from python_scripts.surface_tracking import TrackingConfig, TrackingResult, track_surfaces


class TestTrackSurfaces:
    def test_no_crossing(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        result = track_surfaces(E, D)
        assert isinstance(result, TrackingResult)
        # No crossing: tracked should match input
        np.testing.assert_allclose(result.E, E)
        np.testing.assert_allclose(result.D, D)
        # Identity permutation
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(result.perm[i, j], [0, 1])

    def test_crossing_tracked(self, crossing_2state_3x3):
        E, D, E_A, E_B, D_A, D_B = crossing_2state_3x3
        config = TrackingConfig(ref_point=(0, 0), sigma_E=1.0)
        result = track_surfaces(E, D, config)
        # At crossed points, states should be re-ordered to follow character
        np.testing.assert_array_equal(result.perm[0, 2], [1, 0])
        np.testing.assert_array_equal(result.perm[2, 2], [1, 0])

    def test_random_permutations_recovered(self, random_permuted_5x5):
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(ref_point=(2, 2), sigma_E=0.1)
        result = track_surfaces(E_perm, D_perm, config)

        # Tracked energies should be smooth
        for i in range(4):
            for j in range(5):
                diff = np.max(np.abs(result.E[i + 1, j] - result.E[i, j]))
                assert diff < 0.1
        for i in range(5):
            for j in range(4):
                diff = np.max(np.abs(result.E[i, j + 1] - result.E[i, j]))
                assert diff < 0.1

    def test_sign_flips_fixed(self, sign_flipped_3x3):
        E, D_flipped, D_true, flips = sign_flipped_3x3
        config = TrackingConfig(ref_point=(0, 0))
        result = track_surfaces(E, D_flipped, config)

        # After phase fixing, all dipoles should point in the same direction
        # as the reference point's dipoles (up to numerical precision)
        ref_signs = np.sign(result.D[0, 0, :, :])
        for i in range(3):
            for j in range(3):
                for k in range(2):
                    # Check alignment: dot product should be positive
                    dot = np.dot(result.D[0, 0, k], result.D[i, j, k])
                    assert dot > 0, (
                        f"Dipole at ({i},{j}), state {k} has wrong sign"
                    )

    def test_invalid_E_shape(self):
        with pytest.raises(ValueError, match="E must be 3D"):
            track_surfaces(np.zeros((3, 3)), np.zeros((3, 3, 2, 3)))

    def test_invalid_D_shape(self):
        with pytest.raises(ValueError, match="D must be 4D"):
            track_surfaces(np.zeros((3, 3, 2)), np.zeros((3, 3, 2)))

    def test_incompatible_shapes(self):
        with pytest.raises(ValueError, match="incompatible"):
            track_surfaces(np.zeros((3, 3, 2)), np.zeros((3, 3, 3, 3)))

    def test_output_shapes(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        result = track_surfaces(E, D)
        assert result.E.shape == (3, 3, 2)
        assert result.D.shape == (3, 3, 2, 3)
        assert result.perm.shape == (3, 3, 2)
        assert result.cost_at_assignment.shape == (3, 3)
        assert result.bfs_parent.shape == (3, 3, 2)
        assert result.diagnostics is not None

    def test_default_config(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        # Should work with no config
        result = track_surfaces(E, D)
        assert isinstance(result, TrackingResult)

    def test_diagnostics_populated(self, smooth_2state_3x3):
        E, D = smooth_2state_3x3
        result = track_surfaces(E, D)
        diag = result.diagnostics
        assert diag.neighbor_energy_mismatch >= 0
        assert diag.neighbor_dipole_mismatch >= 0
        assert diag.per_state_energy_roughness.shape == (2,)

    def test_confidence_weighting_pipeline(self, crossing_2state_3x3):
        """Feature A through full pipeline."""
        E, D, E_A, E_B, D_A, D_B = crossing_2state_3x3
        config = TrackingConfig(
            ref_point=(0, 0), sigma_E=1.0, use_confidence=True
        )
        result = track_surfaces(E, D, config)
        assert np.array_equal(result.perm[2, 2], [1, 0])

    def test_priority_queue_pipeline(self, random_permuted_5x5):
        """Feature B through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2), sigma_E=0.1, use_priority_queue=True
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1

    def test_multi_neighbor_pipeline(self, random_permuted_5x5):
        """Feature C through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2), sigma_E=0.1, use_multi_neighbor=True
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1

    def test_all_phase1_features_combined(self, random_permuted_5x5):
        """Features A+B+C combined through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2),
            sigma_E=0.1,
            use_confidence=True,
            use_priority_queue=True,
            use_multi_neighbor=True,
        )
        result = track_surfaces(E_perm, D_perm, config)
        assert result.E.shape == (5, 5, 3)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1

    def test_smart_ref_point_selection(self):
        """Feature F: ref_point='auto' picks smooth region."""
        n_x1, n_x2, n_states = 5, 5, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            for j in range(n_x2):
                E[i, j, 0] = 1.0 + 0.01 * i + 0.02 * j
                E[i, j, 1] = 3.0 + 0.01 * i + 0.02 * j
                D[i, j, 0] = [1.0, 0.0, 0.0]
                D[i, j, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(ref_point="auto")
        result = track_surfaces(E, D, config)
        # Should complete and produce smooth output
        np.testing.assert_allclose(result.E, E)

    def test_subspace_tracking_pipeline(self):
        """Feature H: subspace tracking through full pipeline."""
        n_x1, n_x2, n_states = 3, 3, 2
        E = np.zeros((n_x1, n_x2, n_states))
        D = np.zeros((n_x1, n_x2, n_states, 3))

        for i in range(n_x1):
            for j in range(n_x2):
                E[i, j, 0] = 2.0 + 0.001 * i  # near-degenerate
                E[i, j, 1] = 2.005 + 0.001 * j
                D[i, j, 0] = [1.0, 0.0, 0.0]
                D[i, j, 1] = [0.0, 1.0, 0.0]

        config = TrackingConfig(
            ref_point=(1, 1),
            use_subspace=True,
            degen_threshold=0.1,
            sigma_E=1.0,
        )
        result = track_surfaces(E, D, config)
        assert result.E.shape == (3, 3, 2)

    def test_energy_derivative_cost_pipeline(self, random_permuted_5x5):
        """Feature I: energy-derivative cost through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2),
            sigma_E=0.1,
            w_dE=1.0,
            sigma_dE=0.1,
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1

    def test_dipole_derivative_cost_pipeline(self, random_permuted_5x5):
        """Feature J: dipole-derivative cost through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point=(2, 2),
            sigma_E=0.1,
            w_dD=1.0,
            sigma_dD=0.1,
        )
        result = track_surfaces(E_perm, D_perm, config)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1

    def test_all_features_combined(self, random_permuted_5x5):
        """All features (A-J) combined through full pipeline."""
        E_perm, D_perm, E_true, D_true, perms = random_permuted_5x5
        config = TrackingConfig(
            ref_point="auto",
            sigma_E=0.1,
            use_confidence=True,
            use_priority_queue=True,
            use_multi_neighbor=True,
            n_repair_iter=2,
            adaptive_sigma_E=True,
            w_norm=0.5,
            use_subspace=True,
            degen_threshold=0.01,
            w_dE=0.5,
            sigma_dE=0.1,
            w_dD=0.5,
            sigma_dD=0.1,
        )
        result = track_surfaces(E_perm, D_perm, config)
        assert result.E.shape == (5, 5, 3)
        max_diff = np.max(np.abs(np.diff(result.E, axis=0)))
        assert max_diff < 0.1
