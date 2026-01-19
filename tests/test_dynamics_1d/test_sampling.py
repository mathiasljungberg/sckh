"""Tests for sampling module."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.sampling import (
    sample_even,
    sample_random,
    compute_momentum_distribution,
    sample_momenta,
    create_initial_conditions,
)
from python_scripts.dynamics_1d.constants import CONST


@pytest.fixture
def gaussian_wavefunction(position_grid):
    """Normalized Gaussian wavefunction centered at origin."""
    sigma = 0.1e-10  # 0.1 Angstrom width
    psi = np.exp(-position_grid**2 / (2 * sigma**2))
    dx = position_grid[1] - position_grid[0]
    norm = np.sqrt(np.sum(psi**2) * dx)
    return psi / norm


class TestSampleEven:
    """Tests for even (CDF-based) sampling."""

    def test_correct_number_of_samples(self, position_grid, gaussian_wavefunction):
        """Should return requested number of samples."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2
        n_samples = 20

        x_samples = sample_even(position_grid, psi_squared, n_samples)

        assert len(x_samples) == n_samples

    def test_samples_within_grid_bounds(self, position_grid, gaussian_wavefunction):
        """Samples should be within grid bounds."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2

        x_samples = sample_even(position_grid, psi_squared, 50)

        assert np.all(x_samples >= position_grid[0])
        assert np.all(x_samples <= position_grid[-1])

    def test_samples_concentrated_at_distribution_center(
        self, position_grid, gaussian_wavefunction
    ):
        """Samples should be concentrated where |ψ|² is large."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2

        x_samples = sample_even(position_grid, psi_squared, 100)

        # For Gaussian centered at 0, most samples should be near 0
        # At least 68% within 1 sigma (0.1 Angstrom = 1e-11 m)
        sigma = 0.1e-10
        within_one_sigma = np.sum(np.abs(x_samples) < sigma) / len(x_samples)
        assert within_one_sigma > 0.5  # Should be close to 0.68

    def test_even_distribution_in_cdf_space(self, position_grid, gaussian_wavefunction):
        """Samples should be evenly spaced in CDF space."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2
        n_samples = 10

        x_samples = sample_even(position_grid, psi_squared, n_samples)

        # Compute CDF values at sample points
        dx = position_grid[1] - position_grid[0]
        cdf = np.cumsum(psi_squared) * dx
        cdf_norm = cdf / cdf[-1]

        # Interpolate to get CDF at sample points
        cdf_at_samples = np.interp(x_samples, position_grid, cdf_norm)

        # Should be evenly spaced: 0.05, 0.15, 0.25, ..., 0.95 for n=10
        expected_cdf = (np.arange(n_samples) + 0.5) / n_samples
        np.testing.assert_allclose(cdf_at_samples, expected_cdf, rtol=1e-2)


class TestSampleRandom:
    """Tests for random (inverse CDF) sampling."""

    def test_correct_number_of_samples(self, position_grid, gaussian_wavefunction, rng):
        """Should return requested number of samples."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2
        n_samples = 30

        x_samples = sample_random(position_grid, psi_squared, n_samples, rng)

        assert len(x_samples) == n_samples

    def test_reproducible_with_same_rng(self, position_grid, gaussian_wavefunction):
        """Same RNG seed should give same samples."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2

        rng1 = np.random.default_rng(12345)
        rng2 = np.random.default_rng(12345)

        x1 = sample_random(position_grid, psi_squared, 20, rng1)
        x2 = sample_random(position_grid, psi_squared, 20, rng2)

        np.testing.assert_array_equal(x1, x2)

    def test_samples_follow_distribution(self, position_grid, gaussian_wavefunction, rng):
        """Large number of samples should approximate the distribution."""
        psi_squared = np.abs(gaussian_wavefunction) ** 2
        n_samples = 1000

        x_samples = sample_random(position_grid, psi_squared, n_samples, rng)

        # Compute histogram and compare to expected distribution
        dx = position_grid[1] - position_grid[0]
        bins = position_grid
        hist, _ = np.histogram(x_samples, bins=bins, density=True)

        # Normalize expected distribution
        psi_squared_norm = psi_squared / (np.sum(psi_squared) * dx)

        # Compare histograms (should match approximately)
        # Use only central region where there are enough counts
        mask = psi_squared_norm > psi_squared_norm.max() * 0.1
        mask = mask[:-1]  # Histogram has one fewer bin

        # Large tolerance due to statistical noise
        correlation = np.corrcoef(hist[mask], psi_squared_norm[:-1][mask])[0, 1]
        assert correlation > 0.9


class TestComputeMomentumDistribution:
    """Tests for momentum-space wavefunction via FFT."""

    def test_gaussian_transforms_to_gaussian(self, position_grid, gaussian_wavefunction):
        """Gaussian in position space should give Gaussian in momentum space."""
        p, psi_p_squared = compute_momentum_distribution(
            position_grid, gaussian_wavefunction
        )

        # Find where momentum distribution peaks
        idx_max = np.argmax(psi_p_squared)

        # For Gaussian centered at x=0, momentum should peak at p=0
        np.testing.assert_allclose(p[idx_max], 0.0, atol=p[1] - p[0])

    def test_momentum_distribution_shape(
        self, position_grid, gaussian_wavefunction
    ):
        """Momentum distribution should have correct shape and be positive."""
        p, psi_p_squared = compute_momentum_distribution(
            position_grid, gaussian_wavefunction
        )

        # Should have same length as input
        assert len(p) == len(position_grid)
        assert len(psi_p_squared) == len(position_grid)

        # Should be non-negative (it's |psi|^2)
        assert np.all(psi_p_squared >= 0)


class TestCreateInitialConditions:
    """Tests for creating (x, p) initial condition pairs."""

    def test_mode1_correct_total_samples(
        self, position_grid, gaussian_wavefunction
    ):
        """Mode 1 should give n_x * n_p total samples."""
        n_x, n_p = 5, 7

        x_samples, p_samples = create_initial_conditions(
            position_grid, gaussian_wavefunction, n_x, n_p, mode=1
        )

        assert len(x_samples) == n_x * n_p
        assert len(p_samples) == n_x * n_p

    def test_mode2_correct_total_samples(
        self, position_grid, gaussian_wavefunction, rng
    ):
        """Mode 2 should give n_x total samples (paired x, p)."""
        n_x, n_p = 10, 10  # n_p ignored in mode 2

        x_samples, p_samples = create_initial_conditions(
            position_grid, gaussian_wavefunction, n_x, n_p, mode=2, rng=rng
        )

        assert len(x_samples) == n_x
        assert len(p_samples) == n_x

    def test_mode1_all_combinations(self, position_grid, gaussian_wavefunction):
        """Mode 1 should produce all (x, p) combinations."""
        n_x, n_p = 3, 4

        x_samples, p_samples = create_initial_conditions(
            position_grid, gaussian_wavefunction, n_x, n_p, mode=1
        )

        # Check that we have all combinations
        x_unique = np.unique(x_samples)
        p_unique = np.unique(p_samples)

        assert len(x_unique) == n_x
        assert len(p_unique) == n_p
