"""Tests for PES module."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.pes import PES1D, create_harmonic_pes
from python_scripts.dynamics_1d.constants import CONST


class TestPES1D:
    """Tests for PES1D class."""

    def test_energy_at_grid_points(self, position_grid, harmonic_params):
        """Spline should exactly reproduce energy values at grid points."""
        k, x0 = harmonic_params["k"], harmonic_params["x0"]
        E_expected = 0.5 * k * (position_grid - x0) ** 2

        pes = PES1D(x=position_grid, E=E_expected)
        E_computed = pes.energy(position_grid)

        np.testing.assert_allclose(E_computed, E_expected, rtol=1e-12)

    def test_energy_between_grid_points(self, position_grid, harmonic_params):
        """Spline should accurately interpolate between grid points."""
        k, x0 = harmonic_params["k"], harmonic_params["x0"]
        E_grid = 0.5 * k * (position_grid - x0) ** 2

        pes = PES1D(x=position_grid, E=E_grid)

        # Test at midpoints
        x_mid = (position_grid[:-1] + position_grid[1:]) / 2
        E_expected = 0.5 * k * (x_mid - x0) ** 2
        E_computed = pes.energy(x_mid)

        # Cubic spline interpolation of quadratic function
        # Not exact due to natural boundary conditions
        np.testing.assert_allclose(E_computed, E_expected, rtol=1e-3)

    def test_force_from_harmonic_potential(self, harmonic_pes, harmonic_params):
        """Force should equal -kx for harmonic potential V = 0.5*k*x²."""
        k, x0 = harmonic_params["k"], harmonic_params["x0"]

        # Test at several points (in the interior of the grid)
        x_test = np.array([-0.3e-10, -0.1e-10, 0.0, 0.1e-10, 0.3e-10])
        F_expected = -k * (x_test - x0)
        F_computed = harmonic_pes.force(x_test)

        # Spline derivative accuracy limited by natural boundary conditions
        # Use both rtol and atol to handle values near zero
        np.testing.assert_allclose(F_computed, F_expected, rtol=1e-2, atol=1e-20)

    def test_force_matches_finite_difference(self, position_grid, harmonic_params):
        """Analytical spline derivative should match numerical derivative."""
        k, x0 = harmonic_params["k"], harmonic_params["x0"]
        E_grid = 0.5 * k * (position_grid - x0) ** 2

        pes = PES1D(x=position_grid, E=E_grid)

        # Test at a point away from boundaries
        x_test = 0.1e-10
        h = 1e-14  # Small step for finite difference

        # Finite difference: -dV/dx ≈ -(V(x+h) - V(x-h)) / (2h)
        F_numerical = -(pes.energy(x_test + h) - pes.energy(x_test - h)) / (2 * h)
        F_analytical = pes.force(x_test)

        np.testing.assert_allclose(F_analytical, F_numerical, rtol=1e-6)

    def test_energy_and_force_consistency(self, harmonic_pes):
        """energy_and_force should return same values as individual calls."""
        x_test = np.array([0.0, 0.1e-10, 0.2e-10])

        E1, F1 = harmonic_pes.energy_and_force(x_test)
        E2 = harmonic_pes.energy(x_test)
        F2 = harmonic_pes.force(x_test)

        np.testing.assert_array_equal(E1, E2)
        np.testing.assert_array_equal(F1, F2)

    def test_properties(self, harmonic_pes, position_grid):
        """Test PES property accessors."""
        assert harmonic_pes.x_min == position_grid[0]
        assert harmonic_pes.x_max == position_grid[-1]
        assert harmonic_pes.npoints == len(position_grid)
        np.testing.assert_allclose(
            harmonic_pes.dx, position_grid[1] - position_grid[0], rtol=1e-12
        )


class TestCreateHarmonicPES:
    """Tests for create_harmonic_pes factory function."""

    def test_harmonic_pes_minimum_at_x0(self, position_grid):
        """Harmonic PES should have minimum at x0."""
        x0 = 0.0  # Equilibrium at grid center
        k = 100.0

        pes = create_harmonic_pes(position_grid, x0=x0, k=k)
        x_min, E_min = pes.find_minimum()

        # Minimum location should be close to x0
        # Tolerance ~ 0.2 Angstrom due to optimizer convergence on spline
        np.testing.assert_allclose(x_min, x0, atol=2e-11)
        # Minimum energy should be close to 0
        np.testing.assert_allclose(E_min, 0.0, atol=1e-20)

    def test_harmonic_pes_with_offset(self, position_grid):
        """Harmonic PES should respect energy offset E0."""
        E0 = 1.0e-19  # Energy offset in Joules
        pes = create_harmonic_pes(position_grid, x0=0.0, k=100.0, E0=E0)

        # Minimum energy should be E0
        _, E_min = pes.find_minimum()
        np.testing.assert_allclose(E_min, E0, rtol=0.1)  # 10% tolerance for interpolation
