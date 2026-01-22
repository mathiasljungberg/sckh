"""Tests for dipole module."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from python_scripts.dynamics_1d.dipole import (
    Dipole1D,
    create_constant_dipole,
    create_dipole_from_file,
)


class TestDipole1D:
    """Tests for Dipole1D class."""

    @pytest.fixture
    def position_grid(self):
        """Standard position grid for testing (in SI units)."""
        x_ang = np.linspace(-0.5, 0.5, 101)
        x_m = x_ang * 1e-10
        return x_m

    @pytest.fixture
    def linear_dipole(self, position_grid):
        """Dipole that varies linearly with position."""
        # d_x = x, d_y = 2x, d_z = 3x (in some units)
        d = np.column_stack([
            position_grid * 1e10,  # Scale to reasonable values
            2 * position_grid * 1e10,
            3 * position_grid * 1e10,
        ])
        return Dipole1D(x=position_grid, d=d)

    @pytest.fixture
    def constant_dipole(self, position_grid):
        """Constant dipole moment surface."""
        d_value = np.array([1.0, 2.0, 3.0])
        d = np.tile(d_value, (len(position_grid), 1))
        return Dipole1D(x=position_grid, d=d)

    def test_dipole_at_grid_points(self, linear_dipole, position_grid):
        """Spline should exactly reproduce dipole values at grid points."""
        expected = np.column_stack([
            position_grid * 1e10,
            2 * position_grid * 1e10,
            3 * position_grid * 1e10,
        ])
        computed = linear_dipole.dipole(position_grid)
        np.testing.assert_allclose(computed, expected, rtol=1e-12)

    def test_dipole_between_grid_points(self, linear_dipole, position_grid):
        """Spline should interpolate smoothly between grid points."""
        # Test at midpoints
        x_mid = (position_grid[:-1] + position_grid[1:]) / 2

        # For linear function, spline should be exact
        expected = np.column_stack([
            x_mid * 1e10,
            2 * x_mid * 1e10,
            3 * x_mid * 1e10,
        ])
        computed = linear_dipole.dipole(x_mid)
        np.testing.assert_allclose(computed, expected, rtol=1e-6)

    def test_constant_dipole_interpolation(self, constant_dipole, position_grid):
        """Constant dipole should give same value everywhere."""
        x_test = np.array([position_grid[10], position_grid[50], position_grid[90]])
        expected = np.array([[1.0, 2.0, 3.0]] * 3)
        computed = constant_dipole.dipole(x_test)
        np.testing.assert_allclose(computed, expected, rtol=1e-12)

    def test_dipole_scalar_input(self, linear_dipole):
        """Dipole should handle scalar input correctly."""
        x_test = 0.0  # Equilibrium
        result = linear_dipole.dipole(x_test)

        assert result.shape == (3,)
        np.testing.assert_allclose(result, [0.0, 0.0, 0.0], atol=1e-12)

    def test_dipole_derivative(self, linear_dipole, position_grid):
        """Derivative of linear dipole should be constant."""
        x_test = position_grid[10:90]  # Avoid boundaries
        deriv = linear_dipole.dipole_derivative(x_test)

        # For linear function d = [x, 2x, 3x] * 1e10
        # derivative should be [1, 2, 3] * 1e10 (in units of 1/m)
        expected_deriv = np.array([1e10, 2e10, 3e10])

        for i in range(len(x_test)):
            np.testing.assert_allclose(deriv[i], expected_deriv, rtol=1e-3)

    def test_properties(self, linear_dipole, position_grid):
        """Test property accessors."""
        assert linear_dipole.x_min == position_grid[0]
        assert linear_dipole.x_max == position_grid[-1]
        assert linear_dipole.npoints == len(position_grid)

    def test_invalid_dipole_shape(self, position_grid):
        """Should raise error for invalid dipole array shape."""
        # Wrong number of components
        d_wrong = np.zeros((len(position_grid), 2))
        with pytest.raises(ValueError, match="shape"):
            Dipole1D(x=position_grid, d=d_wrong)

        # 1D array
        d_1d = np.zeros(len(position_grid))
        with pytest.raises(ValueError, match="shape"):
            Dipole1D(x=position_grid, d=d_1d)

    def test_mismatched_lengths(self, position_grid):
        """Should raise error for mismatched x and d lengths."""
        d_wrong = np.zeros((len(position_grid) + 1, 3))
        with pytest.raises(ValueError, match="same number of points"):
            Dipole1D(x=position_grid, d=d_wrong)


class TestCreateConstantDipole:
    """Tests for create_constant_dipole factory function."""

    def test_constant_dipole_creation(self):
        """Create constant dipole with given value."""
        x = np.linspace(-1e-10, 1e-10, 51)
        d_value = np.array([0.5, 1.0, 1.5])

        dipole = create_constant_dipole(x, d_value)

        # Check at various points
        x_test = np.array([x[0], x[25], x[-1]])
        result = dipole.dipole(x_test)

        for i in range(3):
            np.testing.assert_allclose(result[i], d_value, rtol=1e-12)

    def test_constant_dipole_wrong_shape(self):
        """Should raise error for wrong d_value shape."""
        x = np.linspace(-1e-10, 1e-10, 51)

        with pytest.raises(ValueError, match="3 components"):
            create_constant_dipole(x, np.array([1.0, 2.0]))  # Only 2 components


class TestDipoleFileIO:
    """Tests for dipole file I/O."""

    def test_read_write_dipole_file(self):
        """Test reading dipole from file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"

            # Create test file: x d_x d_y d_z
            x_ang = np.linspace(0.5, 2.5, 21)
            d_x = np.sin(x_ang)
            d_y = np.cos(x_ang)
            d_z = np.ones_like(x_ang)

            data = np.column_stack([x_ang, d_x, d_y, d_z])
            np.savetxt(filepath, data)

            # Read back
            dipole = create_dipole_from_file(filepath, units="angstrom")

            # Check properties
            assert dipole.npoints == 21
            np.testing.assert_allclose(dipole.x, x_ang * 1e-10, rtol=1e-10)

            # Check dipole values at grid points
            d_interp = dipole.dipole(dipole.x)
            np.testing.assert_allclose(d_interp[:, 0], d_x, rtol=1e-10)
            np.testing.assert_allclose(d_interp[:, 1], d_y, rtol=1e-10)
            np.testing.assert_allclose(d_interp[:, 2], d_z, rtol=1e-10)

    def test_read_dipole_bohr_units(self):
        """Test reading dipole file with bohr units."""
        from python_scripts.dynamics_1d.constants import CONST

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"

            # Create test file in bohr
            x_bohr = np.linspace(1.0, 5.0, 21)
            d = np.zeros((21, 3))
            d[:, 0] = 1.0  # Constant dipole

            data = np.column_stack([x_bohr, d])
            np.savetxt(filepath, data)

            # Read with bohr units
            dipole = create_dipole_from_file(filepath, units="bohr")

            # Check x conversion
            expected_x = x_bohr * CONST.bohr
            np.testing.assert_allclose(dipole.x, expected_x, rtol=1e-10)
