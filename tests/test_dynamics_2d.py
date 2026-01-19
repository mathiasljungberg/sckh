"""
Pytest test suite for dynamics_2d.py

This test file incorporates the implicit tests from src/python_scripts/dynamics_2d.py
(lines 41-46 and 88-93) into a proper pytest structure.

Note: The dynamics_2d.py module has a bug on line 95 ('ssasda') that will cause
the __main__ block to fail. These tests only import the functions, so they are
not affected by that bug.
"""

import numpy as np
import pytest
from scipy.interpolate import RectBivariateSpline

from python_scripts.dynamics_2d import poly_on_grid, product_poly_2d


# Fixtures for common test data
@pytest.fixture
def coarse_grid():
    """Coarse grid for initial testing (5x5 points)"""
    x = np.linspace(-1, 1, 5)
    y = np.linspace(-1, 1, 5)
    return x, y


@pytest.fixture
def fine_grid():
    """Fine grid for interpolation testing (200x200 points)"""
    xf = np.linspace(-1, 1, 200)
    yf = np.linspace(-1, 1, 200)
    return xf, yf


@pytest.fixture
def poly_coeffs():
    """Polynomial coefficients used in the original tests"""
    return {
        'coeffs_x': [1, 2, 3, 4],
        'coeffs_y': [1, 2, 3]
    }


@pytest.fixture
def product_surface(coarse_grid, poly_coeffs):
    """Generate product polynomial surface on coarse grid"""
    x, y = coarse_grid
    return product_poly_2d(
        x, y, x0=0, y0=0,
        coeffs_x=poly_coeffs['coeffs_x'],
        coeffs_y=poly_coeffs['coeffs_y']
    )


@pytest.fixture
def spline_from_surface(coarse_grid, product_surface):
    """Create spline interpolation from product surface"""
    x, y = coarse_grid
    Z = product_surface['f']
    return RectBivariateSpline(x, y, Z, kx=3, ky=3)


# Test 1: Polynomial evaluation
def test_poly_on_grid():
    """Test that poly_on_grid correctly evaluates polynomials and derivatives"""
    x = np.array([-1.0, 0.0, 1.0])
    coeffs = [1, 2, 3]  # p(x) = 1 + 2x + 3x^2, dp/dx = 2 + 6x

    f, df = poly_on_grid(x, x0=0, coeffs=coeffs)

    # Check polynomial values
    expected_f = np.array([1 - 2 + 3, 1, 1 + 2 + 3])  # [2, 1, 6]
    np.testing.assert_allclose(f, expected_f, rtol=1e-12)

    # Check derivative values
    expected_df = np.array([2 - 6, 2, 2 + 6])  # [-4, 2, 8]
    np.testing.assert_allclose(df, expected_df, rtol=1e-12)


# Test 2: Product polynomial structure
def test_product_poly_2d_structure(product_surface, coarse_grid):
    """Test that product_poly_2d returns correct dictionary structure and shapes"""
    x, y = coarse_grid

    # Check that all expected keys are present
    assert 'f' in product_surface
    assert 'df_dx' in product_surface
    assert 'df_dy' in product_surface

    # Check shapes
    expected_shape = (len(x), len(y))
    assert product_surface['f'].shape == expected_shape
    assert product_surface['df_dx'].shape == expected_shape
    assert product_surface['df_dy'].shape == expected_shape


# Test 3: Spline interpolation at grid points (from lines 41-46)
def test_spline_interpolation_at_grid_points(coarse_grid, product_surface, spline_from_surface):
    """
    Test that spline interpolation matches analytical values at grid points.

    This test extracts the assertions from lines 41-46 of dynamics_2d.py:
    - Function values should match exactly
    - x-derivatives should match exactly
    - y-derivatives should match exactly
    """
    x, y = coarse_grid
    spl = spline_from_surface

    # Test function values (dx=0, dy=0)
    f_spl = spl(x, y, dx=0, dy=0)
    assert np.max(np.abs(f_spl - product_surface['f'])) < 1e-12, \
        "Spline function values don't match at grid points"

    # Test x-derivatives (dx=1, dy=0)
    df_dx_spl = spl(x, y, dx=1, dy=0)
    assert np.max(np.abs(df_dx_spl - product_surface['df_dx'])) < 1e-12, \
        "Spline x-derivatives don't match at grid points"

    # Test y-derivatives (dx=0, dy=1)
    df_dy_spl = spl(x, y, dx=0, dy=1)
    assert np.max(np.abs(df_dy_spl - product_surface['df_dy'])) < 1e-12, \
        "Spline y-derivatives don't match at grid points"


# Test 4: Spline interpolation on fine grid (from lines 88-93)
def test_spline_interpolation_at_fine_grid(fine_grid, poly_coeffs, spline_from_surface):
    """
    Test that spline interpolation matches analytical values on a fine grid.

    This test extracts the assertions from lines 88-93 of dynamics_2d.py:
    - Spline should accurately interpolate between grid points
    - Derivatives should also be accurate on the fine grid
    """
    xf, yf = fine_grid
    spl = spline_from_surface

    # Compute analytical values on fine grid
    retf = product_poly_2d(
        xf, yf, x0=0, y0=0,
        coeffs_x=poly_coeffs['coeffs_x'],
        coeffs_y=poly_coeffs['coeffs_y']
    )

    # Test function values
    f_spl = spl(xf, yf, dx=0, dy=0)
    assert np.max(np.abs(f_spl - retf['f'])) < 1e-12, \
        "Spline function values don't match on fine grid"

    # Test x-derivatives
    df_dx_spl = spl(xf, yf, dx=1, dy=0)
    assert np.max(np.abs(df_dx_spl - retf['df_dx'])) < 1e-12, \
        "Spline x-derivatives don't match on fine grid"

    # Test y-derivatives
    df_dy_spl = spl(xf, yf, dx=0, dy=1)
    assert np.max(np.abs(df_dy_spl - retf['df_dy'])) < 1e-12, \
        "Spline y-derivatives don't match on fine grid"


# Test 5: Known polynomial values
def test_product_poly_2d_values():
    """Test product_poly_2d with known values for verification"""
    # Simple test: use constant polynomials
    x = np.array([0.0, 1.0])
    y = np.array([0.0, 1.0])

    # f(x) = 2, g(y) = 3, so product should be 6 everywhere
    result = product_poly_2d(x, y, x0=0, y0=0, coeffs_x=[2], coeffs_y=[3])

    expected_f = np.full((2, 2), 6.0)
    np.testing.assert_allclose(result['f'], expected_f, rtol=1e-12)

    # Derivatives of constants should be zero
    expected_df = np.zeros((2, 2))
    np.testing.assert_allclose(result['df_dx'], expected_df, rtol=1e-12)
    np.testing.assert_allclose(result['df_dy'], expected_df, rtol=1e-12)


# Test 6: Polynomial shift parameter
def test_poly_on_grid_with_shift():
    """Test that x0 shift parameter works correctly"""
    x = np.array([0.0, 1.0, 2.0])
    coeffs = [1, 2]  # p(x) = 1 + 2x

    # Without shift: p(0)=1, p(1)=3, p(2)=5
    f_no_shift, _ = poly_on_grid(x, x0=0, coeffs=coeffs)
    expected_no_shift = np.array([1, 3, 5])
    np.testing.assert_allclose(f_no_shift, expected_no_shift, rtol=1e-12)

    # With shift x0=1: p(x-1), so p(0)=p(-1)=-1, p(1)=p(0)=1, p(2)=p(1)=3
    f_with_shift, _ = poly_on_grid(x, x0=1, coeffs=coeffs)
    expected_with_shift = np.array([-1, 1, 3])
    np.testing.assert_allclose(f_with_shift, expected_with_shift, rtol=1e-12)
