"""Tests for 2D Potential Energy Surface."""

import numpy as np
import pytest

from python_scripts.dynamics_2d.pes import PES2D, create_harmonic_pes_2d


class TestPES2D:
    """Tests for PES2D class."""

    def test_harmonic_energy(self, harmonic_pes_2d, harmonic_params_2d):
        """Test energy evaluation for harmonic potential."""
        pes = harmonic_pes_2d
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]
        x1_0 = harmonic_params_2d["x1_0"]
        x2_0 = harmonic_params_2d["x2_0"]

        # Test at equilibrium
        E0 = pes.energy(x1_0, x2_0)
        assert E0 == pytest.approx(0.0, abs=1e-20)

        # Test at displaced position
        x1_test = 0.1e-10  # 0.1 Angstrom
        x2_test = 0.2e-10  # 0.2 Angstrom
        E_expected = 0.5 * k1 * (x1_test - x1_0) ** 2 + 0.5 * k2 * (x2_test - x2_0) ** 2
        E_calc = pes.energy(x1_test, x2_test)
        assert E_calc == pytest.approx(E_expected, rel=1e-6)

    def test_harmonic_forces(self, harmonic_pes_2d, harmonic_params_2d):
        """Test force evaluation for harmonic potential."""
        pes = harmonic_pes_2d
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]
        x1_0 = harmonic_params_2d["x1_0"]
        x2_0 = harmonic_params_2d["x2_0"]

        # Test at displaced position
        x1_test = 0.1e-10
        x2_test = 0.2e-10

        # For harmonic: F = -k * (x - x0)
        F1_expected = -k1 * (x1_test - x1_0)
        F2_expected = -k2 * (x2_test - x2_0)

        F1_calc = pes.force_x1(x1_test, x2_test)
        F2_calc = pes.force_x2(x1_test, x2_test)

        assert F1_calc == pytest.approx(F1_expected, rel=1e-4)
        assert F2_calc == pytest.approx(F2_expected, rel=1e-4)

    def test_forces_tuple(self, harmonic_pes_2d, harmonic_params_2d):
        """Test that forces() returns consistent tuple."""
        pes = harmonic_pes_2d

        x1_test = 0.1e-10
        x2_test = 0.2e-10

        F1, F2 = pes.forces(x1_test, x2_test)
        F1_direct = pes.force_x1(x1_test, x2_test)
        F2_direct = pes.force_x2(x1_test, x2_test)

        assert F1 == pytest.approx(F1_direct)
        assert F2 == pytest.approx(F2_direct)

    def test_second_derivatives(self, harmonic_pes_2d, harmonic_params_2d):
        """Test second derivatives for harmonic potential."""
        pes = harmonic_pes_2d
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]

        # For harmonic: d²V/dx² = k
        x1_test = 0.1e-10
        x2_test = 0.2e-10

        d2V_dx1 = pes.second_derivative_x1(x1_test, x2_test)
        d2V_dx2 = pes.second_derivative_x2(x1_test, x2_test)
        d2V_dx1dx2 = pes.mixed_derivative(x1_test, x2_test)

        assert d2V_dx1 == pytest.approx(k1, rel=1e-4)
        assert d2V_dx2 == pytest.approx(k2, rel=1e-4)
        # Mixed derivative should be zero for separable potential
        assert d2V_dx1dx2 == pytest.approx(0.0, abs=1e-10)

    def test_grid_properties(self, harmonic_pes_2d, position_grid_2d):
        """Test grid property accessors."""
        pes = harmonic_pes_2d
        x1, x2 = position_grid_2d

        assert pes.x1_min == pytest.approx(x1[0])
        assert pes.x1_max == pytest.approx(x1[-1])
        assert pes.x2_min == pytest.approx(x2[0])
        assert pes.x2_max == pytest.approx(x2[-1])
        assert pes.npoints_x1 == len(x1)
        assert pes.npoints_x2 == len(x2)

    def test_find_minimum(self, harmonic_pes_2d, harmonic_params_2d):
        """Test minimum finding for harmonic potential."""
        pes = harmonic_pes_2d
        x1_0 = harmonic_params_2d["x1_0"]
        x2_0 = harmonic_params_2d["x2_0"]

        x1_min, x2_min, E_min = pes.find_minimum()

        assert x1_min == pytest.approx(x1_0, abs=1e-12)
        assert x2_min == pytest.approx(x2_0, abs=1e-12)
        assert E_min == pytest.approx(0.0, abs=1e-20)

    def test_vectorized_evaluation(self, harmonic_pes_2d, harmonic_params_2d):
        """Test evaluation at multiple points."""
        pes = harmonic_pes_2d
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]

        # Test with arrays
        x1_arr = np.array([0.0, 0.1e-10, 0.2e-10])
        x2_arr = np.array([0.0, 0.1e-10, 0.1e-10])

        E_arr = pes.energy(x1_arr, x2_arr)
        F1_arr = pes.force_x1(x1_arr, x2_arr)

        # Check shape
        assert E_arr.shape == x1_arr.shape
        assert F1_arr.shape == x1_arr.shape

        # Check values
        for i in range(len(x1_arr)):
            E_expected = 0.5 * k1 * x1_arr[i] ** 2 + 0.5 * k2 * x2_arr[i] ** 2
            assert E_arr[i] == pytest.approx(E_expected, rel=1e-6)


class TestCreateHarmonicPES2D:
    """Tests for create_harmonic_pes_2d factory function."""

    def test_creates_pes_object(self, position_grid_2d):
        """Test that factory creates valid PES2D."""
        x1, x2 = position_grid_2d
        pes = create_harmonic_pes_2d(x1, x2, 0.0, 0.0, 100.0, 200.0)
        assert isinstance(pes, PES2D)

    def test_energy_offset(self, position_grid_2d):
        """Test energy offset parameter."""
        x1, x2 = position_grid_2d
        E0 = 1.0e-19  # 1e-19 Joules offset

        pes = create_harmonic_pes_2d(x1, x2, 0.0, 0.0, 100.0, 100.0, E0=E0)

        # At equilibrium, energy should be E0
        E_min = pes.energy(0.0, 0.0)
        assert E_min == pytest.approx(E0, rel=1e-6)

    def test_different_equilibrium_positions(self, position_grid_2d):
        """Test non-zero equilibrium positions."""
        x1, x2 = position_grid_2d
        x1_0 = 0.1e-10
        x2_0 = -0.1e-10
        k1 = 100.0
        k2 = 200.0

        pes = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k1, k2)

        # At equilibrium, energy should be zero
        E_eq = pes.energy(x1_0, x2_0)
        assert E_eq == pytest.approx(0.0, abs=1e-20)

        # Force at equilibrium should be zero
        F1_eq, F2_eq = pes.forces(x1_0, x2_0)
        assert F1_eq == pytest.approx(0.0, abs=1e-10)
        assert F2_eq == pytest.approx(0.0, abs=1e-10)
