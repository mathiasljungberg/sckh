"""Tests for 2D Potential Energy Surface."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from python_scripts.dynamics_2d.pes import PES2D, create_harmonic_pes_2d, create_pes_from_file_2d
from python_scripts.dynamics_2d.io import write_pes_file_2d, read_pes_file_2d
from python_scripts.dynamics_1d.constants import CONST


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


class TestIndexOrdering:
    """Tests for index ordering in PES file I/O."""

    def test_c_order_roundtrip(self, position_grid_2d):
        """Test writing and reading PES with C-style ordering."""
        x1, x2 = position_grid_2d
        # Create simple test PES
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + 2 * X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test_c_order.dat"

            # Write with C ordering
            write_pes_file_2d(
                filepath,
                x1,
                x2,
                E_orig,
                position_units="angstrom",
                energy_units="hartree",
                index_order="C",
            )

            # Read back with C ordering
            x1_read, x2_read, E_read = read_pes_file_2d(
                filepath,
                position_units="angstrom",
                energy_units="hartree",
                index_order="C",
            )

            # Check values match
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(E_read, E_orig)

    def test_fortran_order_roundtrip(self, position_grid_2d):
        """Test writing and reading PES with Fortran-style ordering."""
        x1, x2 = position_grid_2d
        # Create simple test PES
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + 2 * X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test_f_order.dat"

            # Write with F ordering
            write_pes_file_2d(
                filepath,
                x1,
                x2,
                E_orig,
                position_units="angstrom",
                energy_units="hartree",
                index_order="F",
            )

            # Read back with F ordering
            x1_read, x2_read, E_read = read_pes_file_2d(
                filepath,
                position_units="angstrom",
                energy_units="hartree",
                index_order="F",
            )

            # Check values match
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(E_read, E_orig)

    def test_c_vs_f_ordering_different(self, position_grid_2d):
        """Test that C and F ordering produce different file contents."""
        x1, x2 = position_grid_2d
        # Use non-symmetric PES so ordering matters
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + 3 * X2**3  # Asymmetric

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath_c = Path(tmpdir) / "test_c.dat"
            filepath_f = Path(tmpdir) / "test_f.dat"

            # Write with both orderings
            write_pes_file_2d(filepath_c, x1, x2, E_orig, index_order="C")
            write_pes_file_2d(filepath_f, x1, x2, E_orig, index_order="F")

            # Read raw file contents
            with open(filepath_c) as f:
                lines_c = f.readlines()
            with open(filepath_f) as f:
                lines_f = f.readlines()

            # Files should have same length but different content
            assert len(lines_c) == len(lines_f)
            # Skip header and compare data lines
            data_lines_c = [line for line in lines_c if not line.startswith("#")]
            data_lines_f = [line for line in lines_f if not line.startswith("#")]
            # They should differ (except potentially at some points)
            assert data_lines_c != data_lines_f

    def test_wrong_index_order_detected(self, position_grid_2d):
        """Test that reading with wrong index order raises error."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + 2 * X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"

            # Write with C ordering
            write_pes_file_2d(filepath, x1, x2, E_orig, index_order="C")

            # Try to read with F ordering - should fail validation
            with pytest.raises(ValueError, match="Index ordering mismatch"):
                read_pes_file_2d(filepath, index_order="F")

    def test_invalid_index_order(self, position_grid_2d):
        """Test that invalid index_order raises ValueError."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"

            with pytest.raises(ValueError, match="index_order must be"):
                write_pes_file_2d(filepath, x1, x2, E, index_order="INVALID")

            # Write valid file first
            write_pes_file_2d(filepath, x1, x2, E)

            with pytest.raises(ValueError, match="index_order must be"):
                read_pes_file_2d(filepath, index_order="INVALID")


class TestEnergyUnits:
    """Tests for energy unit conversion in PES file I/O."""

    def test_hartree_roundtrip(self, position_grid_2d):
        """Test Hartree energy units."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test_hartree.dat"

            write_pes_file_2d(filepath, x1, x2, E_orig, energy_units="hartree")
            x1_read, x2_read, E_read = read_pes_file_2d(
                filepath, energy_units="hartree"
            )

            np.testing.assert_allclose(E_read, E_orig)

    def test_ev_roundtrip(self, position_grid_2d):
        """Test eV energy units."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test_ev.dat"

            write_pes_file_2d(filepath, x1, x2, E_orig, energy_units="ev")
            x1_read, x2_read, E_read = read_pes_file_2d(filepath, energy_units="ev")

            np.testing.assert_allclose(E_read, E_orig)

    def test_hartree_vs_ev_conversion(self, position_grid_2d):
        """Test that Hartree and eV conversions are consistent."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath_hartree = Path(tmpdir) / "test_hartree.dat"
            filepath_ev = Path(tmpdir) / "test_ev.dat"

            # Write in Hartree
            write_pes_file_2d(
                filepath_hartree, x1, x2, E_orig, energy_units="hartree"
            )
            # Write in eV
            write_pes_file_2d(filepath_ev, x1, x2, E_orig, energy_units="ev")

            # Read both back
            _, _, E_hartree = read_pes_file_2d(
                filepath_hartree, energy_units="hartree"
            )
            _, _, E_ev = read_pes_file_2d(filepath_ev, energy_units="ev")

            # Both should give same result in SI
            np.testing.assert_allclose(E_hartree, E_ev)
            np.testing.assert_allclose(E_hartree, E_orig)

    def test_mixed_units_conversion(self, position_grid_2d):
        """Test reading/writing with different energy units."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"

            # Write in Hartree
            write_pes_file_2d(filepath, x1, x2, E_orig, energy_units="hartree")

            # Read as eV - should apply correct conversion
            _, _, E_read_as_ev = read_pes_file_2d(filepath, energy_units="ev")

            # File contains Hartree values, but we told reader it's eV
            # So E_read_as_ev = (file_value) * eV_to_J
            #                 = (E_orig / hartree_to_J) * eV_to_J
            #                 = E_orig * (eV_to_J / hartree_to_J)
            expected = E_orig * (CONST.eV / CONST.hartree)
            np.testing.assert_allclose(E_read_as_ev, expected)

    def test_invalid_energy_units(self, position_grid_2d):
        """Test that invalid energy_units raises ValueError."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"

            with pytest.raises(ValueError, match="Unknown energy units"):
                write_pes_file_2d(filepath, x1, x2, E, energy_units="INVALID")

            # Write valid file first
            write_pes_file_2d(filepath, x1, x2, E)

            with pytest.raises(ValueError, match="Unknown energy units"):
                read_pes_file_2d(filepath, energy_units="INVALID")


class TestPESFactoryFunction:
    """Tests for create_pes_from_file_2d with new parameters."""

    def test_create_from_file_with_all_parameters(self, position_grid_2d):
        """Test factory function with all new parameters."""
        x1, x2 = position_grid_2d
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        E_orig = X1**2 + X2**2

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "test.dat"

            # Write with specific parameters
            write_pes_file_2d(
                filepath,
                x1,
                x2,
                E_orig,
                position_units="bohr",
                energy_units="ev",
                index_order="F",
            )

            # Create PES object with matching parameters
            pes = create_pes_from_file_2d(
                filepath,
                position_units="bohr",
                energy_units="ev",
                index_order="F",
            )

            # Test a few points
            assert isinstance(pes, PES2D)
            E_test = pes.energy(x1[50], x2[50])
            assert E_test == pytest.approx(E_orig[50, 50], rel=1e-4)
