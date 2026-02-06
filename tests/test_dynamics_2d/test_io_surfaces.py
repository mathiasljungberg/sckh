"""Tests for generic 2D surface I/O and index-order inference."""

from pathlib import Path
import tempfile

import numpy as np
import pytest

from python_scripts.dynamics_2d.io import (
    read_dipole_file_2d,
    read_dipole_file_2d_raw,
    read_surface_file_2d_raw,
)
from python_scripts.dynamics_1d.constants import CONST


def _write_surface_file(
    filepath: Path,
    x1: np.ndarray,
    x2: np.ndarray,
    values: np.ndarray,
    index_order: str,
) -> None:
    """Write synthetic surface file with one or more value columns."""
    rows = []
    n_values = values.shape[0]

    if index_order == "C":
        for i in range(len(x1)):
            for j in range(len(x2)):
                row = [x1[i], x2[j]]
                row.extend(values[k, i, j] for k in range(n_values))
                rows.append(row)
    elif index_order == "F":
        for j in range(len(x2)):
            for i in range(len(x1)):
                row = [x1[i], x2[j]]
                row.extend(values[k, i, j] for k in range(n_values))
                rows.append(row)
    else:
        raise ValueError("index_order must be C or F")

    np.savetxt(filepath, np.array(rows), fmt="%.10e")


class TestReadSurfaceFile2DRaw:
    def test_auto_infers_c_order(self):
        x1 = np.array([0.5, 0.6, 0.7])
        x2 = np.array([1.0, 1.2, 1.4, 1.6])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        v1 = X1 + 2.0 * X2
        values = np.stack([v1], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "surface_c.dat"
            _write_surface_file(filepath, x1, x2, values, "C")

            x1_read, x2_read, vals_read, order = read_surface_file_2d_raw(
                filepath,
                value_columns=[2],
                index_order="auto",
            )

            assert order == "C"
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(vals_read[0], v1)

    def test_auto_infers_f_order(self):
        x1 = np.array([0.5, 0.6, 0.7])
        x2 = np.array([1.0, 1.2, 1.4, 1.6])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        v1 = X1**2 - X2
        values = np.stack([v1], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "surface_f.dat"
            _write_surface_file(filepath, x1, x2, values, "F")

            x1_read, x2_read, vals_read, order = read_surface_file_2d_raw(
                filepath,
                value_columns=[2],
                index_order="auto",
            )

            assert order == "F"
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(vals_read[0], v1)

    def test_wrong_explicit_order_raises(self):
        x1 = np.array([0.5, 0.6, 0.7])
        x2 = np.array([1.0, 1.2, 1.4, 1.6])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        v1 = X1 + X2
        values = np.stack([v1], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "surface_c.dat"
            _write_surface_file(filepath, x1, x2, values, "C")

            with pytest.raises(ValueError, match="Index ordering mismatch"):
                read_surface_file_2d_raw(filepath, value_columns=[2], index_order="F")

    def test_multi_column_values(self):
        x1 = np.array([0.0, 0.1])
        x2 = np.array([0.2, 0.3, 0.4])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        mux = X1 + X2
        muy = X1 - X2
        muz = X1 * X2
        values = np.stack([mux, muy, muz], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"
            _write_surface_file(filepath, x1, x2, values, "C")

            x1_read, x2_read, vals_read, order = read_surface_file_2d_raw(
                filepath,
                value_columns=[2, 3, 4],
                index_order="auto",
            )

            assert order == "C"
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(vals_read[0], mux)
            np.testing.assert_allclose(vals_read[1], muy)
            np.testing.assert_allclose(vals_read[2], muz)

    def test_degenerate_grid_auto_defaults_c(self):
        x1 = np.array([0.0])
        x2 = np.array([0.2, 0.3, 0.4])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        v1 = X1 + X2
        values = np.stack([v1], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "degenerate.dat"
            _write_surface_file(filepath, x1, x2, values, "F")

            _, _, _, order = read_surface_file_2d_raw(
                filepath,
                value_columns=[2],
                index_order="auto",
            )

            assert order == "C"


class TestReadDipoleFile2D:
    def test_read_dipole_file_2d_raw(self):
        x1 = np.array([0.0, 0.1])
        x2 = np.array([0.2, 0.3, 0.4])
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        mux = X1 + X2
        muy = X1 - X2
        muz = X1 * X2
        values = np.stack([mux, muy, muz], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"
            _write_surface_file(filepath, x1, x2, values, "F")

            x1_read, x2_read, d_read, order = read_dipole_file_2d_raw(
                filepath,
                index_order="auto",
            )

            assert order == "F"
            np.testing.assert_allclose(x1_read, x1)
            np.testing.assert_allclose(x2_read, x2)
            np.testing.assert_allclose(d_read, values)

    def test_read_dipole_file_2d_position_conversion(self):
        x1_ang = np.array([0.5, 0.6])
        x2_ang = np.array([0.7, 0.8])
        X1, X2 = np.meshgrid(x1_ang, x2_ang, indexing="ij")
        values = np.stack([X1, X2, X1 + X2], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"
            _write_surface_file(filepath, x1_ang, x2_ang, values, "C")

            x1_si, x2_si, d_read, order = read_dipole_file_2d(
                filepath,
                position_units="angstrom",
                index_order="auto",
            )

            assert order == "C"
            np.testing.assert_allclose(x1_si, x1_ang * 1e-10)
            np.testing.assert_allclose(x2_si, x2_ang * 1e-10)
            np.testing.assert_allclose(d_read, values)

    def test_read_dipole_file_2d_bohr_conversion(self):
        x1_bohr = np.array([1.0, 1.5])
        x2_bohr = np.array([2.0, 2.5])
        X1, X2 = np.meshgrid(x1_bohr, x2_bohr, indexing="ij")
        values = np.stack([X1, X2, X1 - X2], axis=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole_bohr.dat"
            _write_surface_file(filepath, x1_bohr, x2_bohr, values, "C")

            x1_si, x2_si, d_read, _ = read_dipole_file_2d(
                filepath,
                position_units="bohr",
                index_order="auto",
            )

            np.testing.assert_allclose(x1_si, x1_bohr * CONST.bohr)
            np.testing.assert_allclose(x2_si, x2_bohr * CONST.bohr)
            np.testing.assert_allclose(d_read, values)
