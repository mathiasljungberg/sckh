"""File I/O for PES and trajectory data."""

from pathlib import Path
from typing import Tuple, Optional

import numpy as np

from .constants import CONST


def read_pes_file(
    filepath: Path,
    units: str = "angstrom",
    energy_column: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read PES file with multi-column format (x, E1, E2, ...).

    Args:
        filepath: Path to PES file
        units: "angstrom" or "bohr" for x coordinate
        energy_column: Which column to use for energy (1-indexed, default=1 means first energy column)

    Returns:
        x: Position array in meters (SI)
        E: Energy array in Joules (SI)
    """
    data = np.loadtxt(filepath)
    x_raw = data[:, 0]
    E_raw = data[:, energy_column]  # 1-indexed energy column

    # Convert to SI units
    if units.lower() == "angstrom":
        x = x_raw * 1e-10  # Angstrom to meters
    elif units.lower() == "bohr":
        x = x_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown units: {units}")

    E = E_raw * CONST.hartree  # Hartree to Joules

    return x, E


def read_pes_file_raw(
    filepath: Path,
    energy_column: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read PES file without unit conversion.

    Args:
        filepath: Path to PES file
        energy_column: Which column to use for energy (1-indexed)

    Returns:
        x: Position array in original units (typically Angstrom)
        E: Energy array in original units (typically Hartree)
    """
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, energy_column]


def write_trajectory(
    filepath: Path,
    time: np.ndarray,
    x: np.ndarray,
    v: Optional[np.ndarray] = None,
    units: str = "SI",
) -> None:
    """Write trajectory data to file.

    Args:
        filepath: Output file path
        time: Time array
        x: Position array
        v: Velocity array (optional)
        units: "SI" or "user" (Angstrom/fs)
    """
    if units.lower() == "si":
        header = "# time(s) x(m)"
        if v is not None:
            header += " v(m/s)"
            data = np.column_stack([time, x, v])
        else:
            data = np.column_stack([time, x])
    else:
        # Convert to user-friendly units
        time_fs = time * 1e15  # s to fs
        x_ang = x * 1e10  # m to Angstrom
        header = "# time(fs) x(Angstrom)"
        if v is not None:
            v_ang_fs = v * 1e-5  # m/s to Angstrom/fs
            header += " v(Angstrom/fs)"
            data = np.column_stack([time_fs, x_ang, v_ang_fs])
        else:
            data = np.column_stack([time_fs, x_ang])

    np.savetxt(filepath, data, header=header, fmt="%16.8E")


def write_distribution(
    filepath: Path,
    x: np.ndarray,
    density: np.ndarray,
    units: str = "SI",
) -> None:
    """Write position distribution to file (Fortran-compatible format).

    Args:
        filepath: Output file path
        x: Position array
        density: Probability density array
        units: "SI" or "user" (Angstrom)
    """
    if units.lower() == "si":
        data = np.column_stack([x, density])
    else:
        x_ang = x * 1e10  # m to Angstrom
        # Density needs to be rescaled: P(x_ang) = P(x_m) * 1e-10
        density_ang = density * 1e-10
        data = np.column_stack([x_ang, density_ang])

    np.savetxt(filepath, data, fmt="%16.6E")


def write_eigenstate(
    filepath: Path,
    x: np.ndarray,
    psi: np.ndarray,
    eigenvalue: float,
    units: str = "SI",
) -> None:
    """Write eigenstate wavefunction to file.

    Args:
        filepath: Output file path
        x: Position grid
        psi: Wavefunction on grid
        eigenvalue: Energy eigenvalue
        units: "SI" or "user"
    """
    if units.lower() == "si":
        header = f"# Eigenvalue: {eigenvalue:.10E} J\n# x(m) psi(x)"
        data = np.column_stack([x, psi])
    else:
        x_ang = x * 1e10
        eigenvalue_eV = eigenvalue / CONST.eV
        header = f"# Eigenvalue: {eigenvalue_eV:.10E} eV\n# x(Angstrom) psi(x)"
        data = np.column_stack([x_ang, psi])

    np.savetxt(filepath, data, header=header, fmt="%16.8E")
