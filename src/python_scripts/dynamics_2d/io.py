"""File I/O for 2D PES, trajectory, and ground state data."""

from pathlib import Path
from typing import Tuple, Optional, TYPE_CHECKING

import numpy as np

from python_scripts.dynamics_1d.constants import CONST

if TYPE_CHECKING:
    from .trajectory import TrajectoryResult2D
    from .vibrational import ProductGroundState


def _loadtxt_fortran(filepath: Path) -> np.ndarray:
    """Load text file with Fortran D-format exponent support.

    Handles Fortran double-precision format where 'D' or 'd' is used
    instead of 'E' for exponents (e.g., '-0.345376681424D-03').

    Args:
        filepath: Path to data file

    Returns:
        numpy array of floats
    """
    from io import StringIO

    with open(filepath) as f:
        text = f.read().replace("D", "E").replace("d", "e")

    return np.loadtxt(StringIO(text))


def read_pes_file_2d(
    filepath: Path,
    units: str = "angstrom",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D PES file.

    Expected format:
        x1  x2  E
        ...

    Where x1 and x2 are coordinates and E is energy. The file should
    contain data on a regular grid with x1 varying slowest and x2
    varying fastest.

    Args:
        filepath: Path to PES file
        units: "angstrom" or "bohr" for coordinates

    Returns:
        x1: Unique x1 grid values (SI: meters)
        x2: Unique x2 grid values (SI: meters)
        E: Energy on 2D grid (SI: Joules), shape (n_x1, n_x2)
    """
    data = _loadtxt_fortran(filepath)

    x1_raw = data[:, 0]
    x2_raw = data[:, 1]
    E_raw = data[:, 2]

    # Get unique grid values
    x1_unique = np.unique(x1_raw)
    x2_unique = np.unique(x2_raw)

    n_x1 = len(x1_unique)
    n_x2 = len(x2_unique)

    # Reshape energy to 2D grid
    # Assume data is ordered with x2 varying fastest (row-major for fixed x1)
    E_2d = E_raw.reshape((n_x1, n_x2))

    # Convert to SI units
    if units.lower() == "angstrom":
        x1 = x1_unique * 1e-10  # Angstrom to meters
        x2 = x2_unique * 1e-10
    elif units.lower() == "bohr":
        x1 = x1_unique * CONST.bohr
        x2 = x2_unique * CONST.bohr
    else:
        raise ValueError(f"Unknown units: {units}")

    E = E_2d * CONST.hartree  # Hartree to Joules

    return x1, x2, E


def read_pes_file_2d_raw(
    filepath: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D PES file without unit conversion.

    Args:
        filepath: Path to PES file

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        E: Energy on 2D grid in original units, shape (n_x1, n_x2)
    """
    data = _loadtxt_fortran(filepath)

    x1_raw = data[:, 0]
    x2_raw = data[:, 1]
    E_raw = data[:, 2]

    x1_unique = np.unique(x1_raw)
    x2_unique = np.unique(x2_raw)

    n_x1 = len(x1_unique)
    n_x2 = len(x2_unique)

    E_2d = E_raw.reshape((n_x1, n_x2))

    return x1_unique, x2_unique, E_2d


def write_trajectory_2d(
    filepath: Path,
    traj: "TrajectoryResult2D",
    units: str = "SI",
) -> None:
    """Write 2D trajectory data to file.

    Args:
        filepath: Output file path
        traj: TrajectoryResult2D object
        units: "SI" or "user" (Angstrom/fs)
    """
    if units.lower() == "si":
        header = "# time(s) x1(m) x2(m) v1(m/s) v2(m/s)"
        data = np.column_stack(
            [traj.time, traj.x1, traj.x2, traj.v1, traj.v2]
        )
    else:
        # Convert to user-friendly units
        time_fs = traj.time * 1e15  # s to fs
        x1_ang = traj.x1 * 1e10  # m to Angstrom
        x2_ang = traj.x2 * 1e10
        v1_ang_fs = traj.v1 * 1e-5  # m/s to Angstrom/fs
        v2_ang_fs = traj.v2 * 1e-5
        header = "# time(fs) x1(Angstrom) x2(Angstrom) v1(Angstrom/fs) v2(Angstrom/fs)"
        data = np.column_stack(
            [time_fs, x1_ang, x2_ang, v1_ang_fs, v2_ang_fs]
        )

    np.savetxt(filepath, data, header=header, fmt="%16.8E")


def write_ground_state_2d(
    filepath: Path,
    ground_state: "ProductGroundState",
    units: str = "SI",
) -> None:
    """Write 2D product ground state information to file.

    Writes the 1D wavefunctions and energies that make up the product state.

    Args:
        filepath: Output file path
        ground_state: ProductGroundState object
        units: "SI" or "user"
    """
    if units.lower() == "si":
        header = (
            f"# Product ground state\n"
            f"# E1 = {ground_state.E1:.10E} J\n"
            f"# E2 = {ground_state.E2:.10E} J\n"
            f"# E_total = {ground_state.E_total:.10E} J\n"
            f"# Column 1: x1 (m)\n"
            f"# Column 2: psi1(x1)\n"
            f"# Column 3: x2 (m)\n"
            f"# Column 4: psi2(x2)"
        )
        # Pad shorter array with zeros if grids have different lengths
        n1 = len(ground_state.x1_grid)
        n2 = len(ground_state.x2_grid)
        n_max = max(n1, n2)

        x1_padded = np.zeros(n_max)
        psi1_padded = np.zeros(n_max)
        x2_padded = np.zeros(n_max)
        psi2_padded = np.zeros(n_max)

        x1_padded[:n1] = ground_state.x1_grid
        psi1_padded[:n1] = ground_state.psi1
        x2_padded[:n2] = ground_state.x2_grid
        psi2_padded[:n2] = ground_state.psi2

        data = np.column_stack([x1_padded, psi1_padded, x2_padded, psi2_padded])
    else:
        E1_eV = ground_state.E1 / CONST.eV
        E2_eV = ground_state.E2 / CONST.eV
        E_total_eV = ground_state.E_total / CONST.eV
        header = (
            f"# Product ground state\n"
            f"# E1 = {E1_eV:.10E} eV\n"
            f"# E2 = {E2_eV:.10E} eV\n"
            f"# E_total = {E_total_eV:.10E} eV\n"
            f"# Column 1: x1 (Angstrom)\n"
            f"# Column 2: psi1(x1)\n"
            f"# Column 3: x2 (Angstrom)\n"
            f"# Column 4: psi2(x2)"
        )
        n1 = len(ground_state.x1_grid)
        n2 = len(ground_state.x2_grid)
        n_max = max(n1, n2)

        x1_padded = np.zeros(n_max)
        psi1_padded = np.zeros(n_max)
        x2_padded = np.zeros(n_max)
        psi2_padded = np.zeros(n_max)

        x1_padded[:n1] = ground_state.x1_grid * 1e10  # m to Angstrom
        psi1_padded[:n1] = ground_state.psi1
        x2_padded[:n2] = ground_state.x2_grid * 1e10
        psi2_padded[:n2] = ground_state.psi2

        data = np.column_stack([x1_padded, psi1_padded, x2_padded, psi2_padded])

    np.savetxt(filepath, data, header=header, fmt="%16.8E")


def write_pes_file_2d(
    filepath: Path,
    x1: np.ndarray,
    x2: np.ndarray,
    E: np.ndarray,
    units: str = "angstrom",
) -> None:
    """Write 2D PES to file.

    Args:
        filepath: Output file path
        x1: Grid points for x1 (SI: meters)
        x2: Grid points for x2 (SI: meters)
        E: Energy on 2D grid (SI: Joules), shape (len(x1), len(x2))
        units: "angstrom" or "bohr" for output coordinates
    """
    if units.lower() == "angstrom":
        x1_out = x1 * 1e10  # m to Angstrom
        x2_out = x2 * 1e10
    elif units.lower() == "bohr":
        x1_out = x1 / CONST.bohr
        x2_out = x2 / CONST.bohr
    else:
        raise ValueError(f"Unknown units: {units}")

    E_out = E / CONST.hartree  # J to Hartree

    # Create output array with x1, x2, E columns
    rows = []
    for i, x1_val in enumerate(x1_out):
        for j, x2_val in enumerate(x2_out):
            rows.append([x1_val, x2_val, E_out[i, j]])

    data = np.array(rows)
    header = f"# x1({units}) x2({units}) E(Hartree)"
    np.savetxt(filepath, data, header=header, fmt="%16.8E")
