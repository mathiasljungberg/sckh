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


def _validate_index_ordering(
    x1_raw: np.ndarray,
    x2_raw: np.ndarray,
    x1_unique: np.ndarray,
    x2_unique: np.ndarray,
    index_order: str,
) -> None:
    """Validate that raw data follows the expected index ordering.

    Creates the expected coordinate arrays from the outer product of
    x1_unique and x2_unique, then compares with the actual raw data.

    Args:
        x1_raw: Raw x1 coordinates from file
        x2_raw: Raw x2 coordinates from file
        x1_unique: Unique x1 values (sorted)
        x2_unique: Unique x2 values (sorted)
        index_order: "C" or "F" for index ordering

    Raises:
        ValueError: If data ordering doesn't match the specified index_order
    """
    # Create expected coordinate arrays based on index ordering
    X1_grid, X2_grid = np.meshgrid(x1_unique, x2_unique, indexing="ij")

    if index_order.upper() == "C":
        # C order: x2 varies fastest (row-major)
        x1_expected = X1_grid.flatten(order="C")
        x2_expected = X2_grid.flatten(order="C")
    elif index_order.upper() == "F":
        # F order: x1 varies fastest (column-major)
        x1_expected = X1_grid.flatten(order="F")
        x2_expected = X2_grid.flatten(order="F")
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

    # Compare with actual raw data
    if not (np.allclose(x1_raw, x1_expected) and np.allclose(x2_raw, x2_expected)):
        # Determine the likely correct ordering for helpful error message
        if index_order.upper() == "C":
            suggestion = "F"
        else:
            suggestion = "C"

        raise ValueError(
            f"Index ordering mismatch: Data does not match expected {index_order} "
            f"order (x2 fast for 'C', x1 fast for 'F'). "
            f"Try index_order='{suggestion}' instead."
        )


def read_pes_file_2d(
    filepath: Path,
    position_units: str = "angstrom",
    energy_units: str = "hartree",
    index_order: str = "C",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D PES file with unit conversion.

    Expected format:
        x1  x2  E
        ...

    Where x1 and x2 are coordinates and E is energy.

    Args:
        filepath: Path to PES file
        position_units: "angstrom" or "bohr" for coordinates
        energy_units: "hartree" or "ev" for energy
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering

    Returns:
        x1: Unique x1 grid values (SI: meters)
        x2: Unique x2 grid values (SI: meters)
        E: Energy on 2D grid (SI: Joules), shape (n_x1, n_x2)

    Raises:
        ValueError: If position_units, energy_units, or index_order are invalid
    """
    # Read raw data
    x1_raw, x2_raw, E_raw = read_pes_file_2d_raw(filepath, index_order)

    # Convert position units
    if position_units.lower() == "angstrom":
        x1 = x1_raw * 1e-10
        x2 = x2_raw * 1e-10
    elif position_units.lower() == "bohr":
        x1 = x1_raw * CONST.bohr
        x2 = x2_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    # Convert energy units
    if energy_units.lower() == "hartree":
        E = E_raw * CONST.hartree  # Hartree to Joules
    elif energy_units.lower() == "ev":
        E = E_raw * CONST.eV  # eV to Joules
    else:
        raise ValueError(f"Unknown energy units: {energy_units}")

    return x1, x2, E


def read_pes_file_2d_raw(
    filepath: Path,
    index_order: str = "C",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D PES file without unit conversion.

    Args:
        filepath: Path to PES file
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        E: Energy on 2D grid in original units, shape (n_x1, n_x2)

    Raises:
        ValueError: If index_order is not "C" or "F", or if data layout
                    doesn't match the specified ordering
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

    # Validate that the data ordering matches the specified index_order
    # This must be done BEFORE reshaping
    _validate_index_ordering(x1_raw, x2_raw, x1_unique, x2_unique, index_order)

    # Reshape according to index order
    if index_order.upper() == "C":
        # x2 varies fastest (C-style, row-major)
        E_2d = E_raw.reshape((n_x1, n_x2), order="C")
    elif index_order.upper() == "F":
        # x1 varies fastest (Fortran-style, column-major)
        E_2d = E_raw.reshape((n_x1, n_x2), order="F")
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

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


def read_trajectory_2d(filepath: Path) -> "TrajectoryResult2D":
    """Read 2D trajectory from file.

    Expects format: time x1 x2 v1 v2 (in SI or user units).
    Auto-detects units based on header.

    Args:
        filepath: Path to trajectory file

    Returns:
        TrajectoryResult2D object
    """
    from .trajectory import TrajectoryResult2D

    # Read header to detect units
    with open(filepath) as f:
        header = f.readline()

    data = _loadtxt_fortran(filepath)

    time = data[:, 0]
    x1 = data[:, 1]
    x2 = data[:, 2]
    v1 = data[:, 3] if data.shape[1] > 3 else np.zeros_like(x1)
    v2 = data[:, 4] if data.shape[1] > 4 else np.zeros_like(x1)

    # Detect and convert units
    if "fs" in header.lower() or "angstrom" in header.lower():
        # User units: fs and Angstrom
        time = time * 1e-15  # fs to s
        x1 = x1 * 1e-10  # Angstrom to m
        x2 = x2 * 1e-10
        v1 = v1 * 1e5  # Angstrom/fs to m/s
        v2 = v2 * 1e5
    # else assume SI units

    # Compute acceleration (not available from file, set to zero)
    a1 = np.zeros_like(x1)
    a2 = np.zeros_like(x2)

    return TrajectoryResult2D(
        time=time,
        x1=x1,
        x2=x2,
        v1=v1,
        v2=v2,
        a1=a1,
        a2=a2,
        x1_0=x1[0],
        x2_0=x2[0],
        p1_0=0.0,  # Unknown from file
        p2_0=0.0,
    )


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


def read_dipole_file_2d(
    filepath: Path,
    position_units: str = "angstrom",
    index_order: str = "C",
    dipole_components: int = 3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D dipole file with unit conversion.

    Expected format for 3 components:
        x1  x2  d_x  d_y  d_z
        ...

    Expected format for 1 component:
        x1  x2  d
        ...

    Where x1 and x2 are coordinates and d (or d_x, d_y, d_z) are dipole components.

    Args:
        filepath: Path to dipole file
        position_units: "angstrom" or "bohr" for coordinates
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
        dipole_components: Number of dipole components in file (1 or 3).
            If 1, the file contains |d|^2 (dipole squared). The square root
            is taken and placed in the z-component (index 2), with x and y
            components set to 0.

    Returns:
        x1: Unique x1 grid values (SI: meters)
        x2: Unique x2 grid values (SI: meters)
        d: Dipole on 2D grid in atomic units, shape (n_x1, n_x2, 3)

    Raises:
        ValueError: If position_units or index_order are invalid
    """
    # Read raw data
    x1_raw, x2_raw, d_raw = read_dipole_file_2d_raw(
        filepath, index_order, dipole_components
    )

    # Convert position units
    if position_units.lower() == "angstrom":
        x1 = x1_raw * 1e-10
        x2 = x2_raw * 1e-10
    elif position_units.lower() == "bohr":
        x1 = x1_raw * CONST.bohr
        x2 = x2_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    # Dipole is kept in atomic units (as in Fortran code)
    return x1, x2, d_raw


def read_dipole_file_2d_raw(
    filepath: Path,
    index_order: str = "C",
    dipole_components: int = 3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read 2D dipole file without unit conversion.

    Args:
        filepath: Path to dipole file
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
        dipole_components: Number of dipole components in file (1 or 3).
            If 1, the file contains |d|^2 (dipole squared). The square root
            is taken and placed in the z-component (index 2), with x and y
            components set to 0.

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        d: Dipole on 2D grid in original units, shape (n_x1, n_x2, 3)

    Raises:
        ValueError: If index_order is not "C" or "F", or if data layout
                    doesn't match the specified ordering
    """
    if dipole_components not in (1, 3):
        raise ValueError(f"dipole_components must be 1 or 3, got {dipole_components}")

    data = _loadtxt_fortran(filepath)

    min_cols = 2 + dipole_components  # x1, x2, plus dipole columns

    if data.ndim == 1:
        raise ValueError(
            f"Dipole file must have at least {min_cols} columns, got 1D array"
        )

    if data.shape[1] < min_cols:
        raise ValueError(
            f"Dipole file must have at least {min_cols} columns "
            f"(x1, x2, {'d' if dipole_components == 1 else 'd_x, d_y, d_z'}), "
            f"got {data.shape[1]}"
        )

    x1_raw = data[:, 0]
    x2_raw = data[:, 1]

    if dipole_components == 3:
        d_raw = data[:, 2:5]  # d_x, d_y, d_z columns
    else:
        # Single component: file contains |d|^2, take sqrt and place in z-component
        d_squared = data[:, 2]
        d_raw = np.zeros((len(d_squared), 3))
        d_raw[:, 2] = np.sqrt(d_squared)  # z-component = sqrt(|d|^2)

    # Get unique grid values
    x1_unique = np.unique(x1_raw)
    x2_unique = np.unique(x2_raw)
    n_x1 = len(x1_unique)
    n_x2 = len(x2_unique)

    # Validate that the data ordering matches the specified index_order
    # This must be done BEFORE reshaping
    _validate_index_ordering(x1_raw, x2_raw, x1_unique, x2_unique, index_order)

    # Reshape according to index order
    # d_2d has shape (n_x1, n_x2, 3)
    if index_order.upper() == "C":
        # x2 varies fastest (C-style, row-major)
        d_2d = d_raw.reshape((n_x1, n_x2, 3), order="C")
    elif index_order.upper() == "F":
        # x1 varies fastest (Fortran-style, column-major)
        # Need to handle the 3 components separately
        d_2d = np.zeros((n_x1, n_x2, 3))
        for i in range(3):
            d_2d[:, :, i] = d_raw[:, i].reshape((n_x1, n_x2), order="F")
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

    return x1_unique, x2_unique, d_2d


def write_dipole_file_2d(
    filepath: Path,
    x1: np.ndarray,
    x2: np.ndarray,
    d: np.ndarray,
    position_units: str = "angstrom",
    index_order: str = "C",
    dipole_components: int = 3,
) -> None:
    """Write 2D dipole to file.

    Args:
        filepath: Output file path
        x1: Grid points for x1 (SI: meters)
        x2: Grid points for x2 (SI: meters)
        d: Dipole on 2D grid in atomic units, shape (len(x1), len(x2), 3)
        position_units: "angstrom" or "bohr" for output coordinates
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
        dipole_components: Number of dipole components to write (1 or 3).
            If 1, writes |d|^2 (the z-component squared, index 2).

    Raises:
        ValueError: If position_units, index_order, or dipole_components are invalid
    """
    if dipole_components not in (1, 3):
        raise ValueError(f"dipole_components must be 1 or 3, got {dipole_components}")

    # Convert from SI to output units
    if position_units.lower() == "angstrom":
        x1_out = x1 * 1e10
        x2_out = x2 * 1e10
    elif position_units.lower() == "bohr":
        x1_out = x1 / CONST.bohr
        x2_out = x2 / CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    # Create output array with proper ordering
    rows = []
    if index_order.upper() == "C":
        # x2 varies fastest
        for i in range(len(x1_out)):
            for j in range(len(x2_out)):
                if dipole_components == 3:
                    rows.append([x1_out[i], x2_out[j], d[i, j, 0], d[i, j, 1], d[i, j, 2]])
                else:
                    # Write |d|^2 (z-component squared)
                    rows.append([x1_out[i], x2_out[j], d[i, j, 2] ** 2])
    elif index_order.upper() == "F":
        # x1 varies fastest
        for j in range(len(x2_out)):
            for i in range(len(x1_out)):
                if dipole_components == 3:
                    rows.append([x1_out[i], x2_out[j], d[i, j, 0], d[i, j, 1], d[i, j, 2]])
                else:
                    # Write |d|^2 (z-component squared)
                    rows.append([x1_out[i], x2_out[j], d[i, j, 2] ** 2])
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

    data = np.array(rows)
    if dipole_components == 3:
        header = f"# x1({position_units}) x2({position_units}) d_x d_y d_z (atomic units)"
    else:
        header = f"# x1({position_units}) x2({position_units}) |d|^2 (atomic units)"
    np.savetxt(filepath, data, header=header, fmt="%16.8E")


def write_pes_file_2d(
    filepath: Path,
    x1: np.ndarray,
    x2: np.ndarray,
    E: np.ndarray,
    position_units: str = "angstrom",
    energy_units: str = "hartree",
    index_order: str = "C",
) -> None:
    """Write 2D PES to file.

    Args:
        filepath: Output file path
        x1: Grid points for x1 (SI: meters)
        x2: Grid points for x2 (SI: meters)
        E: Energy on 2D grid (SI: Joules), shape (len(x1), len(x2))
        position_units: "angstrom" or "bohr" for output coordinates
        energy_units: "hartree" or "ev" for output energy
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering

    Raises:
        ValueError: If position_units, energy_units, or index_order are invalid
    """
    # Convert from SI to output units
    if position_units.lower() == "angstrom":
        x1_out = x1 * 1e10
        x2_out = x2 * 1e10
    elif position_units.lower() == "bohr":
        x1_out = x1 / CONST.bohr
        x2_out = x2 / CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    if energy_units.lower() == "hartree":
        E_out = E / CONST.hartree
    elif energy_units.lower() == "ev":
        E_out = E / CONST.eV
    else:
        raise ValueError(f"Unknown energy units: {energy_units}")

    # Create output array with proper ordering
    rows = []
    if index_order.upper() == "C":
        # x2 varies fastest
        for i in range(len(x1_out)):
            for j in range(len(x2_out)):
                rows.append([x1_out[i], x2_out[j], E_out[i, j]])
    elif index_order.upper() == "F":
        # x1 varies fastest
        for j in range(len(x2_out)):
            for i in range(len(x1_out)):
                rows.append([x1_out[i], x2_out[j], E_out[i, j]])
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

    data = np.array(rows)
    header = f"# x1({position_units}) x2({position_units}) E({energy_units})"
    np.savetxt(filepath, data, header=header, fmt="%16.8E")
