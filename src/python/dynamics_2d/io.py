"""File I/O for 2D PES, dipole surfaces, trajectory, and ground state data."""

from pathlib import Path
from typing import List, Sequence, Tuple, Optional, TYPE_CHECKING

import numpy as np

from dynamics_1d.constants import CONST

if TYPE_CHECKING:
    from .density import TrajectoryDensityResult2D
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


def infer_index_order_2d(
    x1_raw: np.ndarray,
    x2_raw: np.ndarray,
    x1_unique: np.ndarray,
    x2_unique: np.ndarray,
) -> str:
    """Infer C/Fortran ordering from the coordinate sequence."""
    X1_grid, X2_grid = np.meshgrid(x1_unique, x2_unique, indexing="ij")

    x1_expected_c = X1_grid.flatten(order="C")
    x2_expected_c = X2_grid.flatten(order="C")
    x1_expected_f = X1_grid.flatten(order="F")
    x2_expected_f = X2_grid.flatten(order="F")

    c_match = np.allclose(x1_raw, x1_expected_c) and np.allclose(x2_raw, x2_expected_c)
    f_match = np.allclose(x1_raw, x1_expected_f) and np.allclose(x2_raw, x2_expected_f)

    if c_match and not f_match:
        return "C"
    if f_match and not c_match:
        return "F"
    if c_match and f_match:
        # Degenerate grids (single row or column) can be valid in both orders.
        return "C"

    raise ValueError(
        "Could not infer index ordering: coordinate sequence matches neither "
        "C order (x2 fast) nor F order (x1 fast)."
    )


def read_surface_file_2d_raw(
    filepath: Path,
    value_columns: Sequence[int],
    index_order: str = "C",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read generic 2D surface file without unit conversion.

    Expected format:
        x1  x2  value_1 [value_2 ...]

    Args:
        filepath: Path to surface file
        value_columns: Zero-based columns to read as values
        index_order: "C", "F", or "auto"

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        values: Values on 2D grid, shape (n_values, n_x1, n_x2)
        resolved_order: Inferred/resolved order ("C" or "F")
    """
    data = np.atleast_2d(_loadtxt_fortran(filepath))
    if data.shape[1] < 3:
        raise ValueError(
            f"Surface file must have at least 3 columns (x1 x2 value), got {data.shape[1]}"
        )
    if len(value_columns) == 0:
        raise ValueError("value_columns must contain at least one column index")

    x1_raw = data[:, 0]
    x2_raw = data[:, 1]

    x1_unique = np.unique(x1_raw)
    x2_unique = np.unique(x2_raw)
    n_x1 = len(x1_unique)
    n_x2 = len(x2_unique)

    expected_points = n_x1 * n_x2
    if len(x1_raw) != expected_points:
        raise ValueError(
            f"Surface file is not a full rectangular grid: got {len(x1_raw)} points, "
            f"expected {expected_points} = {n_x1}*{n_x2}"
        )

    max_col = data.shape[1] - 1
    for col in value_columns:
        if col < 2 or col > max_col:
            raise ValueError(
                f"Requested value column {col} out of range for file with "
                f"{data.shape[1]} columns"
            )

    index_order_upper = index_order.upper()
    if index_order_upper == "AUTO":
        resolved_order = infer_index_order_2d(x1_raw, x2_raw, x1_unique, x2_unique)
    elif index_order_upper in {"C", "F"}:
        _validate_index_ordering(x1_raw, x2_raw, x1_unique, x2_unique, index_order_upper)
        resolved_order = index_order_upper
    else:
        raise ValueError(f"index_order must be 'C', 'F', or 'auto', got {index_order}")

    values = np.stack(
        [
            data[:, col].reshape((n_x1, n_x2), order=resolved_order)
            for col in value_columns
        ],
        axis=0,
    )
    return x1_unique, x2_unique, values, resolved_order


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
        index_order: "C", "F", or "auto" for data ordering

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
        index_order: "C", "F", or "auto" for data ordering

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        E: Energy on 2D grid in original units, shape (n_x1, n_x2)

    Raises:
        ValueError: If index_order is invalid, or if data layout doesn't match
                    the specified ordering
    """
    x1_unique, x2_unique, values, _ = read_surface_file_2d_raw(
        filepath,
        value_columns=[2],
        index_order=index_order,
    )
    return x1_unique, x2_unique, values[0]


def read_dipole_file_2d_raw(
    filepath: Path,
    index_order: str = "C",
    dipole_components: int = 3,
    return_resolved_order: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read 2D dipole file without unit conversion.

    Args:
        filepath: Path to dipole file
        index_order: "C", "F", or "auto" for data ordering
        dipole_components: Number of dipole components in file (1 or 3)
        return_resolved_order: Return inferred/resolved ordering if True

    Returns:
        x1: Unique x1 grid values in original units
        x2: Unique x2 grid values in original units
        d: Dipole on 2D grid in original units, shape (n_x1, n_x2, 3)
        resolved_order: Optional inferred/resolved order ("C" or "F")
    """
    if dipole_components not in (1, 3):
        raise ValueError(f"dipole_components must be 1 or 3, got {dipole_components}")

    if dipole_components == 3:
        x1_unique, x2_unique, values, resolved_order = read_surface_file_2d_raw(
            filepath,
            value_columns=[2, 3, 4],
            index_order=index_order,
        )
        d_2d = np.transpose(values, (1, 2, 0))
    else:
        x1_unique, x2_unique, values, resolved_order = read_surface_file_2d_raw(
            filepath,
            value_columns=[2],
            index_order=index_order,
        )
        d_2d = np.zeros((len(x1_unique), len(x2_unique), 3))
        d_2d[:, :, 2] = np.sqrt(values[0])

    if return_resolved_order:
        return x1_unique, x2_unique, d_2d, resolved_order
    return x1_unique, x2_unique, d_2d


def read_dipole_file_2d(
    filepath: Path,
    position_units: str = "angstrom",
    index_order: str = "C",
    dipole_components: int = 3,
    return_resolved_order: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read 2D dipole file with position unit conversion."""
    raw = read_dipole_file_2d_raw(
        filepath,
        index_order=index_order,
        dipole_components=dipole_components,
        return_resolved_order=return_resolved_order,
    )
    if return_resolved_order:
        x1_raw, x2_raw, d_raw, resolved_order = raw
    else:
        x1_raw, x2_raw, d_raw = raw

    if position_units.lower() == "angstrom":
        x1 = x1_raw * 1e-10
        x2 = x2_raw * 1e-10
    elif position_units.lower() == "bohr":
        x1 = x1_raw * CONST.bohr
        x2 = x2_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    if return_resolved_order:
        return x1, x2, d_raw, resolved_order
    return x1, x2, d_raw


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


def read_trajectories_2d(
    directory: Path,
    basename: str,
) -> List["TrajectoryResult2D"]:
    """Read all trajectory files written by ``DynamicsRunner2D.save_results``.

    Args:
        directory: Directory containing trajectory files.
        basename: File-name prefix (matches ``DynamicsConfig2D.outfile``).
            Files are expected to be named ``{basename}_traj_*.dat``.

    Returns:
        List of TrajectoryResult2D in sorted filename order. Empty list if
        no matching files are found.
    """
    directory = Path(directory)
    files = sorted(directory.glob(f"{basename}_traj_*.dat"))
    return [read_trajectory_2d(f) for f in files]


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


def _write_scalar_field_2d(
    filepath: Path,
    x1: np.ndarray,
    x2: np.ndarray,
    values: np.ndarray,
    *,
    header: str,
    position_units: str,
    index_order: str,
) -> None:
    """Write a 2D scalar field as `x1  x2  value` rows.

    Args:
        filepath: Output file path
        x1: Grid points for x1 (SI: meters)
        x2: Grid points for x2 (SI: meters)
        values: Field on the 2D grid, shape (len(x1), len(x2)), already in
            the desired output units
        header: Full header text (may be multi-line; np.savetxt prepends '# ')
        position_units: "angstrom" or "bohr" for output coordinates
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
    """
    if position_units.lower() == "angstrom":
        x1_out = x1 * 1e10
        x2_out = x2 * 1e10
    elif position_units.lower() == "bohr":
        x1_out = x1 / CONST.bohr
        x2_out = x2 / CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    rows = []
    if index_order.upper() == "C":
        for i in range(len(x1_out)):
            for j in range(len(x2_out)):
                rows.append([x1_out[i], x2_out[j], values[i, j]])
    elif index_order.upper() == "F":
        for j in range(len(x2_out)):
            for i in range(len(x1_out)):
                rows.append([x1_out[i], x2_out[j], values[i, j]])
    else:
        raise ValueError(f"index_order must be 'C' or 'F', got {index_order}")

    data = np.array(rows)
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
    if energy_units.lower() == "hartree":
        E_out = E / CONST.hartree
    elif energy_units.lower() == "ev":
        E_out = E / CONST.eV
    else:
        raise ValueError(f"Unknown energy units: {energy_units}")

    header = f"# x1({position_units}) x2({position_units}) E({energy_units})"
    _write_scalar_field_2d(
        filepath,
        x1,
        x2,
        E_out,
        header=header,
        position_units=position_units,
        index_order=index_order,
    )


def write_density_2d(
    filepath: Path,
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
    density: np.ndarray,
    *,
    time: Optional[float] = None,
    step_index: Optional[int] = None,
    fwhm_SI: Optional[float] = None,
    position_units: str = "angstrom",
    index_order: str = "C",
) -> None:
    """Write a 2D trajectory density (SI: 1/m^2) to file.

    Args:
        filepath: Output file path
        x1_grid: Grid points for x1 (SI: meters)
        x2_grid: Grid points for x2 (SI: meters)
        density: Density on the 2D grid (SI: 1/m^2), shape (len(x1), len(x2))
        time: Time of the snapshot (SI: seconds); written into header as fs
        step_index: Step index for the snapshot; written into header
        fwhm_SI: FWHM used for Gaussian broadening (SI: meters); written into
            header in `position_units`
        position_units: "angstrom" or "bohr" for output coordinates
        index_order: "C" (x2 fast) or "F" (x1 fast) for data ordering
    """
    pu = position_units.lower()
    if pu == "angstrom":
        rho_out = density * 1e-20  # 1/m^2 -> 1/Angstrom^2
        fwhm_user = fwhm_SI * 1e10 if fwhm_SI is not None else None
    elif pu == "bohr":
        rho_out = density * (CONST.bohr ** 2)  # 1/m^2 -> 1/Bohr^2
        fwhm_user = fwhm_SI / CONST.bohr if fwhm_SI is not None else None
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    header_lines = []
    if time is not None:
        header_lines.append(f"time(fs) = {time * 1e15:.10E}")
    if step_index is not None:
        header_lines.append(f"step_index = {int(step_index)}")
    if fwhm_user is not None:
        header_lines.append(f"fwhm({position_units}) = {fwhm_user:.10E}")
    header_lines.append(
        f"x1({position_units}) x2({position_units}) rho(1/{position_units}^2)"
    )
    header = "\n".join(header_lines)

    _write_scalar_field_2d(
        filepath,
        x1_grid,
        x2_grid,
        rho_out,
        header=header,
        position_units=position_units,
        index_order=index_order,
    )


def write_density_timeseries(
    output_dir: Path,
    result: "TrajectoryDensityResult2D",
    *,
    basename: str = "density",
    position_units: str = "angstrom",
    index_order: str = "C",
) -> List[Path]:
    """Write one density file per time step in a TrajectoryDensityResult2D.

    Args:
        output_dir: Directory for the per-step files (created if missing).
        result: TrajectoryDensityResult2D from ensemble_density_timeseries().
        basename: Filename prefix; files are named
            f"{basename}_step_{step_index:06d}.dat".
        position_units: "angstrom" or "bohr" for output coordinates.
        index_order: "C" (x2 fast) or "F" (x1 fast).

    Returns:
        List of written file paths, one per time step (in the same order as
        result.step_indices).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []
    for i, step in enumerate(result.step_indices):
        path = output_dir / f"{basename}_step_{int(step):06d}.dat"
        write_density_2d(
            path,
            result.x1_grid,
            result.x2_grid,
            result.density[i],
            time=float(result.time[i]),
            step_index=int(step),
            fwhm_SI=result.fwhm_SI,
            position_units=position_units,
            index_order=index_order,
        )
        paths.append(path)
    return paths
