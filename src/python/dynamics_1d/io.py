"""File I/O for PES, trajectory, dipole, and spectrum data."""

from pathlib import Path
from typing import Tuple, Optional, TYPE_CHECKING

import numpy as np

from .constants import CONST

if TYPE_CHECKING:
    from .trajectory import TrajectoryResult


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
        text = f.read().replace('D', 'E').replace('d', 'e')

    return np.loadtxt(StringIO(text))


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
    data = _loadtxt_fortran(filepath)
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
    data = _loadtxt_fortran(filepath)
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


def read_dipole_file(
    filepath: Path,
    units: str = "angstrom",
) -> Tuple[np.ndarray, np.ndarray]:
    """Read dipole file with format: x d_x d_y d_z.

    Args:
        filepath: Path to dipole file
        units: "angstrom" or "bohr" for x coordinate

    Returns:
        x: Position array in meters (SI)
        d: Dipole array (n_points, 3) in atomic units
    """
    data = _loadtxt_fortran(filepath)

    if data.ndim == 1:
        raise ValueError(
            f"Dipole file must have at least 4 columns (x, d_x, d_y, d_z), "
            f"got 1D array"
        )

    if data.shape[1] < 4:
        raise ValueError(
            f"Dipole file must have at least 4 columns (x, d_x, d_y, d_z), "
            f"got {data.shape[1]}"
        )

    x_raw = data[:, 0]
    d = data[:, 1:4]  # d_x, d_y, d_z columns

    # Convert x to SI units
    if units.lower() == "angstrom":
        x = x_raw * 1e-10  # Angstrom to meters
    elif units.lower() == "bohr":
        x = x_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown units: {units}")

    # Dipole is kept in atomic units (as in Fortran code)
    return x, d


def read_trajectory_file(filepath: Path) -> "TrajectoryResult":
    """Read trajectory from file.

    Expects format: time(s) x(m) v(m/s) or time(fs) x(Angstrom) v(Angstrom/fs).
    Auto-detects units based on header or magnitude.

    Args:
        filepath: Path to trajectory file

    Returns:
        TrajectoryResult object
    """
    from .trajectory import TrajectoryResult

    # Read header to detect units
    with open(filepath) as f:
        header = f.readline()

    data = _loadtxt_fortran(filepath)

    time = data[:, 0]
    x = data[:, 1]
    v = data[:, 2] if data.shape[1] > 2 else np.zeros_like(x)

    # Detect and convert units
    if "fs" in header.lower() or "angstrom" in header.lower():
        # User units: fs and Angstrom
        time = time * 1e-15  # fs to s
        x = x * 1e-10  # Angstrom to m
        v = v * 1e5  # Angstrom/fs to m/s
    # else assume SI units

    # Compute acceleration (not available from file, set to zero)
    a = np.zeros_like(x)

    return TrajectoryResult(
        time=time,
        x=x,
        v=v,
        a=a,
        x0=x[0],
        p0=0.0,  # Unknown from file
    )


def write_spectrum(
    filepath: Path,
    omega: np.ndarray,
    sigma: np.ndarray,
    header: Optional[str] = None,
) -> None:
    """Write spectrum to file in Fortran-compatible format.

    Args:
        filepath: Output file path
        omega: Frequency array in eV
        sigma: Cross-section array (same length as omega)
        header: Optional header comment
    """
    if header is None:
        header = "omega(eV)  sigma"

    data = np.column_stack([omega, sigma])
    np.savetxt(filepath, data, header=header, fmt="%16.6E", comments="# ")


def write_spectrum_per_final(
    filepath_base: Path,
    omega: np.ndarray,
    sigma_f: np.ndarray,
) -> None:
    """Write per-final-state spectra to files.

    Creates files: {filepath_base}_final_1.dat, {filepath_base}_final_2.dat, etc.

    Args:
        filepath_base: Base path for output files (without extension)
        omega: Frequency array in eV
        sigma_f: Cross-section array (n_final, n_omega)
    """
    n_final = sigma_f.shape[0]

    for j in range(n_final):
        filepath = Path(f"{filepath_base}_final_{j + 1}.dat")
        write_spectrum(
            filepath,
            omega,
            sigma_f[j],
            header=f"omega(eV)  sigma_f (final state {j + 1})",
        )
