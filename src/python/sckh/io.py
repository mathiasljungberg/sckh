"""I/O for SCKH trajectories and spectrum data."""

from pathlib import Path
from typing import Optional, TYPE_CHECKING

import numpy as np

from dynamics_1d.constants import CONST

if TYPE_CHECKING:
    from .trajectory import SCKHTrajectory


def _parse_fortran_float(s: str) -> float:
    """Parse a Fortran-format float (handles D/d exponent notation)."""
    return float(s.replace('D', 'E').replace('d', 'e'))



def read_sckh_trajectory(filepath: Path) -> "SCKHTrajectory":
    """Read SCKH trajectory file in Fortran format.

    File format (per timestep):
        time (fs), E_gs (Hartree), E_n (Hartree), E_IP1s (Hartree),
        "ntrans" nfinal, then nfinal lines of E_trans D_x D_y D_z.

    Applies the Fortran derivation:
        E_f = E_gs - E_trans + E_IP1s * (eV/Hartree)  (all in Hartree)
    Then converts E_n and E_f to eV.

    Args:
        filepath: Path to SCKH trajectory file

    Returns:
        SCKHTrajectory with converted units
    """
    from .trajectory import SCKHTrajectory

    with open(filepath) as f:
        lines = f.readlines()

    # First pass: determine nfinal from the first "ntrans" line
    # Lines per block before transitions: time, E_gs, E_n, E_IP1s, ntrans_line = 5
    ntrans_line = lines[4].split()
    nfinal = int(ntrans_line[1])

    # Lines per timestep: 5 + nfinal
    lines_per_step = 5 + nfinal
    ntsteps = len(lines) // lines_per_step

    # Allocate arrays
    time = np.zeros(ntsteps)
    E_gs = np.zeros(ntsteps)
    E_n_raw = np.zeros(ntsteps)
    E_IP1s = np.zeros(ntsteps)
    E_trans = np.zeros((nfinal, ntsteps))
    D_fn = np.zeros((nfinal, ntsteps, 3))

    # Parse each timestep
    for i in range(ntsteps):
        base = i * lines_per_step
        time[i] = _parse_fortran_float(lines[base].strip())
        E_gs[i] = _parse_fortran_float(lines[base + 1].strip())
        E_n_raw[i] = _parse_fortran_float(lines[base + 2].strip())
        E_IP1s[i] = _parse_fortran_float(lines[base + 3].strip())

        # Parse "ntrans nfinal" line and verify
        parts = lines[base + 4].split()
        ntrans = int(parts[1])
        if ntrans != nfinal:
            raise ValueError(
                f"Inconsistent ntrans at step {i}: {ntrans} != {nfinal}"
            )

        # Parse transition data
        for j in range(nfinal):
            vals = lines[base + 5 + j].split()
            vals = [_parse_fortran_float(v) for v in vals]
            E_trans[j, i] = vals[0]
            D_fn[j, i, 0] = vals[1]
            D_fn[j, i, 1] = vals[2]
            D_fn[j, i, 2] = vals[3]

    # Compute E_f following Fortran: E_f = E_gs - E_trans + E_IP1s * (eV/Hartree)
    # All quantities in Hartree at this point
    eV_over_hartree = CONST.eV / CONST.hartree
    E_f_hartree = np.zeros((nfinal, ntsteps))
    for i in range(ntsteps):
        E_f_hartree[:, i] = E_gs[i] - E_trans[:, i] + E_IP1s[i] * eV_over_hartree

    # Convert to output units
    time_s = time * 1e-15  # fs to seconds
    E_n_eV = E_n_raw * CONST.hartree2eV  # Hartree to eV
    E_f_eV = E_f_hartree * CONST.hartree2eV  # Hartree to eV

    return SCKHTrajectory(
        time=time_s,
        E_gs=E_gs,  # Keep in Hartree
        E_n=E_n_eV,
        E_f=E_f_eV,
        D_fn=D_fn,
        E_IP1s=E_IP1s,
    )


def write_sckh_trajectory(
    filepath: Path,
    sckh_traj: "SCKHTrajectory",
) -> None:
    """Write SCKH trajectory file in Fortran format.

    Converts from internal units (eV, seconds) back to file units
    (Hartree, femtoseconds).

    Args:
        filepath: Output file path
        sckh_traj: SCKHTrajectory to write
    """
    ntsteps = sckh_traj.ntsteps
    nfinal = sckh_traj.n_final

    # Convert units for output
    time_fs = sckh_traj.time * 1e15  # seconds to fs
    E_n_hartree = sckh_traj.E_n / CONST.hartree2eV  # eV to Hartree
    E_f_hartree = sckh_traj.E_f / CONST.hartree2eV  # eV to Hartree

    # Compute E_trans = E_gs - E_f_hartree + E_IP1s * (eV/Hartree)
    # Inverting: E_f = E_gs - E_trans + E_IP1s * (eV/Hartree)
    # => E_trans = E_gs - E_f + E_IP1s * (eV/Hartree)
    eV_over_hartree = CONST.eV / CONST.hartree
    E_trans = np.zeros((nfinal, ntsteps))
    for i in range(ntsteps):
        E_trans[:, i] = (
            sckh_traj.E_gs[i]
            - E_f_hartree[:, i]
            + sckh_traj.E_IP1s[i] * eV_over_hartree
        )

    with open(filepath, 'w') as f:
        for i in range(ntsteps):
            f.write(f"  {time_fs[i]:20.10E}\n")
            f.write(f"  {sckh_traj.E_gs[i]:20.10E}\n")
            f.write(f"  {E_n_hartree[i]:20.10E}\n")
            f.write(f"  {sckh_traj.E_IP1s[i]:20.10E}\n")
            f.write(f"ntrans {nfinal:5d}\n")
            for j in range(nfinal):
                f.write(
                    f"{E_trans[j, i]:15.8f}"
                    f"{sckh_traj.D_fn[j, i, 0]:17.8E}"
                    f"{sckh_traj.D_fn[j, i, 1]:17.8E}"
                    f"{sckh_traj.D_fn[j, i, 2]:17.8E}\n"
                )


def read_sckh_trajectory_list(
    filepath: Path,
) -> list:
    """Read SCKH trajectory file list (trajectories.dat format).

    Each line has: filename weight
    Lines starting with '#' or '!' are ignored.

    Args:
        filepath: Path to trajectory list file

    Returns:
        List of (filepath, weight) tuples. Paths are relative to the
        directory containing the list file.
    """
    base_dir = Path(filepath).parent
    result = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            parts = line.split()
            traj_path = base_dir / parts[0]
            weight = float(parts[1]) if len(parts) > 1 else 1.0
            result.append((traj_path, weight))
    return result


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
