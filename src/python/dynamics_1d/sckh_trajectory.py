"""SCKH trajectory conversion from dynamics results.

The SCKHTrajectory dataclass itself lives in sckh.trajectory and is
re-exported here for backwards compatibility.
"""

from typing import List

import numpy as np

from .constants import CONST
from sckh.trajectory import SCKHTrajectory  # canonical location

# Re-export for backwards compatibility
__all__ = ["SCKHTrajectory", "sckh_trajectory_from_dynamics_1d", "sckh_trajectory_from_dynamics_2d"]


def sckh_trajectory_from_dynamics_1d(
    traj,
    pes_gs,
    pes_n,
    pes_f_list: List,
    dipoles: List,
    dipole_mode: str = "DIPOLE",
) -> SCKHTrajectory:
    """Create SCKHTrajectory from 1D dynamics trajectory and surfaces.

    Evaluates energies and dipoles at trajectory positions. Sets E_IP1s=0.

    Args:
        traj: TrajectoryResult with time and positions
        pes_gs: Ground state PES (PES1D)
        pes_n: Intermediate state PES (PES1D)
        pes_f_list: List of final state PES (PES1D)
        dipoles: List of Dipole1D objects (one per final state)
        dipole_mode: "DIPOLE" (position-dependent) or "FC" (constant=1)

    Returns:
        SCKHTrajectory with E_IP1s=0
    """
    x = traj.x
    nsteps = len(x)
    n_final = len(pes_f_list)

    # Evaluate energies at trajectory positions (in Joules from PES)
    E_gs_J = pes_gs.energy(x)
    E_n_J = pes_n.energy(x)

    E_f_J = np.zeros((n_final, nsteps))
    for i, pes in enumerate(pes_f_list):
        E_f_J[i] = pes.energy(x)

    # Evaluate dipoles
    D_fn = np.zeros((n_final, nsteps, 3))
    if dipole_mode.upper() == "FC":
        D_fn[:, :, :] = 1.0
    else:
        for i, dipole in enumerate(dipoles):
            D_fn[i] = dipole.dipole(x)

    # Convert to SCKH units
    E_gs_hartree = E_gs_J / CONST.hartree
    E_n_eV = E_n_J / CONST.eV
    E_f_eV = E_f_J / CONST.eV

    return SCKHTrajectory(
        time=traj.time.copy(),
        E_gs=E_gs_hartree,
        E_n=E_n_eV,
        E_f=E_f_eV,
        D_fn=D_fn,
        E_IP1s=np.zeros(nsteps),
    )


def sckh_trajectory_from_dynamics_2d(
    traj,
    pes_gs,
    pes_n,
    pes_f_list: List,
    dipoles: List,
    dipole_mode: str = "DIPOLE",
) -> SCKHTrajectory:
    """Create SCKHTrajectory from 2D dynamics trajectory and surfaces.

    Evaluates energies and dipoles at trajectory positions (x1, x2). Sets E_IP1s=0.

    Args:
        traj: TrajectoryResult2D with time, x1, x2 positions
        pes_gs: Ground state PES (PES2D)
        pes_n: Intermediate state PES (PES2D)
        pes_f_list: List of final state PES (PES2D)
        dipoles: List of Dipole2D objects (one per final state)
        dipole_mode: "DIPOLE" (position-dependent) or "FC" (constant=1)

    Returns:
        SCKHTrajectory with E_IP1s=0
    """
    x1 = traj.x1
    x2 = traj.x2
    nsteps = len(x1)
    n_final = len(pes_f_list)

    # Evaluate energies at trajectory positions (in Joules from PES)
    E_gs_J = pes_gs.energy(x1, x2)
    E_n_J = pes_n.energy(x1, x2)

    E_f_J = np.zeros((n_final, nsteps))
    for i, pes in enumerate(pes_f_list):
        E_f_J[i] = pes.energy(x1, x2)

    # Evaluate dipoles
    D_fn = np.zeros((n_final, nsteps, 3))
    if dipole_mode.upper() == "FC":
        D_fn[:, :, :] = 1.0
    else:
        for i, dipole in enumerate(dipoles):
            D_fn[i] = dipole.dipole(x1, x2)

    # Convert to SCKH units
    E_gs_hartree = E_gs_J / CONST.hartree
    E_n_eV = E_n_J / CONST.eV
    E_f_eV = E_f_J / CONST.eV

    return SCKHTrajectory(
        time=traj.time.copy(),
        E_gs=E_gs_hartree,
        E_n=E_n_eV,
        E_f=E_f_eV,
        D_fn=D_fn,
        E_IP1s=np.zeros(nsteps),
    )
