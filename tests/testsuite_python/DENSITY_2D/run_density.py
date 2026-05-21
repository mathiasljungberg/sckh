"""Trajectory density on the dynamics grid — minimal example.

Reads the trajectories committed under ``input/``, deposits each
trajectory point onto the interpolation grid with a Gaussian of given
FWHM (summed over trajectories, normalized to one), and writes one
density file per selected time step.

Run directly as::

    python run_density.py

Outputs land in ``rundir/density_step_NNNNNN.dat``.
"""

import shutil
from pathlib import Path

import numpy as np

from dynamics_2d import (
    read_trajectories_2d,
    trajectories_density_timeseries,
    write_density_timeseries,
)


HERE = Path(__file__).resolve().parent
TRAJ_DIR = HERE / "input"
TRAJ_BASENAME = "dynamics_2d_out"
RUNDIR = HERE / "rundir"

# Grid that the SCKH_PES_2D trajectories were propagated on (Angstrom).
GRID_START_ANG = -0.35
GRID_DX_ANG = 0.05
GRID_NPOINTS = 34

# 5 snapshots across the trajectory.
SELECTED_STEPS = [0, 256, 512, 768, 1023]

# FWHM of the broadening Gaussian, in Angstrom. The grid spacing is
# 0.05 Å, so 0.1 Å (~ 2 grid points) is a reasonable default.
FWHM_ANG = 0.1


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    trajectories = read_trajectories_2d(TRAJ_DIR, TRAJ_BASENAME)
    if not trajectories:
        raise FileNotFoundError(
            f"No trajectory files {TRAJ_BASENAME}_traj_*.dat found in {TRAJ_DIR}."
        )

    # SI grids (meters) matching the SCKH_PES_2D setup.
    x1_grid = (GRID_START_ANG + np.arange(GRID_NPOINTS) * GRID_DX_ANG) * 1e-10
    x2_grid = x1_grid.copy()

    density_result = trajectories_density_timeseries(
        x1_grid,
        x2_grid,
        trajectories,
        fwhm=FWHM_ANG,
        step_indices=SELECTED_STEPS,
        position_units="angstrom",
    )

    paths = write_density_timeseries(
        RUNDIR,
        density_result,
        basename="density",
        position_units="angstrom",
    )
    print(f"Read {len(trajectories)} trajectories from {TRAJ_DIR}")
    print(f"Wrote {len(paths)} density files to {RUNDIR}")


if __name__ == "__main__":
    main()
