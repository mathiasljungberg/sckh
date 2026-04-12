"""Standalone SCKH spectrum example — pure Python API.

Reads pre-computed SCKH trajectory files from ``./input/`` and writes the
computed spectrum to ``./rundir/``.  Run directly as::

    python run_python.py

This file is both a test fixture and a copy-pasteable usage example of
the high-level ``sckh`` API.
"""

import shutil
from pathlib import Path

from sckh import (
    compute_spectrum_from_sckh,
    read_sckh_trajectories,
    write_spectrum_result,
)


HERE = Path(__file__).resolve().parent
INPUT = HERE / "input"
RUNDIR = HERE / "rundir"

GAMMA_FWHM = 0.18  # eV

TRAJECTORY_FILES = [
    INPUT / "mdpentamer_1_1_1_1.combined",
    INPUT / "mdpentamer_1_1_1_2.combined",
    INPUT / "mdpentamer_1_1_1_3.combined",
]


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    trajs = read_sckh_trajectories(TRAJECTORY_FILES)
    result = compute_spectrum_from_sckh(trajs, gamma_hwhm=GAMMA_FWHM / 2)
    write_spectrum_result(RUNDIR / "pentamer_XES_sigma", result)


if __name__ == "__main__":
    main()
