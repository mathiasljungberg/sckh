"""Standalone SCKH_PES spectrum example — YAML-driven.

Loads the full workflow configuration from ``sckh_pes.yaml`` and
produces the same outputs as ``run_python.py``.  The pipeline is kept
in two explicit stages:

    dynamics_1d  →  SCKH trajectories (+ E_mean, D_ni)
    sckh         →  spectrum

Run directly as::

    python run_yaml.py
"""

import shutil
from pathlib import Path

from dynamics_1d import compute_sckh_trajectories_1d
from sckh import SCKHSpectrumCalculator, load_full_config


HERE = Path(__file__).resolve().parent
RUNDIR = HERE / "rundir"


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    config = load_full_config(HERE / "sckh_pes.yaml")

    # Stage 1: dynamics + surface interpolation → SCKH trajectories.
    inputs = compute_sckh_trajectories_1d(config)

    # Stage 2: SCKH spectrum from those trajectories.
    calc = SCKHSpectrumCalculator(config)
    result = calc.compute_spectrum(
        inputs.trajectories,
        E_mean=inputs.E_mean,
        D_ni=inputs.D_ni,
    )

    calc.write_spectrum_result(RUNDIR / "testout_sigma", result)


if __name__ == "__main__":
    main()
