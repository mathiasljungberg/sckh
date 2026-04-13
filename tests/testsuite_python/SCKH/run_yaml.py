"""Standalone SCKH spectrum example — YAML-driven.

Reads the workflow configuration from ``sckh.yaml`` and produces the same
outputs as ``run_python.py``.  Demonstrates the ``load_full_config()``
entry point for pure-SCKH workflows.  Run directly as::

    python run_yaml.py
"""

import shutil
from pathlib import Path

from sckh import (
    SCKHSpectrumCalculator,
    SpectrumConfig,
    load_full_config,
    read_sckh_trajectories,
)


HERE = Path(__file__).resolve().parent
RUNDIR = HERE / "rundir"


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    cfg = load_full_config(HERE / "sckh.yaml").sckh

    trajs = read_sckh_trajectories(cfg.trajectory_files)

    calc = SCKHSpectrumCalculator(SpectrumConfig(gamma_fwhm=cfg.gamma_fwhm))
    result = calc.compute_spectrum(trajs)
    calc.write_spectrum_result(RUNDIR / cfg.outfile, result)


if __name__ == "__main__":
    main()
