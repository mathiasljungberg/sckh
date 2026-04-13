"""Standalone SCKH_PES spectrum example — YAML-driven.

Reads the full workflow configuration from ``sckh_pes.yaml`` and
produces the same outputs as ``run_python.py``.  Demonstrates the
``load_full_config()`` entry point for a dynamics + interpolation +
spectrum pipeline.

Run directly as::

    python run_yaml.py
"""

import shutil
from pathlib import Path

from dynamics_1d import compute_spectrum_1d
from sckh import load_full_config, write_spectrum_result


HERE = Path(__file__).resolve().parent
RUNDIR = HERE / "rundir"


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    config = load_full_config(HERE / "sckh_pes.yaml")
    result = compute_spectrum_1d(config)
    write_spectrum_result(RUNDIR / "testout_sigma", result)


if __name__ == "__main__":
    main()
