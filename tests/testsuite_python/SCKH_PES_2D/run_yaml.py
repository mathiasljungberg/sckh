"""Standalone SCKH_PES_2D spectrum example — YAML-driven.

Loads the full workflow configuration from ``sckh_pes_2d.yaml`` and produces
the same outputs as ``run_python.py``.

Run directly as::

    python run_yaml.py
"""

import shutil
from pathlib import Path

from dynamics_2d import DynamicsRunner2D, SCKHSpectrumCalculator, SurfaceInterpolator2D, load_config


HERE = Path(__file__).resolve().parent
RUNDIR = HERE / "rundir"


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    config = load_config(HERE / "sckh_pes_2d.yaml")

    runner = DynamicsRunner2D(config)
    dynamics_result = runner.run(verbose=False)
    runner.save_results(dynamics_result, RUNDIR)

    interp = SurfaceInterpolator2D(config)
    interp.load_surfaces()

    sckh_trajs = interp.trajectories_to_sckh(dynamics_result.trajectories)
    e_mean = interp.compute_mean_transition_energy(dynamics_result.trajectories)
    d_ni = interp.compute_D_ni()

    calc = SCKHSpectrumCalculator(config)
    spectrum_result = calc.compute_spectrum(
        sckh_trajs,
        E_mean=e_mean,
        D_ni=d_ni,
    )
    calc.write_spectrum_result(
        RUNDIR / f"{config.dynamics2d.outfile}_sigma",
        spectrum_result,
    )


if __name__ == "__main__":
    main()
