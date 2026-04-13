"""Standalone SCKH_PES spectrum example — pure Python API.

Builds a :class:`FullConfig` in Python, runs the 1D dynamics +
surface-interpolation stage to produce SCKH trajectories, then feeds
those into the SCKH spectrum calculator.  Reads PES and dipole
surfaces from ``./input/`` and writes the computed spectrum files to
``./rundir/``.

The two stages are kept explicit so that the script doubles as a
copy-pasteable example of both packages:

    dynamics_1d  →  SCKH trajectories (+ E_mean, D_ni)
    sckh         →  spectrum

Run directly as::

    python run_python.py
"""

import shutil
from pathlib import Path

from dynamics_1d import compute_sckh_trajectories_1d
from dynamics_1d.config import (
    DynamicsConfig,
    GridConfig,
    SamplingConfig,
    TimeConfig,
)
from dynamics_1d.spectrum_config import InterpolationConfig
from sckh import FullConfig, SCKHSpectrumCalculator, SpectrumConfig


HERE = Path(__file__).resolve().parent
INPUT = HERE / "input"
RUNDIR = HERE / "rundir"


def build_config() -> FullConfig:
    """Construct the FullConfig matching the SCKH_PES test case."""
    dynamics = DynamicsConfig(
        mu=1.0078825,
        grid=GridConfig(start=0.5, dx=0.025, npoints=77),
        time=TimeConfig(dt=0.1, nsteps=512),
        sampling=SamplingConfig(
            mode=1,
            npoints_x=10,
            npoints_mom=10,
            compatibility_mode="fortran",
        ),
        pes_initial=INPUT / "test_energy_tot.dat",
        pes_dynamics=INPUT / "test_energy_tot_exc.dat",
        units="angstrom",
        outfile="testout",
    )

    # Final-state surfaces are listed in the same order as the Fortran
    # pes_file_list_f.txt / dipole_file_list_f.txt (highest index first).
    final_indices = [9, 8, 7, 6, 5, 4]
    interpolation = InterpolationConfig(
        dipole_mode="DIPOLE",
        pes_final_list=[INPUT / f"pesfile_{i}.dat" for i in final_indices],
        dipole_final_list=[INPUT / f"dipolefile_{i}.dat" for i in final_indices],
        pes_lp_corr=INPUT / "test_lp_lp_energy.dat",
        compatibility_mode="fortran",
    )

    spectrum = SpectrumConfig(gamma_fwhm=0.18)

    return FullConfig(
        dynamics1d=dynamics,
        interpolation1d=interpolation,
        spectrum=spectrum,
    )


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    config = build_config()

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
