"""Standalone SCKH_PES spectrum example — pure Python API.

Builds a :class:`FullConfig` in Python and runs the full SCKH workflow
(dynamics → surface interpolation → FFT spectrum) in a single call.
Reads PES and dipole surfaces from ``./input/`` and writes the computed
spectrum files to ``./rundir/``.  Usable as a copy-pasteable example of
the high-level 1D dynamics + SCKH API.

Run directly as::

    python run_python.py
"""

import shutil
from pathlib import Path

from dynamics_1d import compute_spectrum_1d
from dynamics_1d.config import (
    DynamicsConfig,
    GridConfig,
    SamplingConfig,
    TimeConfig,
)
from dynamics_1d.spectrum_config import InterpolationConfig
from sckh import FullConfig, SpectrumConfig, write_spectrum_result


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
    result = compute_spectrum_1d(config)
    write_spectrum_result(RUNDIR / "testout_sigma", result)


if __name__ == "__main__":
    main()
