"""Standalone SCKH_PES_2D spectrum example — pure Python API.

Builds a :class:`FullConfig` in Python, runs the 2D dynamics stage to produce
trajectories, then interpolates surfaces along those trajectories and computes
an SCKH spectrum.

The two stages are kept explicit:

    dynamics_2d  ->  trajectories
    sckh         ->  spectrum

Run directly as::

    python run_python.py
"""

import shutil
from pathlib import Path

from dynamics_1d.spectrum_config import SpectrumConfig
from dynamics_2d import (
    DynamicsConfig2D,
    DynamicsRunner2D,
    FullConfig,
    GridConfig2D,
    InterpolationConfig2D,
    SamplingConfig2D,
    SCKHSpectrumCalculator,
    SurfaceInterpolator2D,
    TimeConfig,
)


HERE = Path(__file__).resolve().parent
INPUT = HERE / "input"
RUNDIR = HERE / "rundir"


def build_config() -> FullConfig:
    """Construct the FullConfig matching the SCKH_PES_2D test case."""
    dynamics = DynamicsConfig2D(
        mu1=1.0078825,
        mu2=1.0078825,
        grid_x1=GridConfig2D(start=-0.35, dx=0.05, npoints=34),
        grid_x2=GridConfig2D(start=-0.35, dx=0.05, npoints=34),
        time=TimeConfig(dt=0.1, nsteps=1024),
        sampling=SamplingConfig2D(
            mode=1,
            npoints_x1=2,
            npoints_x2=2,
            npoints_p1=2,
            npoints_p2=2,
        ),
        pes_initial=INPUT / "tmpGS.pes.dat",
        pes_dynamics=INPUT / "tmpCH.pes.dat",
        index_order="F",
        outfile="dynamics_2d_out",
    )

    interpolation = InterpolationConfig2D(
        dipole_mode="DIPOLE",
        pes_final_list=[
            INPUT / "tmpVH001.pes.dat",
            INPUT / "tmpVH002.pes.dat",
            INPUT / "tmpVH003.pes.dat",
        ],
        dipole_final_list=[
            INPUT / "tmpVH001.dip.dat",
            INPUT / "tmpVH002.dip.dat",
            INPUT / "tmpVH003.dip.dat",
        ],
        dipole_components=3,
    )

    spectrum = SpectrumConfig(gamma_fwhm=0.15)

    return FullConfig(
        dynamics2d=dynamics,
        interpolation2d=interpolation,
        spectrum=spectrum,
    )


def main() -> None:
    if RUNDIR.exists():
        shutil.rmtree(RUNDIR)
    RUNDIR.mkdir()

    config = build_config()

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
