"""High-level 1D SCKH workflow helpers.

Wraps the standard ``DynamicsRunner`` + ``SurfaceInterpolator`` +
``SCKHSpectrumCalculator`` pipeline into a single call so that example
scripts and user code do not need to repeat the same three-stage
boilerplate.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sckh.config import FullConfig
    from sckh.spectrum import SpectrumResult


def compute_spectrum_1d(config: "FullConfig") -> "SpectrumResult":
    """Run the full 1D SCKH workflow and return the spectrum.

    Steps:
        1. Run classical trajectory dynamics (``DynamicsRunner``).
        2. Interpolate PES and dipole surfaces along trajectories, and
           compute the mean transition energy + initial-state dipole
           (``SurfaceInterpolator``).
        3. Compute the SCKH spectrum via FFT (``SCKHSpectrumCalculator``).

    The config must have ``dynamics1d``, ``interpolation1d``, and
    ``spectrum`` populated.

    Args:
        config: FullConfig with 1D dynamics + interpolation + spectrum.

    Returns:
        SpectrumResult with omega, sigma_tot, sigma_f, E_mean, etc.
    """
    from sckh.spectrum import SCKHSpectrumCalculator

    from .interpolation import SurfaceInterpolator
    from .trajectory import DynamicsRunner

    runner = DynamicsRunner(config)
    result_dyn = runner.run(verbose=False)

    interp = SurfaceInterpolator(config)
    interp.load_surfaces()
    sckh_trajs = interp.trajectories_to_sckh(result_dyn.trajectories)
    E_mean = interp.compute_mean_transition_energy(result_dyn.trajectories)
    D_ni = interp.compute_D_ni()

    calc = SCKHSpectrumCalculator(config)
    return calc.compute_spectrum(sckh_trajs, E_mean=E_mean, D_ni=D_ni)
