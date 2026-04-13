"""High-level 1D dynamics + surface-interpolation workflow.

Runs classical trajectory dynamics and evaluates PES/dipole surfaces
along those trajectories to produce SCKH trajectories, plus the
auxiliary quantities (mean transition energy, initial-state dipole)
that the SCKH spectrum code needs.

This function deliberately stays inside the ``dynamics_1d`` domain:
it does not compute a spectrum.  Pass its output to
:class:`sckh.SCKHSpectrumCalculator` (or :func:`sckh.compute_spectrum_from_sckh`)
to obtain a spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, List

import numpy as np

if TYPE_CHECKING:
    from sckh.config import FullConfig
    from sckh.trajectory import SCKHTrajectory


@dataclass
class SCKHInputs1D:
    """Output of the 1D dynamics + surface-interpolation stage.

    Bundles the quantities required as inputs to the SCKH spectrum
    calculator.

    Attributes:
        trajectories: SCKH trajectories with energies and dipoles
            evaluated along each dynamics trajectory.
        E_mean: Mean transition energy (eV) evaluated from the
            dynamics and final-state surfaces.
        D_ni: Initial-state dipole vector (shape ``(3,)``).
    """

    trajectories: List["SCKHTrajectory"]
    E_mean: float
    D_ni: np.ndarray


def compute_sckh_trajectories_1d(config: "FullConfig") -> SCKHInputs1D:
    """Run 1D dynamics and build SCKH trajectories from surface interpolation.

    Steps:
        1. Run classical trajectory dynamics (``DynamicsRunner``).
        2. Evaluate final-state PES and dipole surfaces along each
           trajectory, and compute the mean transition energy and
           initial-state dipole (``SurfaceInterpolator``).

    The config must have ``dynamics1d`` and ``interpolation1d``
    populated.  Spectrum computation itself is not performed here —
    pass the returned ``SCKHInputs1D`` to
    :class:`sckh.SCKHSpectrumCalculator`.

    Args:
        config: FullConfig with 1D dynamics + interpolation.

    Returns:
        SCKHInputs1D with ``trajectories``, ``E_mean``, and ``D_ni``.
    """
    from .interpolation import SurfaceInterpolator
    from .trajectory import DynamicsRunner

    runner = DynamicsRunner(config)
    result_dyn = runner.run(verbose=False)

    interp = SurfaceInterpolator(config)
    interp.load_surfaces()
    trajectories = interp.trajectories_to_sckh(result_dyn.trajectories)
    E_mean = interp.compute_mean_transition_energy(result_dyn.trajectories)
    D_ni = interp.compute_D_ni()

    return SCKHInputs1D(
        trajectories=trajectories,
        E_mean=E_mean,
        D_ni=D_ni,
    )
