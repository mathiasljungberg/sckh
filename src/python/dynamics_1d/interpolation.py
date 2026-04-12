"""Surface interpolation for 1D dynamics trajectories.

Evaluates PES and dipole surfaces at 1D trajectory positions to produce
SCKHTrajectory objects. Spectrum computation is handled by
sckh.spectrum.SCKHSpectrumCalculator.
"""

from typing import List, Optional

import numpy as np

from .constants import CONST
from .dipole import Dipole1D, create_dipole_from_file, create_constant_dipole
from .io import read_pes_file, read_pes_file_raw, read_trajectory_file
from .pes import PES1D, create_pes_from_file
from .config import DynamicsConfig
from .trajectory import TrajectoryResult, DynamicsRunner

from sckh.trajectory import SCKHTrajectory
from sckh.config import FullConfig

from .spectrum_config import InterpolationConfig


class SurfaceInterpolator:
    """Evaluate PES and dipole surfaces along 1D trajectory positions.

    Converts dynamics trajectories into SCKHTrajectory objects by
    interpolating loaded surfaces at each trajectory position.

    Accepts either a FullConfig or individual configs::

        # From FullConfig (YAML workflow)
        interp = SurfaceInterpolator(config)

        # From specific configs
        interp = SurfaceInterpolator(
            dynamics_config=dyn_config, interp_config=interp_config
        )
    """

    def __init__(
        self,
        config=None,
        *,
        dynamics_config: Optional[DynamicsConfig] = None,
        interp_config: Optional[InterpolationConfig] = None,
    ):
        if config is not None and isinstance(config, FullConfig):
            self.dynamics_config = config.dynamics1d
            self.interp_config = config.interpolation1d
        else:
            self.dynamics_config = dynamics_config
            self.interp_config = interp_config

        if self.dynamics_config is None:
            raise ValueError("dynamics_config (or FullConfig with dynamics1d) required")
        if self.interp_config is None:
            raise ValueError("interp_config (or FullConfig with interpolation1d) required")

        # PES surfaces
        self.pes_n: Optional[PES1D] = None  # Intermediate state
        self.pes_f: List[PES1D] = []  # Final states

        # Dipole surfaces
        self.dipoles: List[Dipole1D] = []  # One per final state
        self.dipole_initial: Optional[Dipole1D] = None

        # Computed values
        self.E_mean: float = 0.0  # Mean transition energy (eV)
        self._surfaces_loaded = False

    def load_surfaces(self) -> None:
        """Load all PES and dipole surfaces from files."""
        units = self.dynamics_config.units

        # Load intermediate state PES (same as dynamics PES)
        self.pes_n = create_pes_from_file(
            self.dynamics_config.pes_dynamics, units=units
        )

        # Load final state PES with optional energy correction
        if self.interp_config.pes_lp_corr:
            # Apply energy shift so that first final state matches correction PES
            # Following Fortran: shift = E_lp_corr - E_f(1,:)
            _, E_corr = read_pes_file(
                self.interp_config.pes_lp_corr, units=units
            )
            _, E_f0 = read_pes_file(
                self.interp_config.pes_final_list[0], units=units
            )
            shift = E_corr - E_f0

            self.pes_f = []
            for path in self.interp_config.pes_final_list:
                x, E = read_pes_file(path, units=units)
                self.pes_f.append(PES1D(x=x, E=E + shift))
        else:
            # Load without energy correction
            self.pes_f = [
                create_pes_from_file(path, units=units)
                for path in self.interp_config.pes_final_list
            ]

        # Load dipole surfaces based on mode
        if self.interp_config.dipole_mode == "DIPOLE":
            self.dipoles = [
                create_dipole_from_file(path, units=units)
                for path in self.interp_config.dipole_final_list
            ]
        elif self.interp_config.dipole_mode == "FC":
            # Franck-Condon: constant dipole = 1
            x = self.pes_n.x
            self.dipoles = [
                create_constant_dipole(x, np.array([1.0, 1.0, 1.0]))
                for _ in range(len(self.pes_f))
            ]
        elif self.interp_config.dipole_mode == "DIPOLE_X0":
            # Dipole frozen at equilibrium - will be set during interpolation
            self.dipoles = [
                create_dipole_from_file(path, units=units)
                for path in self.interp_config.dipole_final_list
            ]

        # Load initial state dipole if provided
        if self.interp_config.dipole_initial:
            self.dipole_initial = create_dipole_from_file(
                self.interp_config.dipole_initial, units=units
            )

        self._surfaces_loaded = True

    def load_trajectories(self) -> List[TrajectoryResult]:
        """Load trajectories from files or run dynamics.

        Returns:
            List of TrajectoryResult objects
        """
        if self.interp_config.trajectory_files:
            # Load pre-computed trajectories
            return [
                read_trajectory_file(path)
                for path in self.interp_config.trajectory_files
            ]
        else:
            # Run dynamics to generate trajectories
            runner = DynamicsRunner(self.dynamics_config)
            result = runner.run()
            return result.trajectories

    def interpolate_along_trajectory(
        self,
        traj: TrajectoryResult,
    ) -> dict:
        """Interpolate energies and dipoles along trajectory positions.

        Args:
            traj: Trajectory result with positions x

        Returns:
            Dictionary with:
                E_n: (ntsteps,) intermediate state energies in Joules
                E_f: (nfinal, ntsteps) final state energies in Joules
                D_fn: (nfinal, ntsteps, 3) dipole moments
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        x = traj.x
        nsteps = len(x)
        n_final = len(self.pes_f)

        # Intermediate state energy
        E_n = self.pes_n.energy(x)

        # Final state energies
        E_f = np.zeros((n_final, nsteps))
        for i, pes in enumerate(self.pes_f):
            E_f[i] = pes.energy(x)

        # Dipole moments
        D_fn = np.zeros((n_final, nsteps, 3))

        if self.interp_config.dipole_mode == "DIPOLE_X0":
            # Frozen at equilibrium position
            pes_i = create_pes_from_file(
                self.dynamics_config.pes_initial,
                units=self.dynamics_config.units,
            )
            x_eq, _ = pes_i.find_minimum()

            for i, dipole in enumerate(self.dipoles):
                d_eq = dipole.dipole(x_eq)
                D_fn[i, :, :] = d_eq  # Broadcast to all time steps
        else:
            for i, dipole in enumerate(self.dipoles):
                D_fn[i] = dipole.dipole(x)

        return {"E_n": E_n, "E_f": E_f, "D_fn": D_fn}

    def trajectories_to_sckh(
        self,
        trajectories: List[TrajectoryResult],
    ) -> List[SCKHTrajectory]:
        """Convert dynamics trajectories to SCKHTrajectory objects.

        Evaluates loaded PES/dipole surfaces at trajectory positions and
        packages the results as SCKHTrajectory objects.  The resulting
        objects can be inspected, saved with ``write_sckh_trajectory``,
        or passed to ``compute_spectrum_from_sckh``.

        Args:
            trajectories: List of dynamics trajectories

        Returns:
            List of SCKHTrajectory with energies in eV, dipoles in a.u.
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        sckh_trajs = []
        for traj in trajectories:
            interp = self.interpolate_along_trajectory(traj)
            nsteps = len(traj.x)
            sckh_trajs.append(SCKHTrajectory(
                time=traj.time,
                E_gs=np.zeros(nsteps),
                E_n=interp["E_n"] / CONST.eV,
                E_f=interp["E_f"] / CONST.eV,
                D_fn=interp["D_fn"],
                E_IP1s=np.zeros(nsteps),
            ))
        return sckh_trajs

    def compute_D_ni(self) -> np.ndarray:
        """Compute initial state dipole vector.

        Returns dipole_initial evaluated at equilibrium, or [1,1,1] default.
        """
        if self.dipole_initial:
            pes_i = create_pes_from_file(
                self.dynamics_config.pes_initial,
                units=self.dynamics_config.units,
            )
            x_eq, _ = pes_i.find_minimum()
            return self.dipole_initial.dipole(x_eq)
        return np.array([1.0, 1.0, 1.0])

    def compute_mean_transition_energy(
        self,
        trajectories: List[TrajectoryResult],
    ) -> float:
        """Compute mean transition energy.

        In "standard" mode: Finds the true equilibrium position x_eq from the
        initial state PES and evaluates E_mean = E_n(x_eq) - E_f(x_eq).

        In "fortran" mode: Uses minloc(E_i_inp) to find an index in the input
        PES file, then uses that index into the trajectory arrays. This matches
        Fortran's implementation for validation purposes.

        Args:
            trajectories: List of trajectories

        Returns:
            Mean transition energy in eV
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        mode = self.interp_config.compatibility_mode

        if mode == "fortran":
            return self._compute_E_mean_fortran(trajectories)
        else:
            return self._compute_E_mean_standard()

    def _compute_E_mean_standard(self) -> float:
        """Standard E_mean: evaluate at true equilibrium position."""
        pes_i = create_pes_from_file(
            self.dynamics_config.pes_initial,
            units=self.dynamics_config.units,
        )
        x_eq, _ = pes_i.find_minimum()

        E_n_at_eq = self.pes_n.energy(x_eq)
        E_f_at_eq = self.pes_f[-1].energy(x_eq)

        E_mean = (E_n_at_eq - E_f_at_eq) / CONST.eV
        return E_mean

    def _compute_E_mean_fortran(
        self,
        trajectories: List[TrajectoryResult],
    ) -> float:
        """Fortran-compatible E_mean: use minloc index into trajectory."""
        _, E_i_raw = read_pes_file_raw(self.dynamics_config.pes_initial)
        ind = int(np.argmin(E_i_raw))

        x_traj = trajectories[0].x

        E_n_at_ind = self.pes_n.energy(x_traj[ind])
        E_f_at_ind = self.pes_f[-1].energy(x_traj[ind])

        E_mean = (E_n_at_ind - E_f_at_ind) / CONST.eV
        return E_mean


# Backwards compatibility
SpectrumCalculator = SurfaceInterpolator
