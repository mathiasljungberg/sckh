"""Surface interpolation for 2D dynamics trajectories.

Evaluates PES and dipole surfaces at 2D trajectory positions to produce
SCKHTrajectory objects. Spectrum computation is handled by
sckh.spectrum.SCKHSpectrumCalculator.
"""

from typing import List, Optional

import numpy as np

from dynamics_1d.constants import CONST

from sckh.trajectory import SCKHTrajectory
from sckh.config import FullConfig

from .config import DynamicsConfig2D, InterpolationConfig2D
from .dipole import Dipole2D, create_dipole_from_file_2d, create_constant_dipole_2d
from .io import read_pes_file_2d, read_pes_file_2d_raw
from .pes import PES2D, create_pes_from_file_2d
from .trajectory import TrajectoryResult2D, DynamicsRunner2D


class SurfaceInterpolator2D:
    """Evaluate 2D PES and dipole surfaces along trajectory positions.

    Produces SCKHTrajectory objects that can be passed to
    SCKHSpectrumCalculator for spectrum computation.

    Usage::

        interp = SurfaceInterpolator2D(config)
        interp.load_surfaces()
        sckh_trajs = interp.trajectories_to_sckh(trajectories)
    """

    def __init__(
        self,
        config=None,
        *,
        dynamics_config: Optional[DynamicsConfig2D] = None,
        interp_config: Optional[InterpolationConfig2D] = None,
    ):
        """Initialize with FullConfig or individual configs.

        Args:
            config: FullConfig (extracts dynamics2d + interpolation2d)
            dynamics_config: DynamicsConfig2D (if not using FullConfig)
            interp_config: InterpolationConfig2D (if not using FullConfig)
        """
        if config is not None and isinstance(config, FullConfig):
            self.dynamics_config = config.dynamics2d
            self.interp_config = config.interpolation2d
        elif dynamics_config is not None:
            self.dynamics_config = dynamics_config
            self.interp_config = interp_config
        else:
            raise ValueError(
                "Provide a FullConfig or dynamics_config + interp_config"
            )

        # PES surfaces
        self.pes_n: Optional[PES2D] = None  # Intermediate state
        self.pes_f: List[PES2D] = []  # Final states

        # Dipole surfaces
        self.dipoles: List[Dipole2D] = []  # One per final state
        self.dipole_initial: Optional[Dipole2D] = None

        # Computed values
        self.E_mean: float = 0.0  # Mean transition energy (eV)
        self._surfaces_loaded = False

    def load_surfaces(self) -> None:
        """Load all PES and dipole surfaces from files."""
        pos_units = self.dynamics_config.position_units
        energy_units = self.dynamics_config.energy_units
        index_order = self.dynamics_config.index_order

        # Load intermediate state PES
        if self.interp_config.pes_intermediate:
            pes_file = self.interp_config.pes_intermediate
        else:
            pes_file = self.dynamics_config.pes_dynamics

        self.pes_n = create_pes_from_file_2d(
            pes_file,
            position_units=pos_units,
            energy_units=energy_units,
            index_order=index_order,
        )

        # Load final state PES with optional energy correction
        if self.interp_config.pes_lp_corr:
            # Apply energy shift so that first final state matches correction PES
            _, _, E_corr = read_pes_file_2d(
                self.interp_config.pes_lp_corr,
                position_units=pos_units,
                energy_units=energy_units,
                index_order=index_order,
            )
            _, _, E_f0 = read_pes_file_2d(
                self.interp_config.pes_final_list[0],
                position_units=pos_units,
                energy_units=energy_units,
                index_order=index_order,
            )
            shift = E_corr - E_f0

            self.pes_f = []
            for path in self.interp_config.pes_final_list:
                x1, x2, E = read_pes_file_2d(
                    path,
                    position_units=pos_units,
                    energy_units=energy_units,
                    index_order=index_order,
                )
                self.pes_f.append(PES2D(x1=x1, x2=x2, E=E + shift))
        else:
            # Load without energy correction
            self.pes_f = [
                create_pes_from_file_2d(
                    path,
                    position_units=pos_units,
                    energy_units=energy_units,
                    index_order=index_order,
                )
                for path in self.interp_config.pes_final_list
            ]

        # Load dipole surfaces based on mode
        dipole_components = self.interp_config.dipole_components
        if self.interp_config.dipole_mode == "DIPOLE":
            self.dipoles = [
                create_dipole_from_file_2d(
                    path,
                    position_units=pos_units,
                    index_order=index_order,
                    dipole_components=dipole_components,
                )
                for path in self.interp_config.dipole_final_list
            ]
        elif self.interp_config.dipole_mode == "FC":
            # Franck-Condon: constant dipole = 1
            self.dipoles = [
                create_constant_dipole_2d(
                    self.pes_n.x1, self.pes_n.x2, np.array([1.0, 1.0, 1.0])
                )
                for _ in range(len(self.pes_f))
            ]
        elif self.interp_config.dipole_mode == "DIPOLE_X0":
            # Dipole frozen at equilibrium - will be handled during interpolation
            self.dipoles = [
                create_dipole_from_file_2d(
                    path,
                    position_units=pos_units,
                    index_order=index_order,
                    dipole_components=dipole_components,
                )
                for path in self.interp_config.dipole_final_list
            ]

        # Load initial state dipole if provided
        if self.interp_config.dipole_initial:
            self.dipole_initial = create_dipole_from_file_2d(
                self.interp_config.dipole_initial,
                position_units=pos_units,
                index_order=index_order,
                dipole_components=dipole_components,
            )

        self._surfaces_loaded = True

    def load_trajectories(self) -> List[TrajectoryResult2D]:
        """Load trajectories from files or run dynamics.

        Returns:
            List of TrajectoryResult2D objects
        """
        if self.interp_config.trajectory_files:
            # Load pre-computed trajectories
            from .io import read_trajectory_2d
            return [
                read_trajectory_2d(path)
                for path in self.interp_config.trajectory_files
            ]
        else:
            # Run dynamics to generate trajectories
            full_config = FullConfig(dynamics2d=self.dynamics_config)
            runner = DynamicsRunner2D(full_config)
            result = runner.run()
            return result.trajectories

    def interpolate_along_trajectory(
        self,
        traj: TrajectoryResult2D,
    ) -> dict:
        """Interpolate energies and dipoles along 2D trajectory positions.

        Args:
            traj: Trajectory result with positions x1, x2

        Returns:
            Dictionary with:
                E_n: (ntsteps,) intermediate state energies in Joules
                E_f: (nfinal, ntsteps) final state energies in Joules
                D_fn: (nfinal, ntsteps, 3) dipole moments
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        x1 = traj.x1
        x2 = traj.x2
        nsteps = len(x1)
        n_final = len(self.pes_f)

        # Intermediate state energy: E_n(t) = V_n(x1(t), x2(t))
        E_n = self.pes_n.energy(x1, x2)

        # Final state energies
        E_f = np.zeros((n_final, nsteps))
        for i, pes in enumerate(self.pes_f):
            E_f[i] = pes.energy(x1, x2)

        # Dipole moments
        D_fn = np.zeros((n_final, nsteps, 3))

        if self.interp_config.dipole_mode == "DIPOLE_X0":
            # Frozen at equilibrium position
            pes_i = create_pes_from_file_2d(
                self.dynamics_config.pes_initial,
                position_units=self.dynamics_config.position_units,
                energy_units=self.dynamics_config.energy_units,
                index_order=self.dynamics_config.index_order,
            )
            x1_eq, x2_eq, _ = pes_i.find_minimum()

            for i, dipole in enumerate(self.dipoles):
                d_eq = dipole.dipole(x1_eq, x2_eq)
                D_fn[i, :, :] = d_eq  # Broadcast to all time steps
        else:
            for i, dipole in enumerate(self.dipoles):
                D_fn[i] = dipole.dipole(x1, x2)

        return {"E_n": E_n, "E_f": E_f, "D_fn": D_fn}

    def trajectories_to_sckh(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> List[SCKHTrajectory]:
        """Convert 2D dynamics trajectories to SCKHTrajectory objects.

        Evaluates loaded PES/dipole surfaces at 2D trajectory positions
        and packages the results as SCKHTrajectory objects.

        Args:
            trajectories: List of 2D dynamics trajectories

        Returns:
            List of SCKHTrajectory with energies in eV, dipoles in a.u.
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        sckh_trajs = []
        for traj in trajectories:
            interp = self.interpolate_along_trajectory(traj)
            nsteps = len(traj.x1)
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

        Returns dipole_initial evaluated at 2D equilibrium, or [1,1,1] default.
        """
        if self.dipole_initial:
            pes_i = create_pes_from_file_2d(
                self.dynamics_config.pes_initial,
                position_units=self.dynamics_config.position_units,
                energy_units=self.dynamics_config.energy_units,
                index_order=self.dynamics_config.index_order,
            )
            x1_eq, x2_eq, _ = pes_i.find_minimum()
            return self.dipole_initial.dipole(x1_eq, x2_eq)
        return np.array([1.0, 1.0, 1.0])

    def compute_mean_transition_energy(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> float:
        """Compute mean transition energy.

        In "standard" mode: Finds the true 2D equilibrium position (x1_eq, x2_eq)
        from the initial state PES and evaluates E_mean = E_n(eq) - E_f(eq).

        In "fortran" mode: Uses trajectory-based indexing for compatibility.

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
        """Standard E_mean: evaluate at true 2D equilibrium position."""
        pes_i = create_pes_from_file_2d(
            self.dynamics_config.pes_initial,
            position_units=self.dynamics_config.position_units,
            energy_units=self.dynamics_config.energy_units,
            index_order=self.dynamics_config.index_order,
        )
        x1_eq, x2_eq, _ = pes_i.find_minimum()

        E_n_at_eq = self.pes_n.energy(x1_eq, x2_eq)
        E_f_at_eq = self.pes_f[-1].energy(x1_eq, x2_eq)

        E_mean = (E_n_at_eq - E_f_at_eq) / CONST.eV
        return E_mean

    def _compute_E_mean_fortran(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> float:
        """Fortran-compatible E_mean: use trajectory-based indexing."""
        _, _, E_i_raw = read_pes_file_2d_raw(
            self.dynamics_config.pes_initial,
            index_order=self.dynamics_config.index_order,
        )
        ind = np.unravel_index(np.argmin(E_i_raw), E_i_raw.shape)

        x1_traj = trajectories[0].x1
        x2_traj = trajectories[0].x2

        idx = min(ind[0] * E_i_raw.shape[1] + ind[1], len(x1_traj) - 1)
        E_n_at_ind = self.pes_n.energy(x1_traj[idx], x2_traj[idx])
        E_f_at_ind = self.pes_f[-1].energy(x1_traj[idx], x2_traj[idx])

        E_mean = (E_n_at_ind - E_f_at_ind) / CONST.eV
        return E_mean


# Backwards compatibility
SpectrumCalculator2D = SurfaceInterpolator2D
