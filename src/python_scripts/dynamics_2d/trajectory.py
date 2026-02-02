"""High-level trajectory management for 2D dynamics."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from python_scripts.dynamics_1d.constants import CONST

from .config import DynamicsConfig2D, FullConfig2D
from .integrators import run_trajectory_2d
from .pes import PES2D, create_pes_from_file_2d
from .sampling import create_initial_conditions_2d, get_n_trajectories
from .vibrational import ProductGroundState, solve_product_ground_state


@dataclass
class TrajectoryResult2D:
    """Result from a single 2D trajectory."""

    time: np.ndarray  # Time points (SI: seconds)
    x1: np.ndarray  # Positions in x1 (SI: meters)
    x2: np.ndarray  # Positions in x2 (SI: meters)
    v1: np.ndarray  # Velocities in x1 (SI: m/s)
    v2: np.ndarray  # Velocities in x2 (SI: m/s)
    a1: np.ndarray  # Accelerations in x1 (SI: m/s^2)
    a2: np.ndarray  # Accelerations in x2 (SI: m/s^2)
    x1_0: float  # Initial position in x1
    x2_0: float  # Initial position in x2
    p1_0: float  # Initial momentum in x1
    p2_0: float  # Initial momentum in x2


@dataclass
class EnsembleResult2D:
    """Result from ensemble of 2D trajectories."""

    trajectories: List[TrajectoryResult2D]
    config: DynamicsConfig2D  # Configuration used
    ground_state: ProductGroundState  # Ground state information
    x1_grid: np.ndarray  # Position grid for x1
    x2_grid: np.ndarray  # Position grid for x2
    x1_eq: float  # Equilibrium position for x1 (SI: meters)
    x2_eq: float  # Equilibrium position for x2 (SI: meters)

    @property
    def n_trajectories(self) -> int:
        return len(self.trajectories)

    @property
    def n_steps(self) -> int:
        if self.trajectories:
            return len(self.trajectories[0].time)
        return 0

    def get_positions_at_time(
        self, step: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get all (x1, x2) positions at a given time step."""
        x1 = np.array([t.x1[step] for t in self.trajectories])
        x2 = np.array([t.x2[step] for t in self.trajectories])
        return x1, x2

    def get_velocities_at_time(
        self, step: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get all (v1, v2) velocities at a given time step."""
        v1 = np.array([t.v1[step] for t in self.trajectories])
        v2 = np.array([t.v2[step] for t in self.trajectories])
        return v1, v2

    def get_time_array(self) -> np.ndarray:
        """Get time array from first trajectory."""
        if self.trajectories:
            return self.trajectories[0].time
        return np.array([])


def _extract_1d_slices(
    pes_2d: PES2D,
    x1_grid: np.ndarray,
    x2_grid: np.ndarray,
) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """Extract 1D PES slices at equilibrium geometry.

    Finds the equilibrium of the 2D PES and extracts 1D slices:
        V1(x1) = V(x1, x2_eq)
        V2(x2) = V(x1_eq, x2)

    Args:
        pes_2d: 2D potential energy surface
        x1_grid: Grid points for x1 (SI: meters)
        x2_grid: Grid points for x2 (SI: meters)

    Returns:
        x1_eq: Equilibrium x1 position (SI: meters)
        x2_eq: Equilibrium x2 position (SI: meters)
        V1: 1D potential along x1 at x2=x2_eq (SI: Joules)
        V2: 1D potential along x2 at x1=x1_eq (SI: Joules)
    """
    # Find equilibrium geometry
    x1_eq, x2_eq, _ = pes_2d.find_minimum()

    # Extract 1D slices at equilibrium
    V1 = np.array([pes_2d.energy(x1, x2_eq) for x1 in x1_grid])
    V2 = np.array([pes_2d.energy(x1_eq, x2) for x2 in x2_grid])

    return x1_eq, x2_eq, V1, V2


class DynamicsRunner2D:
    """Main class for running 2D classical dynamics simulations.

    Usage:
        config = load_config("dynamics_2d.yaml")
        runner = DynamicsRunner2D(config)
        result = runner.run()
    """

    def __init__(self, config: FullConfig2D):
        self.full_config = config
        self.config = config.dynamics2d

        # Convert masses to SI
        self.m1_SI = self.config.mu1 * CONST.u
        self.m2_SI = self.config.mu2 * CONST.u

        # Convert time step to SI
        self.dt_SI = self.config.time.dt * 1e-15  # fs to seconds

        # Build position grids in SI units
        if self.config.position_units.lower() == "angstrom":
            conv_factor = 1e-10
        elif self.config.position_units.lower() == "bohr":
            conv_factor = CONST.bohr
        else:
            raise ValueError(f"Unknown position units: {self.config.position_units}")

        self.x1_grid = (
            self.config.grid_x1.start
            + np.arange(self.config.grid_x1.npoints) * self.config.grid_x1.dx
        ) * conv_factor
        self.x2_grid = (
            self.config.grid_x2.start
            + np.arange(self.config.grid_x2.npoints) * self.config.grid_x2.dx
        ) * conv_factor

        self.dx1_SI = self.config.grid_x1.dx * conv_factor
        self.dx2_SI = self.config.grid_x2.dx * conv_factor

        # Load 2D PES for dynamics
        self.pes = create_pes_from_file_2d(
            self.config.pes_dynamics,
            position_units=self.config.position_units,
            energy_units=self.config.energy_units,
            index_order=self.config.index_order,
        )

        # Load 2D PES for initial state
        self.pes_initial = create_pes_from_file_2d(
            self.config.pes_initial,
            position_units=self.config.position_units,
            energy_units=self.config.energy_units,
            index_order=self.config.index_order,
        )

        # Placeholder for ground state and equilibrium
        self._ground_state: Optional[ProductGroundState] = None
        self._x1_eq: Optional[float] = None
        self._x2_eq: Optional[float] = None

    def solve_ground_state(self) -> None:
        """Solve product approximation ground state.

        Finds equilibrium of initial PES and computes 1D wavefunctions
        from slices at equilibrium:
            psi1(x1) from V(x1, x2_eq)
            psi2(x2) from V(x1_eq, x2)
        """
        # Extract 1D slices at equilibrium
        self._x1_eq, self._x2_eq, V1, V2 = _extract_1d_slices(
            self.pes_initial, self.x1_grid, self.x2_grid
        )

        # Solve using product approximation with the 1D slices
        self._ground_state = solve_product_ground_state(
            self.x1_grid,
            V1,
            self.m1_SI,
            self.x2_grid,
            V2,
            self.m2_SI,
        )

    def set_ground_state(
        self,
        ground_state: ProductGroundState,
        x1_eq: float = 0.0,
        x2_eq: float = 0.0,
    ) -> None:
        """Manually set the ground state.

        Args:
            ground_state: ProductGroundState object
            x1_eq: Equilibrium x1 position (SI: meters)
            x2_eq: Equilibrium x2 position (SI: meters)
        """
        self._ground_state = ground_state
        self._x1_eq = x1_eq
        self._x2_eq = x2_eq

    @property
    def ground_state(self) -> ProductGroundState:
        if self._ground_state is None:
            raise RuntimeError(
                "Ground state not set. Call solve_ground_state() first."
            )
        return self._ground_state

    @property
    def x1_eq(self) -> float:
        if self._x1_eq is None:
            raise RuntimeError(
                "Equilibrium not set. Call solve_ground_state() first."
            )
        return self._x1_eq

    @property
    def x2_eq(self) -> float:
        if self._x2_eq is None:
            raise RuntimeError(
                "Equilibrium not set. Call solve_ground_state() first."
            )
        return self._x2_eq

    def run(self, verbose: bool = False) -> EnsembleResult2D:
        """Run ensemble of 2D trajectories.

        Args:
            verbose: Print progress information

        Returns:
            EnsembleResult2D with all trajectory data
        """
        # Ensure ground state is available
        if self._ground_state is None:
            if verbose:
                print("Solving product ground state...")
            self.solve_ground_state()

        # Calculate number of trajectories
        n_traj = get_n_trajectories(
            self.config.sampling.npoints_x1,
            self.config.sampling.npoints_x2,
            self.config.sampling.npoints_p1,
            self.config.sampling.npoints_p2,
            self.config.sampling.mode,
        )

        if verbose:
            print(
                f"Sampling initial conditions: "
                f"{self.config.sampling.npoints_x1} x {self.config.sampling.npoints_x2} positions, "
                f"{self.config.sampling.npoints_p1} x {self.config.sampling.npoints_p2} momenta = "
                f"{n_traj} trajectories"
            )

        # Sample initial conditions
        x1_init, x2_init, p1_init, p2_init = create_initial_conditions_2d(
            self.x1_grid,
            self._ground_state.psi1,
            self.x2_grid,
            self._ground_state.psi2,
            self.config.sampling.npoints_x1,
            self.config.sampling.npoints_x2,
            self.config.sampling.npoints_p1,
            self.config.sampling.npoints_p2,
            self.config.sampling.mode,
        )

        # Run trajectories
        trajectories = []
        for i, (x1_0, x2_0, p1_0, p2_0) in enumerate(
            zip(x1_init, x2_init, p1_init, p2_init)
        ):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Running trajectory {i + 1}/{n_traj}")

            # Convert momenta to velocities
            v1_0 = p1_0 / self.m1_SI
            v2_0 = p2_0 / self.m2_SI

            time, x1, x2, v1, v2, a1, a2 = run_trajectory_2d(
                x1_0=x1_0,
                x2_0=x2_0,
                v1_0=v1_0,
                v2_0=v2_0,
                force_func=self.pes.forces,
                m1=self.m1_SI,
                m2=self.m2_SI,
                dt=self.dt_SI,
                nsteps=self.config.time.nsteps,
            )

            trajectories.append(
                TrajectoryResult2D(
                    time=time,
                    x1=x1,
                    x2=x2,
                    v1=v1,
                    v2=v2,
                    a1=a1,
                    a2=a2,
                    x1_0=x1_0,
                    x2_0=x2_0,
                    p1_0=p1_0,
                    p2_0=p2_0,
                )
            )

        if verbose:
            print(f"Completed {n_traj} trajectories")

        return EnsembleResult2D(
            trajectories=trajectories,
            config=self.config,
            ground_state=self._ground_state,
            x1_grid=self.x1_grid,
            x2_grid=self.x2_grid,
            x1_eq=self._x1_eq,
            x2_eq=self._x2_eq,
        )

    def save_results(
        self,
        result: EnsembleResult2D,
        output_dir: Path,
        units: str = "user",
    ) -> None:
        """Save results to files.

        Args:
            result: EnsembleResult2D from run()
            output_dir: Directory for output files
            units: "SI" or "user" (Angstrom/fs)
        """
        from .io import write_trajectory_2d, write_ground_state_2d

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        basename = self.config.outfile

        # Save ground state info
        write_ground_state_2d(
            output_dir / f"{basename}_ground_state.dat",
            result.ground_state,
            units=units,
        )

        # Save each trajectory
        for i, traj in enumerate(result.trajectories):
            write_trajectory_2d(
                output_dir / f"{basename}_traj_{i:04d}.dat",
                traj,
                units=units,
            )
