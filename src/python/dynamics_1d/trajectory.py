"""High-level trajectory management."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from .config import DynamicsConfig
from .spectrum_config import FullConfig
from .constants import CONST
from .integrators import run_trajectory
from .io import write_trajectory, write_distribution, write_eigenstate
from .pes import PES1D, create_pes_from_file
from .sampling import create_initial_conditions
from .vibrational import solve_vibrational


@dataclass
class TrajectoryResult:
    """Result from a single trajectory."""

    time: np.ndarray  # Time points (SI: seconds)
    x: np.ndarray  # Positions (SI: meters)
    v: np.ndarray  # Velocities (SI: m/s)
    a: np.ndarray  # Accelerations (SI: m/s^2)
    x0: float  # Initial position
    p0: float  # Initial momentum


@dataclass
class EnsembleResult:
    """Result from ensemble of trajectories."""

    trajectories: List[TrajectoryResult]
    config: DynamicsConfig  # Dynamics config (convenience access)
    ground_state_energy: float  # Ground state eigenvalue (J)
    x_grid: np.ndarray  # Position grid used
    psi_ground: np.ndarray  # Ground state wavefunction
    full_config: Optional[FullConfig] = None  # Full config

    @property
    def n_trajectories(self) -> int:
        return len(self.trajectories)

    @property
    def n_steps(self) -> int:
        if self.trajectories:
            return len(self.trajectories[0].time)
        return 0

    def get_positions_at_time(self, step: int) -> np.ndarray:
        """Get all positions at a given time step."""
        return np.array([t.x[step] for t in self.trajectories])

    def get_velocities_at_time(self, step: int) -> np.ndarray:
        """Get all velocities at a given time step."""
        return np.array([t.v[step] for t in self.trajectories])

    def get_time_array(self) -> np.ndarray:
        """Get time array from first trajectory."""
        if self.trajectories:
            return self.trajectories[0].time
        return np.array([])

    def compute_position_histogram(
        self,
        step: int,
        bins: np.ndarray,
    ) -> np.ndarray:
        """Compute position histogram at a given time step.

        Args:
            step: Time step index
            bins: Bin edges for histogram

        Returns:
            Histogram counts (normalized to density)
        """
        positions = self.get_positions_at_time(step)
        counts, _ = np.histogram(positions, bins=bins, density=True)
        return counts


class DynamicsRunner:
    """Main class for running classical dynamics simulations.

    Usage:
        config = load_full_config("dynamics.yaml")
        runner = DynamicsRunner(config)
        result = runner.run()
    """

    def __init__(self, config: FullConfig):
        self.full_config = config
        self.config = config.dynamics

        # Use dynamics config for setup
        dyn = self.config

        # Convert mass to SI
        self.mass_SI = dyn.mu * CONST.u

        # Convert time step to SI
        self.dt_SI = dyn.time.dt * 1e-15  # fs to seconds

        # Build position grid in SI units
        if dyn.units.lower() == "angstrom":
            grid_start_SI = dyn.grid.start * 1e-10
            dx_SI = dyn.grid.dx * 1e-10
        elif dyn.units.lower() == "bohr":
            grid_start_SI = dyn.grid.start * CONST.bohr
            dx_SI = dyn.grid.dx * CONST.bohr
        else:
            raise ValueError(f"Unknown units: {dyn.units}")

        self.x_grid = grid_start_SI + np.arange(dyn.grid.npoints) * dx_SI
        self.dx_SI = dx_SI

        # Load PES for dynamics
        self.pes = create_pes_from_file(
            dyn.pes_dynamics,
            units=dyn.units,
        )

        # Load initial state PES for vibrational problem
        self.pes_initial = create_pes_from_file(
            dyn.pes_initial,
            units=dyn.units,
        )

        # Placeholder for ground state
        self._psi_ground = None
        self._eigenvalue_ground = None

    def solve_ground_state(self) -> None:
        """Solve vibrational problem to get ground state wavefunction."""
        # Get potential on grid
        V = self.pes_initial.energy(self.x_grid)

        # Solve Schrödinger equation
        eigenvalues, eigenvectors = solve_vibrational(
            self.x_grid,
            V,
            self.mass_SI,
            n_states=1,  # Only need ground state
        )

        self._eigenvalue_ground = eigenvalues[0]
        self._psi_ground = eigenvectors[:, 0]

    def set_ground_state(self, psi: np.ndarray, eigenvalue: float = 0.0):
        """Manually set the ground state wavefunction.

        Args:
            psi: Wavefunction on grid (normalized)
            eigenvalue: Ground state energy (J)
        """
        self._psi_ground = psi
        self._eigenvalue_ground = eigenvalue

    @property
    def psi_ground(self) -> np.ndarray:
        if self._psi_ground is None:
            raise RuntimeError(
                "Ground state not set. Call solve_ground_state() first."
            )
        return self._psi_ground

    @property
    def ground_state_energy(self) -> float:
        if self._eigenvalue_ground is None:
            return 0.0
        return self._eigenvalue_ground

    def run(self, verbose: bool = False) -> EnsembleResult:
        """Run ensemble of trajectories.

        Args:
            verbose: Print progress information

        Returns:
            EnsembleResult with all trajectory data
        """
        # Ensure ground state is available
        if self._psi_ground is None:
            if verbose:
                print("Solving vibrational problem for ground state...")
            self.solve_ground_state()

        # Sample initial conditions
        if verbose:
            print(
                f"Sampling initial conditions: {self.config.sampling.npoints_x} x "
                f"{self.config.sampling.npoints_mom} = "
                f"{self.config.sampling.npoints_x * self.config.sampling.npoints_mom} trajectories"
            )

        x_init, p_init = create_initial_conditions(
            self.x_grid,
            self._psi_ground,
            self.config.sampling.npoints_x,
            self.config.sampling.npoints_mom,
            self.config.sampling.mode,
            compatibility_mode=self.config.sampling.compatibility_mode,
        )

        if verbose:
            print(
                "Sampled positions and momenta"
            )
            for x,p in zip(x_init, p_init):
                print(f"{x}, {p}") 

        n_traj = len(x_init)

        # Run trajectories
        trajectories = []
        for i, (x0, p0) in enumerate(zip(x_init, p_init)):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Running trajectory {i + 1}/{n_traj}")

            v0 = p0 / self.mass_SI

            time, x, v, a = run_trajectory(
                x0=x0,
                v0=v0,
                force_func=self.pes.force,
                mass=self.mass_SI,
                dt=self.dt_SI,
                nsteps=self.config.time.nsteps,
            )

            trajectories.append(
                TrajectoryResult(time=time, x=x, v=v, a=a, x0=x0, p0=p0)
            )

        if verbose:
            print(f"Completed {n_traj} trajectories")

        return EnsembleResult(
            trajectories=trajectories,
            config=self.config,
            ground_state_energy=self._eigenvalue_ground,
            x_grid=self.x_grid,
            psi_ground=self._psi_ground,
            full_config=self.full_config,
        )

    def save_results(
        self,
        result: EnsembleResult,
        output_dir: Path,
        units: str = "user",
    ) -> None:
        """Save results to files.

        Args:
            result: EnsembleResult from run()
            output_dir: Directory for output files
            units: "SI" or "user" (Angstrom/fs)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        basename = self.config.outfile

        # Save ground state
        write_eigenstate(
            output_dir / f"{basename}_ground_state.dat",
            result.x_grid,
            result.psi_ground,
            result.ground_state_energy,
            units=units,
        )

        # Save initial position distribution
        psi_squared = np.abs(result.psi_ground) ** 2
        write_distribution(
            output_dir / f"{basename}_initial_distribution.dat",
            result.x_grid,
            psi_squared,
            units=units,
        )

        # Save each trajectory
        for i, traj in enumerate(result.trajectories):
            write_trajectory(
                output_dir / f"{basename}_traj_{i:04d}.dat",
                traj.time,
                traj.x,
                traj.v,
                units=units,
            )
