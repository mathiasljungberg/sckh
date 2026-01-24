"""Configuration management for 2D dynamics with YAML support."""

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class GridConfig2D:
    """DVR grid configuration for one dimension."""

    start: float  # Starting point in Angstroms
    dx: float  # Grid spacing in Angstroms
    npoints: int  # Number of grid points


@dataclass
class TimeConfig:
    """Time stepping configuration.

    Reuses the same structure as dynamics_1d.
    """

    dt: float  # Time step in femtoseconds
    nsteps: int  # Number of time steps


@dataclass
class SamplingConfig2D:
    """Initial condition sampling configuration for 2D.

    Attributes:
        mode: Sampling mode (1=even, 2=random)
        npoints_x1: Number of position samples for x1
        npoints_x2: Number of position samples for x2
        npoints_p1: Number of momentum samples for p1
        npoints_p2: Number of momentum samples for p2
    """

    mode: int = 1  # 1=even, 2=random
    npoints_x1: int = 10  # Number of position samples for x1
    npoints_x2: int = 10  # Number of position samples for x2
    npoints_p1: int = 10  # Number of momentum samples for p1
    npoints_p2: int = 10  # Number of momentum samples for p2


@dataclass
class DynamicsConfig2D:
    """Complete 2D dynamics configuration.

    Supports two degrees of freedom with potentially different masses.
    Uses product wavefunction approximation for initial state:
    psi(x1, x2) ~ psi1(x1) * psi2(x2)

    The initial state wavefunctions are computed by cutting 1D slices
    from the 2D initial PES at the equilibrium geometry:
        V1(x1) = V(x1, x2_eq)
        V2(x2) = V(x1_eq, x2)
    """

    # Physical parameters
    mu1: float  # Reduced mass for x1 in amu
    mu2: float  # Reduced mass for x2 in amu

    # Grids for each coordinate
    grid_x1: GridConfig2D
    grid_x2: GridConfig2D

    # Time stepping
    time: TimeConfig

    # Sampling
    sampling: SamplingConfig2D

    # PES files
    pes_dynamics: Path  # 2D PES for dynamics (intermediate state)
    pes_initial: Path  # 2D PES for initial state (sliced at equilibrium)

    # Optional
    units: str = "angstrom"  # "angstrom" or "bohr"
    outfile: str = "dynamics_2d_out"  # Output file basename


@dataclass
class FullConfig2D:
    """Combined configuration for 2D dynamics.

    This class wraps DynamicsConfig2D and allows for future extension
    with spectrum calculation (similar to 1D FullConfig).
    """

    dynamics2d: DynamicsConfig2D


def load_config(yaml_path: Path) -> FullConfig2D:
    """Load configuration from YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        FullConfig2D object with dynamics2d parameters

    Example YAML format:
        dynamics2d:
          mu1: 1.0078825
          mu2: 15.999
          grid_x1: {start: 0.5, dx: 0.025, npoints: 77}
          grid_x2: {start: 1.0, dx: 0.025, npoints: 77}
          time: {dt: 0.1, nsteps: 512}
          sampling: {mode: 1, npoints_x1: 10, npoints_x2: 10, npoints_p1: 10, npoints_p2: 10}
          pes_initial: "pes_initial_2d.dat"
          pes_dynamics: "pes_intermediate_2d.dat"
    """
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    dyn_data = data["dynamics2d"]

    dynamics_config = DynamicsConfig2D(
        mu1=dyn_data["mu1"],
        mu2=dyn_data["mu2"],
        grid_x1=GridConfig2D(**dyn_data["grid_x1"]),
        grid_x2=GridConfig2D(**dyn_data["grid_x2"]),
        time=TimeConfig(**dyn_data["time"]),
        sampling=SamplingConfig2D(**dyn_data.get("sampling", {})),
        pes_dynamics=Path(dyn_data["pes_dynamics"]),
        pes_initial=Path(dyn_data["pes_initial"]),
        units=dyn_data.get("units", "angstrom"),
        outfile=dyn_data.get("outfile", "dynamics_2d_out"),
    )

    return FullConfig2D(dynamics2d=dynamics_config)


def save_config(config: FullConfig2D, yaml_path: Path) -> None:
    """Save configuration to YAML file.

    Args:
        config: FullConfig2D object
        yaml_path: Path to output YAML file
    """
    dyn = config.dynamics2d
    data = {
        "dynamics2d": {
            "mu1": dyn.mu1,
            "mu2": dyn.mu2,
            "grid_x1": {
                "start": dyn.grid_x1.start,
                "dx": dyn.grid_x1.dx,
                "npoints": dyn.grid_x1.npoints,
            },
            "grid_x2": {
                "start": dyn.grid_x2.start,
                "dx": dyn.grid_x2.dx,
                "npoints": dyn.grid_x2.npoints,
            },
            "time": {
                "dt": dyn.time.dt,
                "nsteps": dyn.time.nsteps,
            },
            "sampling": {
                "mode": dyn.sampling.mode,
                "npoints_x1": dyn.sampling.npoints_x1,
                "npoints_x2": dyn.sampling.npoints_x2,
                "npoints_p1": dyn.sampling.npoints_p1,
                "npoints_p2": dyn.sampling.npoints_p2,
            },
            "pes_dynamics": str(dyn.pes_dynamics),
            "pes_initial": str(dyn.pes_initial),
            "units": dyn.units,
            "outfile": dyn.outfile,
        },
    }

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
