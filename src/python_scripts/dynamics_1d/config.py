"""Configuration management with YAML support."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass
class GridConfig:
    """DVR grid configuration."""

    start: float  # Starting point in Angstroms
    dx: float  # Grid spacing in Angstroms
    npoints: int  # Number of grid points


@dataclass
class TimeConfig:
    """Time stepping configuration."""

    dt: float  # Time step in femtoseconds
    nsteps: int  # Number of time steps


@dataclass
class SamplingConfig:
    """Initial condition sampling configuration.

    Attributes:
        mode: Sampling mode (1=even, 2=random)
        npoints_x: Number of position samples
        npoints_mom: Number of momentum samples
        compatibility_mode: Algorithm mode for sampling
            - "standard": Proper CDF integration with linear interpolation
            - "fortran": Match Fortran implementation (spline, cumsum, step-lookup)
    """

    mode: int = 1  # 1=even, 2=random
    npoints_x: int = 10  # Number of position samples
    npoints_mom: int = 10  # Number of momentum samples
    compatibility_mode: str = "standard"


@dataclass
class DynamicsConfig:
    """Complete dynamics configuration."""

    # Physical parameters
    mu: float  # Reduced mass in amu

    # Grid
    grid: GridConfig

    # Time stepping
    time: TimeConfig

    # Sampling
    sampling: SamplingConfig

    # PES files
    pes_initial: Path  # Initial state PES
    pes_dynamics: Path  # PES for dynamics (intermediate state)

    # Optional
    units: str = "angstrom"  # "angstrom" or "bohr"
    outfile: str = "dynamics_out"  # Output file basename


def load_config(yaml_path: Path) -> DynamicsConfig:
    """Load configuration from YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        DynamicsConfig object with all parameters
    """
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    return DynamicsConfig(
        mu=data["mu"],
        grid=GridConfig(**data["grid"]),
        time=TimeConfig(**data["time"]),
        sampling=SamplingConfig(**data.get("sampling", {})),
        pes_initial=Path(data["pes_initial"]),
        pes_dynamics=Path(data["pes_dynamics"]),
        units=data.get("units", "angstrom"),
        outfile=data.get("outfile", "dynamics_out"),
    )


def save_config(config: DynamicsConfig, yaml_path: Path) -> None:
    """Save configuration to YAML file.

    Args:
        config: DynamicsConfig object
        yaml_path: Path to output YAML file
    """
    data = {
        "mu": config.mu,
        "grid": {
            "start": config.grid.start,
            "dx": config.grid.dx,
            "npoints": config.grid.npoints,
        },
        "time": {
            "dt": config.time.dt,
            "nsteps": config.time.nsteps,
        },
        "sampling": {
            "mode": config.sampling.mode,
            "npoints_x": config.sampling.npoints_x,
            "npoints_mom": config.sampling.npoints_mom,
            "compatibility_mode": config.sampling.compatibility_mode,
        },
        "pes_initial": str(config.pes_initial),
        "pes_dynamics": str(config.pes_dynamics),
        "units": config.units,
        "outfile": config.outfile,
    }

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
