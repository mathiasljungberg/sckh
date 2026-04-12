"""Configuration management for 2D dynamics with YAML support."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

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
    position_units: str = "angstrom"  # "angstrom" or "bohr" for positions
    energy_units: str = "hartree"  # "hartree" or "ev" for energy
    index_order: str = "C"  # "C" (x2 fast) or "F" (x1 fast) for data ordering
    outfile: str = "dynamics_2d_out"  # Output file basename


@dataclass
class InterpolationConfig2D:
    """Configuration for 2D surface interpolation (PES + dipole evaluation).

    Attributes:
        dipole_mode: "DIPOLE", "FC", or "DIPOLE_X0"
        pes_final_list: List of paths to final state PES files
        dipole_final_list: List of paths to final state dipole files
        pes_intermediate: Optional path to intermediate PES (if different from dynamics)
        pes_lp_corr: Optional path to energy correction PES
        dipole_initial: Optional path to initial state dipole
        trajectory_files: Optional list of pre-computed trajectory files
        compatibility_mode: "standard" or "fortran"
        dipole_components: Number of dipole components in file (1 or 3)
    """

    dipole_mode: str = "DIPOLE"  # "DIPOLE", "FC", or "DIPOLE_X0"
    pes_final_list: List[Path] = field(default_factory=list)
    dipole_final_list: List[Path] = field(default_factory=list)
    pes_intermediate: Optional[Path] = None
    pes_lp_corr: Optional[Path] = None
    dipole_initial: Optional[Path] = None
    trajectory_files: Optional[List[Path]] = None
    compatibility_mode: str = "standard"
    dipole_components: int = 3  # 1 or 3


# Backwards compatibility
SpectrumConfig2D = InterpolationConfig2D


from sckh.config import FullConfig



def load_config(yaml_path: Path) -> FullConfig:
    """Load configuration from YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        FullConfig object with dynamics2d and optional spectrum parameters

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
          position_units: "angstrom"  # Optional, default "angstrom"
          energy_units: "hartree"     # Optional, default "hartree"
          index_order: "C"            # Optional, default "C"

        spectrum:  # Optional
          gamma_fwhm: 0.18
          dipole_mode: "DIPOLE"
          pes_final_list:
            - "pes_final_1_2d.dat"
            - "pes_final_2_2d.dat"
          dipole_final_list:
            - "dipole_1_2d.dat"
            - "dipole_2_2d.dat"
          compatibility_mode: "standard"
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
        position_units=dyn_data.get("position_units", "angstrom"),
        energy_units=dyn_data.get("energy_units", "hartree"),
        index_order=dyn_data.get("index_order", "C"),
        outfile=dyn_data.get("outfile", "dynamics_2d_out"),
    )

    # Parse optional interpolation config (legacy key: "spectrum")
    interp_config = None
    spec_config = None
    interp_data = data.get("interpolation2d") or data.get("spectrum")
    if interp_data:
        interp_config = InterpolationConfig2D(
            dipole_mode=interp_data.get("dipole_mode", "DIPOLE"),
            pes_final_list=[
                Path(p) for p in interp_data.get("pes_final_list", [])
            ],
            dipole_final_list=[
                Path(p) for p in interp_data.get("dipole_final_list", [])
            ],
            pes_intermediate=Path(interp_data["pes_intermediate"]) if interp_data.get("pes_intermediate") else None,
            pes_lp_corr=Path(interp_data["pes_lp_corr"]) if interp_data.get("pes_lp_corr") else None,
            dipole_initial=Path(interp_data["dipole_initial"]) if interp_data.get("dipole_initial") else None,
            trajectory_files=[
                Path(p) for p in interp_data.get("trajectory_files", [])
            ] if interp_data.get("trajectory_files") else None,
            compatibility_mode=interp_data.get("compatibility_mode", "standard"),
            dipole_components=interp_data.get("dipole_components", 3),
        )

    # Parse spectrum config (gamma_fwhm can be in "spectrum" section or legacy combined)
    spec_data = data.get("spectrum", {})
    if spec_data.get("gamma_fwhm") is not None:
        from sckh.config import SpectrumConfig
        spec_config = SpectrumConfig(
            gamma_fwhm=spec_data["gamma_fwhm"],
            dt=spec_data.get("dt"),
        )

    return FullConfig(
        dynamics2d=dynamics_config,
        interpolation2d=interp_config,
        spectrum=spec_config,
    )


def save_config(config: FullConfig, yaml_path: Path) -> None:
    """Save configuration to YAML file.

    Args:
        config: FullConfig object
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
            "position_units": dyn.position_units,
            "energy_units": dyn.energy_units,
            "index_order": dyn.index_order,
            "outfile": dyn.outfile,
        },
    }

    # Add interpolation config if present
    if config.interpolation2d is not None:
        interp = config.interpolation2d
        interp_data = {
            "dipole_mode": interp.dipole_mode,
            "pes_final_list": [str(p) for p in interp.pes_final_list],
            "dipole_final_list": [str(p) for p in interp.dipole_final_list],
            "compatibility_mode": interp.compatibility_mode,
            "dipole_components": interp.dipole_components,
        }
        if interp.pes_intermediate:
            interp_data["pes_intermediate"] = str(interp.pes_intermediate)
        if interp.pes_lp_corr:
            interp_data["pes_lp_corr"] = str(interp.pes_lp_corr)
        if interp.dipole_initial:
            interp_data["dipole_initial"] = str(interp.dipole_initial)
        if interp.trajectory_files:
            interp_data["trajectory_files"] = [str(p) for p in interp.trajectory_files]
        data["interpolation2d"] = interp_data

    # Add spectrum config if present
    if config.spectrum is not None:
        spec_data = {"gamma_fwhm": config.spectrum.gamma_fwhm}
        if config.spectrum.dt is not None:
            spec_data["dt"] = config.spectrum.dt
        data["spectrum"] = spec_data

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
