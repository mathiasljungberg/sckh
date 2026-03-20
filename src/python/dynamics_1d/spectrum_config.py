"""Configuration for spectrum calculation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml

from .config import DynamicsConfig, GridConfig, TimeConfig, SamplingConfig


@dataclass
class SpectrumConfig:
    """Spectrum calculation configuration.

    Attributes:
        gamma_fwhm: Lifetime broadening FWHM in eV
        dipole_mode: "DIPOLE" (full dipole surface), "FC" (Franck-Condon, constant=1),
                     or "DIPOLE_X0" (dipole frozen at equilibrium)
        pes_final_list: List of final state PES file paths
        dipole_final_list: List of dipole file paths (one per final state)
        dipole_initial: Optional initial state dipole file (for D_ni)
        trajectory_files: Optional list of pre-computed trajectory files
        pes_lp_corr: Optional PES file for energy correction. When provided, all
                     final state PES energies are shifted so that the first final
                     state matches this PES (e.g., lone-pair corrected energies).
        compatibility_mode: Algorithm mode for numerical methods.
            - "standard": Use best-practice algorithms (proper CDF integration,
              direct equilibrium position lookup for E_mean)
            - "fortran": Match Fortran implementation exactly (for validation)
    """

    gamma_fwhm: float  # eV
    dipole_mode: str = "DIPOLE"
    pes_final_list: List[Path] = field(default_factory=list)
    dipole_final_list: List[Path] = field(default_factory=list)
    dipole_initial: Optional[Path] = None
    trajectory_files: Optional[List[Path]] = None
    pes_lp_corr: Optional[Path] = None
    compatibility_mode: str = "standard"

    def __post_init__(self):
        """Validate configuration."""
        valid_modes = ["DIPOLE", "FC", "DIPOLE_X0"]
        if self.dipole_mode.upper() not in valid_modes:
            raise ValueError(
                f"dipole_mode must be one of {valid_modes}, got {self.dipole_mode}"
            )
        self.dipole_mode = self.dipole_mode.upper()

        valid_compat_modes = ["standard", "fortran"]
        if self.compatibility_mode.lower() not in valid_compat_modes:
            raise ValueError(
                f"compatibility_mode must be one of {valid_compat_modes}, "
                f"got {self.compatibility_mode}"
            )
        self.compatibility_mode = self.compatibility_mode.lower()

        if self.dipole_mode == "DIPOLE" and len(self.dipole_final_list) == 0:
            raise ValueError(
                "dipole_final_list required when dipole_mode is 'DIPOLE'"
            )

        if (
            self.dipole_mode == "DIPOLE"
            and len(self.dipole_final_list) != len(self.pes_final_list)
        ):
            raise ValueError(
                f"dipole_final_list ({len(self.dipole_final_list)}) must match "
                f"pes_final_list ({len(self.pes_final_list)})"
            )

    @property
    def gamma_hwhm(self) -> float:
        """Half-width at half maximum in eV."""
        return self.gamma_fwhm / 2.0

    @property
    def n_final_states(self) -> int:
        """Number of final states."""
        return len(self.pes_final_list)


@dataclass
class FullConfig:
    """Combined configuration for dynamics and spectrum calculation.

    This class holds both dynamics and spectrum configurations,
    allowing the full workflow from trajectory generation to spectrum output.
    """

    dynamics: DynamicsConfig
    spectrum: SpectrumConfig


def load_full_config(yaml_path: Path) -> FullConfig:
    """Load combined dynamics and spectrum configuration from YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        FullConfig object with both dynamics and spectrum parameters

    Example YAML format:
        dynamics:
          mu: 1.0078825
          grid: {start: 0.5, dx: 0.025, npoints: 77}
          time: {dt: 0.1, nsteps: 512}
          sampling: {mode: 1, npoints_x: 10, npoints_mom: 10}
          pes_initial: "pes_initial.dat"
          pes_dynamics: "pes_intermediate.dat"

        spectrum:
          gamma_fwhm: 0.18
          dipole_mode: "DIPOLE"
          pes_final_list:
            - "pes_final_1.dat"
            - "pes_final_2.dat"
          dipole_final_list:
            - "dipole_1.dat"
            - "dipole_2.dat"
    """
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    # Parse dynamics config
    dyn_data = data["dynamics"]
    spec_data = data["spectrum"]

    # Get compatibility_mode - spectrum setting takes precedence, then dynamics sampling
    compat_mode = spec_data.get(
        "compatibility_mode",
        dyn_data.get("sampling", {}).get("compatibility_mode", "standard")
    )

    # Build sampling config with unified compatibility_mode
    sampling_data = dyn_data.get("sampling", {})
    sampling_data["compatibility_mode"] = compat_mode

    dynamics_config = DynamicsConfig(
        mu=dyn_data["mu"],
        grid=GridConfig(**dyn_data["grid"]),
        time=TimeConfig(**dyn_data["time"]),
        sampling=SamplingConfig(**sampling_data),
        pes_initial=Path(dyn_data["pes_initial"]),
        pes_dynamics=Path(dyn_data["pes_dynamics"]),
        units=dyn_data.get("units", "angstrom"),
        outfile=dyn_data.get("outfile", "dynamics_out"),
    )

    # Parse spectrum config (spec_data already defined above)
    spectrum_config = SpectrumConfig(
        gamma_fwhm=spec_data["gamma_fwhm"],
        dipole_mode=spec_data.get("dipole_mode", "DIPOLE"),
        pes_final_list=[Path(p) for p in spec_data.get("pes_final_list", [])],
        dipole_final_list=[Path(p) for p in spec_data.get("dipole_final_list", [])],
        dipole_initial=(
            Path(spec_data["dipole_initial"])
            if spec_data.get("dipole_initial")
            else None
        ),
        trajectory_files=(
            [Path(p) for p in spec_data["trajectory_files"]]
            if spec_data.get("trajectory_files")
            else None
        ),
        pes_lp_corr=(
            Path(spec_data["pes_lp_corr"])
            if spec_data.get("pes_lp_corr")
            else None
        ),
        compatibility_mode=compat_mode,
    )

    return FullConfig(dynamics=dynamics_config, spectrum=spectrum_config)


def save_full_config(config: FullConfig, yaml_path: Path) -> None:
    """Save combined configuration to YAML file.

    Args:
        config: FullConfig object
        yaml_path: Path to output YAML file
    """
    data = {
        "dynamics": {
            "mu": config.dynamics.mu,
            "grid": {
                "start": config.dynamics.grid.start,
                "dx": config.dynamics.grid.dx,
                "npoints": config.dynamics.grid.npoints,
            },
            "time": {
                "dt": config.dynamics.time.dt,
                "nsteps": config.dynamics.time.nsteps,
            },
            "sampling": {
                "mode": config.dynamics.sampling.mode,
                "npoints_x": config.dynamics.sampling.npoints_x,
                "npoints_mom": config.dynamics.sampling.npoints_mom,
            },
            "pes_initial": str(config.dynamics.pes_initial),
            "pes_dynamics": str(config.dynamics.pes_dynamics),
            "units": config.dynamics.units,
            "outfile": config.dynamics.outfile,
        },
        "spectrum": {
            "gamma_fwhm": config.spectrum.gamma_fwhm,
            "dipole_mode": config.spectrum.dipole_mode,
            "pes_final_list": [str(p) for p in config.spectrum.pes_final_list],
            "dipole_final_list": [str(p) for p in config.spectrum.dipole_final_list],
            "compatibility_mode": config.spectrum.compatibility_mode,
        },
    }

    if config.spectrum.dipole_initial:
        data["spectrum"]["dipole_initial"] = str(config.spectrum.dipole_initial)

    if config.spectrum.trajectory_files:
        data["spectrum"]["trajectory_files"] = [
            str(p) for p in config.spectrum.trajectory_files
        ]

    if config.spectrum.pes_lp_corr:
        data["spectrum"]["pes_lp_corr"] = str(config.spectrum.pes_lp_corr)

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
