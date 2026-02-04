"""Configuration dataclasses for KH 1D spectrum calculation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml


@dataclass
class GridConfig:
    """Configuration for spatial grid.

    Attributes:
        start: Starting position in Angstrom
        dx: Grid spacing in Angstrom
        npoints: Number of grid points (should be odd for Fourier grid)
    """

    start: float
    dx: float
    npoints: int

    def __post_init__(self):
        if self.npoints % 2 == 0:
            raise ValueError(
                f"npoints must be odd for Fourier grid method, got {self.npoints}"
            )


@dataclass
class BroadeningConfig:
    """Configuration for spectral broadening.

    Attributes:
        gamma_fwhm: Full width at half maximum (FWHM) in eV
    """

    gamma_fwhm: float

    @property
    def gamma_hwhm(self) -> float:
        """Half width at half maximum (HWHM) in eV."""
        return self.gamma_fwhm / 2.0


@dataclass
class FrequencyGridConfig:
    """Configuration for output frequency grid.

    Attributes:
        omega_start: Starting frequency in eV
        omega_end: Ending frequency in eV
        n_omega: Number of frequency points
    """

    omega_start: float
    omega_end: float
    n_omega: int

    def get_omega_array(self):
        """Generate the frequency array."""
        import numpy as np

        return np.linspace(self.omega_start, self.omega_end, self.n_omega)


@dataclass
class KHConfig:
    """Complete configuration for Kramers-Heisenberg calculation.

    Attributes:
        mu: Reduced mass in atomic mass units (amu)
        grid: Spatial grid configuration
        broadening: Spectral broadening configuration
        frequency: Output frequency grid configuration
        pes_initial: Path to initial state PES file
        pes_intermediate: Path to intermediate (core-excited) state PES file
        pes_final_list: List of paths to final state PES files
        dipole_final_list: List of paths to dipole files for each final state
        dipole_mode: Dipole calculation mode ("FC", "DIPOLE", or "DIPOLE_X0")
        n_vib_states: Number of vibrational states to compute (default: all)
        energy_column_initial: Column index for initial state energy (1-indexed)
        energy_column_intermediate: Column index for intermediate state energy
        energy_column_final: Column index for final state energy
    """

    mu: float  # amu
    grid: GridConfig
    broadening: BroadeningConfig
    frequency: FrequencyGridConfig
    pes_initial: Path
    pes_intermediate: Path
    pes_final_list: List[Path]
    dipole_final_list: List[Path]
    dipole_mode: str = "FC"
    n_vib_states: Optional[int] = None
    pes_lp_corr: Optional[Path] = None
    energy_column_initial: int = 1
    energy_column_intermediate: int = 1
    energy_column_final: int = 1

    def __post_init__(self):
        # Convert paths to Path objects
        self.pes_initial = Path(self.pes_initial)
        self.pes_intermediate = Path(self.pes_intermediate)
        self.pes_final_list = [Path(p) for p in self.pes_final_list]
        self.dipole_final_list = [Path(p) for p in self.dipole_final_list]
        if self.pes_lp_corr is not None:
            self.pes_lp_corr = Path(self.pes_lp_corr)

        # Validate dipole mode
        valid_modes = ("FC", "DIPOLE", "DIPOLE_X0")
        if self.dipole_mode.upper() not in valid_modes:
            raise ValueError(
                f"dipole_mode must be one of {valid_modes}, got '{self.dipole_mode}'"
            )
        self.dipole_mode = self.dipole_mode.upper()

        # Validate list lengths match
        if len(self.pes_final_list) != len(self.dipole_final_list):
            raise ValueError(
                f"pes_final_list and dipole_final_list must have same length: "
                f"{len(self.pes_final_list)} vs {len(self.dipole_final_list)}"
            )

    @property
    def n_final_states(self) -> int:
        """Number of final electronic states."""
        return len(self.pes_final_list)


def load_config(yaml_path: Path) -> KHConfig:
    """Load KH configuration from YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        KHConfig object with all parameters
    """
    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    data = raw.get("kh_1d", raw)

    return KHConfig(
        mu=data["mu"],
        grid=GridConfig(**data["grid"]),
        broadening=BroadeningConfig(**data["broadening"]),
        frequency=FrequencyGridConfig(**data["frequency"]),
        pes_initial=Path(data["pes_initial"]),
        pes_intermediate=Path(data["pes_intermediate"]),
        pes_final_list=[Path(p) for p in data["pes_final_list"]],
        dipole_final_list=[Path(p) for p in data["dipole_final_list"]],
        dipole_mode=data.get("dipole_mode", "FC"),
        n_vib_states=data.get("n_vib_states"),
        pes_lp_corr=Path(data["pes_lp_corr"]) if data.get("pes_lp_corr") else None,
        energy_column_initial=data.get("energy_column_initial", 1),
        energy_column_intermediate=data.get("energy_column_intermediate", 1),
        energy_column_final=data.get("energy_column_final", 1),
    )


def save_config(config: KHConfig, yaml_path: Path) -> None:
    """Save KH configuration to YAML file.

    Args:
        config: KHConfig object
        yaml_path: Path to output YAML file
    """
    data = {
        "kh_1d": {
            "mu": config.mu,
            "grid": {
                "start": config.grid.start,
                "dx": config.grid.dx,
                "npoints": config.grid.npoints,
            },
            "broadening": {
                "gamma_fwhm": config.broadening.gamma_fwhm,
            },
            "frequency": {
                "omega_start": config.frequency.omega_start,
                "omega_end": config.frequency.omega_end,
                "n_omega": config.frequency.n_omega,
            },
            "pes_initial": str(config.pes_initial),
            "pes_intermediate": str(config.pes_intermediate),
            "pes_final_list": [str(p) for p in config.pes_final_list],
            "dipole_final_list": [str(p) for p in config.dipole_final_list],
            "dipole_mode": config.dipole_mode,
            "n_vib_states": config.n_vib_states,
            "pes_lp_corr": str(config.pes_lp_corr) if config.pes_lp_corr else None,
            "energy_column_initial": config.energy_column_initial,
            "energy_column_intermediate": config.energy_column_intermediate,
            "energy_column_final": config.energy_column_final,
        }
    }

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
