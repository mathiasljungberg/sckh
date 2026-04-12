"""Configuration for surface interpolation.

InterpolationConfig lives here (1D-specific).
SpectrumConfig, FullConfig, and load/save functions have moved to
sckh.config and are re-exported for backwards compatibility.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from .config import DynamicsConfig, GridConfig, TimeConfig, SamplingConfig


@dataclass
class InterpolationConfig:
    """Configuration for surface interpolation (PES + dipole evaluation).

    Attributes:
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
    def n_final_states(self) -> int:
        """Number of final states."""
        return len(self.pes_final_list)


# Re-exports from sckh.config for backwards compatibility
from sckh.config import SpectrumConfig, FullConfig, load_full_config, save_full_config
