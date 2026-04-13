"""Configuration for SCKH spectrum computation and unified FullConfig.

FullConfig ties together dynamics, interpolation, and spectrum configs.
Dimension-specific configs (DynamicsConfig, InterpolationConfig, etc.)
are imported lazily to avoid circular dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from dynamics_1d.config import DynamicsConfig
    from dynamics_1d.spectrum_config import InterpolationConfig


@dataclass
class SpectrumConfig:
    """Configuration for SCKH spectrum computation.

    Attributes:
        gamma_fwhm: Lifetime broadening FWHM in eV
        dt: Optional time step for trajectory interpolation (seconds).
            If None, the trajectory's own time step is used.
            Smaller dt gives wider spectral range (FFT range ~ 1/(2*dt)).
    """

    gamma_fwhm: float  # eV
    dt: Optional[float] = None  # seconds

    @property
    def gamma_hwhm(self) -> float:
        """Half-width at half maximum in eV."""
        return self.gamma_fwhm / 2.0


@dataclass
class SCKHRunConfig:
    """Configuration for a pure-SCKH workflow (no dynamics, no interpolation).

    Used when the spectrum is computed directly from pre-computed SCKH
    trajectory files (e.g. produced by the Fortran ``sckh_main`` binary
    or by prior Python dynamics + interpolation runs).

    Attributes:
        gamma_fwhm: Lifetime broadening FWHM in eV.
        trajectory_files: List of SCKH trajectory file paths.
        outfile: Base name for output spectrum files (no extension).
            Produces ``{outfile}.dat`` plus ``{outfile}_final_*.dat``.
    """

    gamma_fwhm: float
    trajectory_files: List[Path] = field(default_factory=list)
    outfile: str = "sckh_spectrum"


@dataclass
class FullConfig:
    """Unified configuration for dynamics, interpolation, and spectrum.

    Supports both 1D and 2D dynamics via dimension-tagged fields.
    Exactly one of dynamics1d/dynamics2d should be set, and the
    matching interpolation config if surface interpolation is needed.

    For pure-SCKH workflows (spectrum from pre-computed trajectories),
    use the ``sckh`` field instead.
    """

    dynamics1d: Optional[DynamicsConfig] = None
    dynamics2d: Optional[object] = None  # DynamicsConfig2D
    interpolation1d: Optional[InterpolationConfig] = None
    interpolation2d: Optional[object] = None  # InterpolationConfig2D
    spectrum: Optional[SpectrumConfig] = None
    sckh: Optional[SCKHRunConfig] = None


def load_full_config(yaml_path: Path) -> FullConfig:
    """Load configuration from YAML file.

    Supports two YAML formats:

    New format (recommended)::

        dynamics1d:
          mu: 1.0078825
          ...
        interpolation1d:
          dipole_mode: "DIPOLE"
          pes_final_list: [...]
          ...
        spectrum:
          gamma_fwhm: 0.18

    Legacy format (backwards compatible)::

        dynamics:
          mu: 1.0078825
          ...
        spectrum:
          gamma_fwhm: 0.18
          dipole_mode: "DIPOLE"
          pes_final_list: [...]
          ...

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        FullConfig object
    """
    from dynamics_1d.config import DynamicsConfig, GridConfig, TimeConfig, SamplingConfig
    from dynamics_1d.spectrum_config import InterpolationConfig

    yaml_path = Path(yaml_path).resolve()
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    base_dir = yaml_path.parent

    # Pure-SCKH workflow: no dynamics needed, just trajectory files + gamma.
    if "sckh" in data:
        return FullConfig(sckh=_parse_sckh(data["sckh"], base_dir))

    # Detect format: new (dynamics1d/dynamics2d) vs legacy (dynamics)
    if "dynamics1d" in data:
        dyn_data = data["dynamics1d"]
        interp_data = data.get("interpolation1d", {})
        spec_data = data.get("spectrum", {})

        compat_mode = interp_data.get(
            "compatibility_mode",
            dyn_data.get("sampling", {}).get("compatibility_mode", "standard"),
        )

        return FullConfig(
            dynamics1d=_parse_dynamics1d(dyn_data, compat_mode, base_dir),
            interpolation1d=_parse_interpolation(interp_data, compat_mode, base_dir) if interp_data else None,
            spectrum=_parse_spectrum(spec_data) if spec_data.get("gamma_fwhm") else None,
        )

    elif "dynamics" in data:
        # Legacy format: dynamics + spectrum (combined interpolation + spectrum)
        dyn_data = data["dynamics"]
        combined = data.get("spectrum", {})

        compat_mode = combined.get(
            "compatibility_mode",
            dyn_data.get("sampling", {}).get("compatibility_mode", "standard"),
        )

        return FullConfig(
            dynamics1d=_parse_dynamics1d(dyn_data, compat_mode, base_dir),
            interpolation1d=_parse_interpolation(combined, compat_mode, base_dir),
            spectrum=_parse_spectrum(combined) if combined.get("gamma_fwhm") else None,
        )

    else:
        raise ValueError(
            "YAML must contain one of: 'sckh', 'dynamics1d', 'dynamics2d', "
            "or 'dynamics' section"
        )


def save_full_config(config: FullConfig, yaml_path: Path) -> None:
    """Save configuration to YAML file (new format).

    Args:
        config: FullConfig object
        yaml_path: Path to output YAML file
    """
    data = {}

    if config.dynamics1d is not None:
        dyn = config.dynamics1d
        data["dynamics1d"] = {
            "mu": dyn.mu,
            "grid": {
                "start": dyn.grid.start,
                "dx": dyn.grid.dx,
                "npoints": dyn.grid.npoints,
            },
            "time": {
                "dt": dyn.time.dt,
                "nsteps": dyn.time.nsteps,
            },
            "sampling": {
                "mode": dyn.sampling.mode,
                "npoints_x": dyn.sampling.npoints_x,
                "npoints_mom": dyn.sampling.npoints_mom,
            },
            "pes_initial": str(dyn.pes_initial),
            "pes_dynamics": str(dyn.pes_dynamics),
            "units": dyn.units,
            "outfile": dyn.outfile,
        }

    if config.interpolation1d is not None:
        interp = config.interpolation1d
        interp_data = {
            "dipole_mode": interp.dipole_mode,
            "pes_final_list": [str(p) for p in interp.pes_final_list],
            "dipole_final_list": [str(p) for p in interp.dipole_final_list],
            "compatibility_mode": interp.compatibility_mode,
        }
        if interp.dipole_initial:
            interp_data["dipole_initial"] = str(interp.dipole_initial)
        if interp.trajectory_files:
            interp_data["trajectory_files"] = [
                str(p) for p in interp.trajectory_files
            ]
        if interp.pes_lp_corr:
            interp_data["pes_lp_corr"] = str(interp.pes_lp_corr)
        data["interpolation1d"] = interp_data

    if config.spectrum is not None:
        spec_data = {"gamma_fwhm": config.spectrum.gamma_fwhm}
        if config.spectrum.dt is not None:
            spec_data["dt"] = config.spectrum.dt
        data["spectrum"] = spec_data

    if config.sckh is not None:
        data["sckh"] = {
            "gamma_fwhm": config.sckh.gamma_fwhm,
            "trajectory_files": [str(p) for p in config.sckh.trajectory_files],
            "outfile": config.sckh.outfile,
        }

    with open(yaml_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def _resolve(p, base_dir: Path) -> Path:
    """Resolve a YAML path string against ``base_dir`` if relative."""
    path = Path(p)
    return path if path.is_absolute() else (base_dir / path)


def _parse_dynamics1d(dyn_data: dict, compat_mode: str, base_dir: Path):
    """Parse 1D dynamics config from YAML data.

    Relative file paths are resolved against ``base_dir`` so the YAML is
    portable regardless of the caller's working directory.
    """
    from dynamics_1d.config import DynamicsConfig, GridConfig, TimeConfig, SamplingConfig

    sampling_data = dyn_data.get("sampling", {})
    sampling_data["compatibility_mode"] = compat_mode

    return DynamicsConfig(
        mu=dyn_data["mu"],
        grid=GridConfig(**dyn_data["grid"]),
        time=TimeConfig(**dyn_data["time"]),
        sampling=SamplingConfig(**sampling_data),
        pes_initial=_resolve(dyn_data["pes_initial"], base_dir),
        pes_dynamics=_resolve(dyn_data["pes_dynamics"], base_dir),
        units=dyn_data.get("units", "angstrom"),
        outfile=dyn_data.get("outfile", "dynamics_out"),
    )


def _parse_interpolation(interp_data: dict, compat_mode: str, base_dir: Path):
    """Parse interpolation config from YAML data.

    Relative file paths are resolved against ``base_dir``.
    """
    from dynamics_1d.spectrum_config import InterpolationConfig

    return InterpolationConfig(
        dipole_mode=interp_data.get("dipole_mode", "DIPOLE"),
        pes_final_list=[
            _resolve(p, base_dir) for p in interp_data.get("pes_final_list", [])
        ],
        dipole_final_list=[
            _resolve(p, base_dir) for p in interp_data.get("dipole_final_list", [])
        ],
        dipole_initial=(
            _resolve(interp_data["dipole_initial"], base_dir)
            if interp_data.get("dipole_initial")
            else None
        ),
        trajectory_files=(
            [_resolve(p, base_dir) for p in interp_data["trajectory_files"]]
            if interp_data.get("trajectory_files")
            else None
        ),
        pes_lp_corr=(
            _resolve(interp_data["pes_lp_corr"], base_dir)
            if interp_data.get("pes_lp_corr")
            else None
        ),
        compatibility_mode=compat_mode,
    )


def _parse_spectrum(spec_data: dict) -> SpectrumConfig:
    """Parse spectrum config from YAML data."""
    return SpectrumConfig(
        gamma_fwhm=spec_data["gamma_fwhm"],
        dt=spec_data.get("dt"),
    )


def _parse_sckh(sckh_data: dict, base_dir: Path) -> SCKHRunConfig:
    """Parse pure-SCKH config from YAML data.

    Trajectory file paths are resolved relative to ``base_dir`` (the
    directory containing the YAML file) so that absolute paths are
    available to callers regardless of their working directory.
    """
    return SCKHRunConfig(
        gamma_fwhm=sckh_data["gamma_fwhm"],
        trajectory_files=[
            _resolve(p, base_dir) for p in sckh_data.get("trajectory_files", [])
        ],
        outfile=sckh_data.get("outfile", "sckh_spectrum"),
    )
