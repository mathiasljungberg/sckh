"""
kh_1d: Kramers-Heisenberg 1D X-ray emission spectrum calculation.

This package provides tools for computing non-resonant X-ray emission
spectra (XES) using the Kramers-Heisenberg formalism on 1D potential
energy surfaces.

Example:
    >>> from kh_1d import KHConfig, KHCalculator, GridConfig, BroadeningConfig, FrequencyGridConfig
    >>> config = KHConfig(
    ...     mu=1.0078,
    ...     grid=GridConfig(start=0.5, dx=0.025, npoints=77),
    ...     broadening=BroadeningConfig(gamma_fwhm=0.3),
    ...     frequency=FrequencyGridConfig(omega_start=517, omega_end=528, n_omega=1000),
    ...     pes_initial="pes_initial.dat",
    ...     pes_intermediate="pes_intermediate.dat",
    ...     pes_final_list=["pes_final_1.dat"],
    ...     dipole_final_list=["dipole_final_1.dat"],
    ... )
    >>> calc = KHCalculator(config)
    >>> result = calc.run()
    >>> result.write_spectrum("spectrum.dat")
"""

from .config import (
    GridConfig,
    BroadeningConfig,
    FrequencyGridConfig,
    KHConfig,
    load_config,
    save_config,
)
from .dipole_matrix import (
    compute_fc_overlap,
    compute_dipole_overlap,
    compute_transition_dipoles_fc,
    compute_transition_dipoles_fc_loop,
    compute_transition_dipoles_full,
    compute_transition_dipoles_full_loop,
    compute_dipole_matrix_elements,
)
from .amplitude import (
    compute_amplitude_F_res,
    compute_amplitude_F_res_loop,
    compute_amplitude_F_nonres,
    compute_amplitude_F_nonres_loop,
    compute_amplitude_F,
    compute_cross_section_from_F,
    compute_XES_nonres,
    compute_XES_nonres_vectorized,
    compute_XES_per_final_state,
    compute_XES_per_final_state_loop,
)
from .spectrum import (
    VibrationalState,
    KHResult,
    KHCalculator,
)

__all__ = [
    # Config
    "GridConfig",
    "BroadeningConfig",
    "FrequencyGridConfig",
    "KHConfig",
    "load_config",
    "save_config",
    # Dipole matrix
    "compute_fc_overlap",
    "compute_dipole_overlap",
    "compute_transition_dipoles_fc",
    "compute_transition_dipoles_fc_loop",
    "compute_transition_dipoles_full",
    "compute_transition_dipoles_full_loop",
    "compute_dipole_matrix_elements",
    # Amplitude
    "compute_amplitude_F_res",
    "compute_amplitude_F_res_loop",
    "compute_amplitude_F_nonres",
    "compute_amplitude_F_nonres_loop",
    "compute_amplitude_F",
    "compute_cross_section_from_F",
    "compute_XES_nonres",
    "compute_XES_nonres_vectorized",
    "compute_XES_per_final_state",
    "compute_XES_per_final_state_loop",
    # Spectrum
    "VibrationalState",
    "KHResult",
    "KHCalculator",
]
