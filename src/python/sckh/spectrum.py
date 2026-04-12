"""SCKH spectrum computation via FFT.

Implements the Semi-Classical Kramers-Heisenberg approach for computing
X-ray emission spectra from SCKHTrajectory objects.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from dynamics_1d.constants import CONST

from .config import FullConfig, SpectrumConfig
from .trajectory import SCKHTrajectory
from .io import write_spectrum, write_spectrum_per_final


@dataclass
class SpectrumResult:
    """Result from spectrum calculation.

    Attributes:
        omega: Frequency grid in eV
        sigma_tot: Total cross-section (n_omega,)
        sigma_f: Per-final-state cross-section (n_final, n_omega)
        n_trajectories: Number of trajectories used
        E_mean: Mean transition energy used for frequency shift (eV)
    """

    omega: np.ndarray
    sigma_tot: np.ndarray
    sigma_f: np.ndarray
    n_trajectories: int
    E_mean: float


def compute_energy_phase(
    E_n: np.ndarray,
    E_f: np.ndarray,
    E_mean: float,
    time: np.ndarray,
) -> np.ndarray:
    """Compute energy phase factor exp(i * integral(E_n - E_f - E_mean) dt / hbar).

    Uses cumulative sum (rectangular rule) for the integral, matching Fortran.

    Args:
        E_n: Intermediate state energies along trajectory (ntsteps,) in Joules
        E_f: Final state energies along trajectory (ntsteps,) in Joules
        E_mean: Mean transition energy in eV
        time: Time array (ntsteps,) in seconds

    Returns:
        Complex phase factor array (ntsteps,)
    """
    dt = time[1] - time[0]

    # Energy difference: E_n - E_f - E_mean (convert E_mean from eV to J)
    delta_E = E_n - E_f - E_mean * CONST.eV

    # Cumulative integral: integral(delta_E) dt using rectangular rule
    # Matches Fortran: int_W_I(i) = int_W_I(i-1) + delta_E(i)
    int_W = np.cumsum(delta_E) * dt

    # Phase factor: exp(-i * integral / hbar)
    # Negative sign matches Fortran's factor = -1.0 for "negative" case
    phase = np.exp(-1j * int_W / CONST.hbar)

    return phase


def compute_F_if(
    e_factor: np.ndarray,
    D_fn: np.ndarray,
    D_ni: np.ndarray,
    time: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Compute transition amplitude via FFT.

    F_if(omega) = FFT[ D_fn(t) * e_factor(t) * exp(-gamma*t/hbar) ]

    Args:
        e_factor: Energy phase factor (ntsteps,) complex
        D_fn: Final state dipole along trajectory (ntsteps, 3)
        D_ni: Initial state dipole (3,) - constant or at equilibrium
        time: Time array (ntsteps,) in seconds
        gamma: HWHM broadening in eV

    Returns:
        F_if: Transition amplitude (n_omega, 3, 3) complex
    """
    nsteps = len(time)

    # Decay factor: exp(-gamma*t/hbar)
    decay = np.exp(-gamma * CONST.eV * time / CONST.hbar)

    F_if = np.zeros((nsteps, 3, 3), dtype=complex)

    for m1 in range(3):  # Final state polarization
        # Integrand: D(t) * phase(t) * decay(t)
        integrand = D_fn[:, m1] * e_factor * decay

        # FFT (backward = inverse, matches Fortran's fft_c2c_1d_backward)
        # Fortran FFTW backward is unnormalized IFFT * n
        # numpy.fft.ifft is normalized (divides by n), so we multiply by n
        # Note: Fortran doesn't multiply by dt, and it cancels in normalization
        fft_result = np.fft.ifft(integrand) * nsteps

        # Reorder: FFT output is [0, 1, ..., n/2-1, -n/2, ..., -1]
        # We want [-n/2, ..., -1, 0, 1, ..., n/2-1] for physical frequencies
        fft_result = np.fft.fftshift(fft_result)

        for m2 in range(3):  # Initial state polarization
            F_if[:, m1, m2] = D_ni[m2] * fft_result

    return F_if


def get_frequency_grid(time: np.ndarray, E_mean: float) -> np.ndarray:
    """Generate frequency grid matching FFT output.

    Reorders FFT frequencies and shifts by mean transition energy.

    Args:
        time: Time array (ntsteps,) in seconds
        E_mean: Mean transition energy in eV

    Returns:
        Frequencies in eV, matching fftshift ordering
    """
    nsteps = len(time)
    # Match Fortran: time_l = (ntsteps-1) * delta_t
    T = time[-1] - time[0]

    # FFT frequencies in Hz
    # After fftshift: [-n/2, ..., -1, 0, 1, ..., n/2-1] * (1/T)
    freq_indices = np.arange(nsteps) - nsteps // 2
    freq = freq_indices / T  # in Hz

    # Convert to angular frequency then to energy
    # omega = 2*pi*f, E = hbar*omega
    omega = 2 * np.pi * CONST.hbar * freq / CONST.eV

    # Shift by mean transition energy
    omega += E_mean

    return omega


class SCKHSpectrumCalculator:
    """Compute SCKH spectrum from SCKHTrajectory objects.

    This is the common spectrum computation class used by all workflows.
    Accepts either a FullConfig or a SpectrumConfig directly::

        # From FullConfig (YAML workflow)
        calc = SCKHSpectrumCalculator(config)

        # From specific config
        calc = SCKHSpectrumCalculator(
            spectrum_config=SpectrumConfig(gamma_fwhm=0.18)
        )
    """

    def __init__(
        self,
        config=None,
        *,
        spectrum_config: Optional[SpectrumConfig] = None,
    ):
        if config is not None and isinstance(config, SpectrumConfig):
            self.config = config
        elif config is not None and hasattr(config, "spectrum"):
            self.config = config.spectrum
        elif spectrum_config is not None:
            self.config = spectrum_config
        else:
            raise ValueError(
                "Provide a FullConfig, a SpectrumConfig positionally, "
                "or spectrum_config=SpectrumConfig(...)"
            )

        if self.config is None:
            raise ValueError("No spectrum config provided (FullConfig.spectrum is None)")

    def compute_spectrum(
        self,
        sckh_trajs: List[SCKHTrajectory],
        E_mean: Optional[float] = None,
        D_ni: Optional[np.ndarray] = None,
    ) -> SpectrumResult:
        """Compute spectrum from SCKH trajectory data.

        Args:
            sckh_trajs: List of SCKHTrajectory objects
            E_mean: Mean transition energy in eV. If None, auto-computed.
            D_ni: Initial state dipole (3,). Defaults to [1, 1, 1].

        Returns:
            SpectrumResult with omega, sigma_tot, sigma_f arrays
        """
        if self.config.dt is not None:
            sckh_trajs = [t.interpolate_dt(self.config.dt) for t in sckh_trajs]
        return compute_spectrum_from_sckh(
            sckh_trajs,
            gamma_hwhm=self.config.gamma_hwhm,
            E_mean=E_mean,
            D_ni=D_ni,
        )

    def save_results(
        self,
        result: SpectrumResult,
        output_dir: Path,
        basename: str = "spectrum",
    ) -> None:
        """Save spectrum results to files.

        Args:
            result: SpectrumResult from compute_spectrum()
            output_dir: Directory for output files
            basename: Prefix for output filenames
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        write_spectrum(
            output_dir / f"{basename}_sigma.dat",
            result.omega,
            result.sigma_tot,
        )

        write_spectrum_per_final(
            output_dir / f"{basename}_sigma",
            result.omega,
            result.sigma_f,
        )


def compute_spectrum_from_sckh(
    sckh_trajs: List,
    gamma_hwhm: float,
    E_mean: Optional[float] = None,
    D_ni: Optional[np.ndarray] = None,
) -> SpectrumResult:
    """Compute spectrum directly from SCKH trajectory data.

    This bypasses PES/dipole surface loading and interpolation, using
    pre-computed energies and dipoles from SCKH trajectory files.

    Args:
        sckh_trajs: List of SCKHTrajectory objects
        gamma_hwhm: Half-width at half maximum broadening in eV
        E_mean: Mean transition energy in eV. If None, computed as
                time-average of E_n - E_f[last] from first trajectory.
        D_ni: Initial state dipole (3,). Defaults to [1, 1, 1].

    Returns:
        SpectrumResult with omega, sigma_tot, sigma_f arrays
    """
    if len(sckh_trajs) == 0:
        raise ValueError("No trajectories provided")

    if D_ni is None:
        D_ni = np.array([1.0, 1.0, 1.0])

    first = sckh_trajs[0]
    time = first.time
    nsteps = first.ntsteps
    n_final = first.n_final

    # Compute E_mean if not provided (matching Fortran: average over trajectory)
    if E_mean is None:
        E_mean = np.mean(first.E_n - first.E_f[-1])

    # Get frequency grid
    omega = get_frequency_grid(time, E_mean)

    # Accumulate |F_if|^2 over all trajectories
    sigma_mm = np.zeros((n_final, nsteps, 3, 3))

    for sckh in sckh_trajs:
        # Convert energies from eV to Joules for phase computation
        E_n_J = sckh.E_n * CONST.eV
        E_f_J = sckh.E_f * CONST.eV

        for f_idx in range(n_final):
            e_factor = compute_energy_phase(
                E_n_J, E_f_J[f_idx], E_mean, sckh.time
            )
            F_if = compute_F_if(
                e_factor, sckh.D_fn[f_idx], D_ni, sckh.time, gamma_hwhm
            )
            sigma_mm[f_idx] += np.abs(F_if) ** 2

    # Sum over polarizations (diagonal elements)
    sigma_f = np.zeros((n_final, nsteps))
    for m in range(3):
        sigma_f += sigma_mm[:, :, m, m]

    # Total cross-section
    sigma_tot = np.sum(sigma_f, axis=0)

    # Normalize
    T = time[-1] - time[0]
    norm = np.sum(sigma_tot) * 2 * np.pi * CONST.hbar / (T * CONST.eV)
    if norm > 0:
        sigma_tot = sigma_tot / norm
        sigma_f = sigma_f / norm

    return SpectrumResult(
        omega=omega,
        sigma_tot=sigma_tot,
        sigma_f=sigma_f,
        n_trajectories=len(sckh_trajs),
        E_mean=E_mean,
    )
