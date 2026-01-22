"""Spectrum calculation via FFT of classical trajectories.

Implements the SCKH (Semi-Classical Kramers-Heisenberg) approach for computing
X-ray emission spectra from classical molecular dynamics trajectories.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from .constants import CONST
from .dipole import Dipole1D, create_dipole_from_file, create_constant_dipole
from .io import (
    read_trajectory_file,
    write_spectrum,
    write_spectrum_per_final,
)
from .pes import PES1D, create_pes_from_file
from .spectrum_config import FullConfig
from .trajectory import TrajectoryResult, DynamicsRunner


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
    """Compute energy phase factor exp(i * ∫(E_n - E_f - E_mean) dt / ℏ).

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

    # Cumulative integral: ∫(delta_E) dt using rectangular rule
    # Matches Fortran: int_W_I(i) = int_W_I(i-1) + delta_E(i)
    int_W = np.cumsum(delta_E) * dt

    # Phase factor: exp(-i * integral / ℏ)
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

    F_if(ω) = FFT[ D_fn(t) * e_factor(t) * exp(-γt/ℏ) ]

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
    dt = time[1] - time[0]

    # Decay factor: exp(-γt/ℏ)
    decay = np.exp(-gamma * CONST.eV * time / CONST.hbar)

    F_if = np.zeros((nsteps, 3, 3), dtype=complex)

    for m1 in range(3):  # Final state polarization
        # Integrand: D(t) * phase(t) * decay(t)
        integrand = D_fn[:, m1] * e_factor * decay

        # FFT (backward = inverse, matches Fortran's fft_c2c_1d_backward)
        # Fortran FFTW backward is unnormalized IFFT * n
        # numpy.fft.ifft is normalized (divides by n), so we multiply by n
        # Also multiply by dt to get the integral approximation
        fft_result = np.fft.ifft(integrand) * nsteps * dt

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
    T = time[-1] - time[0] + (time[1] - time[0])  # Total time span

    # FFT frequencies in Hz
    # After fftshift: [-n/2, ..., -1, 0, 1, ..., n/2-1] * (1/T)
    freq_indices = np.arange(nsteps) - nsteps // 2
    freq = freq_indices / T  # in Hz

    # Convert to angular frequency then to energy
    # ω = 2πf, E = ℏω
    omega = 2 * np.pi * CONST.hbar * freq / CONST.eV

    # Shift by mean transition energy
    omega += E_mean

    return omega


class SpectrumCalculator:
    """Calculator for SCKH spectrum from classical trajectories.

    Usage:
        config = load_full_config("config.yaml")
        calculator = SpectrumCalculator(config)
        result = calculator.run()
    """

    def __init__(self, config: FullConfig):
        """Initialize calculator with configuration.

        Args:
            config: Combined dynamics and spectrum configuration
        """
        self.config = config

        # PES surfaces
        self.pes_n: Optional[PES1D] = None  # Intermediate state
        self.pes_f: List[PES1D] = []  # Final states

        # Dipole surfaces
        self.dipoles: List[Dipole1D] = []  # One per final state
        self.dipole_initial: Optional[Dipole1D] = None

        # Computed values
        self.E_mean: float = 0.0  # Mean transition energy (eV)
        self._surfaces_loaded = False

    def load_surfaces(self) -> None:
        """Load all PES and dipole surfaces from files."""
        units = self.config.dynamics.units

        # Load intermediate state PES (same as dynamics PES)
        self.pes_n = create_pes_from_file(
            self.config.dynamics.pes_dynamics, units=units
        )

        # Load final state PES
        self.pes_f = [
            create_pes_from_file(path, units=units)
            for path in self.config.spectrum.pes_final_list
        ]

        # Load dipole surfaces based on mode
        if self.config.spectrum.dipole_mode == "DIPOLE":
            self.dipoles = [
                create_dipole_from_file(path, units=units)
                for path in self.config.spectrum.dipole_final_list
            ]
        elif self.config.spectrum.dipole_mode == "FC":
            # Franck-Condon: constant dipole = 1
            x = self.pes_n.x
            self.dipoles = [
                create_constant_dipole(x, np.array([1.0, 1.0, 1.0]))
                for _ in range(len(self.pes_f))
            ]
        elif self.config.spectrum.dipole_mode == "DIPOLE_X0":
            # Dipole frozen at equilibrium - will be set during interpolation
            self.dipoles = [
                create_dipole_from_file(path, units=units)
                for path in self.config.spectrum.dipole_final_list
            ]

        # Load initial state dipole if provided
        if self.config.spectrum.dipole_initial:
            self.dipole_initial = create_dipole_from_file(
                self.config.spectrum.dipole_initial, units=units
            )

        self._surfaces_loaded = True

    def load_trajectories(self) -> List[TrajectoryResult]:
        """Load trajectories from files or run dynamics.

        Returns:
            List of TrajectoryResult objects
        """
        if self.config.spectrum.trajectory_files:
            # Load pre-computed trajectories
            return [
                read_trajectory_file(path)
                for path in self.config.spectrum.trajectory_files
            ]
        else:
            # Run dynamics to generate trajectories
            runner = DynamicsRunner(self.config.dynamics)
            result = runner.run()
            return result.trajectories

    def interpolate_along_trajectory(
        self,
        traj: TrajectoryResult,
    ) -> dict:
        """Interpolate energies and dipoles along trajectory positions.

        Args:
            traj: Trajectory result with positions x

        Returns:
            Dictionary with:
                E_n: (ntsteps,) intermediate state energies in Joules
                E_f: (nfinal, ntsteps) final state energies in Joules
                D_fn: (nfinal, ntsteps, 3) dipole moments
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        x = traj.x
        nsteps = len(x)
        n_final = len(self.pes_f)

        # Intermediate state energy
        E_n = self.pes_n.energy(x)

        # Final state energies
        E_f = np.zeros((n_final, nsteps))
        for i, pes in enumerate(self.pes_f):
            E_f[i] = pes.energy(x)

        # Dipole moments
        D_fn = np.zeros((n_final, nsteps, 3))

        if self.config.spectrum.dipole_mode == "DIPOLE_X0":
            # Frozen at equilibrium position
            # Find equilibrium from initial state PES
            pes_i = create_pes_from_file(
                self.config.dynamics.pes_initial,
                units=self.config.dynamics.units,
            )
            x_eq, _ = pes_i.find_minimum()

            for i, dipole in enumerate(self.dipoles):
                d_eq = dipole.dipole(x_eq)
                D_fn[i, :, :] = d_eq  # Broadcast to all time steps
        else:
            for i, dipole in enumerate(self.dipoles):
                D_fn[i] = dipole.dipole(x)

        return {"E_n": E_n, "E_f": E_f, "D_fn": D_fn}

    def compute_mean_transition_energy(
        self,
        trajectories: List[TrajectoryResult],
    ) -> float:
        """Compute mean transition energy from first trajectory.

        Following Fortran: E_mean = E_n(x_eq) - E_f[last](x_eq)

        Args:
            trajectories: List of trajectories

        Returns:
            Mean transition energy in eV
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        # Find equilibrium position from initial state PES
        pes_i = create_pes_from_file(
            self.config.dynamics.pes_initial,
            units=self.config.dynamics.units,
        )
        x_eq, _ = pes_i.find_minimum()

        # Compute transition energy at equilibrium
        E_n_eq = self.pes_n.energy(x_eq)
        E_f_eq = self.pes_f[-1].energy(x_eq)  # Last final state

        E_mean = (E_n_eq - E_f_eq) / CONST.eV  # Convert to eV

        return E_mean

    def compute_spectrum(
        self,
        trajectories: List[TrajectoryResult],
    ) -> SpectrumResult:
        """Main spectrum calculation.

        For each trajectory:
            1. Interpolate energies and dipoles
            2. Compute energy phase factors
            3. Compute transition amplitudes via FFT
            4. Accumulate |F_if|² into cross-section

        Args:
            trajectories: List of trajectory results

        Returns:
            SpectrumResult with omega, sigma_tot, sigma_f arrays
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        if len(trajectories) == 0:
            raise ValueError("No trajectories provided")

        # Get time array from first trajectory
        time = trajectories[0].time
        nsteps = len(time)
        n_final = len(self.pes_f)
        gamma = self.config.spectrum.gamma_hwhm

        # Compute mean transition energy
        self.E_mean = self.compute_mean_transition_energy(trajectories)

        # Get frequency grid
        omega = get_frequency_grid(time, self.E_mean)

        # Get initial state dipole (D_ni)
        # Using constant value or from equilibrium position
        if self.dipole_initial:
            pes_i = create_pes_from_file(
                self.config.dynamics.pes_initial,
                units=self.config.dynamics.units,
            )
            x_eq, _ = pes_i.find_minimum()
            D_ni = self.dipole_initial.dipole(x_eq)
        else:
            # Default: unit dipole for all components
            D_ni = np.array([1.0, 1.0, 1.0])

        # Accumulate |F_if|² over all trajectories
        # Shape: (n_final, n_omega, 3, 3) for polarization components
        sigma_mm = np.zeros((n_final, nsteps, 3, 3))

        for traj in trajectories:
            # Interpolate along trajectory
            interp = self.interpolate_along_trajectory(traj)
            E_n = interp["E_n"]
            E_f = interp["E_f"]
            D_fn = interp["D_fn"]

            # For each final state
            for f_idx in range(n_final):
                # Compute energy phase factor
                e_factor = compute_energy_phase(
                    E_n, E_f[f_idx], self.E_mean, traj.time
                )

                # Compute transition amplitude via FFT
                F_if = compute_F_if(
                    e_factor, D_fn[f_idx], D_ni, traj.time, gamma
                )

                # Accumulate |F_if|²
                sigma_mm[f_idx] += np.abs(F_if) ** 2

        # Sum over polarizations to get per-final-state cross-section
        # sigma_f[f, omega] = sum_{m1} sigma_mm[f, omega, m1, m1]
        sigma_f = np.zeros((n_final, nsteps))
        for m in range(3):
            sigma_f += sigma_mm[:, :, m, m]

        # Total cross-section
        sigma_tot = np.sum(sigma_f, axis=0)

        # Normalize
        T = time[-1] - time[0] + (time[1] - time[0])
        norm = np.sum(sigma_tot) * 2 * np.pi * CONST.hbar / (T * CONST.eV)
        if norm > 0:
            sigma_tot = sigma_tot / norm
            sigma_f = sigma_f / norm

        return SpectrumResult(
            omega=omega,
            sigma_tot=sigma_tot,
            sigma_f=sigma_f,
            n_trajectories=len(trajectories),
            E_mean=self.E_mean,
        )

    def run(self, verbose: bool = False) -> SpectrumResult:
        """Full workflow: load data, get trajectories, compute spectrum.

        Args:
            verbose: Print progress information

        Returns:
            SpectrumResult with computed spectrum
        """
        if verbose:
            print("Loading PES and dipole surfaces...")
        self.load_surfaces()

        if verbose:
            print("Loading/computing trajectories...")
        trajectories = self.load_trajectories()

        if verbose:
            print(f"Computing spectrum from {len(trajectories)} trajectories...")
        result = self.compute_spectrum(trajectories)

        if verbose:
            print(f"Mean transition energy: {result.E_mean:.2f} eV")
            print(f"Frequency range: {result.omega[0]:.2f} to {result.omega[-1]:.2f} eV")

        return result

    def save_results(
        self,
        result: SpectrumResult,
        output_dir: Path,
    ) -> None:
        """Save spectrum results to files.

        Args:
            result: SpectrumResult from compute_spectrum()
            output_dir: Directory for output files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        basename = self.config.dynamics.outfile

        # Write total spectrum
        write_spectrum(
            output_dir / f"{basename}_sigma.dat",
            result.omega,
            result.sigma_tot,
        )

        # Write per-final-state spectra
        write_spectrum_per_final(
            output_dir / f"{basename}_sigma",
            result.omega,
            result.sigma_f,
        )
