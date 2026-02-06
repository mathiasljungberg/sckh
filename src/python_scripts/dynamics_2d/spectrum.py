"""Spectrum calculation via FFT of classical trajectories for 2D dynamics.

Implements the SCKH (Semi-Classical Kramers-Heisenberg) approach for computing
X-ray emission spectra from 2D classical molecular dynamics trajectories.

Reuses core functions from dynamics_1d.spectrum that operate on 1D time series.
"""

from pathlib import Path
from typing import List, Optional

import numpy as np

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_1d.spectrum import (
    compute_energy_phase,
    compute_F_if,
    get_frequency_grid,
    SpectrumResult,
)
from python_scripts.dynamics_1d.io import write_spectrum, write_spectrum_per_final

from .config import DynamicsConfig2D, SpectrumConfig2D, FullConfig2D
from .dipole import Dipole2D, create_dipole_from_file_2d, create_constant_dipole_2d
from .io import read_pes_file_2d, read_pes_file_2d_raw
from .pes import PES2D, create_pes_from_file_2d
from .trajectory import TrajectoryResult2D, DynamicsRunner2D


class SpectrumCalculator2D:
    """Calculator for SCKH spectrum from 2D classical trajectories.

    Usage:
        dynamics_config = load_config("config.yaml").dynamics2d
        spectrum_config = SpectrumConfig2D(...)
        calculator = SpectrumCalculator2D(dynamics_config, spectrum_config)
        result = calculator.run()
    """

    def __init__(
        self,
        config: FullConfig2D,
    ):
        """Initialize calculator with configuration.

        Args:
            dynamics_config: Dynamics configuration (for PES files and units)
            spectrum_config: Spectrum calculation configuration
        """
        self.dynamics_config = config.dynamics2d
        self.spectrum_config = config.spectrum

        # PES surfaces
        self.pes_n: Optional[PES2D] = None  # Intermediate state
        self.pes_f: List[PES2D] = []  # Final states

        # Dipole surfaces
        self.dipoles: List[Dipole2D] = []  # One per final state
        self.dipole_initial: Optional[Dipole2D] = None

        # Computed values
        self.E_mean: float = 0.0  # Mean transition energy (eV)
        self._surfaces_loaded = False

    def load_surfaces(self) -> None:
        """Load all PES and dipole surfaces from files."""
        pos_units = self.dynamics_config.position_units
        energy_units = self.dynamics_config.energy_units
        index_order = self.dynamics_config.index_order

        # Load intermediate state PES
        if self.spectrum_config.pes_intermediate:
            pes_file = self.spectrum_config.pes_intermediate
        else:
            pes_file = self.dynamics_config.pes_dynamics

        self.pes_n = create_pes_from_file_2d(
            pes_file,
            position_units=pos_units,
            energy_units=energy_units,
            index_order=index_order,
        )

        # Load final state PES with optional energy correction
        if self.spectrum_config.pes_lp_corr:
            # Apply energy shift so that first final state matches correction PES
            x1_corr, x2_corr, E_corr = read_pes_file_2d(
                self.spectrum_config.pes_lp_corr,
                position_units=pos_units,
                energy_units=energy_units,
                index_order=index_order,
            )
            x1_f0, x2_f0, E_f0 = read_pes_file_2d(
                self.spectrum_config.pes_final_list[0],
                position_units=pos_units,
                energy_units=energy_units,
                index_order=index_order,
            )
            shift = E_corr - E_f0

            self.pes_f = []
            for path in self.spectrum_config.pes_final_list:
                x1, x2, E = read_pes_file_2d(
                    path,
                    position_units=pos_units,
                    energy_units=energy_units,
                    index_order=index_order,
                )
                self.pes_f.append(PES2D(x1=x1, x2=x2, E=E + shift))
        else:
            # Load without energy correction
            self.pes_f = [
                create_pes_from_file_2d(
                    path,
                    position_units=pos_units,
                    energy_units=energy_units,
                    index_order=index_order,
                )
                for path in self.spectrum_config.pes_final_list
            ]

        # Load dipole surfaces based on mode
        dipole_components = self.spectrum_config.dipole_components
        if self.spectrum_config.dipole_mode == "DIPOLE":
            self.dipoles = [
                create_dipole_from_file_2d(
                    path,
                    position_units=pos_units,
                    index_order=index_order,
                    dipole_components=dipole_components,
                )
                for path in self.spectrum_config.dipole_final_list
            ]
        elif self.spectrum_config.dipole_mode == "FC":
            # Franck-Condon: constant dipole = 1
            self.dipoles = [
                create_constant_dipole_2d(
                    self.pes_n.x1, self.pes_n.x2, np.array([1.0, 1.0, 1.0])
                )
                for _ in range(len(self.pes_f))
            ]
        elif self.spectrum_config.dipole_mode == "DIPOLE_X0":
            # Dipole frozen at equilibrium - will be handled during interpolation
            self.dipoles = [
                create_dipole_from_file_2d(
                    path,
                    position_units=pos_units,
                    index_order=index_order,
                    dipole_components=dipole_components,
                )
                for path in self.spectrum_config.dipole_final_list
            ]

        # Load initial state dipole if provided
        if self.spectrum_config.dipole_initial:
            self.dipole_initial = create_dipole_from_file_2d(
                self.spectrum_config.dipole_initial,
                position_units=pos_units,
                index_order=index_order,
                dipole_components=dipole_components,
            )

        self._surfaces_loaded = True

    def load_trajectories(self) -> List[TrajectoryResult2D]:
        """Load trajectories from files or run dynamics.

        Returns:
            List of TrajectoryResult2D objects
        """
        if self.spectrum_config.trajectory_files:
            # Load pre-computed trajectories
            from .io import read_trajectory_2d
            return [
                read_trajectory_2d(path)
                for path in self.spectrum_config.trajectory_files
            ]
        else:
            # Run dynamics to generate trajectories
            from .config import FullConfig2D
            full_config = FullConfig2D(dynamics2d=self.dynamics_config)
            runner = DynamicsRunner2D(full_config)
            result = runner.run()
            return result.trajectories

    def interpolate_along_trajectory(
        self,
        traj: TrajectoryResult2D,
    ) -> dict:
        """Interpolate energies and dipoles along 2D trajectory positions.

        KEY DIFFERENCE from 1D: evaluates 2D surfaces at (x1(t), x2(t)).

        Args:
            traj: Trajectory result with positions x1, x2

        Returns:
            Dictionary with:
                E_n: (ntsteps,) intermediate state energies in Joules
                E_f: (nfinal, ntsteps) final state energies in Joules
                D_fn: (nfinal, ntsteps, 3) dipole moments
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        x1 = traj.x1
        x2 = traj.x2
        nsteps = len(x1)
        n_final = len(self.pes_f)

        # Intermediate state energy: E_n(t) = V_n(x1(t), x2(t))
        E_n = self.pes_n.energy(x1, x2)

        # Final state energies
        E_f = np.zeros((n_final, nsteps))
        for i, pes in enumerate(self.pes_f):
            E_f[i] = pes.energy(x1, x2)

        # Dipole moments
        D_fn = np.zeros((n_final, nsteps, 3))

        if self.spectrum_config.dipole_mode == "DIPOLE_X0":
            # Frozen at equilibrium position
            # Find equilibrium from initial state PES
            pes_i = create_pes_from_file_2d(
                self.dynamics_config.pes_initial,
                position_units=self.dynamics_config.position_units,
                energy_units=self.dynamics_config.energy_units,
                index_order=self.dynamics_config.index_order,
            )
            x1_eq, x2_eq, _ = pes_i.find_minimum()

            for i, dipole in enumerate(self.dipoles):
                d_eq = dipole.dipole(x1_eq, x2_eq)
                D_fn[i, :, :] = d_eq  # Broadcast to all time steps
        else:
            for i, dipole in enumerate(self.dipoles):
                D_fn[i] = dipole.dipole(x1, x2)

        return {"E_n": E_n, "E_f": E_f, "D_fn": D_fn}

    def compute_mean_transition_energy(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> float:
        """Compute mean transition energy.

        In "standard" mode: Finds the true 2D equilibrium position (x1_eq, x2_eq)
        from the initial state PES and evaluates E_mean = E_n(eq) - E_f(eq).

        In "fortran" mode: Uses trajectory-based indexing for compatibility.

        Args:
            trajectories: List of trajectories

        Returns:
            Mean transition energy in eV
        """
        if not self._surfaces_loaded:
            self.load_surfaces()

        mode = self.spectrum_config.compatibility_mode

        if mode == "fortran":
            return self._compute_E_mean_fortran(trajectories)
        else:
            return self._compute_E_mean_standard()

    def _compute_E_mean_standard(self) -> float:
        """Standard E_mean: evaluate at true 2D equilibrium position."""
        # Find 2D equilibrium position from initial state PES
        pes_i = create_pes_from_file_2d(
            self.dynamics_config.pes_initial,
            position_units=self.dynamics_config.position_units,
            energy_units=self.dynamics_config.energy_units,
            index_order=self.dynamics_config.index_order,
        )
        x1_eq, x2_eq, _ = pes_i.find_minimum()

        # Evaluate transition energy at equilibrium
        E_n_at_eq = self.pes_n.energy(x1_eq, x2_eq)
        E_f_at_eq = self.pes_f[-1].energy(x1_eq, x2_eq)  # Last final state

        E_mean = (E_n_at_eq - E_f_at_eq) / CONST.eV  # Convert to eV
        return E_mean

    def _compute_E_mean_fortran(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> float:
        """Fortran-compatible E_mean: use trajectory-based indexing.

        Following Fortran pattern for 2D case.
        """
        # Find minimum energy index in initial state PES
        _, _, E_i_raw = read_pes_file_2d_raw(
            self.dynamics_config.pes_initial,
            index_order=self.dynamics_config.index_order,
        )
        ind = np.unravel_index(np.argmin(E_i_raw), E_i_raw.shape)

        # Use trajectory positions near the minimum
        x1_traj = trajectories[0].x1
        x2_traj = trajectories[0].x2

        # Get energies along trajectory at a representative index
        idx = min(ind[0] * E_i_raw.shape[1] + ind[1], len(x1_traj) - 1)
        E_n_at_ind = self.pes_n.energy(x1_traj[idx], x2_traj[idx])
        E_f_at_ind = self.pes_f[-1].energy(x1_traj[idx], x2_traj[idx])

        E_mean = (E_n_at_ind - E_f_at_ind) / CONST.eV
        return E_mean

    def compute_spectrum(
        self,
        trajectories: List[TrajectoryResult2D],
    ) -> SpectrumResult:
        """Main spectrum calculation.

        For each trajectory:
            1. Interpolate energies and dipoles along 2D trajectory
            2. Compute energy phase factors (reuses 1D function)
            3. Compute transition amplitudes via FFT (reuses 1D function)
            4. Accumulate |F_if|² into cross-section

        Args:
            trajectories: List of 2D trajectory results

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
        gamma = self.spectrum_config.gamma_hwhm

        # Compute mean transition energy
        self.E_mean = self.compute_mean_transition_energy(trajectories)

        # Get frequency grid (reuses 1D function)
        omega = get_frequency_grid(time, self.E_mean)

        # Get initial state dipole (D_ni)
        if self.dipole_initial:
            pes_i = create_pes_from_file_2d(
                self.dynamics_config.pes_initial,
                position_units=self.dynamics_config.position_units,
                energy_units=self.dynamics_config.energy_units,
                index_order=self.dynamics_config.index_order,
            )
            x1_eq, x2_eq, _ = pes_i.find_minimum()
            D_ni = self.dipole_initial.dipole(x1_eq, x2_eq)
        else:
            # Default: unit dipole for all components
            D_ni = np.array([1.0, 1.0, 1.0])

        # Accumulate |F_if|² over all trajectories
        # Shape: (n_final, n_omega, 3, 3) for polarization components
        sigma_mm = np.zeros((n_final, nsteps, 3, 3))

        for traj in trajectories:
            # Interpolate along 2D trajectory
            interp = self.interpolate_along_trajectory(traj)
            E_n = interp["E_n"]
            E_f = interp["E_f"]
            D_fn = interp["D_fn"]

            # For each final state
            for f_idx in range(n_final):
                # Compute energy phase factor (reuses 1D function)
                e_factor = compute_energy_phase(
                    E_n, E_f[f_idx], self.E_mean, traj.time
                )

                # Compute transition amplitude via FFT (reuses 1D function)
                F_if = compute_F_if(
                    e_factor, D_fn[f_idx], D_ni, traj.time, gamma
                )

                # Accumulate |F_if|²
                sigma_mm[f_idx] += np.abs(F_if) ** 2

        # Sum over polarizations to get per-final-state cross-section
        sigma_f = np.zeros((n_final, nsteps))
        for m in range(3):
            sigma_f += sigma_mm[:, :, m, m]

        # Total cross-section
        sigma_tot = np.sum(sigma_f, axis=0)

        # Normalize (match Fortran: time_l = (nsteps-1) * dt)
        T = time[-1] - time[0]
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
            print("Loading 2D PES and dipole surfaces...")
        self.load_surfaces()

        if verbose:
            print("Loading/computing 2D trajectories...")
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

        basename = self.dynamics_config.outfile

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
