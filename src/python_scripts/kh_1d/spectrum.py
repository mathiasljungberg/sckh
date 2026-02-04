"""Main KH spectrum calculator and result classes.

Provides KHCalculator class for computing Kramers-Heisenberg XES spectra
from potential energy surfaces and dipole moment files.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_1d.pes import PES1D, create_pes_from_file
from python_scripts.dynamics_1d.dipole import Dipole1D, create_dipole_from_file
from python_scripts.dynamics_1d.vibrational import solve_vibrational
from python_scripts.dynamics_1d.io import read_pes_file

from .config import KHConfig
from .dipole_matrix import compute_dipole_matrix_elements
from .amplitude import compute_XES_per_final_state


@dataclass
class VibrationalState:
    """Container for vibrational eigenstates of one electronic state.

    Attributes:
        eigenvalues: Energy eigenvalues in Joules, shape (n_states,)
        eigenvectors: Wavefunctions on grid (columns), shape (npoints, n_states)
        x: Position grid in meters, shape (npoints,)
    """

    eigenvalues: np.ndarray  # (n_states,) in Joules
    eigenvectors: np.ndarray  # (npoints, n_states)
    x: np.ndarray  # (npoints,) in meters

    @property
    def eigenvalues_eV(self) -> np.ndarray:
        """Eigenvalues in eV."""
        return self.eigenvalues / CONST.eV

    @property
    def n_states(self) -> int:
        """Number of vibrational states."""
        return len(self.eigenvalues)


@dataclass
class KHResult:
    """Result container for KH spectrum calculation.

    Attributes:
        omega: Frequency grid in eV, shape (n_omega,)
        sigma_tot: Total cross-section (normalized), shape (n_omega,)
        sigma_per_final_elec: Cross-section per final electronic state,
            shape (n_final_elec, n_omega)
        eigenvalues_initial: Initial state vibrational eigenvalues (eV)
        eigenvalues_intermediate: Intermediate state vibrational eigenvalues (eV)
        eigenvalues_final: Final state eigenvalues per electronic state (list of arrays)
        D_ni: Transition dipoles initial->intermediate
        D_fn_list: Transition dipoles intermediate->final (per electronic state)
    """

    omega: np.ndarray
    sigma_tot: np.ndarray
    sigma_per_final_elec: np.ndarray
    eigenvalues_initial: np.ndarray
    eigenvalues_intermediate: np.ndarray
    eigenvalues_final: List[np.ndarray]
    D_ni: np.ndarray
    D_fn_list: List[np.ndarray]

    def write_spectrum(self, filepath: Path) -> None:
        """Write total spectrum to file in Fortran-compatible format."""
        from python_scripts.dynamics_1d.io import write_spectrum

        write_spectrum(filepath, self.omega, self.sigma_tot)

    def write_spectrum_per_final(self, filepath_base: Path) -> None:
        """Write per-final-electronic-state spectra."""
        from python_scripts.dynamics_1d.io import write_spectrum

        for j, sigma in enumerate(self.sigma_per_final_elec):
            filepath = Path(f"{filepath_base}_states_{j + 1}.dat")
            write_spectrum(
                filepath,
                self.omega,
                sigma,
                header=f"omega(eV) sigma (final elec. state {j + 1})",
            )


class KHCalculator:
    """Kramers-Heisenberg spectrum calculator.

    Computes XES spectra using the Kramers-Heisenberg formula for
    non-resonant scattering.

    Example:
        >>> config = KHConfig(...)
        >>> calc = KHCalculator(config)
        >>> result = calc.run()
        >>> result.write_spectrum("spectrum.dat")
    """

    def __init__(self, config: KHConfig):
        """Initialize calculator with configuration.

        Args:
            config: KHConfig object with all calculation parameters
        """
        self.config = config

        # Internal storage (populated by load_surfaces and solve_vibrational_problem)
        self._pes_i: Optional[PES1D] = None
        self._pes_n: Optional[PES1D] = None
        self._pes_f_list: List[PES1D] = []
        self._dipole_f_list: List[Dipole1D] = []

        self._vib_i: Optional[VibrationalState] = None
        self._vib_n: Optional[VibrationalState] = None
        self._vib_f_list: List[VibrationalState] = []

        # Grid
        self._x: Optional[np.ndarray] = None
        self._dx: Optional[float] = None

    def _setup_grid(self) -> Tuple[np.ndarray, float]:
        """Create spatial grid from configuration.

        Returns:
            x: Position grid in meters
            dx: Grid spacing in meters
        """
        cfg = self.config.grid
        # Convert from Angstrom to meters
        start_m = cfg.start * 1e-10
        dx_m = cfg.dx * 1e-10

        x = np.array([start_m + i * dx_m for i in range(cfg.npoints)])
        return x, dx_m

    def load_surfaces(self) -> None:
        """Load PES and dipole surfaces from files."""
        cfg = self.config

        # Setup grid
        self._x, self._dx = self._setup_grid()

        # Load initial state PES
        self._pes_i = create_pes_from_file(
            cfg.pes_initial,
            units="angstrom",
            energy_column=cfg.energy_column_initial,
        )

        # Load intermediate state PES
        self._pes_n = create_pes_from_file(
            cfg.pes_intermediate,
            units="angstrom",
            energy_column=cfg.energy_column_intermediate,
        )

        # Load final state PES with optional LP energy correction
        self._pes_f_list = []
        self._dipole_f_list = []

        if cfg.pes_lp_corr:
            # Apply energy shift so that first final state matches correction PES
            x_corr, E_corr = read_pes_file(
                cfg.pes_lp_corr, units="angstrom",
                energy_column=cfg.energy_column_final,
            )
            x_f0, E_f0 = read_pes_file(
                cfg.pes_final_list[0], units="angstrom",
                energy_column=cfg.energy_column_final,
            )
            shift = E_corr - E_f0

            for pes_path in cfg.pes_final_list:
                x, E = read_pes_file(
                    pes_path, units="angstrom",
                    energy_column=cfg.energy_column_final,
                )
                self._pes_f_list.append(PES1D(x=x, E=E + shift))
        else:
            for pes_path in cfg.pes_final_list:
                pes_f = create_pes_from_file(
                    pes_path, units="angstrom",
                    energy_column=cfg.energy_column_final,
                )
                self._pes_f_list.append(pes_f)

        # Load dipole files
        for dipole_path in cfg.dipole_final_list:
            dipole_f = create_dipole_from_file(dipole_path, units="angstrom")
            self._dipole_f_list.append(dipole_f)

    def solve_vibrational_problem(self) -> None:
        """Solve vibrational Schrödinger equation for all electronic states."""
        if self._x is None:
            self._setup_grid()
            self._x, self._dx = self._setup_grid()

        cfg = self.config
        mass = cfg.mu * CONST.u  # Convert amu to kg

        n_states = cfg.n_vib_states  # None means all states

        # Get potential on grid for each electronic state
        V_i = self._pes_i.energy(self._x)
        V_n = self._pes_n.energy(self._x)

        # Solve for initial state
        eig_i, vec_i = solve_vibrational(self._x, V_i, mass, n_states)
        self._vib_i = VibrationalState(
            eigenvalues=eig_i, eigenvectors=vec_i, x=self._x
        )

        # Solve for intermediate state
        eig_n, vec_n = solve_vibrational(self._x, V_n, mass, n_states)
        self._vib_n = VibrationalState(
            eigenvalues=eig_n, eigenvectors=vec_n, x=self._x
        )

        # Solve for each final electronic state
        self._vib_f_list = []
        for pes_f in self._pes_f_list:
            V_f = pes_f.energy(self._x)
            eig_f, vec_f = solve_vibrational(self._x, V_f, mass, n_states)
            vib_f = VibrationalState(
                eigenvalues=eig_f, eigenvectors=vec_f, x=self._x
            )
            self._vib_f_list.append(vib_f)

    def compute_spectrum_nonres(self, verbose: bool = False) -> KHResult:
        """Compute non-resonant XES spectrum.

        Args:
            verbose: Print progress information

        Returns:
            KHResult containing spectrum and intermediate data
        """
        cfg = self.config
        gamma = cfg.broadening.gamma_hwhm  # HWHM in eV
        omega_grid = cfg.frequency.get_omega_array()

        # Get eigenvalues in eV
        E_i = self._vib_i.eigenvalues_eV
        E_n = self._vib_n.eigenvalues_eV

        if verbose:
            print(f"Initial state: {len(E_i)} vibrational states")
            print(f"  Ground state energy: {E_i[0]:.6f} eV")
            print(f"Intermediate state: {len(E_n)} vibrational states")
            print(f"  Ground state energy: {E_n[0]:.6f} eV")

        # Compute spectrum for each final electronic state
        n_omega = len(omega_grid)
        n_final_elec = cfg.n_final_states
        sigma_per_final_elec = np.zeros((n_final_elec, n_omega))

        D_fn_list = []
        D_ni = None  # Will be computed once (same for all final states)

        for j, (vib_f, dipole_f) in enumerate(
            zip(self._vib_f_list, self._dipole_f_list)
        ):
            if verbose:
                print(f"Final electronic state {j + 1}: {vib_f.n_states} vib. states")

            E_f = vib_f.eigenvalues_eV

            # Get dipole on grid
            dipole_on_grid = dipole_f.dipole(self._x)

            # Compute transition dipoles
            # For initial state, we only use the ground vibrational state (v=0)
            vib_i_ground = self._vib_i.eigenvectors[:, 0:1]  # Keep 2D shape

            D_ni_j, D_fn_j = compute_dipole_matrix_elements(
                vib_i=vib_i_ground,
                vib_n=self._vib_n.eigenvectors,
                vib_f=vib_f.eigenvectors,
                dipole_on_grid=dipole_on_grid,
                mode=cfg.dipole_mode,
                dx=self._dx,
                x=self._x,
                x0=self._pes_i.find_minimum()[0] if cfg.dipole_mode == "DIPOLE_X0" else None,
            )

            if D_ni is None:
                # D_ni has shape (n_states_n, 1, 3) -> squeeze to (n_states_n, 3)
                D_ni = D_ni_j[:, 0, :]

            D_fn_list.append(D_fn_j)

            # Compute spectrum for this final electronic state
            sigma_j = compute_XES_per_final_state(
                E_i=E_i[0],  # Ground vibrational state
                E_n=E_n,
                E_f=E_f,
                D_ni=D_ni,
                D_fn=D_fn_j,
                omega_grid=omega_grid,
                gamma=gamma,
            )

            # Sum over final vibrational states
            sigma_per_final_elec[j, :] = np.sum(sigma_j, axis=0)

        # Total spectrum: sum over final electronic states
        sigma_tot = np.sum(sigma_per_final_elec, axis=0)

        # Normalize to unit integral
        d_omega = omega_grid[1] - omega_grid[0]
        norm = np.sum(sigma_tot) * d_omega
        if norm > 0:
            sigma_tot = sigma_tot / norm
            sigma_per_final_elec = sigma_per_final_elec / norm

        # Collect eigenvalues for output
        eigenvalues_final = [vib_f.eigenvalues_eV for vib_f in self._vib_f_list]

        return KHResult(
            omega=omega_grid,
            sigma_tot=sigma_tot,
            sigma_per_final_elec=sigma_per_final_elec,
            eigenvalues_initial=E_i,
            eigenvalues_intermediate=E_n,
            eigenvalues_final=eigenvalues_final,
            D_ni=D_ni,
            D_fn_list=D_fn_list,
        )

    def run(self, verbose: bool = False) -> KHResult:
        """Run full KH calculation workflow.

        1. Load PES and dipole surfaces
        2. Solve vibrational problems
        3. Compute spectrum

        Args:
            verbose: Print progress information

        Returns:
            KHResult containing spectrum and intermediate data
        """
        if verbose:
            print("Loading potential energy and dipole surfaces...")
        self.load_surfaces()

        if verbose:
            print("Solving vibrational Schrödinger equations...")
        self.solve_vibrational_problem()

        if verbose:
            print("Computing KH spectrum...")
        result = self.compute_spectrum_nonres(verbose=verbose)

        if verbose:
            print("Done.")

        return result
