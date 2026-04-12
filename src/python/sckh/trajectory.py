"""SCKH trajectory data structure.

An SCKH trajectory stores pre-computed spectroscopic data (energies and
transition dipole moments) at each timestep, as used by the SCKH
(Semi-Classical Kramers-Heisenberg) spectrum computation.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class SCKHTrajectory:
    """Pre-computed spectroscopic data along a trajectory (SCKH format).

    Attributes:
        time: Time points (ntsteps,) in seconds
        E_gs: Ground state energy (ntsteps,) in Hartree
        E_n: Intermediate state energy (ntsteps,) in eV
        E_f: Final state energies (nfinal, ntsteps) in eV
        D_fn: Transition dipole moments (nfinal, ntsteps, 3) in atomic units
        E_IP1s: Ionization potential (ntsteps,) in Hartree (zero for dynamics-generated)
    """

    time: np.ndarray
    E_gs: np.ndarray
    E_n: np.ndarray
    E_f: np.ndarray
    D_fn: np.ndarray
    E_IP1s: np.ndarray

    @property
    def n_final(self) -> int:
        """Number of final states."""
        return self.E_f.shape[0]

    @property
    def ntsteps(self) -> int:
        """Number of time steps."""
        return len(self.time)

    def interpolate_dt(self, dt_new: float) -> "SCKHTrajectory":
        """Interpolate trajectory to a new time step.

        Useful when the dynamics time step is too large for the desired
        spectral range.  The FFT frequency range is proportional to
        1/(2*dt), so a smaller dt gives wider spectral coverage.

        Args:
            dt_new: New time step in seconds

        Returns:
            New SCKHTrajectory with interpolated arrays
        """
        from scipy.interpolate import CubicSpline

        time_new = np.arange(self.time[0], self.time[-1], dt_new)
        # Ensure we include a point at or very near the end
        if time_new[-1] < self.time[-1] - dt_new * 0.5:
            time_new = np.append(time_new, self.time[-1])

        # Interpolate 1D arrays
        E_gs_new = CubicSpline(self.time, self.E_gs)(time_new)
        E_n_new = CubicSpline(self.time, self.E_n)(time_new)
        E_IP1s_new = CubicSpline(self.time, self.E_IP1s)(time_new)

        # Interpolate per-final-state arrays
        nsteps_new = len(time_new)
        E_f_new = np.zeros((self.n_final, nsteps_new))
        D_fn_new = np.zeros((self.n_final, nsteps_new, 3))

        for f in range(self.n_final):
            E_f_new[f] = CubicSpline(self.time, self.E_f[f])(time_new)
            for c in range(3):
                D_fn_new[f, :, c] = CubicSpline(
                    self.time, self.D_fn[f, :, c]
                )(time_new)

        return SCKHTrajectory(
            time=time_new,
            E_gs=E_gs_new,
            E_n=E_n_new,
            E_f=E_f_new,
            D_fn=D_fn_new,
            E_IP1s=E_IP1s_new,
        )
