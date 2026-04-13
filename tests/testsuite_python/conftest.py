"""Shared fixtures for the Python testsuite.

These tests validate the Python implementation against Fortran reference output,
reusing the same input files and reference data from the Fortran testsuite.
"""

from pathlib import Path

import numpy as np
import pytest


# Root of the Fortran testsuite (contains SCKH/, SCKH_PES/, etc.)
TESTSUITE_FORTRAN_DIR = Path(__file__).parent.parent / "testsuite"


def assert_spectrum_matches_ref(
    omega: np.ndarray,
    sigma: np.ndarray,
    ref_path: Path,
    rel_rms_threshold: float = 0.005,
    peak_energy_tol: float = 0.15,
) -> None:
    """Compare computed spectrum against Fortran reference file.

    Uses interpolation-based comparison to handle small differences in
    frequency grids between Python (standard numpy FFT conventions) and
    Fortran.  This is more robust than point-by-point comparison when
    the grids differ slightly (e.g. odd-n fftshift conventions).

    Args:
        omega: Computed frequency grid (eV)
        sigma: Computed cross-section
        ref_path: Path to Fortran reference file (two columns: omega, sigma)
        rel_rms_threshold: Max allowed RMS difference relative to peak
        peak_energy_tol: Max allowed peak position difference (eV)
    """
    ref = np.loadtxt(ref_path)
    ref_omega = ref[:, 0]
    ref_sigma = ref[:, 1]

    assert len(omega) == len(ref_omega), (
        f"Length mismatch: computed {len(omega)} vs reference {len(ref_omega)}"
    )

    # Peak position comparison
    peak_computed = omega[np.argmax(sigma)]
    peak_ref = ref_omega[np.argmax(ref_sigma)]
    assert abs(peak_computed - peak_ref) < peak_energy_tol, (
        f"Peak mismatch: computed {peak_computed:.3f} vs ref {peak_ref:.3f} eV"
    )

    # Interpolation-based comparison on common grid
    omega_common = np.linspace(
        max(omega[0], ref_omega[0]),
        min(omega[-1], ref_omega[-1]),
        500,
    )
    py_interp = np.interp(omega_common, omega, sigma)
    ref_interp = np.interp(omega_common, ref_omega, ref_sigma)

    # Only compare where signal is significant (> 1% of max)
    mask = ref_interp > ref_interp.max() * 0.01
    if mask.any():
        rms_diff = np.sqrt(np.mean((py_interp[mask] - ref_interp[mask]) ** 2))
        rel_rms = rms_diff / np.max(ref_interp)

        assert rel_rms < rel_rms_threshold, (
            f"RMS difference {rel_rms*100:.2f}% of max > "
            f"threshold {rel_rms_threshold*100:.2f}%"
        )
