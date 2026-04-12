"""SCKH testsuite: spectrum from pre-computed trajectory files.

Reads SCKH trajectory files from the Fortran testsuite, computes the
spectrum using compute_spectrum_from_sckh, and compares against the
Fortran reference output.

This test mirrors tests/testsuite/SCKH/ and also serves as an example
of computing spectra from SCKH trajectory files.

Example usage (as a script)::

    from dynamics_1d import (
        read_sckh_trajectory,
        read_sckh_trajectory_list,
        compute_spectrum_from_sckh,
    )

    # Read trajectory list
    traj_entries = read_sckh_trajectory_list("trajectories.dat")

    # Read each trajectory file
    trajs = [read_sckh_trajectory(path) for path, weight in traj_entries]

    # Compute spectrum
    result = compute_spectrum_from_sckh(trajs, gamma_hwhm=0.09)

    # result.omega  -> frequency grid (eV)
    # result.sigma_tot -> total cross-section
    # result.sigma_f -> per-final-state cross-section
"""

from pathlib import Path

import numpy as np
import pytest

from dynamics_1d import (
    read_sckh_trajectory,
    read_sckh_trajectory_list,
    compute_spectrum_from_sckh,
)

from .conftest import TESTSUITE_FORTRAN_DIR, assert_spectrum_matches_ref


SCKH_DIR = TESTSUITE_FORTRAN_DIR / "SCKH"
INPUT_DIR = SCKH_DIR / "input"
REF_DIR = SCKH_DIR / "ref"

# Parameters from SCKH/input/input.txt
GAMMA_FWHM = 0.18  # eV
GAMMA_HWHM = GAMMA_FWHM / 2.0
NTSTEPS = 161
NFINAL = 20
NTRAJ = 3


@pytest.fixture
def sckh_trajectories():
    """Load all SCKH trajectory files from the Fortran testsuite."""
    traj_entries = read_sckh_trajectory_list(INPUT_DIR / "trajectories.dat")
    assert len(traj_entries) == NTRAJ

    trajs = [read_sckh_trajectory(path) for path, _weight in traj_entries]
    return trajs


class TestSCKHTrajectoryReading:
    """Test that SCKH trajectory files are read correctly."""

    def test_trajectory_count(self, sckh_trajectories):
        """Should read the correct number of trajectories."""
        assert len(sckh_trajectories) == NTRAJ

    def test_trajectory_dimensions(self, sckh_trajectories):
        """Each trajectory should have the expected shape."""
        for traj in sckh_trajectories:
            assert traj.ntsteps == NTSTEPS
            assert traj.n_final == NFINAL
            assert traj.D_fn.shape == (NFINAL, NTSTEPS, 3)

    def test_time_grid(self, sckh_trajectories):
        """Time grid should start at 0 and have correct spacing."""
        traj = sckh_trajectories[0]
        assert traj.time[0] == 0.0
        dt = traj.time[1] - traj.time[0]
        np.testing.assert_allclose(dt, 0.25e-15, rtol=1e-6)

    def test_energies_in_ev(self, sckh_trajectories):
        """E_n and E_f should be in eV (large magnitude for core states)."""
        traj = sckh_trajectories[0]
        # Core-excited state energy ~-125 Hartree -> ~-3400 eV
        assert abs(traj.E_n[0]) > 1000
        # Final states also in eV
        assert abs(traj.E_f[0, 0]) > 1000

    def test_trajectory_list_weights(self):
        """Trajectory list should contain weights."""
        entries = read_sckh_trajectory_list(INPUT_DIR / "trajectories.dat")
        for _path, weight in entries:
            assert weight == 1.0


class TestSCKHSpectrum:
    """Test spectrum computation from SCKH trajectory files."""

    def test_spectrum_matches_fortran_reference(self, sckh_trajectories):
        """Computed spectrum should match Fortran reference output.

        This is the main validation test, using the same comparison
        metric as the Fortran testsuite (my_diff_sum_ref).
        """
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )

        assert_spectrum_matches_ref(
            result.omega,
            result.sigma_tot,
            REF_DIR / "pentamer_XES_sigma.dat",
        )

    def test_spectrum_shape(self, sckh_trajectories):
        """Spectrum arrays should have correct dimensions."""
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )

        assert len(result.omega) == NTSTEPS
        assert len(result.sigma_tot) == NTSTEPS
        assert result.sigma_f.shape == (NFINAL, NTSTEPS)
        assert result.n_trajectories == NTRAJ

    def test_spectrum_is_normalized(self, sckh_trajectories):
        """Spectrum integral should be close to 1."""
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )

        integral = np.trapezoid(result.sigma_tot, result.omega)
        assert 0.9 < integral < 1.1, f"Integral={integral:.4f}, expected ~1.0"

    def test_peak_in_expected_range(self, sckh_trajectories):
        """Spectrum peak should be in the X-ray emission energy range."""
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )

        peak_idx = np.argmax(result.sigma_tot)
        peak_energy = result.omega[peak_idx]
        # XES for this system: peaks around 520-530 eV
        assert 518 < peak_energy < 530, f"Peak at {peak_energy:.2f} eV"

    def test_e_mean_reasonable(self, sckh_trajectories):
        """E_mean should be in the expected energy range."""
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )
        # E_mean is computed from E_n - E_f[last], should be ~520 eV
        assert 515 < result.E_mean < 530, f"E_mean={result.E_mean:.2f} eV"

    def test_per_final_state_spectra_match_reference(self, sckh_trajectories):
        """Per-final-state spectra should match Fortran references."""
        result = compute_spectrum_from_sckh(
            sckh_trajectories, gamma_hwhm=GAMMA_HWHM
        )

        for j in range(NFINAL):
            ref_file = REF_DIR / f"pentamer_XES_sigma_final_{j + 1}.dat"
            if ref_file.exists():
                assert_spectrum_matches_ref(
                    result.omega,
                    result.sigma_f[j],
                    ref_file,
                )
