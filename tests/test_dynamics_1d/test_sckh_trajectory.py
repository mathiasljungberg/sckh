"""Tests for SCKH trajectory I/O and spectrum computation."""

import numpy as np
import pytest
from pathlib import Path

from dynamics_1d.constants import CONST
from dynamics_1d.sckh_trajectory import (
    SCKHTrajectory,
    sckh_trajectory_from_dynamics_1d,
)
from dynamics_1d.io import read_sckh_trajectory, write_sckh_trajectory
from dynamics_1d.spectrum import compute_spectrum_from_sckh
from dynamics_1d.pes import PES1D, create_harmonic_pes
from dynamics_1d.dipole import create_constant_dipole
from dynamics_1d.trajectory import TrajectoryResult


def _make_sample_sckh_file(filepath, ntsteps=10, nfinal=2):
    """Write a sample SCKH trajectory file in Fortran format."""
    with open(filepath, 'w') as f:
        for i in range(ntsteps):
            time_fs = i * 0.5  # 0.5 fs steps
            E_gs = -76.0 + i * 0.001  # Hartree
            E_n = -75.0 + i * 0.002  # Hartree
            E_IP1s = 0.0  # Hartree
            f.write(f"  {time_fs:20.10E}\n")
            f.write(f"  {E_gs:20.10E}\n")
            f.write(f"  {E_n:20.10E}\n")
            f.write(f"  {E_IP1s:20.10E}\n")
            f.write(f"ntrans {nfinal:5d}\n")
            for j in range(nfinal):
                E_trans = 0.5 + j * 0.1 + i * 0.001  # Hartree
                D_x = 0.1 * (j + 1)
                D_y = 0.2 * (j + 1)
                D_z = 0.3 * (j + 1)
                f.write(
                    f"{E_trans:15.8f}{D_x:17.8E}{D_y:17.8E}{D_z:17.8E}\n"
                )


class TestSCKHTrajectory:
    """Tests for SCKHTrajectory dataclass."""

    def test_properties(self):
        """Test n_final and ntsteps properties."""
        ntsteps = 20
        nfinal = 3
        sckh = SCKHTrajectory(
            time=np.zeros(ntsteps),
            E_gs=np.zeros(ntsteps),
            E_n=np.zeros(ntsteps),
            E_f=np.zeros((nfinal, ntsteps)),
            D_fn=np.zeros((nfinal, ntsteps, 3)),
            E_IP1s=np.zeros(ntsteps),
        )
        assert sckh.n_final == nfinal
        assert sckh.ntsteps == ntsteps


class TestReadSCKHTrajectory:
    """Tests for reading SCKH trajectory files."""

    def test_read_basic(self, tmp_path):
        """Read a sample file and verify basic structure."""
        filepath = tmp_path / "traj.dat"
        ntsteps = 10
        nfinal = 2
        _make_sample_sckh_file(filepath, ntsteps=ntsteps, nfinal=nfinal)

        sckh = read_sckh_trajectory(filepath)

        assert sckh.ntsteps == ntsteps
        assert sckh.n_final == nfinal
        assert sckh.time.shape == (ntsteps,)
        assert sckh.E_gs.shape == (ntsteps,)
        assert sckh.E_n.shape == (ntsteps,)
        assert sckh.E_f.shape == (nfinal, ntsteps)
        assert sckh.D_fn.shape == (nfinal, ntsteps, 3)
        assert sckh.E_IP1s.shape == (ntsteps,)

    def test_time_conversion(self, tmp_path):
        """Time should be converted from fs to seconds."""
        filepath = tmp_path / "traj.dat"
        _make_sample_sckh_file(filepath, ntsteps=5)
        sckh = read_sckh_trajectory(filepath)

        # First timestep: 0.0 fs = 0.0 s
        assert sckh.time[0] == 0.0
        # Second timestep: 0.5 fs = 0.5e-15 s
        np.testing.assert_allclose(sckh.time[1], 0.5e-15, rtol=1e-10)

    def test_energy_conversion(self, tmp_path):
        """E_n and E_f should be in eV, E_gs in Hartree."""
        filepath = tmp_path / "traj.dat"
        _make_sample_sckh_file(filepath, ntsteps=5, nfinal=2)
        sckh = read_sckh_trajectory(filepath)

        # E_gs should remain in Hartree (around -76)
        assert abs(sckh.E_gs[0] - (-76.0)) < 0.01

        # E_n should be in eV (around -75 Hartree * 27.2 eV/Hartree ≈ -2040 eV)
        expected_E_n_eV = -75.0 * CONST.hartree2eV
        np.testing.assert_allclose(sckh.E_n[0], expected_E_n_eV, rtol=1e-6)

    def test_e_f_derivation(self, tmp_path):
        """E_f should be derived correctly: E_gs - E_trans + E_IP1s*(eV/Hartree)."""
        filepath = tmp_path / "traj.dat"
        _make_sample_sckh_file(filepath, ntsteps=3, nfinal=1)
        sckh = read_sckh_trajectory(filepath)

        # At step 0: E_gs=-76.0, E_trans=0.5, E_IP1s=0.0
        # E_f_hartree = -76.0 - 0.5 + 0.0 = -76.5
        # E_f_eV = -76.5 * hartree2eV
        expected_eV = -76.5 * CONST.hartree2eV
        np.testing.assert_allclose(sckh.E_f[0, 0], expected_eV, rtol=1e-6)

    def test_dipole_values(self, tmp_path):
        """Dipole moments should be read correctly."""
        filepath = tmp_path / "traj.dat"
        _make_sample_sckh_file(filepath, ntsteps=3, nfinal=2)
        sckh = read_sckh_trajectory(filepath)

        # First final state at step 0: D = (0.1, 0.2, 0.3)
        np.testing.assert_allclose(sckh.D_fn[0, 0, :], [0.1, 0.2, 0.3], rtol=1e-6)
        # Second final state at step 0: D = (0.2, 0.4, 0.6)
        np.testing.assert_allclose(sckh.D_fn[1, 0, :], [0.2, 0.4, 0.6], rtol=1e-6)

    def test_fortran_d_exponent(self, tmp_path):
        """Should handle Fortran D-exponent notation."""
        filepath = tmp_path / "traj.dat"
        with open(filepath, 'w') as f:
            f.write("  0.0000000000D+00\n")  # time
            f.write(" -7.6000000000D+01\n")  # E_gs
            f.write(" -7.5000000000D+01\n")  # E_n
            f.write("  0.0000000000D+00\n")  # E_IP1s
            f.write("ntrans     1\n")
            f.write("     0.50000000  1.00000000D-01  2.00000000D-01  3.00000000D-01\n")

        sckh = read_sckh_trajectory(filepath)
        assert sckh.ntsteps == 1
        np.testing.assert_allclose(sckh.E_gs[0], -76.0, rtol=1e-10)

    def test_nonzero_e_ip1s(self, tmp_path):
        """E_IP1s contribution should be included in E_f derivation."""
        filepath = tmp_path / "traj.dat"
        # Write a file with nonzero E_IP1s
        with open(filepath, 'w') as f:
            E_IP1s = 540.0 / CONST.hartree2eV  # ~540 eV in Hartree ≈ 19.84
            f.write("  0.0000000000E+00\n")
            f.write(f"  {-76.0:20.10E}\n")
            f.write(f"  {-75.0:20.10E}\n")
            f.write(f"  {E_IP1s:20.10E}\n")
            f.write("ntrans     1\n")
            f.write("     0.50000000  1.00000000E-01  2.00000000E-01  3.00000000E-01\n")

        sckh = read_sckh_trajectory(filepath)

        # E_f_hartree = E_gs - E_trans + E_IP1s * (eV/Hartree)
        eV_over_hartree = CONST.eV / CONST.hartree
        E_f_hartree = -76.0 - 0.5 + E_IP1s * eV_over_hartree
        expected_eV = E_f_hartree * CONST.hartree2eV
        np.testing.assert_allclose(sckh.E_f[0, 0], expected_eV, rtol=1e-6)


class TestWriteReadRoundtrip:
    """Tests for write -> read roundtrip."""

    def test_roundtrip(self, tmp_path):
        """Write and read back should preserve data."""
        ntsteps = 8
        nfinal = 3

        original = SCKHTrajectory(
            time=np.linspace(0, 5e-15, ntsteps),
            E_gs=np.linspace(-76.0, -75.9, ntsteps),
            E_n=np.linspace(-2040.0, -2038.0, ntsteps),  # eV
            E_f=np.random.default_rng(42).uniform(-2080, -2070, (nfinal, ntsteps)),
            D_fn=np.random.default_rng(43).uniform(-1, 1, (nfinal, ntsteps, 3)),
            E_IP1s=np.zeros(ntsteps),
        )

        filepath = tmp_path / "roundtrip.dat"
        write_sckh_trajectory(filepath, original)
        recovered = read_sckh_trajectory(filepath)

        np.testing.assert_allclose(recovered.time, original.time, rtol=1e-8)
        np.testing.assert_allclose(recovered.E_gs, original.E_gs, rtol=1e-6)
        np.testing.assert_allclose(recovered.E_n, original.E_n, rtol=1e-5)
        np.testing.assert_allclose(recovered.E_f, original.E_f, rtol=1e-5)
        np.testing.assert_allclose(recovered.D_fn, original.D_fn, rtol=1e-5)
        np.testing.assert_allclose(recovered.E_IP1s, original.E_IP1s, atol=1e-10)

    def test_roundtrip_single_final_state(self, tmp_path):
        """Roundtrip with 1 final state."""
        ntsteps = 5
        original = SCKHTrajectory(
            time=np.arange(ntsteps) * 1e-15,
            E_gs=np.full(ntsteps, -76.0),
            E_n=np.full(ntsteps, -2040.0),
            E_f=np.full((1, ntsteps), -2070.0),
            D_fn=np.ones((1, ntsteps, 3)) * 0.5,
            E_IP1s=np.zeros(ntsteps),
        )

        filepath = tmp_path / "roundtrip_1.dat"
        write_sckh_trajectory(filepath, original)
        recovered = read_sckh_trajectory(filepath)

        assert recovered.n_final == 1
        np.testing.assert_allclose(recovered.E_n, original.E_n, rtol=1e-5)


class TestSCKHFromDynamics1D:
    """Tests for creating SCKH trajectories from 1D dynamics."""

    def test_basic_conversion(self, position_grid, harmonic_params):
        """Convert a 1D trajectory to SCKH format."""
        x_grid = position_grid
        k = harmonic_params["k"]
        x0 = harmonic_params["x0"]

        # Create PES surfaces (harmonic with different offsets)
        pes_gs = create_harmonic_pes(x_grid, x0=x0, k=k, E0=0.0)
        pes_n = create_harmonic_pes(x_grid, x0=x0, k=k, E0=10.0 * CONST.eV)
        pes_f1 = create_harmonic_pes(x_grid, x0=x0, k=k, E0=5.0 * CONST.eV)

        # Create constant dipole
        dipole = create_constant_dipole(x_grid, np.array([1.0, 0.0, 0.0]))

        # Create a simple trajectory
        nsteps = 50
        dt = 0.1e-15
        time = np.arange(nsteps) * dt
        x = x0 + 0.01e-10 * np.sin(harmonic_params["omega"] * time)
        traj = TrajectoryResult(
            time=time, x=x,
            v=np.zeros(nsteps), a=np.zeros(nsteps),
            x0=x[0], p0=0.0,
        )

        sckh = sckh_trajectory_from_dynamics_1d(
            traj, pes_gs, pes_n, [pes_f1], [dipole]
        )

        assert sckh.ntsteps == nsteps
        assert sckh.n_final == 1
        np.testing.assert_array_equal(sckh.time, time)

        # E_n should be ~10 eV (from PES offset)
        np.testing.assert_allclose(sckh.E_n, 10.0, atol=0.1)

    def test_e_ip1s_is_zero(self, position_grid, harmonic_params):
        """E_IP1s should be zero for dynamics-generated trajectories."""
        x_grid = position_grid
        k = harmonic_params["k"]
        x0 = harmonic_params["x0"]

        pes_gs = create_harmonic_pes(x_grid, x0=x0, k=k, E0=0.0)
        pes_n = create_harmonic_pes(x_grid, x0=x0, k=k, E0=10.0 * CONST.eV)
        pes_f = create_harmonic_pes(x_grid, x0=x0, k=k, E0=5.0 * CONST.eV)
        dipole = create_constant_dipole(x_grid, np.array([1.0, 1.0, 1.0]))

        nsteps = 20
        traj = TrajectoryResult(
            time=np.arange(nsteps) * 1e-15,
            x=np.full(nsteps, x0),
            v=np.zeros(nsteps), a=np.zeros(nsteps),
            x0=x0, p0=0.0,
        )

        sckh = sckh_trajectory_from_dynamics_1d(
            traj, pes_gs, pes_n, [pes_f], [dipole]
        )
        np.testing.assert_array_equal(sckh.E_IP1s, 0.0)

    def test_energies_match_pes_evaluation(self, position_grid, harmonic_params):
        """SCKH energies should match direct PES evaluation."""
        x_grid = position_grid
        k = harmonic_params["k"]
        x0 = harmonic_params["x0"]

        E0_gs = 0.0
        E0_n = 10.0 * CONST.eV
        E0_f = 5.0 * CONST.eV

        pes_gs = create_harmonic_pes(x_grid, x0=x0, k=k, E0=E0_gs)
        pes_n = create_harmonic_pes(x_grid, x0=x0, k=k, E0=E0_n)
        pes_f = create_harmonic_pes(x_grid, x0=x0, k=k, E0=E0_f)
        dipole = create_constant_dipole(x_grid, np.array([1.0, 0.0, 0.0]))

        nsteps = 10
        x = np.linspace(-0.1e-10, 0.1e-10, nsteps)
        traj = TrajectoryResult(
            time=np.arange(nsteps) * 1e-15,
            x=x, v=np.zeros(nsteps), a=np.zeros(nsteps),
            x0=x[0], p0=0.0,
        )

        sckh = sckh_trajectory_from_dynamics_1d(
            traj, pes_gs, pes_n, [pes_f], [dipole]
        )

        # Check E_n matches PES evaluation converted to eV
        expected_E_n_eV = pes_n.energy(x) / CONST.eV
        np.testing.assert_allclose(sckh.E_n, expected_E_n_eV, rtol=1e-10)

        # Check E_f matches
        expected_E_f_eV = pes_f.energy(x) / CONST.eV
        np.testing.assert_allclose(sckh.E_f[0], expected_E_f_eV, rtol=1e-10)

    def test_fc_dipole_mode(self, position_grid, harmonic_params):
        """FC mode should set all dipoles to 1."""
        x_grid = position_grid
        k = harmonic_params["k"]
        x0 = harmonic_params["x0"]

        pes_gs = create_harmonic_pes(x_grid, x0=x0, k=k, E0=0.0)
        pes_n = create_harmonic_pes(x_grid, x0=x0, k=k, E0=10.0 * CONST.eV)
        pes_f = create_harmonic_pes(x_grid, x0=x0, k=k, E0=5.0 * CONST.eV)

        nsteps = 5
        traj = TrajectoryResult(
            time=np.arange(nsteps) * 1e-15,
            x=np.full(nsteps, x0),
            v=np.zeros(nsteps), a=np.zeros(nsteps),
            x0=x0, p0=0.0,
        )

        sckh = sckh_trajectory_from_dynamics_1d(
            traj, pes_gs, pes_n, [pes_f], [], dipole_mode="FC"
        )
        np.testing.assert_array_equal(sckh.D_fn, 1.0)


class TestComputeSpectrumFromSCKH:
    """Tests for computing spectra from SCKH trajectories."""

    def _make_constant_sckh(self, E_n_eV=530.0, E_f_eV=520.0, nsteps=512):
        """Create SCKH trajectory with constant energies."""
        dt = 0.1e-15
        time = np.arange(nsteps) * dt
        return SCKHTrajectory(
            time=time,
            E_gs=np.full(nsteps, -76.0),
            E_n=np.full(nsteps, E_n_eV),
            E_f=np.full((1, nsteps), E_f_eV),
            D_fn=np.ones((1, nsteps, 3)),
            E_IP1s=np.zeros(nsteps),
        )

    def test_basic_spectrum(self):
        """Spectrum from constant energies should produce a peak."""
        sckh = self._make_constant_sckh()
        result = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)

        assert result.omega is not None
        assert result.sigma_tot is not None
        assert len(result.omega) == sckh.ntsteps
        assert result.n_trajectories == 1
        # Spectrum should be non-negative
        assert np.all(result.sigma_tot >= -1e-15)

    def test_peak_at_transition_energy(self):
        """Peak should be at E_n - E_f."""
        E_n = 530.0
        E_f = 520.0
        sckh = self._make_constant_sckh(E_n_eV=E_n, E_f_eV=E_f, nsteps=1024)
        result = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)

        # Find peak position
        peak_idx = np.argmax(result.sigma_tot)
        peak_energy = result.omega[peak_idx]

        # Peak should be near E_n - E_f = 10 eV
        np.testing.assert_allclose(peak_energy, E_n - E_f, atol=0.5)

    def test_e_mean_auto_computation(self):
        """Auto-computed E_mean should match time average of E_n - E_f[last]."""
        sckh = self._make_constant_sckh(E_n_eV=530.0, E_f_eV=520.0)
        result = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)

        expected_E_mean = np.mean(sckh.E_n - sckh.E_f[-1])
        np.testing.assert_allclose(result.E_mean, expected_E_mean, rtol=1e-10)

    def test_e_mean_explicit(self):
        """Explicit E_mean should be used when provided."""
        sckh = self._make_constant_sckh()
        E_mean_explicit = 15.0
        result = compute_spectrum_from_sckh(
            [sckh], gamma_hwhm=0.5, E_mean=E_mean_explicit
        )
        assert result.E_mean == E_mean_explicit

    def test_multiple_trajectories(self):
        """Multiple trajectories should be accumulated."""
        sckh1 = self._make_constant_sckh()
        sckh2 = self._make_constant_sckh()
        result = compute_spectrum_from_sckh([sckh1, sckh2], gamma_hwhm=0.5)

        assert result.n_trajectories == 2

    def test_per_final_state(self):
        """sigma_f should have correct shape."""
        nsteps = 128
        nfinal = 3
        dt = 0.1e-15
        time = np.arange(nsteps) * dt
        sckh = SCKHTrajectory(
            time=time,
            E_gs=np.full(nsteps, -76.0),
            E_n=np.full(nsteps, 530.0),
            E_f=np.full((nfinal, nsteps), 520.0),
            D_fn=np.ones((nfinal, nsteps, 3)),
            E_IP1s=np.zeros(nsteps),
        )

        result = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)
        assert result.sigma_f.shape == (nfinal, nsteps)

    def test_normalization(self):
        """Spectrum should be normalized."""
        sckh = self._make_constant_sckh(nsteps=512)
        result = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)

        # Total integral should be close to 1 (within normalization convention)
        dw = result.omega[1] - result.omega[0]
        integral = np.sum(result.sigma_tot) * dw
        # Should be order 1 (exact value depends on convention)
        assert 0.1 < integral < 10.0

    def test_empty_trajectories_raises(self):
        """Should raise ValueError for empty trajectory list."""
        with pytest.raises(ValueError, match="No trajectories"):
            compute_spectrum_from_sckh([], gamma_hwhm=0.5)


class TestEndToEnd:
    """End-to-end tests: dynamics -> SCKH -> write -> read -> spectrum."""

    def test_dynamics_to_spectrum(self, tmp_path, position_grid, harmonic_params):
        """Full pipeline: create trajectory, convert to SCKH, compute spectrum."""
        x_grid = position_grid
        k = harmonic_params["k"]
        x0 = harmonic_params["x0"]

        # Create surfaces
        pes_gs = create_harmonic_pes(x_grid, x0=x0, k=k, E0=0.0)
        pes_n = create_harmonic_pes(x_grid, x0=x0, k=k, E0=10.0 * CONST.eV)
        pes_f = create_harmonic_pes(x_grid, x0=x0, k=k, E0=5.0 * CONST.eV)
        dipole = create_constant_dipole(x_grid, np.array([1.0, 1.0, 1.0]))

        # Create trajectory at equilibrium
        nsteps = 128
        traj = TrajectoryResult(
            time=np.arange(nsteps) * 0.1e-15,
            x=np.full(nsteps, x0),
            v=np.zeros(nsteps), a=np.zeros(nsteps),
            x0=x0, p0=0.0,
        )

        # Convert to SCKH
        sckh = sckh_trajectory_from_dynamics_1d(
            traj, pes_gs, pes_n, [pes_f], [dipole]
        )

        # Write and read back
        filepath = tmp_path / "traj_sckh.dat"
        write_sckh_trajectory(filepath, sckh)
        sckh_read = read_sckh_trajectory(filepath)

        # Compute spectrum from both (should match)
        result1 = compute_spectrum_from_sckh([sckh], gamma_hwhm=0.5)
        result2 = compute_spectrum_from_sckh([sckh_read], gamma_hwhm=0.5)

        np.testing.assert_allclose(
            result1.sigma_tot, result2.sigma_tot, rtol=1e-4
        )
        np.testing.assert_allclose(result1.omega, result2.omega, rtol=1e-4)
