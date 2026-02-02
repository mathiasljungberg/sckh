"""Integration tests for trajectory module."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from python_scripts.dynamics_1d.trajectory import (
    TrajectoryResult,
    EnsembleResult,
    DynamicsRunner,
)
from python_scripts.dynamics_1d.config import (
    DynamicsConfig,
    GridConfig,
    TimeConfig,
    SamplingConfig,
)
from python_scripts.dynamics_1d.spectrum_config import FullConfig, SpectrumConfig
from python_scripts.dynamics_1d.pes import create_harmonic_pes
from python_scripts.dynamics_1d.vibrational import solve_vibrational
from python_scripts.dynamics_1d.constants import CONST


@pytest.fixture
def harmonic_config(tmp_path):
    """Create a configuration for harmonic oscillator dynamics."""
    # Create temporary PES files
    x_ang = np.linspace(0.5, 1.5, 77)
    x_m = x_ang * 1e-10

    # Harmonic potential centered at 1.0 Angstrom
    x0 = 1.0e-10
    k = 500.0
    E_hartree = 0.5 * k * (x_m - x0) ** 2 / CONST.hartree

    # Write PES file
    pes_file = tmp_path / "harmonic_pes.dat"
    np.savetxt(pes_file, np.column_stack([x_ang, E_hartree]))

    dynamics_config = DynamicsConfig(
        mu=1.0,  # 1 amu
        grid=GridConfig(start=0.5, dx=0.025 / 1.5, npoints=77),
        time=TimeConfig(dt=0.5, nsteps=100),
        sampling=SamplingConfig(mode=1, npoints_x=3, npoints_mom=3),
        pes_initial=pes_file,
        pes_dynamics=pes_file,
        units="angstrom",
        outfile="test_dynamics",
    )

    # Minimal spectrum config (not used in dynamics tests)
    spectrum_config = SpectrumConfig(
        gamma_fwhm=0.1,
        dipole_mode="FC",
        pes_final_list=[pes_file],
    )

    return FullConfig(dynamics=dynamics_config, spectrum=spectrum_config)


class TestTrajectoryResult:
    """Tests for TrajectoryResult dataclass."""

    def test_trajectory_result_creation(self):
        """Should create TrajectoryResult with all fields."""
        nsteps = 50
        result = TrajectoryResult(
            time=np.linspace(0, 1e-12, nsteps),
            x=np.zeros(nsteps),
            v=np.ones(nsteps) * 1e-5,
            a=np.zeros(nsteps),
            x0=0.0,
            p0=1e-5 * CONST.u,
        )

        assert len(result.time) == nsteps
        assert result.x0 == 0.0


class TestEnsembleResult:
    """Tests for EnsembleResult dataclass."""

    def test_ensemble_properties(self):
        """Test EnsembleResult property accessors."""
        nsteps = 20
        traj1 = TrajectoryResult(
            time=np.linspace(0, 1e-12, nsteps),
            x=np.linspace(0, 1e-10, nsteps),
            v=np.ones(nsteps) * 1e-5,
            a=np.zeros(nsteps),
            x0=0.0,
            p0=0.0,
        )
        traj2 = TrajectoryResult(
            time=np.linspace(0, 1e-12, nsteps),
            x=np.linspace(0.5e-10, 1.5e-10, nsteps),
            v=np.ones(nsteps) * 2e-5,
            a=np.zeros(nsteps),
            x0=0.5e-10,
            p0=0.0,
        )

        config = DynamicsConfig(
            mu=1.0,
            grid=GridConfig(start=0.0, dx=0.01, npoints=100),
            time=TimeConfig(dt=0.1, nsteps=nsteps),
            sampling=SamplingConfig(),
            pes_initial=Path("dummy.dat"),
            pes_dynamics=Path("dummy.dat"),
        )

        result = EnsembleResult(
            trajectories=[traj1, traj2],
            config=config,
            ground_state_energy=1e-20,
            x_grid=np.linspace(0, 1e-9, 100),
            psi_ground=np.ones(100),
        )

        assert result.n_trajectories == 2
        assert result.n_steps == nsteps

    def test_get_positions_at_time(self):
        """Should return positions from all trajectories at given time step."""
        nsteps = 10
        traj1 = TrajectoryResult(
            time=np.linspace(0, 1e-12, nsteps),
            x=np.full(nsteps, 1e-10),
            v=np.zeros(nsteps),
            a=np.zeros(nsteps),
            x0=1e-10,
            p0=0.0,
        )
        traj2 = TrajectoryResult(
            time=np.linspace(0, 1e-12, nsteps),
            x=np.full(nsteps, 2e-10),
            v=np.zeros(nsteps),
            a=np.zeros(nsteps),
            x0=2e-10,
            p0=0.0,
        )

        result = EnsembleResult(
            trajectories=[traj1, traj2],
            config=None,
            ground_state_energy=0.0,
            x_grid=np.array([]),
            psi_ground=np.array([]),
        )

        positions = result.get_positions_at_time(5)

        np.testing.assert_array_equal(positions, [1e-10, 2e-10])


class TestDynamicsRunner:
    """Integration tests for DynamicsRunner."""

    def test_runner_initialization(self, harmonic_config):
        """Runner should initialize from config."""
        runner = DynamicsRunner(harmonic_config)

        assert runner.mass_SI == 1.0 * CONST.u
        assert runner.dt_SI == 0.5e-15
        assert len(runner.x_grid) == 77

    def test_solve_ground_state(self, harmonic_config):
        """Runner should solve for ground state wavefunction."""
        runner = DynamicsRunner(harmonic_config)
        runner.solve_ground_state()

        # Ground state should be set
        assert runner._psi_ground is not None
        assert runner._eigenvalue_ground is not None

        # Ground state should be normalized
        dx = runner.x_grid[1] - runner.x_grid[0]
        norm = np.sum(runner._psi_ground**2) * dx
        np.testing.assert_allclose(norm, 1.0, rtol=1e-6)

    def test_run_returns_ensemble_result(self, harmonic_config):
        """Run should return EnsembleResult with trajectories."""
        runner = DynamicsRunner(harmonic_config)
        result = runner.run()

        assert isinstance(result, EnsembleResult)
        # 3 x 3 = 9 trajectories for mode=1
        assert result.n_trajectories == 9
        assert result.n_steps == 100

    def test_trajectories_start_at_sampled_positions(self, harmonic_config):
        """Trajectories should start at sampled initial conditions."""
        runner = DynamicsRunner(harmonic_config)
        result = runner.run()

        for traj in result.trajectories:
            # First position should match x0
            np.testing.assert_allclose(traj.x[0], traj.x0, rtol=1e-12)

            # First velocity should match p0/m
            v0_expected = traj.p0 / runner.mass_SI
            np.testing.assert_allclose(traj.v[0], v0_expected, rtol=1e-12)

    def test_energy_conservation_in_trajectories(self, harmonic_config):
        """Energy should be approximately conserved in each trajectory."""
        runner = DynamicsRunner(harmonic_config)
        result = runner.run()

        for traj in result.trajectories:
            # Compute kinetic energy
            E_kin = 0.5 * runner.mass_SI * traj.v**2

            # Compute potential energy
            E_pot = runner.pes.energy(traj.x)

            # Total energy
            E_total = E_kin + E_pot

            # Should be conserved (allow for 5% drift due to spline force approximation)
            E_drift = (E_total[-1] - E_total[0]) / E_total[0]
            assert abs(E_drift) < 0.05

    def test_save_results(self, harmonic_config, tmp_path):
        """Should save results to files."""
        runner = DynamicsRunner(harmonic_config)
        result = runner.run()

        output_dir = tmp_path / "output"
        runner.save_results(result, output_dir)

        # Check files exist
        assert (output_dir / "test_dynamics_ground_state.dat").exists()
        assert (output_dir / "test_dynamics_initial_distribution.dat").exists()
        assert (output_dir / "test_dynamics_traj_0000.dat").exists()

    def test_manual_ground_state(self, harmonic_config):
        """Should accept manually set ground state."""
        runner = DynamicsRunner(harmonic_config)

        # Create a simple Gaussian ground state
        sigma = 0.05e-10
        x0 = 1.0e-10
        psi = np.exp(-((runner.x_grid - x0) ** 2) / (2 * sigma**2))
        dx = runner.x_grid[1] - runner.x_grid[0]
        psi /= np.sqrt(np.sum(psi**2) * dx)

        runner.set_ground_state(psi, eigenvalue=-1e-20)

        # Should be able to run
        result = runner.run()
        assert result.n_trajectories > 0
