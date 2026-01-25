"""Tests for 2D trajectory management and runner."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_2d.trajectory import (
    TrajectoryResult2D,
    EnsembleResult2D,
    DynamicsRunner2D,
)
from python_scripts.dynamics_2d.config import (
    DynamicsConfig2D,
    FullConfig2D,
    GridConfig2D,
    TimeConfig,
    SamplingConfig2D,
)
from python_scripts.dynamics_2d.vibrational import ProductGroundState
from python_scripts.dynamics_2d.io import write_pes_file_2d


class TestTrajectoryResult2D:
    """Tests for TrajectoryResult2D dataclass."""

    def test_create_result(self):
        """Test creating a TrajectoryResult2D."""
        nsteps = 100
        time = np.linspace(0, 1e-12, nsteps)
        x1 = np.sin(time * 1e15)
        x2 = np.cos(time * 1e15)
        v1 = np.cos(time * 1e15)
        v2 = -np.sin(time * 1e15)
        a1 = -np.sin(time * 1e15)
        a2 = -np.cos(time * 1e15)

        result = TrajectoryResult2D(
            time=time,
            x1=x1,
            x2=x2,
            v1=v1,
            v2=v2,
            a1=a1,
            a2=a2,
            x1_0=0.0,
            x2_0=1.0,
            p1_0=1.0,
            p2_0=0.0,
        )

        assert len(result.time) == nsteps
        assert result.x1_0 == 0.0
        assert result.x2_0 == 1.0


class TestEnsembleResult2D:
    """Tests for EnsembleResult2D dataclass."""

    def test_ensemble_properties(self):
        """Test EnsembleResult2D property accessors."""
        nsteps = 50
        time = np.linspace(0, 1e-12, nsteps)

        # Create mock trajectories
        trajectories = []
        for i in range(5):
            traj = TrajectoryResult2D(
                time=time,
                x1=np.random.randn(nsteps) * 1e-10,
                x2=np.random.randn(nsteps) * 1e-10,
                v1=np.random.randn(nsteps) * 100,
                v2=np.random.randn(nsteps) * 100,
                a1=np.random.randn(nsteps) * 1e15,
                a2=np.random.randn(nsteps) * 1e15,
                x1_0=float(i) * 1e-11,
                x2_0=float(i) * 1e-11,
                p1_0=0.0,
                p2_0=0.0,
            )
            trajectories.append(traj)

        # Create mock ground state
        x1_grid = np.linspace(-0.5, 0.5, 101) * 1e-10
        x2_grid = np.linspace(-0.5, 0.5, 101) * 1e-10
        ground_state = ProductGroundState(
            E1=1e-20,
            E2=1e-20,
            E_total=2e-20,
            psi1=np.exp(-x1_grid**2 / (2 * (0.1e-10) ** 2)),
            psi2=np.exp(-x2_grid**2 / (2 * (0.1e-10) ** 2)),
            x1_grid=x1_grid,
            x2_grid=x2_grid,
        )

        # Create mock config (minimal)
        config = DynamicsConfig2D(
            mu1=1.0,
            mu2=2.0,
            grid_x1=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
            grid_x2=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
            time=TimeConfig(dt=1.0, nsteps=nsteps),
            sampling=SamplingConfig2D(),
            pes_dynamics=Path("dummy.dat"),
            pes_initial=Path("dummy_initial.dat"),
        )

        ensemble = EnsembleResult2D(
            trajectories=trajectories,
            config=config,
            ground_state=ground_state,
            x1_grid=x1_grid,
            x2_grid=x2_grid,
            x1_eq=0.0,
            x2_eq=0.0,
        )

        assert ensemble.n_trajectories == 5
        assert ensemble.n_steps == nsteps

        # Test position extraction
        x1_at_0, x2_at_0 = ensemble.get_positions_at_time(0)
        assert len(x1_at_0) == 5
        assert len(x2_at_0) == 5


class TestDynamicsRunner2D:
    """Tests for DynamicsRunner2D class."""

    @pytest.fixture
    def temp_pes_files(self, position_grid_2d, harmonic_params_2d):
        """Create temporary PES files for testing."""
        x1, x2 = position_grid_2d
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create 2D PES file for dynamics
            X1, X2 = np.meshgrid(x1, x2, indexing="ij")
            E_2d = 0.5 * k1 * X1**2 + 0.5 * k2 * X2**2
            pes_2d_path = tmpdir / "pes_2d.dat"
            write_pes_file_2d(
                pes_2d_path, x1, x2, E_2d, position_units="angstrom"
            )

            # Create 2D PES file for initial state (same as dynamics for test)
            pes_initial_path = tmpdir / "pes_initial_2d.dat"
            write_pes_file_2d(
                pes_initial_path, x1, x2, E_2d, position_units="angstrom"
            )

            yield {
                "dir": tmpdir,
                "pes_2d": pes_2d_path,
                "pes_initial": pes_initial_path,
            }

    def test_runner_initialization(self, temp_pes_files):
        """Test DynamicsRunner2D initialization."""
        config = FullConfig2D(
            dynamics2d=DynamicsConfig2D(
                mu1=1.0,
                mu2=2.0,
                grid_x1=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                grid_x2=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                time=TimeConfig(dt=0.1, nsteps=100),
                sampling=SamplingConfig2D(
                    npoints_x1=2, npoints_x2=2, npoints_p1=2, npoints_p2=2
                ),
                pes_dynamics=temp_pes_files["pes_2d"],
                pes_initial=temp_pes_files["pes_initial"],
                position_units="angstrom",
            )
        )

        runner = DynamicsRunner2D(config)

        # Check mass conversion
        assert runner.m1_SI == pytest.approx(1.0 * CONST.u)
        assert runner.m2_SI == pytest.approx(2.0 * CONST.u)

        # Check time step conversion
        assert runner.dt_SI == pytest.approx(0.1e-15)

    def test_solve_ground_state(self, temp_pes_files):
        """Test ground state solver with PES slicing."""
        config = FullConfig2D(
            dynamics2d=DynamicsConfig2D(
                mu1=1.0,
                mu2=2.0,
                grid_x1=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                grid_x2=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                time=TimeConfig(dt=0.1, nsteps=100),
                sampling=SamplingConfig2D(
                    npoints_x1=2, npoints_x2=2, npoints_p1=2, npoints_p2=2
                ),
                pes_dynamics=temp_pes_files["pes_2d"],
                pes_initial=temp_pes_files["pes_initial"],
                position_units="angstrom",
            )
        )

        runner = DynamicsRunner2D(config)
        runner.solve_ground_state()

        gs = runner.ground_state
        assert gs.E_total > 0  # Zero-point energy
        assert len(gs.psi1) == 101
        assert len(gs.psi2) == 101

        # Check equilibrium was found (should be at origin for harmonic)
        assert runner.x1_eq == pytest.approx(0.0, abs=1e-12)
        assert runner.x2_eq == pytest.approx(0.0, abs=1e-12)

    def test_run_trajectories(self, temp_pes_files):
        """Test running ensemble of trajectories."""
        config = FullConfig2D(
            dynamics2d=DynamicsConfig2D(
                mu1=1.0,
                mu2=2.0,
                grid_x1=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                grid_x2=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                time=TimeConfig(dt=0.1, nsteps=50),
                sampling=SamplingConfig2D(
                    npoints_x1=2,
                    npoints_x2=2,
                    npoints_p1=2,
                    npoints_p2=2,
                    mode=1,
                ),
                pes_dynamics=temp_pes_files["pes_2d"],
                pes_initial=temp_pes_files["pes_initial"],
                position_units="angstrom",
            )
        )

        runner = DynamicsRunner2D(config)
        result = runner.run()

        # Should have 2*2*2*2 = 16 trajectories
        assert result.n_trajectories == 16
        assert result.n_steps == 50

        # Check trajectory data
        for traj in result.trajectories:
            assert len(traj.time) == 50
            assert len(traj.x1) == 50
            assert len(traj.x2) == 50

    def test_energy_conservation_in_run(self, temp_pes_files):
        """Test that trajectories conserve energy."""
        config = FullConfig2D(
            dynamics2d=DynamicsConfig2D(
                mu1=1.0,
                mu2=2.0,
                grid_x1=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                grid_x2=GridConfig2D(start=-0.5, dx=0.01, npoints=101),
                time=TimeConfig(dt=0.05, nsteps=200),  # Smaller time step
                sampling=SamplingConfig2D(
                    npoints_x1=1,
                    npoints_x2=1,
                    npoints_p1=1,
                    npoints_p2=1,
                    mode=1,
                ),
                pes_dynamics=temp_pes_files["pes_2d"],
                pes_initial=temp_pes_files["pes_initial"],
                position_units="angstrom",
            )
        )

        runner = DynamicsRunner2D(config)
        result = runner.run()

        # Check energy conservation for first trajectory
        traj = result.trajectories[0]
        from python_scripts.dynamics_2d.integrators import compute_total_energy_2d

        E_total = compute_total_energy_2d(
            traj.x1,
            traj.x2,
            traj.v1,
            traj.v2,
            runner.m1_SI,
            runner.m2_SI,
            runner.pes.energy,
        )

        E_initial = E_total[0]
        E_drift = np.abs(E_total - E_initial) / E_initial

        # Should conserve energy to within 0.1%
        assert np.max(E_drift) < 0.001
