"""Tests for 2D spectrum calculation."""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from python_scripts.dynamics_2d import (
    CONST,
    Dipole2D,
    create_dipole_from_file_2d,
    create_constant_dipole_2d,
    PES2D,
    create_harmonic_pes_2d,
    SpectrumConfig2D,
    SpectrumCalculator2D,
    DynamicsConfig2D,
    FullConfig2D,
    GridConfig2D,
    TimeConfig,
    SamplingConfig2D,
    TrajectoryResult2D,
)
from python_scripts.dynamics_2d.io import (
    read_dipole_file_2d,
    write_dipole_file_2d,
    read_dipole_file_2d_raw,
)


class TestDipole2D:
    """Test Dipole2D class."""

    def test_constant_dipole_creation(self):
        """Test creating a constant dipole surface."""
        x1 = np.linspace(0, 1, 11) * 1e-10  # meters
        x2 = np.linspace(0, 1, 11) * 1e-10
        d_value = np.array([1.0, 2.0, 3.0])

        dipole = create_constant_dipole_2d(x1, x2, d_value)

        assert dipole.npoints_x1 == 11
        assert dipole.npoints_x2 == 11
        assert dipole.d.shape == (11, 11, 3)

    def test_constant_dipole_interpolation(self):
        """Test that constant dipole returns constant values."""
        x1 = np.linspace(0, 1, 11) * 1e-10
        x2 = np.linspace(0, 1, 11) * 1e-10
        d_value = np.array([1.0, 2.0, 3.0])

        dipole = create_constant_dipole_2d(x1, x2, d_value)

        # Test scalar interpolation
        d_interp = dipole.dipole(0.5e-10, 0.5e-10)
        np.testing.assert_allclose(d_interp, d_value, rtol=1e-6)

        # Test array interpolation
        x1_test = np.array([0.2e-10, 0.5e-10, 0.8e-10])
        x2_test = np.array([0.2e-10, 0.5e-10, 0.8e-10])
        d_interp = dipole.dipole(x1_test, x2_test)
        for i in range(3):
            np.testing.assert_allclose(d_interp[i], d_value, rtol=1e-5)

    def test_linear_dipole_interpolation(self):
        """Test interpolation of linear dipole surface."""
        x1 = np.linspace(0, 1, 21) * 1e-10
        x2 = np.linspace(0, 1, 21) * 1e-10

        # Create linear dipole d = [x1, x2, x1+x2] in atomic units
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        d = np.zeros((21, 21, 3))
        d[:, :, 0] = X1 / 1e-10  # scale for atomic units
        d[:, :, 1] = X2 / 1e-10
        d[:, :, 2] = (X1 + X2) / 1e-10

        dipole = Dipole2D(x1=x1, x2=x2, d=d)

        # Test interpolation at midpoint
        d_mid = dipole.dipole(0.5e-10, 0.5e-10)
        np.testing.assert_allclose(d_mid[0], 0.5, rtol=1e-3)
        np.testing.assert_allclose(d_mid[1], 0.5, rtol=1e-3)
        np.testing.assert_allclose(d_mid[2], 1.0, rtol=1e-3)

    def test_dipole_derivative_x1(self):
        """Test derivative with respect to x1."""
        x1 = np.linspace(0, 1, 21) * 1e-10
        x2 = np.linspace(0, 1, 21) * 1e-10

        # Create d = [x1^2, 0, 0]
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        d = np.zeros((21, 21, 3))
        d[:, :, 0] = (X1 / 1e-10) ** 2

        dipole = Dipole2D(x1=x1, x2=x2, d=d)

        # Derivative of x1^2 w.r.t. x1 is 2*x1
        dd_dx1 = dipole.dipole_derivative_x1(0.5e-10, 0.5e-10)
        expected = 2 * 0.5 / 1e-10  # in atomic units / meter
        np.testing.assert_allclose(dd_dx1[0], expected, rtol=0.05)

    def test_dipole_derivative_x2(self):
        """Test derivative with respect to x2."""
        x1 = np.linspace(0, 1, 21) * 1e-10
        x2 = np.linspace(0, 1, 21) * 1e-10

        # Create d = [0, x2^2, 0]
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        d = np.zeros((21, 21, 3))
        d[:, :, 1] = (X2 / 1e-10) ** 2

        dipole = Dipole2D(x1=x1, x2=x2, d=d)

        # Derivative of x2^2 w.r.t. x2 is 2*x2
        dd_dx2 = dipole.dipole_derivative_x2(0.5e-10, 0.5e-10)
        expected = 2 * 0.5 / 1e-10
        np.testing.assert_allclose(dd_dx2[1], expected, rtol=0.05)

    def test_dipole_file_roundtrip(self):
        """Test writing and reading dipole file."""
        x1 = np.linspace(0, 1, 11) * 1e-10
        x2 = np.linspace(0, 1, 11) * 1e-10

        d = np.zeros((11, 11, 3))
        d[:, :, 0] = 1.0
        d[:, :, 1] = 2.0
        d[:, :, 2] = 3.0

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole_2d.dat"

            # Write
            write_dipole_file_2d(filepath, x1, x2, d, position_units="angstrom")

            # Read
            x1_read, x2_read, d_read = read_dipole_file_2d(
                filepath, position_units="angstrom"
            )

            np.testing.assert_allclose(x1_read, x1, rtol=1e-6)
            np.testing.assert_allclose(x2_read, x2, rtol=1e-6)
            np.testing.assert_allclose(d_read, d, rtol=1e-6)


class TestDipoleFileIO:
    """Test dipole file I/O functions."""

    def test_read_dipole_file_2d_c_order(self):
        """Test reading 2D dipole file with C-style ordering."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"

            # Create test file with C ordering (x2 varies fastest)
            with open(filepath, "w") as f:
                for i, x1 in enumerate([0.5, 1.0]):
                    for j, x2 in enumerate([1.0, 2.0, 3.0]):
                        d_x = i + 1.0
                        d_y = j + 1.0
                        d_z = i + j + 1.0
                        f.write(f"{x1} {x2} {d_x} {d_y} {d_z}\n")

            x1, x2, d = read_dipole_file_2d(filepath, index_order="C")

            assert len(x1) == 2
            assert len(x2) == 3
            assert d.shape == (2, 3, 3)
            np.testing.assert_allclose(d[0, 0, 0], 1.0)  # d_x at (0.5, 1.0)
            np.testing.assert_allclose(d[1, 2, 1], 3.0)  # d_y at (1.0, 3.0)

    def test_read_dipole_file_2d_f_order(self):
        """Test reading 2D dipole file with Fortran-style ordering."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"

            # Create test file with F ordering (x1 varies fastest)
            with open(filepath, "w") as f:
                for j, x2 in enumerate([1.0, 2.0, 3.0]):
                    for i, x1 in enumerate([0.5, 1.0]):
                        d_x = i + 1.0
                        d_y = j + 1.0
                        d_z = i + j + 1.0
                        f.write(f"{x1} {x2} {d_x} {d_y} {d_z}\n")

            x1, x2, d = read_dipole_file_2d(filepath, index_order="F")

            assert len(x1) == 2
            assert len(x2) == 3
            assert d.shape == (2, 3, 3)
            np.testing.assert_allclose(d[0, 0, 0], 1.0)
            np.testing.assert_allclose(d[1, 2, 1], 3.0)

    def test_read_dipole_file_2d_single_component(self):
        """Test reading 2D dipole file with single component (|d|^2 format)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole.dat"

            # Create test file with |d|^2 values (C ordering)
            # File contains squared values, reader takes sqrt
            with open(filepath, "w") as f:
                for i, x1 in enumerate([0.5, 1.0]):
                    for j, x2 in enumerate([1.0, 2.0, 3.0]):
                        d_squared = (i + j + 1.0) ** 2  # Write |d|^2
                        f.write(f"{x1} {x2} {d_squared}\n")

            x1, x2, d = read_dipole_file_2d(
                filepath, index_order="C", dipole_components=1
            )

            assert len(x1) == 2
            assert len(x2) == 3
            assert d.shape == (2, 3, 3)
            # x and y components should be 0
            np.testing.assert_allclose(d[:, :, 0], 0.0)
            np.testing.assert_allclose(d[:, :, 1], 0.0)
            # z component should have sqrt of file values
            np.testing.assert_allclose(d[0, 0, 2], 1.0)  # sqrt(1^2) at (0.5, 1.0)
            np.testing.assert_allclose(d[1, 2, 2], 4.0)  # sqrt(16) at (1.0, 3.0)

    def test_write_read_single_component_roundtrip(self):
        """Test writing and reading single-component dipole file.

        Write writes |d|^2, read takes sqrt, so roundtrip preserves values.
        """
        x1 = np.linspace(0, 1, 5) * 1e-10
        x2 = np.linspace(0, 1, 5) * 1e-10

        # Create dipole with only z-component
        d = np.zeros((5, 5, 3))
        d[:, :, 2] = 2.5  # z-component

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "dipole_1comp.dat"

            # Write with single component (writes z^2 = 6.25)
            write_dipole_file_2d(
                filepath, x1, x2, d,
                position_units="angstrom",
                dipole_components=1,
            )

            # Read back with single component (reads sqrt(6.25) = 2.5)
            x1_read, x2_read, d_read = read_dipole_file_2d(
                filepath,
                position_units="angstrom",
                dipole_components=1,
            )

            np.testing.assert_allclose(x1_read, x1, rtol=1e-6)
            np.testing.assert_allclose(x2_read, x2, rtol=1e-6)
            # x and y should be 0, z should match original
            np.testing.assert_allclose(d_read[:, :, 0], 0.0)
            np.testing.assert_allclose(d_read[:, :, 1], 0.0)
            np.testing.assert_allclose(d_read[:, :, 2], 2.5, rtol=1e-6)


class TestSpectrumConfig2D:
    """Test SpectrumConfig2D dataclass."""

    def test_gamma_hwhm_conversion(self):
        """Test FWHM to HWHM conversion."""
        config = SpectrumConfig2D(gamma_fwhm=0.36)
        assert config.gamma_hwhm == pytest.approx(0.18)

    def test_default_values(self):
        """Test default configuration values."""
        config = SpectrumConfig2D(gamma_fwhm=0.18)
        assert config.dipole_mode == "DIPOLE"
        assert config.compatibility_mode == "standard"
        assert config.pes_final_list == []
        assert config.dipole_final_list == []


class TestSpectrumCalculator2D:
    """Test SpectrumCalculator2D class."""

    @pytest.fixture
    def harmonic_setup(self, tmp_path):
        """Create test files for harmonic oscillator case."""
        # Grid parameters
        n = 21
        x1 = np.linspace(0.5, 2.5, n) * 1e-10  # meters
        x2 = np.linspace(0.5, 2.5, n) * 1e-10

        # Equilibrium positions
        x1_0 = 1.5e-10
        x2_0 = 1.5e-10

        # Force constants (reasonable values for molecules)
        k = 50.0  # N/m

        # Energy offsets for different states (in eV, converted to J)
        E_initial = 0.0
        E_intermediate = 520.0 * CONST.eV  # ~520 eV above ground
        E_final = 510.0 * CONST.eV  # ~510 eV above ground

        # Create PES files
        pes_initial = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k, k, E_initial)
        pes_intermediate = create_harmonic_pes_2d(
            x1, x2, x1_0, x2_0, k, k, E_intermediate
        )
        pes_final = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k, k, E_final)

        # Write PES files
        from python_scripts.dynamics_2d.io import write_pes_file_2d

        write_pes_file_2d(
            tmp_path / "pes_initial_2d.dat", x1, x2, pes_initial.E
        )
        write_pes_file_2d(
            tmp_path / "pes_intermediate_2d.dat", x1, x2, pes_intermediate.E
        )
        write_pes_file_2d(
            tmp_path / "pes_final_2d.dat", x1, x2, pes_final.E
        )

        # Create constant dipole files
        d = np.ones((n, n, 3))
        write_dipole_file_2d(tmp_path / "dipole_2d.dat", x1, x2, d)

        # Create dynamics config
        dynamics_config = DynamicsConfig2D(
            mu1=1.0,  # amu
            mu2=1.0,
            grid_x1=GridConfig2D(start=0.5, dx=0.1, npoints=n),
            grid_x2=GridConfig2D(start=0.5, dx=0.1, npoints=n),
            time=TimeConfig(dt=0.1, nsteps=128),
            sampling=SamplingConfig2D(mode=1, npoints_x1=2, npoints_x2=2, npoints_p1=2, npoints_p2=2),
            pes_dynamics=tmp_path / "pes_intermediate_2d.dat",
            pes_initial=tmp_path / "pes_initial_2d.dat",
        )

        # Create spectrum config
        spectrum_config = SpectrumConfig2D(
            gamma_fwhm=0.36,
            dipole_mode="FC",
            pes_final_list=[tmp_path / "pes_final_2d.dat"],
            dipole_final_list=[tmp_path / "dipole_2d.dat"],
        )

        # Create full config
        full_config = FullConfig2D(
            dynamics2d=dynamics_config,
            spectrum=spectrum_config,
        )

        return full_config, tmp_path

    def test_load_surfaces(self, harmonic_setup):
        """Test loading PES and dipole surfaces."""
        full_config, _ = harmonic_setup

        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        assert calc.pes_n is not None
        assert len(calc.pes_f) == 1
        assert len(calc.dipoles) == 1

    def test_interpolate_along_trajectory(self, harmonic_setup):
        """Test interpolation of energies and dipoles along trajectory."""
        full_config, _ = harmonic_setup

        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        # Create a simple trajectory
        nsteps = 64
        time = np.linspace(0, 10, nsteps) * 1e-15  # seconds
        x1 = 1.5e-10 + 0.1e-10 * np.sin(time * 1e14)
        x2 = 1.5e-10 + 0.1e-10 * np.cos(time * 1e14)
        v1 = np.zeros(nsteps)
        v2 = np.zeros(nsteps)
        a1 = np.zeros(nsteps)
        a2 = np.zeros(nsteps)

        traj = TrajectoryResult2D(
            time=time, x1=x1, x2=x2, v1=v1, v2=v2, a1=a1, a2=a2,
            x1_0=x1[0], x2_0=x2[0], p1_0=0.0, p2_0=0.0
        )

        interp = calc.interpolate_along_trajectory(traj)

        assert "E_n" in interp
        assert "E_f" in interp
        assert "D_fn" in interp
        assert interp["E_n"].shape == (nsteps,)
        assert interp["E_f"].shape == (1, nsteps)
        assert interp["D_fn"].shape == (1, nsteps, 3)

    def test_compute_mean_transition_energy(self, harmonic_setup):
        """Test mean transition energy calculation."""
        full_config, _ = harmonic_setup

        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        # Create dummy trajectory
        nsteps = 64
        time = np.linspace(0, 10, nsteps) * 1e-15
        x1 = np.ones(nsteps) * 1.5e-10
        x2 = np.ones(nsteps) * 1.5e-10
        traj = TrajectoryResult2D(
            time=time, x1=x1, x2=x2,
            v1=np.zeros(nsteps), v2=np.zeros(nsteps),
            a1=np.zeros(nsteps), a2=np.zeros(nsteps),
            x1_0=x1[0], x2_0=x2[0], p1_0=0.0, p2_0=0.0
        )

        E_mean = calc.compute_mean_transition_energy([traj])

        # Should be approximately E_intermediate - E_final = 520 - 510 = 10 eV
        assert E_mean == pytest.approx(10.0, rel=0.01)

    def test_compute_spectrum(self, harmonic_setup):
        """Test full spectrum calculation."""
        full_config, tmp_path = harmonic_setup

        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        # Create simple trajectories
        nsteps = 128
        time = np.linspace(0, 10, nsteps) * 1e-15
        trajectories = []
        for _ in range(4):
            x1 = 1.5e-10 + np.random.randn(nsteps) * 0.01e-10
            x2 = 1.5e-10 + np.random.randn(nsteps) * 0.01e-10
            traj = TrajectoryResult2D(
                time=time, x1=x1, x2=x2,
                v1=np.zeros(nsteps), v2=np.zeros(nsteps),
                a1=np.zeros(nsteps), a2=np.zeros(nsteps),
                x1_0=x1[0], x2_0=x2[0], p1_0=0.0, p2_0=0.0
            )
            trajectories.append(traj)

        result = calc.compute_spectrum(trajectories)

        assert result.omega is not None
        assert result.sigma_tot is not None
        assert result.sigma_f is not None
        assert len(result.omega) == nsteps
        assert len(result.sigma_tot) == nsteps
        assert result.sigma_f.shape == (1, nsteps)
        assert result.n_trajectories == 4

    def test_save_results(self, harmonic_setup):
        """Test saving spectrum results."""
        full_config, tmp_path = harmonic_setup

        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        # Create simple trajectory
        nsteps = 64
        time = np.linspace(0, 10, nsteps) * 1e-15
        x1 = np.ones(nsteps) * 1.5e-10
        x2 = np.ones(nsteps) * 1.5e-10
        traj = TrajectoryResult2D(
            time=time, x1=x1, x2=x2,
            v1=np.zeros(nsteps), v2=np.zeros(nsteps),
            a1=np.zeros(nsteps), a2=np.zeros(nsteps),
            x1_0=x1[0], x2_0=x2[0], p1_0=0.0, p2_0=0.0
        )

        result = calc.compute_spectrum([traj])
        calc.save_results(result, tmp_path / "output")

        # Check files exist
        assert (tmp_path / "output" / f"{full_config.dynamics2d.outfile}_sigma.dat").exists()


class TestDipoleModes:
    """Test different dipole modes."""

    @pytest.fixture
    def basic_setup(self, tmp_path):
        """Create basic test setup."""
        n = 11
        x1 = np.linspace(0.5, 2.5, n) * 1e-10
        x2 = np.linspace(0.5, 2.5, n) * 1e-10
        x1_0, x2_0 = 1.5e-10, 1.5e-10
        k = 50.0

        from python_scripts.dynamics_2d.io import write_pes_file_2d

        pes_initial = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k, k, 0.0)
        pes_intermediate = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k, k, 520 * CONST.eV)
        pes_final = create_harmonic_pes_2d(x1, x2, x1_0, x2_0, k, k, 510 * CONST.eV)

        write_pes_file_2d(tmp_path / "pes_i.dat", x1, x2, pes_initial.E)
        write_pes_file_2d(tmp_path / "pes_n.dat", x1, x2, pes_intermediate.E)
        write_pes_file_2d(tmp_path / "pes_f.dat", x1, x2, pes_final.E)

        # Non-constant dipole for testing
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        d = np.zeros((n, n, 3))
        d[:, :, 0] = 1.0 + 0.1 * ((X1 - x1_0) / 1e-10)
        d[:, :, 1] = 1.0 + 0.1 * ((X2 - x2_0) / 1e-10)
        d[:, :, 2] = 1.0
        write_dipole_file_2d(tmp_path / "dipole.dat", x1, x2, d)

        dynamics_config = DynamicsConfig2D(
            mu1=1.0, mu2=1.0,
            grid_x1=GridConfig2D(start=0.5, dx=0.2, npoints=n),
            grid_x2=GridConfig2D(start=0.5, dx=0.2, npoints=n),
            time=TimeConfig(dt=0.1, nsteps=64),
            sampling=SamplingConfig2D(),
            pes_dynamics=tmp_path / "pes_n.dat",
            pes_initial=tmp_path / "pes_i.dat",
        )

        return dynamics_config, tmp_path

    def test_fc_mode(self, basic_setup):
        """Test Franck-Condon mode (constant dipole = 1)."""
        dynamics_config, tmp_path = basic_setup

        spectrum_config = SpectrumConfig2D(
            gamma_fwhm=0.36,
            dipole_mode="FC",
            pes_final_list=[tmp_path / "pes_f.dat"],
            dipole_final_list=[],  # Not needed for FC mode
        )

        full_config = FullConfig2D(dynamics2d=dynamics_config, spectrum=spectrum_config)
        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        assert len(calc.dipoles) == 1
        # FC mode should give constant dipole = 1
        d = calc.dipoles[0].dipole(1.5e-10, 1.5e-10)
        np.testing.assert_allclose(d, [1.0, 1.0, 1.0])

    def test_dipole_mode(self, basic_setup):
        """Test full position-dependent dipole mode."""
        dynamics_config, tmp_path = basic_setup

        spectrum_config = SpectrumConfig2D(
            gamma_fwhm=0.36,
            dipole_mode="DIPOLE",
            pes_final_list=[tmp_path / "pes_f.dat"],
            dipole_final_list=[tmp_path / "dipole.dat"],
        )

        full_config = FullConfig2D(dynamics2d=dynamics_config, spectrum=spectrum_config)
        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        assert len(calc.dipoles) == 1
        # Dipole should vary with position
        d_center = calc.dipoles[0].dipole(1.5e-10, 1.5e-10)
        d_off = calc.dipoles[0].dipole(1.7e-10, 1.7e-10)
        # Should not be exactly equal
        assert not np.allclose(d_center, d_off)

    def test_dipole_x0_mode(self, basic_setup):
        """Test dipole frozen at equilibrium mode."""
        dynamics_config, tmp_path = basic_setup

        spectrum_config = SpectrumConfig2D(
            gamma_fwhm=0.36,
            dipole_mode="DIPOLE_X0",
            pes_final_list=[tmp_path / "pes_f.dat"],
            dipole_final_list=[tmp_path / "dipole.dat"],
        )

        full_config = FullConfig2D(dynamics2d=dynamics_config, spectrum=spectrum_config)
        calc = SpectrumCalculator2D(full_config)
        calc.load_surfaces()

        # Create trajectory away from equilibrium
        nsteps = 32
        time = np.linspace(0, 5, nsteps) * 1e-15
        x1 = np.linspace(1.5, 2.0, nsteps) * 1e-10  # Moves away from equilibrium
        x2 = np.linspace(1.5, 2.0, nsteps) * 1e-10

        traj = TrajectoryResult2D(
            time=time, x1=x1, x2=x2,
            v1=np.zeros(nsteps), v2=np.zeros(nsteps),
            a1=np.zeros(nsteps), a2=np.zeros(nsteps),
            x1_0=x1[0], x2_0=x2[0], p1_0=0.0, p2_0=0.0
        )

        interp = calc.interpolate_along_trajectory(traj)

        # In DIPOLE_X0 mode, dipole should be constant along trajectory
        # (frozen at equilibrium value)
        D_fn = interp["D_fn"][0]  # First (only) final state
        # All time steps should have same dipole
        for i in range(1, nsteps):
            np.testing.assert_allclose(D_fn[i], D_fn[0], rtol=1e-10)
