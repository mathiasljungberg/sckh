"""Tests for full KH spectrum calculation."""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from python_scripts.kh_1d.config import (
    GridConfig, BroadeningConfig, FrequencyGridConfig, KHConfig,
    load_config, save_config,
)
from python_scripts.kh_1d.spectrum import KHCalculator, KHResult, VibrationalState
from python_scripts.dynamics_1d.constants import CONST


class TestConfig:
    """Tests for configuration classes."""

    def test_grid_config_odd_npoints(self):
        """GridConfig should require odd npoints."""
        with pytest.raises(ValueError, match="odd"):
            GridConfig(start=0.5, dx=0.025, npoints=78)

        # Should work with odd npoints
        cfg = GridConfig(start=0.5, dx=0.025, npoints=77)
        assert cfg.npoints == 77

    def test_broadening_config_hwhm(self):
        """BroadeningConfig should compute HWHM from FWHM."""
        cfg = BroadeningConfig(gamma_fwhm=0.4)
        assert cfg.gamma_hwhm == 0.2

    def test_frequency_grid(self):
        """FrequencyGridConfig should generate correct omega array."""
        cfg = FrequencyGridConfig(omega_start=517, omega_end=528, n_omega=101)
        omega = cfg.get_omega_array()

        assert len(omega) == 101
        assert omega[0] == 517
        assert omega[-1] == 528
        assert np.allclose(np.diff(omega), np.diff(omega)[0])  # Uniform

    def test_kh_config_validation(self):
        """KHConfig should validate dipole mode and list lengths."""
        # Invalid dipole mode
        with pytest.raises(ValueError, match="dipole_mode"):
            KHConfig(
                mu=1.0,
                grid=GridConfig(0.5, 0.025, 77),
                broadening=BroadeningConfig(0.3),
                frequency=FrequencyGridConfig(517, 528, 100),
                pes_initial=Path("a.dat"),
                pes_intermediate=Path("b.dat"),
                pes_final_list=[Path("c.dat")],
                dipole_final_list=[Path("d.dat")],
                dipole_mode="INVALID",
            )

        # Mismatched list lengths
        with pytest.raises(ValueError, match="same length"):
            KHConfig(
                mu=1.0,
                grid=GridConfig(0.5, 0.025, 77),
                broadening=BroadeningConfig(0.3),
                frequency=FrequencyGridConfig(517, 528, 100),
                pes_initial=Path("a.dat"),
                pes_intermediate=Path("b.dat"),
                pes_final_list=[Path("c.dat"), Path("e.dat")],
                dipole_final_list=[Path("d.dat")],
            )


class TestVibrationalState:
    """Tests for VibrationalState dataclass."""

    def test_eigenvalues_ev(self):
        """eigenvalues_eV should convert correctly."""
        state = VibrationalState(
            eigenvalues=np.array([CONST.eV, 2 * CONST.eV, 3 * CONST.eV]),
            eigenvectors=np.eye(3),
            x=np.linspace(0, 1, 3),
        )

        np.testing.assert_allclose(state.eigenvalues_eV, [1.0, 2.0, 3.0], rtol=1e-10)

    def test_n_states(self):
        """n_states should return correct count."""
        state = VibrationalState(
            eigenvalues=np.zeros(5),
            eigenvectors=np.zeros((10, 5)),
            x=np.linspace(0, 1, 10),
        )
        assert state.n_states == 5


class TestKHCalculatorHarmonic:
    """Tests using harmonic oscillator PES files."""

    @pytest.fixture
    def harmonic_files(self, tmp_path):
        """Create temporary PES and dipole files for harmonic oscillators."""
        # Grid parameters
        npoints = 77
        x_start = 0.5  # Angstrom
        dx = 0.025  # Angstrom
        x = np.array([x_start + i * dx for i in range(npoints)])

        # Harmonic PES parameters (in Hartree for files)
        k_au = 0.5  # Force constant in au
        x0_initial = 0.96  # Angstrom
        x0_intermediate = 1.00  # Shifted
        x0_final = 0.98

        # Convert to SI for force constant
        # k_au is in Hartree/bohr^2
        # k_SI = k_au * Hartree / bohr^2

        # Initial state PES
        E_initial = 0.5 * k_au * ((x - x0_initial) / (CONST.bohr * 1e10)) ** 2
        pes_i_file = tmp_path / "pes_initial.dat"
        np.savetxt(pes_i_file, np.column_stack([x, E_initial]))

        # Intermediate state PES (core-excited, shifted up)
        E_offset = 530.0 / CONST.hartree2eV  # 530 eV in Hartree
        E_intermediate = E_offset + 0.5 * k_au * ((x - x0_intermediate) / (CONST.bohr * 1e10)) ** 2
        pes_n_file = tmp_path / "pes_intermediate.dat"
        np.savetxt(pes_n_file, np.column_stack([x, E_intermediate]))

        # Final state PES
        E_final_offset = 5.0 / CONST.hartree2eV  # 5 eV in Hartree
        E_final = E_final_offset + 0.5 * k_au * ((x - x0_final) / (CONST.bohr * 1e10)) ** 2
        pes_f_file = tmp_path / "pes_final.dat"
        np.savetxt(pes_f_file, np.column_stack([x, E_final]))

        # Dipole file (constant z-component)
        dipole = np.zeros((npoints, 4))
        dipole[:, 0] = x
        dipole[:, 3] = 0.05  # Small z-component in atomic units
        dipole_file = tmp_path / "dipole.dat"
        np.savetxt(dipole_file, dipole)

        return {
            "pes_initial": pes_i_file,
            "pes_intermediate": pes_n_file,
            "pes_final": pes_f_file,
            "dipole": dipole_file,
            "grid": (x_start, dx, npoints),
        }

    def test_calculator_workflow(self, harmonic_files):
        """Test full calculator workflow with harmonic potentials."""
        x_start, dx, npoints = harmonic_files["grid"]

        config = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=520, omega_end=530, n_omega=101),
            pes_initial=harmonic_files["pes_initial"],
            pes_intermediate=harmonic_files["pes_intermediate"],
            pes_final_list=[harmonic_files["pes_final"]],
            dipole_final_list=[harmonic_files["dipole"]],
            dipole_mode="FC",
            n_vib_states=10,
        )

        calc = KHCalculator(config)
        result = calc.run(verbose=False)

        # Check result structure
        assert isinstance(result, KHResult)
        assert len(result.omega) == 101
        assert len(result.sigma_tot) == 101
        assert result.sigma_per_final_elec.shape == (1, 101)

        # Check eigenvalues are reasonable
        assert len(result.eigenvalues_initial) == 10
        assert len(result.eigenvalues_intermediate) == 10
        assert len(result.eigenvalues_final) == 1
        assert len(result.eigenvalues_final[0]) == 10

        # Check spectrum is normalized
        d_omega = result.omega[1] - result.omega[0]
        integral = np.sum(result.sigma_tot) * d_omega
        assert abs(integral - 1.0) < 0.01  # Within 1%

        # Check spectrum has a peak
        assert np.max(result.sigma_tot) > np.mean(result.sigma_tot) * 2

    def test_harmonic_eigenvalues(self, harmonic_files):
        """Harmonic oscillator eigenvalues should follow E_n = hw(n+1/2)."""
        x_start, dx, npoints = harmonic_files["grid"]

        config = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=520, omega_end=530, n_omega=11),
            pes_initial=harmonic_files["pes_initial"],
            pes_intermediate=harmonic_files["pes_intermediate"],
            pes_final_list=[harmonic_files["pes_final"]],
            dipole_final_list=[harmonic_files["dipole"]],
            n_vib_states=5,
        )

        calc = KHCalculator(config)
        calc.load_surfaces()
        calc.solve_vibrational_problem()

        # Get eigenvalues
        E = calc._vib_i.eigenvalues_eV

        # Energy spacings should be approximately constant for harmonic oscillator
        spacings = np.diff(E)
        assert np.std(spacings) / np.mean(spacings) < 0.1  # Within 10%


class TestKHCalculatorWithTestsuite:
    """Integration tests using Fortran testsuite files."""

    @pytest.fixture
    def testsuite_available(self, testsuite_input_path):
        """Check if testsuite files are available."""
        if not testsuite_input_path.exists():
            pytest.skip("Testsuite not found")
        return testsuite_input_path

    def test_load_testsuite_files(self, testsuite_available):
        """Test loading actual testsuite PES and dipole files."""
        input_path = testsuite_available

        # Check required files exist
        pes_initial = input_path / "test_energy_tot.dat"
        pes_intermediate = input_path / "test_energy_tot_exc.dat"
        pes_final = input_path / "pesfile_4.dat"
        dipole_final = input_path / "dipolefile_4.dat"

        if not all(p.exists() for p in [pes_initial, pes_intermediate, pes_final, dipole_final]):
            pytest.skip("Some testsuite files missing")

        config = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=0.5, dx=0.025, npoints=77),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=517, omega_end=528, n_omega=100),
            pes_initial=pes_initial,
            pes_intermediate=pes_intermediate,
            pes_final_list=[pes_final],
            dipole_final_list=[dipole_final],
            dipole_mode="FC",
        )

        calc = KHCalculator(config)
        calc.load_surfaces()

        # Check surfaces loaded correctly
        assert calc._pes_i is not None
        assert calc._pes_n is not None
        assert len(calc._pes_f_list) == 1
        assert len(calc._dipole_f_list) == 1

    def test_full_testsuite_calculation(self, testsuite_available, testsuite_ref_path):
        """Compare full calculation with Fortran reference output."""
        input_path = testsuite_available

        # Collect all PES and dipole files
        pes_final_list = []
        dipole_final_list = []
        for i in [9, 8, 7, 6, 5, 4]:
            pes_file = input_path / f"pesfile_{i}.dat"
            dipole_file = input_path / f"dipolefile_{i}.dat"
            if pes_file.exists() and dipole_file.exists():
                pes_final_list.append(pes_file)
                dipole_final_list.append(dipole_file)

        if len(pes_final_list) < 6:
            pytest.skip("Not all testsuite final state files available")

        pes_initial = input_path / "test_energy_tot.dat"
        pes_intermediate = input_path / "test_energy_tot_exc.dat"

        if not pes_initial.exists() or not pes_intermediate.exists():
            pytest.skip("Initial/intermediate PES files missing")

        config = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=0.5, dx=0.025, npoints=77),
            broadening=BroadeningConfig(gamma_fwhm=0.3),  # Need to check Fortran value
            frequency=FrequencyGridConfig(omega_start=517, omega_end=528, n_omega=100),
            pes_initial=pes_initial,
            pes_intermediate=pes_intermediate,
            pes_final_list=pes_final_list,
            dipole_final_list=dipole_final_list,
            dipole_mode="FC",
            n_vib_states=20,  # Limit to 20 vibrational states for faster test
        )

        calc = KHCalculator(config)
        result = calc.run(verbose=False)

        # Basic sanity checks
        assert result.sigma_tot.min() >= 0
        assert not np.any(np.isnan(result.sigma_tot))
        assert not np.any(np.isinf(result.sigma_tot))

        # Check normalization
        d_omega = result.omega[1] - result.omega[0]
        integral = np.sum(result.sigma_tot) * d_omega
        assert abs(integral - 1.0) < 0.05

        # Check we have the right number of final states
        assert result.sigma_per_final_elec.shape[0] == 6


class TestKHResultIO:
    """Tests for result I/O."""

    def test_write_spectrum(self, tmp_path):
        """Test writing spectrum to file."""
        result = KHResult(
            omega=np.linspace(500, 510, 11),
            sigma_tot=np.ones(11) / 10,
            sigma_per_final_elec=np.ones((2, 11)) / 20,
            eigenvalues_initial=np.array([0.0, 0.5]),
            eigenvalues_intermediate=np.array([530.0, 530.5]),
            eigenvalues_final=[np.array([5.0, 5.5])],
            D_ni=np.zeros((2, 3)),
            D_fn_list=[np.zeros((2, 2, 3))],
        )

        outfile = tmp_path / "spectrum.dat"
        result.write_spectrum(outfile)

        assert outfile.exists()

        # Read back and check
        data = np.loadtxt(outfile)
        assert data.shape == (11, 2)
        np.testing.assert_allclose(data[:, 0], result.omega)

    def test_write_per_final_spectrum(self, tmp_path):
        """Test writing per-final-state spectra."""
        result = KHResult(
            omega=np.linspace(500, 510, 11),
            sigma_tot=np.ones(11) / 10,
            sigma_per_final_elec=np.ones((3, 11)) / 30,
            eigenvalues_initial=np.array([0.0]),
            eigenvalues_intermediate=np.array([530.0]),
            eigenvalues_final=[np.array([5.0])],
            D_ni=np.zeros((1, 3)),
            D_fn_list=[np.zeros((1, 1, 3))],
        )

        base_path = tmp_path / "test_spectrum"
        result.write_spectrum_per_final(base_path)

        # Should create 3 files
        for j in range(3):
            filepath = tmp_path / f"test_spectrum_states_{j + 1}.dat"
            assert filepath.exists()


class TestYAMLConfig:
    """Tests for YAML config loading and saving."""

    def _make_config(self, tmp_path):
        """Create a KHConfig for testing."""
        return KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=0.5, dx=0.025, npoints=77),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=517, omega_end=528, n_omega=1000),
            pes_initial=Path("pes_initial.dat"),
            pes_intermediate=Path("pes_intermediate.dat"),
            pes_final_list=[Path("pesfile_9.dat"), Path("pesfile_8.dat")],
            dipole_final_list=[Path("dipolefile_9.dat"), Path("dipolefile_8.dat")],
            dipole_mode="FC",
            n_vib_states=20,
            pes_lp_corr=None,
            energy_column_initial=1,
            energy_column_intermediate=1,
            energy_column_final=1,
        )

    def test_save_and_load_roundtrip(self, tmp_path):
        """save_config then load_config should reproduce the original config."""
        original = self._make_config(tmp_path)
        yaml_path = tmp_path / "config.yaml"

        save_config(original, yaml_path)
        loaded = load_config(yaml_path)

        assert loaded.mu == original.mu
        assert loaded.grid.start == original.grid.start
        assert loaded.grid.dx == original.grid.dx
        assert loaded.grid.npoints == original.grid.npoints
        assert loaded.broadening.gamma_fwhm == original.broadening.gamma_fwhm
        assert loaded.frequency.omega_start == original.frequency.omega_start
        assert loaded.frequency.omega_end == original.frequency.omega_end
        assert loaded.frequency.n_omega == original.frequency.n_omega
        assert loaded.pes_initial == original.pes_initial
        assert loaded.pes_intermediate == original.pes_intermediate
        assert loaded.pes_final_list == original.pes_final_list
        assert loaded.dipole_final_list == original.dipole_final_list
        assert loaded.dipole_mode == original.dipole_mode
        assert loaded.n_vib_states == original.n_vib_states
        assert loaded.pes_lp_corr == original.pes_lp_corr
        assert loaded.energy_column_initial == original.energy_column_initial
        assert loaded.energy_column_intermediate == original.energy_column_intermediate
        assert loaded.energy_column_final == original.energy_column_final

    def test_roundtrip_with_pes_lp_corr(self, tmp_path):
        """Round-trip with pes_lp_corr set."""
        config = self._make_config(tmp_path)
        config.pes_lp_corr = Path("corr_pes.dat")

        yaml_path = tmp_path / "config_lp.yaml"
        save_config(config, yaml_path)
        loaded = load_config(yaml_path)

        assert loaded.pes_lp_corr == Path("corr_pes.dat")

    def test_load_without_optional_fields(self, tmp_path):
        """Loading YAML without optional fields should use defaults."""
        yaml_path = tmp_path / "minimal.yaml"
        yaml_path.write_text("""\
kh_1d:
  mu: 1.0
  grid:
    start: 0.5
    dx: 0.025
    npoints: 77
  broadening:
    gamma_fwhm: 0.3
  frequency:
    omega_start: 517
    omega_end: 528
    n_omega: 100
  pes_initial: "a.dat"
  pes_intermediate: "b.dat"
  pes_final_list:
    - "c.dat"
  dipole_final_list:
    - "d.dat"
""")
        config = load_config(yaml_path)
        assert config.dipole_mode == "FC"
        assert config.n_vib_states is None
        assert config.pes_lp_corr is None
        assert config.energy_column_initial == 1


class TestLPCorrection:
    """Tests for pes_lp_corr energy correction."""

    @pytest.fixture
    def harmonic_files_with_corr(self, tmp_path):
        """Create PES/dipole files plus a correction PES."""
        npoints = 77
        x_start = 0.5
        dx = 0.025
        x = np.array([x_start + i * dx for i in range(npoints)])

        k_au = 0.5
        x0_initial = 0.96
        x0_intermediate = 1.00
        x0_final = 0.98

        # Initial
        E_initial = 0.5 * k_au * ((x - x0_initial) / (CONST.bohr * 1e10)) ** 2
        pes_i_file = tmp_path / "pes_initial.dat"
        np.savetxt(pes_i_file, np.column_stack([x, E_initial]))

        # Intermediate
        E_offset = 530.0 / CONST.hartree2eV
        E_intermediate = E_offset + 0.5 * k_au * ((x - x0_intermediate) / (CONST.bohr * 1e10)) ** 2
        pes_n_file = tmp_path / "pes_intermediate.dat"
        np.savetxt(pes_n_file, np.column_stack([x, E_intermediate]))

        # Final state
        E_final_offset = 5.0 / CONST.hartree2eV
        E_final = E_final_offset + 0.5 * k_au * ((x - x0_final) / (CONST.bohr * 1e10)) ** 2
        pes_f_file = tmp_path / "pes_final.dat"
        np.savetxt(pes_f_file, np.column_stack([x, E_final]))

        # Correction PES: same shape but shifted by 0.1 Hartree
        energy_shift = 0.1  # Hartree
        E_corr = E_final + energy_shift
        pes_corr_file = tmp_path / "pes_lp_corr.dat"
        np.savetxt(pes_corr_file, np.column_stack([x, E_corr]))

        # Dipole
        dipole = np.zeros((npoints, 4))
        dipole[:, 0] = x
        dipole[:, 3] = 0.05
        dipole_file = tmp_path / "dipole.dat"
        np.savetxt(dipole_file, dipole)

        return {
            "pes_initial": pes_i_file,
            "pes_intermediate": pes_n_file,
            "pes_final": pes_f_file,
            "pes_corr": pes_corr_file,
            "dipole": dipole_file,
            "grid": (x_start, dx, npoints),
            "energy_shift_hartree": energy_shift,
        }

    def test_lp_correction_shifts_eigenvalues(self, harmonic_files_with_corr):
        """LP correction should shift final state eigenvalues."""
        files = harmonic_files_with_corr
        x_start, dx, npoints = files["grid"]

        # Without correction
        config_no_corr = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=520, omega_end=530, n_omega=11),
            pes_initial=files["pes_initial"],
            pes_intermediate=files["pes_intermediate"],
            pes_final_list=[files["pes_final"]],
            dipole_final_list=[files["dipole"]],
            n_vib_states=5,
        )

        calc_no_corr = KHCalculator(config_no_corr)
        calc_no_corr.load_surfaces()
        calc_no_corr.solve_vibrational_problem()
        E_no_corr = calc_no_corr._vib_f_list[0].eigenvalues_eV

        # With correction
        config_corr = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=520, omega_end=530, n_omega=11),
            pes_initial=files["pes_initial"],
            pes_intermediate=files["pes_intermediate"],
            pes_final_list=[files["pes_final"]],
            dipole_final_list=[files["dipole"]],
            pes_lp_corr=files["pes_corr"],
            n_vib_states=5,
        )

        calc_corr = KHCalculator(config_corr)
        calc_corr.load_surfaces()
        calc_corr.solve_vibrational_problem()
        E_corr = calc_corr._vib_f_list[0].eigenvalues_eV

        # Eigenvalues should be shifted by the energy correction
        shift_eV = files["energy_shift_hartree"] * CONST.hartree2eV
        np.testing.assert_allclose(E_corr - E_no_corr, shift_eV, rtol=1e-6)

    def test_lp_correction_shifts_spectrum(self, harmonic_files_with_corr):
        """LP correction should shift the spectrum peak position."""
        files = harmonic_files_with_corr
        x_start, dx, npoints = files["grid"]

        # Without correction
        config_no_corr = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=518, omega_end=530, n_omega=1001),
            pes_initial=files["pes_initial"],
            pes_intermediate=files["pes_intermediate"],
            pes_final_list=[files["pes_final"]],
            dipole_final_list=[files["dipole"]],
            n_vib_states=10,
        )
        result_no_corr = KHCalculator(config_no_corr).run()
        peak_no_corr = result_no_corr.omega[np.argmax(result_no_corr.sigma_tot)]

        # With correction
        config_corr = KHConfig(
            mu=1.0078825,
            grid=GridConfig(start=x_start, dx=dx, npoints=npoints),
            broadening=BroadeningConfig(gamma_fwhm=0.3),
            frequency=FrequencyGridConfig(omega_start=518, omega_end=530, n_omega=1001),
            pes_initial=files["pes_initial"],
            pes_intermediate=files["pes_intermediate"],
            pes_final_list=[files["pes_final"]],
            dipole_final_list=[files["dipole"]],
            pes_lp_corr=files["pes_corr"],
            n_vib_states=10,
        )
        result_corr = KHCalculator(config_corr).run()
        peak_corr = result_corr.omega[np.argmax(result_corr.sigma_tot)]

        # Peak should shift to lower energy (final state shifted up ->
        # emission energy = E_n - E_f decreases)
        shift_eV = files["energy_shift_hartree"] * CONST.hartree2eV
        assert abs((peak_no_corr - peak_corr) - shift_eV) < 0.2
