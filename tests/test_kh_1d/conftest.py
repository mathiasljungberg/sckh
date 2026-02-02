"""Shared fixtures for kh_1d tests."""

import numpy as np
import pytest
from pathlib import Path

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_1d.pes import PES1D, create_harmonic_pes
from python_scripts.dynamics_1d.dipole import Dipole1D, create_constant_dipole
from python_scripts.dynamics_1d.vibrational import solve_vibrational


@pytest.fixture
def harmonic_grid():
    """Create a standard grid for harmonic oscillator tests."""
    # Grid in meters (like Angstrom but in SI)
    npoints = 77  # Odd number for Fourier grid
    x_start = 0.5e-10  # 0.5 Angstrom in meters
    dx = 0.025e-10  # 0.025 Angstrom in meters
    x = np.array([x_start + i * dx for i in range(npoints)])
    return x, dx


@pytest.fixture
def hydrogen_mass():
    """Hydrogen mass in kg."""
    return 1.0078825 * CONST.u


@pytest.fixture
def harmonic_pes_params():
    """Parameters for a harmonic oscillator PES."""
    return {
        "x0": 0.96e-10,  # Equilibrium at ~0.96 Angstrom
        "k": 500.0,  # Force constant in N/m (reasonable for OH stretch)
        "E0": 0.0,  # Energy offset
    }


@pytest.fixture
def harmonic_pes_initial(harmonic_grid, harmonic_pes_params):
    """Create harmonic PES for initial state."""
    x, dx = harmonic_grid
    return create_harmonic_pes(x, **harmonic_pes_params)


@pytest.fixture
def harmonic_pes_intermediate(harmonic_grid, harmonic_pes_params):
    """Create harmonic PES for intermediate state (slightly shifted)."""
    x, dx = harmonic_grid
    params = harmonic_pes_params.copy()
    params["x0"] = 1.0e-10  # Shifted equilibrium
    params["E0"] = 530.0 * CONST.eV  # Core-excited state energy
    return create_harmonic_pes(x, **params)


@pytest.fixture
def harmonic_pes_final(harmonic_grid, harmonic_pes_params):
    """Create harmonic PES for final state."""
    x, dx = harmonic_grid
    params = harmonic_pes_params.copy()
    params["x0"] = 0.98e-10  # Slightly different equilibrium
    params["E0"] = 5.0 * CONST.eV  # Final state energy
    return create_harmonic_pes(x, **params)


@pytest.fixture
def constant_dipole(harmonic_grid):
    """Create constant dipole moment surface."""
    x, dx = harmonic_grid
    d_value = np.array([0.0, 0.0, 1.0])  # Unit dipole in z-direction
    return create_constant_dipole(x, d_value)


@pytest.fixture
def harmonic_vibrational_states(harmonic_grid, harmonic_pes_initial, hydrogen_mass):
    """Solve vibrational problem for harmonic PES."""
    x, dx = harmonic_grid
    V = harmonic_pes_initial.energy(x)
    eigenvalues, eigenvectors = solve_vibrational(x, V, hydrogen_mass)
    return eigenvalues, eigenvectors, x


@pytest.fixture
def testsuite_path():
    """Path to the KH testsuite."""
    return Path(__file__).parent.parent / "testsuite" / "KH"


@pytest.fixture
def testsuite_input_path(testsuite_path):
    """Path to KH testsuite input files."""
    return testsuite_path / "input"


@pytest.fixture
def testsuite_ref_path(testsuite_path):
    """Path to KH testsuite reference files."""
    return testsuite_path / "ref"
