"""Pytest fixtures for dynamics_1d tests."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.constants import CONST


@pytest.fixture
def harmonic_params():
    """Parameters for a harmonic oscillator test case.

    Uses hydrogen mass and a reasonable force constant.
    """
    mass = 1.0 * CONST.u  # 1 amu in kg (roughly hydrogen)
    k = 500.0  # Force constant in N/m (reasonable for molecular vibration)
    x0 = 0.0  # Equilibrium position in meters
    omega = np.sqrt(k / mass)  # Angular frequency

    return {
        "mass": mass,
        "k": k,
        "x0": x0,
        "omega": omega,
        "period": 2 * np.pi / omega,
    }


@pytest.fixture
def position_grid():
    """Standard position grid for testing (in SI units)."""
    # Grid from -0.5 to 0.5 Angstrom, converted to meters
    x_ang = np.linspace(-0.5, 0.5, 101)
    x_m = x_ang * 1e-10
    return x_m


@pytest.fixture
def fine_position_grid():
    """Fine position grid for high-accuracy tests."""
    x_ang = np.linspace(-1.0, 1.0, 401)
    x_m = x_ang * 1e-10
    return x_m


@pytest.fixture
def harmonic_pes(position_grid, harmonic_params):
    """Harmonic potential energy surface on the standard grid."""
    from python_scripts.dynamics_1d.pes import create_harmonic_pes

    return create_harmonic_pes(
        position_grid,
        x0=harmonic_params["x0"],
        k=harmonic_params["k"],
        E0=0.0,
    )


@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return np.random.default_rng(42)
