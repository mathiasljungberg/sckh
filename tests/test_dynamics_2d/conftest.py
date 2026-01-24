"""Pytest fixtures for dynamics_2d tests."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.constants import CONST


@pytest.fixture
def harmonic_params_2d():
    """Parameters for a 2D separable harmonic oscillator test case.

    Uses different masses and force constants for each coordinate.
    """
    mass1 = 1.0 * CONST.u  # 1 amu in kg
    mass2 = 2.0 * CONST.u  # 2 amu in kg
    k1 = 500.0  # Force constant for x1 (N/m)
    k2 = 300.0  # Force constant for x2 (N/m)
    x1_0 = 0.0  # Equilibrium position for x1 (m)
    x2_0 = 0.0  # Equilibrium position for x2 (m)
    omega1 = np.sqrt(k1 / mass1)
    omega2 = np.sqrt(k2 / mass2)

    return {
        "mass1": mass1,
        "mass2": mass2,
        "k1": k1,
        "k2": k2,
        "x1_0": x1_0,
        "x2_0": x2_0,
        "omega1": omega1,
        "omega2": omega2,
        "period1": 2 * np.pi / omega1,
        "period2": 2 * np.pi / omega2,
    }


@pytest.fixture
def position_grid_2d():
    """Standard position grids for 2D testing (in SI units).

    Returns tuple of (x1_grid, x2_grid).
    """
    # Grid from -0.5 to 0.5 Angstrom, converted to meters
    x1_ang = np.linspace(-0.5, 0.5, 101)
    x2_ang = np.linspace(-0.5, 0.5, 101)
    x1_m = x1_ang * 1e-10
    x2_m = x2_ang * 1e-10
    return x1_m, x2_m


@pytest.fixture
def fine_position_grid_2d():
    """Fine position grids for high-accuracy 2D tests."""
    x1_ang = np.linspace(-1.0, 1.0, 201)
    x2_ang = np.linspace(-1.0, 1.0, 201)
    x1_m = x1_ang * 1e-10
    x2_m = x2_ang * 1e-10
    return x1_m, x2_m


@pytest.fixture
def harmonic_pes_2d(position_grid_2d, harmonic_params_2d):
    """2D separable harmonic potential energy surface."""
    from python_scripts.dynamics_2d.pes import create_harmonic_pes_2d

    x1, x2 = position_grid_2d
    return create_harmonic_pes_2d(
        x1,
        x2,
        x1_0=harmonic_params_2d["x1_0"],
        x2_0=harmonic_params_2d["x2_0"],
        k1=harmonic_params_2d["k1"],
        k2=harmonic_params_2d["k2"],
        E0=0.0,
    )


@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return np.random.default_rng(42)
