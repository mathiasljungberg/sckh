"""Tests for integrators module."""

import numpy as np
import pytest

from python_scripts.dynamics_1d.integrators import (
    velocity_verlet_step,
    run_trajectory,
    compute_kinetic_energy,
    compute_total_energy,
)
from python_scripts.dynamics_1d.constants import CONST


class TestVelocityVerletStep:
    """Tests for velocity_verlet_step function."""

    def test_free_particle_position_update(self):
        """Free particle should move linearly: x = x0 + v*t."""
        mass = 1.0 * CONST.u
        x0, v0, a0 = 0.0, 1.0e-5, 0.0  # Free particle, no force
        dt = 1e-15  # 1 fs

        def zero_force(x):
            return 0.0

        x_new, v_new, a_new = velocity_verlet_step(
            x0, v0, a0, zero_force, mass, dt
        )

        # Position should increase by v*dt
        np.testing.assert_allclose(x_new, x0 + v0 * dt, rtol=1e-12)
        # Velocity should be unchanged
        np.testing.assert_allclose(v_new, v0, rtol=1e-12)
        # Acceleration should be zero
        np.testing.assert_allclose(a_new, 0.0, atol=1e-30)

    def test_constant_force_motion(self):
        """Constant force should give parabolic motion."""
        mass = 1.0 * CONST.u
        F_const = 1e-10  # Constant force in N
        x0, v0 = 0.0, 0.0
        a0 = F_const / mass
        dt = 1e-15

        def const_force(x):
            return F_const

        x_new, v_new, a_new = velocity_verlet_step(
            x0, v0, a0, const_force, mass, dt
        )

        # x = 0.5 * a * t^2
        x_expected = 0.5 * a0 * dt**2
        np.testing.assert_allclose(x_new, x_expected, rtol=1e-12)

        # v = a * t
        v_expected = a0 * dt
        np.testing.assert_allclose(v_new, v_expected, rtol=1e-12)


class TestRunTrajectory:
    """Tests for run_trajectory function."""

    def test_harmonic_oscillator_energy_conservation(self, harmonic_params):
        """Total energy should be conserved for harmonic oscillator."""
        mass = harmonic_params["mass"]
        k = harmonic_params["k"]
        x0_eq = harmonic_params["x0"]

        # Initial displacement from equilibrium
        x0 = 0.1e-10  # 0.1 Angstrom displacement
        v0 = 0.0

        def harmonic_force(x):
            return -k * (x - x0_eq)

        def harmonic_potential(x):
            return 0.5 * k * (x - x0_eq) ** 2

        # Run for several periods with fine time step for good energy conservation
        dt = harmonic_params["period"] / 2000  # 2000 steps per period
        nsteps = 6000

        time, x, v, a = run_trajectory(x0, v0, harmonic_force, mass, dt, nsteps)

        # Compute total energy
        E_total = compute_total_energy(x, v, mass, harmonic_potential)

        # Initial energy
        E0 = E_total[0]

        # Energy should be conserved (Velocity Verlet is symplectic)
        # Allow for small numerical drift proportional to time step squared
        np.testing.assert_allclose(E_total, E0, rtol=1e-4)

    def test_harmonic_oscillator_period(self, harmonic_params):
        """Oscillation period should match T = 2π√(m/k)."""
        mass = harmonic_params["mass"]
        k = harmonic_params["k"]
        x0_eq = harmonic_params["x0"]
        T_expected = harmonic_params["period"]

        # Initial displacement
        x0 = 0.1e-10
        v0 = 0.0

        def harmonic_force(x):
            return -k * (x - x0_eq)

        # Run for ~2 periods with fine time resolution
        dt = T_expected / 500
        nsteps = 1000

        time, x, v, a = run_trajectory(x0, v0, harmonic_force, mass, dt, nsteps)

        # Find zero crossings (x passing through x0_eq)
        # The period is twice the time between first two zero crossings
        crossings = []
        for i in range(1, len(x)):
            if (x[i - 1] - x0_eq) * (x[i] - x0_eq) < 0:
                # Linear interpolation to find exact crossing time
                t_cross = time[i - 1] + (time[i] - time[i - 1]) * abs(
                    x[i - 1] - x0_eq
                ) / abs(x[i] - x[i - 1])
                crossings.append(t_cross)

        # Period is 2x the half-period (time between first two crossings)
        if len(crossings) >= 2:
            T_measured = 2 * (crossings[1] - crossings[0])
            np.testing.assert_allclose(T_measured, T_expected, rtol=1e-3)

    def test_free_particle_linear_motion(self):
        """Free particle should move at constant velocity."""
        mass = 1.0 * CONST.u
        x0 = 0.0
        v0 = 1e-5  # m/s
        dt = 1e-15
        nsteps = 100

        def zero_force(x):
            return 0.0

        time, x, v, a = run_trajectory(x0, v0, zero_force, mass, dt, nsteps)

        # Position should be linear: x = x0 + v0*t
        x_expected = x0 + v0 * time
        np.testing.assert_allclose(x, x_expected, rtol=1e-12)

        # Velocity should be constant
        np.testing.assert_allclose(v, v0, rtol=1e-12)

    def test_output_array_shapes(self):
        """Output arrays should have correct shapes."""
        mass = 1.0 * CONST.u
        nsteps = 50

        def zero_force(x):
            return 0.0

        time, x, v, a = run_trajectory(0.0, 0.0, zero_force, mass, 1e-15, nsteps)

        assert time.shape == (nsteps,)
        assert x.shape == (nsteps,)
        assert v.shape == (nsteps,)
        assert a.shape == (nsteps,)

    def test_initial_conditions_preserved(self):
        """First element should be the initial conditions."""
        mass = 1.0 * CONST.u
        x0, v0 = 1.5e-10, 2.3e-5

        def zero_force(x):
            return 0.0

        time, x, v, a = run_trajectory(x0, v0, zero_force, mass, 1e-15, 10)

        assert time[0] == 0.0
        assert x[0] == x0
        assert v[0] == v0


class TestEnergyFunctions:
    """Tests for energy computation functions."""

    def test_kinetic_energy(self):
        """Kinetic energy should be 0.5*m*v^2."""
        mass = 2.0 * CONST.u
        v = np.array([0.0, 1e-5, 2e-5])

        E_kin = compute_kinetic_energy(v, mass)
        E_expected = 0.5 * mass * v**2

        np.testing.assert_allclose(E_kin, E_expected, rtol=1e-12)

    def test_total_energy(self):
        """Total energy should be kinetic + potential."""
        mass = 1.0 * CONST.u
        x = np.array([0.0, 1e-10, 2e-10])
        v = np.array([1e-5, 2e-5, 3e-5])

        def potential(x):
            return 0.5 * 100 * x**2  # Simple harmonic

        E_total = compute_total_energy(x, v, mass, potential)
        E_expected = 0.5 * mass * v**2 + potential(x)

        np.testing.assert_allclose(E_total, E_expected, rtol=1e-12)
