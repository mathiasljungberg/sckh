"""Tests for 2D numerical integrators."""

import numpy as np
import pytest

from python_scripts.dynamics_2d.integrators import (
    velocity_verlet_step_2d,
    run_trajectory_2d,
    compute_kinetic_energy_2d,
    compute_total_energy_2d,
)


class TestVelocityVerletStep2D:
    """Tests for single Velocity Verlet step."""

    def test_free_particle(self):
        """Test free particle (zero force) propagation."""
        # Initial conditions
        x1_0, x2_0 = 0.0, 0.0
        v1_0, v2_0 = 100.0, 200.0  # m/s
        a1_0, a2_0 = 0.0, 0.0
        m1, m2 = 1e-26, 1e-26  # kg
        dt = 1e-15  # 1 fs

        def zero_force(x1, x2):
            return 0.0, 0.0

        x1, x2, v1, v2, a1, a2 = velocity_verlet_step_2d(
            x1_0, x2_0, v1_0, v2_0, a1_0, a2_0, zero_force, m1, m2, dt
        )

        # For free particle: x = x0 + v*t
        assert x1 == pytest.approx(v1_0 * dt, rel=1e-10)
        assert x2 == pytest.approx(v2_0 * dt, rel=1e-10)
        # Velocity should be unchanged
        assert v1 == pytest.approx(v1_0, rel=1e-10)
        assert v2 == pytest.approx(v2_0, rel=1e-10)

    def test_constant_force(self):
        """Test propagation under constant force."""
        x1_0, x2_0 = 0.0, 0.0
        v1_0, v2_0 = 0.0, 0.0
        m1, m2 = 1e-26, 2e-26
        F1, F2 = 1e-12, 2e-12  # Constant forces
        a1_0 = F1 / m1
        a2_0 = F2 / m2
        dt = 1e-15

        def const_force(x1, x2):
            return F1, F2

        x1, x2, v1, v2, a1, a2 = velocity_verlet_step_2d(
            x1_0, x2_0, v1_0, v2_0, a1_0, a2_0, const_force, m1, m2, dt
        )

        # For constant acceleration: x = 0.5*a*t^2, v = a*t
        assert x1 == pytest.approx(0.5 * a1_0 * dt**2, rel=1e-10)
        assert v1 == pytest.approx(a1_0 * dt, rel=1e-10)
        assert x2 == pytest.approx(0.5 * a2_0 * dt**2, rel=1e-10)
        assert v2 == pytest.approx(a2_0 * dt, rel=1e-10)


class TestRunTrajectory2D:
    """Tests for complete 2D trajectory integration."""

    def test_harmonic_oscillator_energy_conservation(
        self, harmonic_pes_2d, harmonic_params_2d
    ):
        """Test energy conservation for 2D harmonic oscillator."""
        pes = harmonic_pes_2d
        m1 = harmonic_params_2d["mass1"]
        m2 = harmonic_params_2d["mass2"]

        # Initial conditions: displaced from equilibrium
        x1_0 = 0.1e-10  # 0.1 Angstrom
        x2_0 = 0.05e-10
        v1_0 = 0.0
        v2_0 = 0.0

        # Time step and duration
        dt = 0.1e-15  # 0.1 fs
        nsteps = 1000

        time, x1, x2, v1, v2, a1, a2 = run_trajectory_2d(
            x1_0, x2_0, v1_0, v2_0, pes.forces, m1, m2, dt, nsteps
        )

        # Compute total energy along trajectory
        E_total = compute_total_energy_2d(x1, x2, v1, v2, m1, m2, pes.energy)

        # Energy should be conserved (constant within numerical precision)
        E_initial = E_total[0]
        E_drift = np.abs(E_total - E_initial) / E_initial

        # Relative energy drift should be small (< 1e-3 for this time step)
        assert np.max(E_drift) < 1e-3, f"Max energy drift: {np.max(E_drift)}"

    def test_harmonic_oscillator_period(self, harmonic_pes_2d, harmonic_params_2d):
        """Test that oscillation period matches theory."""
        pes = harmonic_pes_2d
        m1 = harmonic_params_2d["mass1"]
        m2 = harmonic_params_2d["mass2"]
        period1 = harmonic_params_2d["period1"]
        period2 = harmonic_params_2d["period2"]

        # Initial conditions: displaced in x1 only
        x1_0 = 0.1e-10
        x2_0 = 0.0
        v1_0 = 0.0
        v2_0 = 0.0

        # Run for ~5 periods of the slower oscillator
        dt = 0.01e-15  # Fine time step
        nsteps = int(5 * max(period1, period2) / dt)

        time, x1, x2, v1, v2, a1, a2 = run_trajectory_2d(
            x1_0, x2_0, v1_0, v2_0, pes.forces, m1, m2, dt, nsteps
        )

        # Find zero crossings to determine period
        # x1 should cross zero at t = T1/4, 3*T1/4, 5*T1/4, etc.
        crossings = []
        for i in range(1, len(x1)):
            if x1[i - 1] > 0 and x1[i] <= 0:
                # Linear interpolation for crossing time
                t_cross = time[i - 1] - x1[i - 1] * (time[i] - time[i - 1]) / (
                    x1[i] - x1[i - 1]
                )
                crossings.append(t_cross)

        # Period is twice the time between consecutive crossings
        if len(crossings) >= 2:
            measured_period = 2 * (crossings[1] - crossings[0])
            assert measured_period == pytest.approx(period1, rel=0.01)

    def test_separable_dynamics(self, harmonic_pes_2d, harmonic_params_2d):
        """Test that separable potential gives independent dynamics."""
        pes = harmonic_pes_2d
        m1 = harmonic_params_2d["mass1"]
        m2 = harmonic_params_2d["mass2"]
        omega1 = harmonic_params_2d["omega1"]
        omega2 = harmonic_params_2d["omega2"]

        # Initial conditions
        x1_0 = 0.1e-10
        x2_0 = 0.05e-10
        v1_0 = 0.0
        v2_0 = 0.0

        dt = 0.05e-15
        nsteps = 500

        time, x1, x2, v1, v2, a1, a2 = run_trajectory_2d(
            x1_0, x2_0, v1_0, v2_0, pes.forces, m1, m2, dt, nsteps
        )

        # For harmonic oscillator: x(t) = A*cos(omega*t)
        # Compare with analytical solution
        x1_analytical = x1_0 * np.cos(omega1 * time)
        x2_analytical = x2_0 * np.cos(omega2 * time)

        # Should match within numerical precision
        # Use absolute tolerance since relative error can be large near zero crossings
        np.testing.assert_allclose(x1, x1_analytical, atol=1e-14, rtol=1e-2)
        np.testing.assert_allclose(x2, x2_analytical, atol=1e-14, rtol=1e-2)


class TestEnergyFunctions:
    """Tests for energy computation functions."""

    def test_kinetic_energy(self, harmonic_params_2d):
        """Test kinetic energy computation."""
        m1 = harmonic_params_2d["mass1"]
        m2 = harmonic_params_2d["mass2"]

        v1 = np.array([100.0, 200.0, 0.0])
        v2 = np.array([50.0, 0.0, 100.0])

        KE = compute_kinetic_energy_2d(v1, v2, m1, m2)

        # Expected: KE = 0.5*m1*v1^2 + 0.5*m2*v2^2
        KE_expected = 0.5 * m1 * v1**2 + 0.5 * m2 * v2**2
        np.testing.assert_allclose(KE, KE_expected)

    def test_total_energy(self, harmonic_pes_2d, harmonic_params_2d):
        """Test total energy computation."""
        pes = harmonic_pes_2d
        m1 = harmonic_params_2d["mass1"]
        m2 = harmonic_params_2d["mass2"]
        k1 = harmonic_params_2d["k1"]
        k2 = harmonic_params_2d["k2"]

        x1 = np.array([0.1e-10, 0.2e-10])
        x2 = np.array([0.05e-10, 0.1e-10])
        v1 = np.array([100.0, 0.0])
        v2 = np.array([50.0, 100.0])

        E_total = compute_total_energy_2d(x1, x2, v1, v2, m1, m2, pes.energy)

        # Expected
        KE = 0.5 * m1 * v1**2 + 0.5 * m2 * v2**2
        PE = 0.5 * k1 * x1**2 + 0.5 * k2 * x2**2
        E_expected = KE + PE

        np.testing.assert_allclose(E_total, E_expected, rtol=1e-6)
