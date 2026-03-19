"""Numerical integrators for classical dynamics."""

from typing import Callable, Tuple

import numpy as np


def velocity_verlet_step(
    x: float,
    v: float,
    a: float,
    force_func: Callable[[float], float],
    mass: float,
    dt: float,
) -> Tuple[float, float, float]:
    """Single step of Velocity Verlet integration.

    Algorithm:
        1. x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
        2. a(t+dt) = F(x(t+dt)) / m
        3. v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt

    Args:
        x: Current position (m)
        v: Current velocity (m/s)
        a: Current acceleration (m/s^2)
        force_func: Function that returns force given position (N)
        mass: Particle mass (kg)
        dt: Time step (s)

    Returns:
        x_new, v_new, a_new: Updated position, velocity, acceleration
    """
    # Update position
    x_new = x + v * dt + 0.5 * a * dt**2

    # Compute new acceleration
    a_new = force_func(x_new) / mass

    # Update velocity
    v_new = v + 0.5 * (a + a_new) * dt

    return x_new, v_new, a_new


def run_trajectory(
    x0: float,
    v0: float,
    force_func: Callable[[float], float],
    mass: float,
    dt: float,
    nsteps: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run a complete trajectory using Velocity Verlet.

    Args:
        x0: Initial position (m)
        v0: Initial velocity (m/s)
        force_func: Function that returns force given position (N)
        mass: Particle mass (kg)
        dt: Time step (s)
        nsteps: Number of time steps

    Returns:
        time: Time array (s)
        x: Position array (m)
        v: Velocity array (m/s)
        a: Acceleration array (m/s^2)
    """
    # Allocate arrays
    time = np.zeros(nsteps)
    x = np.zeros(nsteps)
    v = np.zeros(nsteps)
    a = np.zeros(nsteps)

    # Initial conditions
    x[0] = x0
    v[0] = v0
    a[0] = force_func(x0) / mass

    # Integration loop
    for i in range(1, nsteps):
        time[i] = i * dt
        x[i], v[i], a[i] = velocity_verlet_step(
            x[i - 1], v[i - 1], a[i - 1], force_func, mass, dt
        )

    return time, x, v, a


def compute_kinetic_energy(v: np.ndarray, mass: float) -> np.ndarray:
    """Compute kinetic energy along trajectory.

    Args:
        v: Velocity array (m/s)
        mass: Particle mass (kg)

    Returns:
        Kinetic energy array (J)
    """
    return 0.5 * mass * v**2


def compute_total_energy(
    x: np.ndarray,
    v: np.ndarray,
    mass: float,
    potential_func: Callable[[np.ndarray], np.ndarray],
) -> np.ndarray:
    """Compute total energy along trajectory.

    Args:
        x: Position array (m)
        v: Velocity array (m/s)
        mass: Particle mass (kg)
        potential_func: Function that returns potential energy given position (J)

    Returns:
        Total energy array (J)
    """
    kinetic = compute_kinetic_energy(v, mass)
    potential = potential_func(x)
    return kinetic + potential
