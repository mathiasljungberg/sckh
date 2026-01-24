"""Numerical integrators for 2D classical dynamics."""

from typing import Callable, Tuple

import numpy as np


def velocity_verlet_step_2d(
    x1: float,
    x2: float,
    v1: float,
    v2: float,
    a1: float,
    a2: float,
    force_func: Callable[[float, float], Tuple[float, float]],
    m1: float,
    m2: float,
    dt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Single step of Velocity Verlet integration for 2D system.

    Algorithm:
        1. x1(t+dt) = x1(t) + v1(t)*dt + 0.5*a1(t)*dt^2
           x2(t+dt) = x2(t) + v2(t)*dt + 0.5*a2(t)*dt^2
        2. F1(t+dt), F2(t+dt) = force_func(x1(t+dt), x2(t+dt))
           a1(t+dt) = F1(t+dt) / m1
           a2(t+dt) = F2(t+dt) / m2
        3. v1(t+dt) = v1(t) + 0.5*(a1(t) + a1(t+dt))*dt
           v2(t+dt) = v2(t) + 0.5*(a2(t) + a2(t+dt))*dt

    Args:
        x1: Current position in x1 (m)
        x2: Current position in x2 (m)
        v1: Current velocity in x1 (m/s)
        v2: Current velocity in x2 (m/s)
        a1: Current acceleration in x1 (m/s^2)
        a2: Current acceleration in x2 (m/s^2)
        force_func: Function returning (F1, F2) forces given (x1, x2)
        m1: Mass for x1 coordinate (kg)
        m2: Mass for x2 coordinate (kg)
        dt: Time step (s)

    Returns:
        x1_new, x2_new, v1_new, v2_new, a1_new, a2_new
    """
    # Update positions
    x1_new = x1 + v1 * dt + 0.5 * a1 * dt**2
    x2_new = x2 + v2 * dt + 0.5 * a2 * dt**2

    # Compute new forces and accelerations
    F1_new, F2_new = force_func(x1_new, x2_new)
    a1_new = F1_new / m1
    a2_new = F2_new / m2

    # Update velocities
    v1_new = v1 + 0.5 * (a1 + a1_new) * dt
    v2_new = v2 + 0.5 * (a2 + a2_new) * dt

    return x1_new, x2_new, v1_new, v2_new, a1_new, a2_new


def run_trajectory_2d(
    x1_0: float,
    x2_0: float,
    v1_0: float,
    v2_0: float,
    force_func: Callable[[float, float], Tuple[float, float]],
    m1: float,
    m2: float,
    dt: float,
    nsteps: int,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Run a complete 2D trajectory using Velocity Verlet.

    Args:
        x1_0: Initial position in x1 (m)
        x2_0: Initial position in x2 (m)
        v1_0: Initial velocity in x1 (m/s)
        v2_0: Initial velocity in x2 (m/s)
        force_func: Function returning (F1, F2) forces given (x1, x2)
        m1: Mass for x1 coordinate (kg)
        m2: Mass for x2 coordinate (kg)
        dt: Time step (s)
        nsteps: Number of time steps

    Returns:
        time: Time array (s)
        x1: Position array for x1 (m)
        x2: Position array for x2 (m)
        v1: Velocity array for x1 (m/s)
        v2: Velocity array for x2 (m/s)
        a1: Acceleration array for x1 (m/s^2)
        a2: Acceleration array for x2 (m/s^2)
    """
    # Allocate arrays
    time = np.zeros(nsteps)
    x1 = np.zeros(nsteps)
    x2 = np.zeros(nsteps)
    v1 = np.zeros(nsteps)
    v2 = np.zeros(nsteps)
    a1 = np.zeros(nsteps)
    a2 = np.zeros(nsteps)

    # Initial conditions
    x1[0] = x1_0
    x2[0] = x2_0
    v1[0] = v1_0
    v2[0] = v2_0

    # Initial accelerations
    F1_0, F2_0 = force_func(x1_0, x2_0)
    a1[0] = F1_0 / m1
    a2[0] = F2_0 / m2

    # Integration loop
    for i in range(1, nsteps):
        time[i] = i * dt
        x1[i], x2[i], v1[i], v2[i], a1[i], a2[i] = velocity_verlet_step_2d(
            x1[i - 1],
            x2[i - 1],
            v1[i - 1],
            v2[i - 1],
            a1[i - 1],
            a2[i - 1],
            force_func,
            m1,
            m2,
            dt,
        )

    return time, x1, x2, v1, v2, a1, a2


def compute_kinetic_energy_2d(
    v1: np.ndarray,
    v2: np.ndarray,
    m1: float,
    m2: float,
) -> np.ndarray:
    """Compute kinetic energy along 2D trajectory.

    T = 0.5 * m1 * v1^2 + 0.5 * m2 * v2^2

    Args:
        v1: Velocity array for x1 (m/s)
        v2: Velocity array for x2 (m/s)
        m1: Mass for x1 (kg)
        m2: Mass for x2 (kg)

    Returns:
        Kinetic energy array (J)
    """
    return 0.5 * m1 * v1**2 + 0.5 * m2 * v2**2


def compute_total_energy_2d(
    x1: np.ndarray,
    x2: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    m1: float,
    m2: float,
    potential_func: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> np.ndarray:
    """Compute total energy along 2D trajectory.

    E = T + V = 0.5*m1*v1^2 + 0.5*m2*v2^2 + V(x1, x2)

    Args:
        x1: Position array for x1 (m)
        x2: Position array for x2 (m)
        v1: Velocity array for x1 (m/s)
        v2: Velocity array for x2 (m/s)
        m1: Mass for x1 (kg)
        m2: Mass for x2 (kg)
        potential_func: Function returning potential energy V(x1, x2)

    Returns:
        Total energy array (J)
    """
    kinetic = compute_kinetic_energy_2d(v1, v2, m1, m2)

    # Compute potential energy point by point
    potential = np.array([potential_func(x1[i], x2[i]) for i in range(len(x1))])

    return kinetic + potential
