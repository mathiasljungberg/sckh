"""Initial condition sampling for 2D dynamics.

Samples initial conditions from product wavefunction approximation:
    |psi(x1, x2)|^2 ~ |psi1(x1)|^2 * |psi2(x2)|^2

Each coordinate is sampled independently using the 1D sampling methods.
"""

from typing import Optional, Tuple

import numpy as np

from python_scripts.dynamics_1d.sampling import (
    sample_even,
    sample_random,
    sample_momenta,
)


def create_initial_conditions_2d(
    x1_grid: np.ndarray,
    psi1: np.ndarray,
    x2_grid: np.ndarray,
    psi2: np.ndarray,
    n_x1: int,
    n_x2: int,
    n_p1: int,
    n_p2: int,
    mode: int = 1,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample initial conditions for 2D dynamics from product wavefunction.

    For mode=1 (even sampling): Creates n_x1 * n_x2 * n_p1 * n_p2 total
    trajectories by taking all combinations of position and momentum samples.

    For mode=2 (random sampling): Creates n_x1 pairs where all coordinates
    are sampled randomly.

    Args:
        x1_grid: Position grid for x1 (m)
        psi1: Wavefunction for x1
        x2_grid: Position grid for x2 (m)
        psi2: Wavefunction for x2
        n_x1: Number of position samples for x1
        n_x2: Number of position samples for x2
        n_p1: Number of momentum samples for p1
        n_p2: Number of momentum samples for p2
        mode: Sampling mode (1=even, 2=random)
        rng: Random number generator (for mode=2)

    Returns:
        x1_init: Array of initial x1 positions (m)
        x2_init: Array of initial x2 positions (m)
        p1_init: Array of initial p1 momenta (kg*m/s)
        p2_init: Array of initial p2 momenta (kg*m/s)
    """
    psi1_squared = np.abs(psi1) ** 2
    psi2_squared = np.abs(psi2) ** 2

    if mode == 1:
        # Even sampling: all combinations
        x1_samples = sample_even(x1_grid, psi1_squared, n_x1)
        x2_samples = sample_even(x2_grid, psi2_squared, n_x2)
        p1_samples = sample_momenta(x1_grid, psi1, n_p1, mode=1)
        p2_samples = sample_momenta(x2_grid, psi2, n_p2, mode=1)

        # Create all combinations using meshgrid
        # Order: iterate over p2, then p1, then x2, then x1 (innermost to outermost)
        x1_mesh, x2_mesh, p1_mesh, p2_mesh = np.meshgrid(
            x1_samples, x2_samples, p1_samples, p2_samples, indexing="ij"
        )

        return (
            x1_mesh.flatten(),
            x2_mesh.flatten(),
            p1_mesh.flatten(),
            p2_mesh.flatten(),
        )

    else:
        # Random sampling: paired samples
        if rng is None:
            rng = np.random.default_rng()

        # For random mode, use n_x1 as the total number of samples
        n_total = n_x1

        x1_samples = sample_random(x1_grid, psi1_squared, n_total, rng)
        x2_samples = sample_random(x2_grid, psi2_squared, n_total, rng)
        p1_samples = sample_momenta(x1_grid, psi1, n_total, mode=2, rng=rng)
        p2_samples = sample_momenta(x2_grid, psi2, n_total, mode=2, rng=rng)

        return x1_samples, x2_samples, p1_samples, p2_samples


def get_n_trajectories(
    n_x1: int,
    n_x2: int,
    n_p1: int,
    n_p2: int,
    mode: int = 1,
) -> int:
    """Calculate total number of trajectories for given sampling parameters.

    Args:
        n_x1: Number of position samples for x1
        n_x2: Number of position samples for x2
        n_p1: Number of momentum samples for p1
        n_p2: Number of momentum samples for p2
        mode: Sampling mode (1=even, 2=random)

    Returns:
        Total number of trajectories
    """
    if mode == 1:
        # Even sampling: all combinations
        return n_x1 * n_x2 * n_p1 * n_p2
    else:
        # Random sampling: paired samples
        return n_x1


def sample_positions_2d(
    x1_grid: np.ndarray,
    psi1: np.ndarray,
    x2_grid: np.ndarray,
    psi2: np.ndarray,
    n_x1: int,
    n_x2: int,
    mode: int = 1,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Sample only positions (not momenta) from product wavefunction.

    Useful for diagnostics and visualization.

    Args:
        x1_grid: Position grid for x1 (m)
        psi1: Wavefunction for x1
        x2_grid: Position grid for x2 (m)
        psi2: Wavefunction for x2
        n_x1: Number of position samples for x1
        n_x2: Number of position samples for x2
        mode: Sampling mode (1=even, 2=random)
        rng: Random number generator (for mode=2)

    Returns:
        x1_samples, x2_samples: Arrays of sampled positions
    """
    psi1_squared = np.abs(psi1) ** 2
    psi2_squared = np.abs(psi2) ** 2

    if mode == 1:
        x1_samples = sample_even(x1_grid, psi1_squared, n_x1)
        x2_samples = sample_even(x2_grid, psi2_squared, n_x2)

        # All combinations
        x1_mesh, x2_mesh = np.meshgrid(x1_samples, x2_samples, indexing="ij")
        return x1_mesh.flatten(), x2_mesh.flatten()
    else:
        if rng is None:
            rng = np.random.default_rng()

        n_total = n_x1
        x1_samples = sample_random(x1_grid, psi1_squared, n_total, rng)
        x2_samples = sample_random(x2_grid, psi2_squared, n_total, rng)
        return x1_samples, x2_samples
