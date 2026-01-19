"""Initial condition sampling from ground state distribution."""

from typing import Tuple, Optional

import numpy as np
from scipy.interpolate import CubicSpline

from .constants import CONST


def sample_even(
    x: np.ndarray,
    psi_squared: np.ndarray,
    n_samples: int,
) -> np.ndarray:
    """Sample positions evenly from |ψ|² distribution using inverse CDF.

    Divides the cumulative distribution into n_samples equal parts
    and places sample points at the center of each interval.

    Args:
        x: Position grid (m)
        psi_squared: |ψ|² probability density
        n_samples: Number of samples to generate

    Returns:
        x_samples: Sampled positions (m)
    """
    dx = x[1] - x[0]

    # Compute CDF
    cdf = np.cumsum(psi_squared) * dx
    cdf_norm = cdf / cdf[-1]  # Normalize to [0, 1]

    # Target CDF values at center of each interval
    target_cdf = (np.arange(n_samples) + 0.5) / n_samples

    # Interpolate to find x values (inverse CDF)
    x_samples = np.interp(target_cdf, cdf_norm, x)

    return x_samples


def sample_random(
    x: np.ndarray,
    psi_squared: np.ndarray,
    n_samples: int,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Sample positions randomly from |ψ|² using inverse CDF transform.

    Args:
        x: Position grid (m)
        psi_squared: |ψ|² probability density
        n_samples: Number of samples to generate
        rng: Random number generator (for reproducibility)

    Returns:
        x_samples: Sampled positions (m)
    """
    if rng is None:
        rng = np.random.default_rng()

    dx = x[1] - x[0]

    # Compute CDF
    cdf = np.cumsum(psi_squared) * dx
    cdf_norm = cdf / cdf[-1]

    # Generate uniform random numbers and invert CDF
    u = rng.uniform(0, 1, n_samples)
    x_samples = np.interp(u, cdf_norm, x)

    return x_samples


def compute_momentum_distribution(
    x: np.ndarray,
    psi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute momentum-space wavefunction via FFT.

    Args:
        x: Position grid (m)
        psi: Wavefunction on grid

    Returns:
        p: Momentum grid (kg*m/s)
        psi_p_squared: |ψ(p)|² in momentum space
    """
    dx = x[1] - x[0]
    n = len(x)

    # FFT to momentum space
    psi_p = np.fft.fft(psi) * dx  # Include dx for proper normalization

    # Momentum grid
    # FFT frequencies are in cycles per sample, convert to momentum
    freq = np.fft.fftfreq(n, dx)  # cycles per meter
    p = freq * 2 * CONST.pi * CONST.hbar  # momentum in kg*m/s

    # Reorder for proper momentum axis (negative to positive)
    p = np.fft.fftshift(p)
    psi_p = np.fft.fftshift(psi_p)

    psi_p_squared = np.abs(psi_p) ** 2

    return p, psi_p_squared


def sample_momenta(
    x: np.ndarray,
    psi: np.ndarray,
    n_samples: int,
    mode: int = 1,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Sample momenta from momentum-space wavefunction via FFT.

    Args:
        x: Position grid (m)
        psi: Wavefunction (not |ψ|²)
        n_samples: Number of momentum samples
        mode: 1=even, 2=random
        rng: Random number generator (for mode=2)

    Returns:
        p_samples: Sampled momenta (kg*m/s)
    """
    p, psi_p_squared = compute_momentum_distribution(x, psi)

    # Sample using chosen method
    if mode == 1:
        p_samples = sample_even(p, psi_p_squared, n_samples)
    else:
        p_samples = sample_random(p, psi_p_squared, n_samples, rng)

    return p_samples


def create_initial_conditions(
    x: np.ndarray,
    psi: np.ndarray,
    n_x: int,
    n_p: int,
    mode: int = 1,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Create all (x, p) initial condition pairs.

    For mode=1 (even sampling): creates n_x * n_p total trajectories
    by taking all combinations of position and momentum samples.

    For mode=2 (random sampling): creates n_x pairs where both position
    and momentum are sampled randomly.

    Args:
        x: Position grid (m)
        psi: Ground state wavefunction
        n_x: Number of position samples
        n_p: Number of momentum samples
        mode: Sampling mode (1=even, 2=random)
        rng: Random number generator (for mode=2)

    Returns:
        x_samples: Array of initial positions
        p_samples: Array of initial momenta
    """
    psi_squared = np.abs(psi) ** 2

    if mode == 1:
        # Even sampling: all combinations of x and p
        x_samp = sample_even(x, psi_squared, n_x)
        p_samp = sample_momenta(x, psi, n_p, mode=1)

        # Create all combinations
        x_mesh, p_mesh = np.meshgrid(x_samp, p_samp)
        return x_mesh.flatten(), p_mesh.flatten()

    else:
        # Random sampling: paired (x, p) samples
        if rng is None:
            rng = np.random.default_rng()

        x_samp = sample_random(x, psi_squared, n_x, rng)
        p_samp = sample_momenta(x, psi, n_x, mode=2, rng=rng)

        return x_samp, p_samp


def get_sample_weights(
    x_samples: np.ndarray,
    p_samples: np.ndarray,
    x_grid: np.ndarray,
    psi: np.ndarray,
    mode: int = 1,
) -> np.ndarray:
    """Compute weights for each sample for weighted averaging.

    For even sampling, weights are uniform (1/n_samples).
    This function can be extended for importance sampling.

    Args:
        x_samples: Sampled positions
        p_samples: Sampled momenta
        x_grid: Original position grid
        psi: Wavefunction
        mode: Sampling mode

    Returns:
        weights: Normalized weights summing to 1
    """
    n_samples = len(x_samples)
    return np.ones(n_samples) / n_samples
