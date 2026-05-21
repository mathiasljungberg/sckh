"""Tests for trajectory density on the 2D dynamics grid."""

from pathlib import Path

import numpy as np
import pytest

from dynamics_1d.constants import CONST
from dynamics_2d.config import (
    DynamicsConfig2D,
    GridConfig2D,
    SamplingConfig2D,
    TimeConfig,
)
from dynamics_2d.density import (
    TrajectoryDensityResult2D,
    ensemble_density_at_step,
    ensemble_density_timeseries,
    fwhm_to_sigma,
    gaussian_2d_on_grid,
    positions_to_grid_density,
    trajectories_density_at_step,
    trajectories_density_timeseries,
)
from dynamics_2d.io import (
    read_trajectories_2d,
    write_density_2d,
    write_density_timeseries,
    write_trajectory_2d,
)
from dynamics_2d.trajectory import EnsembleResult2D, TrajectoryResult2D
from dynamics_2d.vibrational import ProductGroundState


def _make_ensemble(
    nsteps: int = 5,
    n_traj: int = 4,
    nx1: int = 41,
    nx2: int = 41,
    grid_half_width_ang: float = 0.5,
):
    """Build a small synthetic EnsembleResult2D for tests.

    Trajectories sit at fixed positions on the (-0.1, +0.1)Å square so they
    are well inside the grid; FWHM in tests is chosen large enough that
    discretization error is small.
    """
    x1_grid = np.linspace(-grid_half_width_ang, grid_half_width_ang, nx1) * 1e-10
    x2_grid = np.linspace(-grid_half_width_ang, grid_half_width_ang, nx2) * 1e-10
    time = np.linspace(0.0, 1e-15 * nsteps, nsteps)

    rng = np.random.default_rng(0)
    fixed_positions = (rng.uniform(-0.1, 0.1, n_traj) * 1e-10,
                       rng.uniform(-0.1, 0.1, n_traj) * 1e-10)

    trajectories = []
    for i in range(n_traj):
        x1 = np.full(nsteps, fixed_positions[0][i])
        x2 = np.full(nsteps, fixed_positions[1][i])
        trajectories.append(
            TrajectoryResult2D(
                time=time,
                x1=x1,
                x2=x2,
                v1=np.zeros(nsteps),
                v2=np.zeros(nsteps),
                a1=np.zeros(nsteps),
                a2=np.zeros(nsteps),
                x1_0=float(x1[0]),
                x2_0=float(x2[0]),
                p1_0=0.0,
                p2_0=0.0,
            )
        )

    ground_state = ProductGroundState(
        E1=0.0,
        E2=0.0,
        E_total=0.0,
        psi1=np.zeros(nx1),
        psi2=np.zeros(nx2),
        x1_grid=x1_grid,
        x2_grid=x2_grid,
    )

    config = DynamicsConfig2D(
        mu1=1.0,
        mu2=2.0,
        grid_x1=GridConfig2D(
            start=-grid_half_width_ang, dx=2 * grid_half_width_ang / (nx1 - 1),
            npoints=nx1
        ),
        grid_x2=GridConfig2D(
            start=-grid_half_width_ang, dx=2 * grid_half_width_ang / (nx2 - 1),
            npoints=nx2
        ),
        time=TimeConfig(dt=1.0, nsteps=nsteps),
        sampling=SamplingConfig2D(),
        pes_dynamics=Path("dummy.dat"),
        pes_initial=Path("dummy_initial.dat"),
    )

    return EnsembleResult2D(
        trajectories=trajectories,
        config=config,
        ground_state=ground_state,
        x1_grid=x1_grid,
        x2_grid=x2_grid,
        x1_eq=0.0,
        x2_eq=0.0,
    )


def _integral(density: np.ndarray, x1: np.ndarray, x2: np.ndarray) -> float:
    dx1 = x1[1] - x1[0]
    dx2 = x2[1] - x2[0]
    return float(density.sum() * dx1 * dx2)


def test_fwhm_to_sigma_roundtrip():
    fwhm = 0.123
    sigma = fwhm_to_sigma(fwhm)
    assert sigma == pytest.approx(fwhm / (2 * np.sqrt(2 * np.log(2))))


def test_gaussian_2d_integrates_to_one():
    x1 = np.linspace(-1.0, 1.0, 401) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 401) * 1e-10
    sigma = 0.1e-10  # fully contained
    g = gaussian_2d_on_grid(x1, x2, 0.0, 0.0, sigma)
    assert _integral(g, x1, x2) == pytest.approx(1.0, abs=1e-4)


def test_gaussian_2d_peak_location():
    x1 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x1_0 = 0.3e-10
    x2_0 = -0.2e-10
    g = gaussian_2d_on_grid(x1, x2, x1_0, x2_0, sigma=0.05e-10)
    i, j = np.unravel_index(np.argmax(g), g.shape)
    assert x1[i] == pytest.approx(x1_0, abs=(x1[1] - x1[0]))
    assert x2[j] == pytest.approx(x2_0, abs=(x2[1] - x2[0]))


def test_positions_to_grid_single_point_matches_gaussian():
    x1 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 201) * 1e-10
    sigma = 0.1e-10
    g = gaussian_2d_on_grid(x1, x2, 0.0, 0.0, sigma)
    g_renorm = g / _integral(g, x1, x2)

    d = positions_to_grid_density(x1, x2, np.array([0.0]), np.array([0.0]), sigma)

    np.testing.assert_allclose(d, g_renorm, rtol=1e-10, atol=1e-14)
    assert _integral(d, x1, x2) == pytest.approx(1.0, abs=1e-10)


def test_positions_to_grid_normalization():
    x1 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 201) * 1e-10
    rng = np.random.default_rng(1)
    x1_pts = rng.uniform(-0.3, 0.3, 25) * 1e-10
    x2_pts = rng.uniform(-0.3, 0.3, 25) * 1e-10
    d = positions_to_grid_density(x1, x2, x1_pts, x2_pts, sigma=0.05e-10)
    assert _integral(d, x1, x2) == pytest.approx(1.0, abs=1e-10)


def test_positions_to_grid_two_points_symmetry():
    x1 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 201) * 1e-10
    x1_pts = np.array([0.2e-10, -0.2e-10])
    x2_pts = np.array([0.0, 0.0])
    d = positions_to_grid_density(x1, x2, x1_pts, x2_pts, sigma=0.05e-10)
    np.testing.assert_allclose(d, d[::-1, :], atol=1e-14)


def test_ensemble_density_timeseries_shape_and_selection():
    ensemble = _make_ensemble(nsteps=6, n_traj=3, nx1=41, nx2=41)

    res_all = ensemble_density_timeseries(ensemble, fwhm=0.1, position_units="angstrom")
    assert res_all.density.shape == (6, 41, 41)
    assert res_all.step_indices.tolist() == list(range(6))
    np.testing.assert_allclose(res_all.time, ensemble.get_time_array())

    res_sub = ensemble_density_timeseries(
        ensemble, fwhm=0.1, step_indices=[0, 3, 5], position_units="angstrom"
    )
    assert res_sub.density.shape == (3, 41, 41)
    assert res_sub.step_indices.tolist() == [0, 3, 5]
    for slice_ in res_sub.density:
        assert _integral(slice_, ensemble.x1_grid, ensemble.x2_grid) == pytest.approx(
            1.0, abs=1e-10
        )


def test_ensemble_density_unit_conversion():
    ensemble = _make_ensemble(nsteps=2, n_traj=3, nx1=41, nx2=41)
    fwhm_ang = 0.12
    fwhm_bohr = fwhm_ang * 1e-10 / CONST.bohr

    d_ang = ensemble_density_at_step(ensemble, step=0, fwhm=fwhm_ang,
                                     position_units="angstrom")
    d_bohr = ensemble_density_at_step(ensemble, step=0, fwhm=fwhm_bohr,
                                      position_units="bohr")
    np.testing.assert_allclose(d_ang, d_bohr, rtol=1e-12, atol=1e-14)


def test_write_density_2d_roundtrip(tmp_path):
    x1 = np.linspace(-1.0, 1.0, 41) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 41) * 1e-10
    sigma = 0.15e-10
    rho_SI = gaussian_2d_on_grid(x1, x2, 0.0, 0.0, sigma)  # 1/m^2
    rho_SI = rho_SI / _integral(rho_SI, x1, x2)

    path = tmp_path / "rho.dat"
    write_density_2d(
        path, x1, x2, rho_SI, time=1e-15, step_index=7, fwhm_SI=sigma,
        position_units="angstrom", index_order="C",
    )

    data = np.loadtxt(path)
    assert data.shape == (41 * 41, 3)

    # C ordering: x2 varies fastest. After reshape (41, 41) rows are x1 blocks.
    x1_col = data[:, 0].reshape(41, 41)
    x2_col = data[:, 1].reshape(41, 41)
    rho_col = data[:, 2].reshape(41, 41)

    np.testing.assert_allclose(x1_col[:, 0], x1 * 1e10)
    np.testing.assert_allclose(x2_col[0, :], x2 * 1e10)

    # density was converted to 1/Angstrom^2 (multiply SI value by 1e-20)
    np.testing.assert_allclose(rho_col, rho_SI * 1e-20, rtol=1e-6, atol=1e-30)

    # Integral in printed units is also 1
    dx1_ang = (x1[1] - x1[0]) * 1e10
    dx2_ang = (x2[1] - x2[0]) * 1e10
    assert rho_col.sum() * dx1_ang * dx2_ang == pytest.approx(1.0, abs=1e-4)


def test_write_density_2d_header_contents(tmp_path):
    x1 = np.linspace(-0.5, 0.5, 11) * 1e-10
    x2 = np.linspace(-0.5, 0.5, 11) * 1e-10
    rho = np.zeros((11, 11))
    rho[5, 5] = 1.0  # arbitrary placeholder

    path = tmp_path / "rho_header.dat"
    write_density_2d(
        path, x1, x2, rho, time=2.5e-15, step_index=42, fwhm_SI=0.1e-10,
        position_units="angstrom",
    )

    header_lines = [l for l in path.read_text().splitlines() if l.startswith("#")]
    joined = "\n".join(header_lines)
    assert "time(fs)" in joined
    assert "step_index = 42" in joined
    assert "fwhm(angstrom)" in joined
    assert "rho(1/angstrom^2)" in joined


def test_write_density_timeseries_files(tmp_path):
    ensemble = _make_ensemble(nsteps=5, n_traj=3, nx1=21, nx2=21)
    result = ensemble_density_timeseries(
        ensemble, fwhm=0.1, step_indices=[0, 2, 4], position_units="angstrom"
    )

    paths = write_density_timeseries(tmp_path, result, basename="rho")
    assert len(paths) == 3
    expected = [
        tmp_path / "rho_step_000000.dat",
        tmp_path / "rho_step_000002.dat",
        tmp_path / "rho_step_000004.dat",
    ]
    assert paths == expected
    for p in expected:
        assert p.exists()

    # Spot-check first file: header has the right step_index, integral ~ 1 in Å^-2.
    text = paths[0].read_text()
    assert "step_index = 0" in text

    data = np.loadtxt(paths[0])
    rho_col = data[:, 2].reshape(21, 21)
    dx_ang = (ensemble.x1_grid[1] - ensemble.x1_grid[0]) * 1e10
    dy_ang = (ensemble.x2_grid[1] - ensemble.x2_grid[0]) * 1e10
    assert rho_col.sum() * dx_ang * dy_ang == pytest.approx(1.0, abs=1e-3)


def test_positions_to_grid_density_validates_shapes():
    x1 = np.linspace(-1.0, 1.0, 21) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 21) * 1e-10
    with pytest.raises(ValueError, match="same shape"):
        positions_to_grid_density(x1, x2, np.array([0.0, 0.1e-10]), np.array([0.0]),
                                  sigma=0.05e-10)


def test_trajectories_density_matches_ensemble():
    ensemble = _make_ensemble(nsteps=4, n_traj=3, nx1=31, nx2=31)
    step_indices = [0, 2, 3]

    res_ens = ensemble_density_timeseries(
        ensemble, fwhm=0.1, step_indices=step_indices, position_units="angstrom"
    )
    res_list = trajectories_density_timeseries(
        ensemble.x1_grid, ensemble.x2_grid, ensemble.trajectories,
        fwhm=0.1, step_indices=step_indices, position_units="angstrom",
    )
    np.testing.assert_allclose(res_ens.density, res_list.density, rtol=1e-12, atol=1e-14)
    np.testing.assert_array_equal(res_ens.step_indices, res_list.step_indices)
    np.testing.assert_allclose(res_ens.time, res_list.time)

    d_ens = ensemble_density_at_step(ensemble, 2, fwhm=0.1, position_units="angstrom")
    d_list = trajectories_density_at_step(
        ensemble.x1_grid, ensemble.x2_grid, ensemble.trajectories,
        step=2, fwhm=0.1, position_units="angstrom",
    )
    np.testing.assert_allclose(d_ens, d_list, rtol=1e-12, atol=1e-14)


def test_trajectories_density_timeseries_empty_list_raises():
    x1 = np.linspace(-1.0, 1.0, 11) * 1e-10
    x2 = np.linspace(-1.0, 1.0, 11) * 1e-10
    with pytest.raises(ValueError, match="empty"):
        trajectories_density_timeseries(x1, x2, [], fwhm=0.1)


def test_read_trajectories_2d_roundtrip(tmp_path):
    ensemble = _make_ensemble(nsteps=3, n_traj=4, nx1=11, nx2=11)
    basename = "ens"
    for i, traj in enumerate(ensemble.trajectories):
        write_trajectory_2d(tmp_path / f"{basename}_traj_{i:04d}.dat", traj, units="user")

    read_back = read_trajectories_2d(tmp_path, basename)
    assert len(read_back) == ensemble.n_trajectories
    for original, loaded in zip(ensemble.trajectories, read_back):
        np.testing.assert_allclose(loaded.time, original.time, atol=1e-25)
        np.testing.assert_allclose(loaded.x1, original.x1, rtol=1e-6, atol=1e-20)
        np.testing.assert_allclose(loaded.x2, original.x2, rtol=1e-6, atol=1e-20)


def test_read_trajectories_2d_missing_dir_returns_empty(tmp_path):
    assert read_trajectories_2d(tmp_path, "nothing") == []


def test_density_from_saved_trajectories_matches_runtime(tmp_path):
    """Save → load → density equals the runtime ensemble density."""
    ensemble = _make_ensemble(nsteps=4, n_traj=3, nx1=21, nx2=21)
    basename = "ens"
    for i, traj in enumerate(ensemble.trajectories):
        write_trajectory_2d(tmp_path / f"{basename}_traj_{i:04d}.dat", traj, units="user")

    loaded = read_trajectories_2d(tmp_path, basename)
    res_loaded = trajectories_density_timeseries(
        ensemble.x1_grid, ensemble.x2_grid, loaded,
        fwhm=0.1, step_indices=[0, 2, 3], position_units="angstrom",
    )
    res_runtime = ensemble_density_timeseries(
        ensemble, fwhm=0.1, step_indices=[0, 2, 3], position_units="angstrom"
    )
    # Tolerance reflects the %16.8E roundtrip in write_trajectory_2d.
    np.testing.assert_allclose(res_loaded.density, res_runtime.density, rtol=1e-5, atol=1e-12)
