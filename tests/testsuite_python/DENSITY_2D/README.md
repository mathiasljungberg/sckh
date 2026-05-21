# DENSITY_2D Example Case

Minimal example of placing a 2D trajectory swarm onto the dynamics grid
with a Gaussian-broadened density, normalized to one per time step.

Reads trajectories from `input/` (a snapshot of the trajectories produced
by the `SCKH_PES_2D` case) so this example is fully self-contained.

Files:

- `run_density.py`: builds the grid, reads trajectories, computes the
  density for a few selected time steps, writes one file per step.
- `input/`: committed trajectory files (`dynamics_2d_out_traj_*.dat`).
- `ref/`: committed reference density output used by the regression test.
- `test_density.py`: runs the example and compares `rundir/` to `ref/`.

Run manually:

```bash
python run_density.py
```

Output is written to `rundir/density_step_NNNNNN.dat`.
