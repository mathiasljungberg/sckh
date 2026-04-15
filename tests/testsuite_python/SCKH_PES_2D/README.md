# SCKH_PES_2D Example Case

This directory is both:

- a runnable example of the Python `dynamics_2d` + `sckh` workflow
- a regression test case

It demonstrates the two-stage 2D workflow:

1. run 2D dynamics to produce trajectories and save them
2. interpolate PES and dipole surfaces along those trajectories and compute the spectrum

Files:

- `run_python.py`: pure Python example that builds the configuration in code
- `run_yaml.py`: YAML-driven example using `sckh_pes_2d.yaml`
- `input/`: self-contained 2D PES and dipole input files
- `ref/`: committed reference output used by the regression test
- `test_sckh_pes_2d.py`: test wrapper that runs both examples and compares output
  to `ref/`

Run manually:

```bash
python run_python.py
python run_yaml.py
```

Both scripts write their output to `rundir/`.
