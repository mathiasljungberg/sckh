# SCKH_PES Example Case

This directory is both:

- a runnable example of the Python `dynamics_1d` + `sckh` workflow
- a regression test case

It demonstrates the two-stage 1D workflow:

1. run dynamics and surface interpolation to produce SCKH trajectories
2. compute the spectrum from those trajectories

Files:

- `run_python.py`: pure Python example that builds the configuration in code
- `run_yaml.py`: YAML-driven example using `sckh_pes.yaml`
- `input/`: self-contained PES and dipole input files
- `ref/`: committed reference output used by the regression test
- `test_sckh_pes.py`: test wrapper that runs both examples and compares output
  to `ref/`

Run manually:

```bash
python run_python.py
python run_yaml.py
```

Both scripts write their output to `rundir/`.
