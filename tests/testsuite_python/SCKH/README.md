# SCKH Example Case

This directory is both:

- a runnable example of the Python `sckh` API
- a regression test case

It demonstrates the pure SCKH workflow starting from precomputed trajectory
files.

Files:

- `run_python.py`: pure Python example using `read_sckh_trajectories(...)` and
  `SCKHSpectrumCalculator`
- `run_yaml.py`: YAML-driven example using `sckh.yaml`
- `input/`: self-contained input trajectory files
- `ref/`: committed reference output used by the regression test
- `test_sckh.py`: test wrapper that runs both examples and compares output to
  `ref/`

Run manually:

```bash
python run_python.py
python run_yaml.py
```

Both scripts write their output to `rundir/`.
