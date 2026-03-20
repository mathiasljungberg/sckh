sckh
====

Programs for calculation of XAS and XES including vibrational effects

Fortran programs:

- **XAS** - quantum mechanical Fermi Golden Rule for vibrational effects in XAS, 1D PES
- **KH** - quantum mechanical Kramers-Heisenberg equations for vibrational effects in XES, 1D PES
- **KH_resonant** - resonant KH, 1D PES
- **KH_resonant_el** - resonant KH, electronic part only (no vibrations)
- **SCKH** - semiclassical KH for XES, using external trajectories
- **SCKH_PES** - semiclassical KH for XES, using 1D PES

Python packages:

- **dynamics_1d** - classical dynamics on 1D PES with Fourier grid vibrational solver
- **dynamics_2d** - 2D polynomial PES with spline interpolation
- **kh_1d** - 1D Kramers-Heisenberg calculations


Directory structure
===================

```
src/
├── fortran/            # Fortran source code and CMakeLists.txt
├── python/             # Python packages (dynamics_1d, dynamics_2d, kh_1d, tools)
tests/
├── test_dynamics_1d/   # Python tests
├── test_dynamics_2d/
├── test_kh_1d/
├── testsuite/          # Fortran integration tests
scripts/                # Streamlit app for 2D surface visualization
examples/               # Example configurations and data
utils/                  # Legacy standalone programs (archived, see utils/README.md)
```


Installation
============

```bash
# Install Python package with dev dependencies
uv sync --extra dev
```


Building Fortran code
=====================

Requires a Fortran compiler (gfortran or ifort), BLAS, LAPACK, and FFTW3.

```bash
cmake -B build src/fortran
cmake --build build
```

This places the compiled binaries (`sckh_main`, `vib_finite_diff`, etc.) in the `build/` directory.

To build a single executable:

```bash
cmake --build build --target sckh_main
```

CMake options can be passed to customize the build, e.g.:

```bash
cmake -B build src/fortran -DCMAKE_Fortran_COMPILER=gfortran -DFFTW_LIB="-L/usr/lib -lfftw3 -lm"
```


Testing
=======

All tests are run via pytest. By default, only Python tests run.

Python tests
------------

```bash
# Run Python tests (default)
uv run pytest

# Run with coverage report
uv run pytest --cov=dynamics_1d --cov-report=term
```

Fortran tests
-------------

Fortran tests require compiled binaries in `build/` (see above).

```bash
# Run Fortran integration tests only
uv run pytest -m fortran

# Run all tests (Python + Fortran)
uv run pytest -o "addopts="
```

The Fortran tests automatically look for binaries in `build/` at the project root.
This can be overridden by setting the `SCKH_PATH` environment variable.


