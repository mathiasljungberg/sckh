sckh
====

Programs for calculation of XAS and XES including vibrational effects

The main development is done in the src directory, and the code pieces will be moved there evenutally. 
Right now the following routines are avaiable in the main program:


XAS - quantum mechanical solution of the Fermi Golden Rule for vibrational effects in XAS, 1d PES.
KH - quantum mechanical solution of the Kramers-Heisenberg equations for vibrational effects in XES, 1d PES.
KH_resonant - resonant KH, 1d PES.
KH_resonant_el - resoant KH, only electronic part, no vibrations.
SCKH  - semiclassical solution of the KH equations for XES, using external trajectories.
SCKH_PES  - semiclassical solution of the KH equations for XES, using 1d PES. 

Next is to get the non-adiabatic effects working, as well as the resonant semiclassical case.
For the time being specialized working programs can be found in the util directory.

Working codes (in util directory):

KH    - quantum mechanical solution of the Kramers-Heisenberg equations for vibrational effects in XES
SCKH  - semiclassical solution of the KH equations for XES
sinc_DVR  - solution of 1d vibrational problem on a grid with a sinc DVR
DVR_PO  - potential optimized DV
Kramers-Heisenberg_resonant - like KH but for the resonant case (in the SCKH_resonant folder)

these codes lack a logical build sequence right now...

Codes under development (in util directory):

XAS_nonadiabatic  - XAS and XAS including nonadiabatic effects. Also attempt to unify codes into one
SCKH_resonant  - like SCKH but for the reoannt case


Testing
=======

Python tests
------------

The Python scripts in `src/python_scripts/` have pytest-based tests.

To run the tests:

```bash
# Install test dependencies (first time only)
uv pip install pytest pytest-cov

# Run all Python tests
uv run pytest tests/test_dynamics_2d.py -v

# Run with coverage report
uv run pytest tests/test_dynamics_2d.py --cov=python_scripts.dynamics_2d --cov-report=term
```

Fortran tests
-------------

The Fortran test suite uses a custom test runner. See `tests/testsuite/README` for details.


