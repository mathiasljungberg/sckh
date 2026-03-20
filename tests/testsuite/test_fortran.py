"""Pytest wrapper for Fortran integration tests.

Each test directory contains a tests.py with a run_test() function that
returns 0 on success. This module parametrizes them as individual pytest cases.
"""

import importlib.util
import os
import sys

import pytest

# Tests from run_tests.py (known to run fast enough)
FORTRAN_TESTS = [
    "KH",
    "KH_resonant",
    "KH_resonant_el",
    "SCKH",
    "SCKH_PES",
    "SCKH_resonant_PES",
    "XAS",
    "SCXAS_PES",
    "vib_finite_diff",
]

TESTSUITE_DIR = os.path.dirname(__file__)


@pytest.mark.fortran
@pytest.mark.parametrize("test_name", FORTRAN_TESTS)
def test_fortran_integration(test_name, sckh_path):
    """Run a Fortran integration test by invoking its tests.py run_test()."""
    test_dir = os.path.join(TESTSUITE_DIR, test_name)
    test_file = os.path.join(test_dir, "tests.py")

    if not os.path.exists(test_file):
        pytest.skip(f"Test file not found: {test_file}")

    # Ensure test_mod is importable from testsuite directory
    if TESTSUITE_DIR not in sys.path:
        sys.path.insert(0, TESTSUITE_DIR)

    # Load the test module
    spec = importlib.util.spec_from_file_location(
        f"fortran_test_{test_name}", test_file
    )
    module = importlib.util.module_from_spec(spec)

    # Run from the test directory (tests.py uses relative paths)
    old_dir = os.getcwd()
    try:
        os.chdir(test_dir)
        spec.loader.exec_module(module)
        err_code = module.run_test()
    finally:
        os.chdir(old_dir)

    assert err_code == 0, f"Fortran test '{test_name}' failed with error code {err_code}"
