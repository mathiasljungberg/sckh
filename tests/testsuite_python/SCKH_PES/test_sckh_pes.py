"""Regression test for the SCKH_PES standalone examples.

Invokes each runner script (``run_python.py`` and ``run_yaml.py``),
then diffs the files produced in ``rundir/`` against the committed
``ref/`` outputs.

Both ``ref/`` and the runners' ``rundir/`` are produced by the same
Python code path, so a strict numerical comparison is appropriate.
"""

from pathlib import Path

import pytest

from testsuite_python.helpers import assert_runner_matches_reference


HERE = Path(__file__).resolve().parent


@pytest.mark.parametrize("runner", ["run_python", "run_yaml"])
def test_runner_matches_reference(runner: str) -> None:
    assert_runner_matches_reference(HERE, runner)
