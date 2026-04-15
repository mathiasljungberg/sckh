"""Regression test for the SCKH_PES_2D standalone examples."""

from pathlib import Path

import pytest

from testsuite_python.helpers import assert_runner_matches_reference


HERE = Path(__file__).resolve().parent


@pytest.mark.parametrize("runner", ["run_python", "run_yaml"])
def test_runner_matches_reference(runner: str) -> None:
    assert_runner_matches_reference(HERE, runner)
