"""Regression test for the DENSITY_2D standalone example."""

from pathlib import Path

from testsuite_python.helpers import assert_runner_matches_reference


HERE = Path(__file__).resolve().parent


def test_runner_matches_reference() -> None:
    assert_runner_matches_reference(HERE, "run_density")
