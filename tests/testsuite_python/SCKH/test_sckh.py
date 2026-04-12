"""Regression test for the SCKH standalone examples.

Invokes each runner script (``run_python.py`` and ``run_yaml.py``),
then diffs the files produced in ``rundir/`` against the committed
``ref/`` outputs.

Both ``ref/`` and the runners' ``rundir/`` are produced by the same
Python code path, so a strict numerical comparison is appropriate.
"""

import importlib.util
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
REF = HERE / "ref"
RUNDIR = HERE / "rundir"


def _load_runner(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("runner", ["run_python", "run_yaml"])
def test_runner_matches_reference(runner: str) -> None:
    _load_runner(runner).main()

    ref_files = sorted(REF.glob("*.dat"))
    assert ref_files, "No reference files found in ref/"

    for ref_path in ref_files:
        out_path = RUNDIR / ref_path.name
        assert out_path.exists(), f"Missing output: {out_path.name}"

        ref = np.loadtxt(ref_path)
        out = np.loadtxt(out_path)
        np.testing.assert_allclose(
            out, ref, atol=1e-12, rtol=1e-10,
            err_msg=f"Mismatch in {ref_path.name}",
        )
