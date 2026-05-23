"""Shared helpers for tests/testsuite_python case tests."""

import importlib.util
from pathlib import Path
from types import ModuleType

import numpy as np


def load_runner(case_dir: Path, name: str) -> ModuleType:
    """Load a case-local runner module from ``run_*.py``."""
    spec = importlib.util.spec_from_file_location(name, case_dir / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def assert_runner_matches_reference(case_dir: Path, runner: str) -> None:
    """Run a case-local example script and compare ``rundir`` to ``ref``."""
    ref_dir = case_dir / "ref"
    run_dir = case_dir / "rundir"

    load_runner(case_dir, runner).main()

    ref_files = sorted(ref_dir.glob("*.dat"))
    assert ref_files, "No reference files found in ref/"

    for ref_path in ref_files:
        out_path = run_dir / ref_path.name
        assert out_path.exists(), f"Missing output: {out_path.name}"

        ref = np.loadtxt(ref_path)
        out = np.loadtxt(out_path)
        # Reference files are stored with ~9 significant figures (%.8E) and the
        # ground-state eigenvector solve is sensitive to the LAPACK/BLAS backend,
        # so results vary at the ~1e-6 level across platforms. Use tolerances
        # that absorb this cross-platform noise while still catching regressions.
        np.testing.assert_allclose(
            out, ref, atol=1e-8, rtol=1e-5,
            err_msg=f"Mismatch in {ref_path.name}",
        )
