import os
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parents[2]
DEFAULT_BUILD_DIR = PROJECT_ROOT / "build"


@pytest.fixture
def sckh_path():
    """Return path to Fortran binaries.

    Checks SCKH_PATH env var first, then falls back to <project_root>/build/.
    Skips if neither contains the binaries.
    """
    # Explicit env var takes precedence
    path = os.environ.get("SCKH_PATH")
    if path and Path(path).is_dir():
        return path

    # Fall back to build/ in project root
    if (DEFAULT_BUILD_DIR / "sckh_main").exists():
        path = str(DEFAULT_BUILD_DIR)
        os.environ["SCKH_PATH"] = path
        return path

    pytest.skip(
        "Fortran binaries not found. "
        "Build with 'cmake -B build src/fortran && cmake --build build' "
        "or set SCKH_PATH."
    )
