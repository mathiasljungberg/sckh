import pytest


def pytest_collection_modifyitems(config, items):
    """Auto-mark tests based on their location."""
    for item in items:
        path = str(item.fspath)
        if "/tests/testsuite/" in path:
            item.add_marker(pytest.mark.fortran)
