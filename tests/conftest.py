import pytest


def pytest_collection_modifyitems(config, items):
    """Auto-mark tests based on their location."""
    for item in items:
        if "testsuite" in str(item.fspath):
            item.add_marker(pytest.mark.fortran)
