"""Options for tests that use an existing native build."""

import sys

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--backend-package-dir",
        help="Package root containing polyfempy and the native backend to test.",
    )


@pytest.fixture(autouse=True)
def legacy_example_import(monkeypatch):
    # The pinned external example repository predates the single package name.
    import polyfempy
    monkeypatch.setitem(sys.modules, "polyfem", polyfempy)
