"""Smoke tests for the current public import surface."""

import os
from pathlib import Path
import subprocess
import sys


def test_import_does_not_change_host_streams_or_environment():
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "python")
    env.pop("POLYFEMPY_WINDOWS_AUTOCONFIG", None)
    env.pop("POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG", None)
    env.pop("KMP_DUPLICATE_LIB_OK", None)
    env["PYTHONIOENCODING"] = "cp1252"
    code = '''
import os, sys
from unittest.mock import patch
before = (sys.stdout, sys.stderr)
with patch("os.system") as shell:
    import polyfempy
assert sys.stdout is before[0]
assert sys.stderr is before[1]
assert "KMP_DUPLICATE_LIB_OK" not in os.environ
shell.assert_not_called()
'''
    subprocess.run([sys.executable, "-c", code], cwd=root, env=env, check=True)


def test_unsupported_callback_fails_before_loading_backend():
    import pytest
    from polyfempy import solve

    with pytest.raises(TypeError, match="sidesets_func"):
        solve(sidesets_func=lambda *_: 1)


def test_public_api_imports():
    from polyfempy import Result, solve

    assert solve is not None
    assert Result is not None


def test_polyfem_is_an_automatic_alias():
    import polyfem
    import polyfempy

    assert polyfem.__all__ == polyfempy.__all__
    for name in polyfempy.__all__:
        assert getattr(polyfem, name) is getattr(polyfempy, name)


def test_generated_api_package_is_the_preferred_generated_surface():
    from polyfempy.generated_api import generated_api as polyfem
    from polyfempy.generated_api import generated_class

    assert polyfem.Root is generated_class.Root


def test_legacy_api_and_generated_packages_are_removed():
    from pathlib import Path

    package_dir = Path(__file__).resolve().parents[1] / "python" / "polyfempy"

    assert not (package_dir / "api").exists()
    assert not (package_dir / "generated").exists()
    assert not (package_dir / "runtime" / "_solve_pipeline.py").exists()
    assert not (package_dir / "runtime" / "_solve_outputs.py").exists()


def test_solve_module_all_only_recommends_solve():
    import importlib
    import inspect

    solve_module = importlib.import_module("polyfempy.solve")

    assert solve_module.__all__ == ["Result", "solve"]
    assert "sampled_vtu_fallback" not in inspect.signature(solve_module.solve).parameters
    assert not hasattr(solve_module, "populate_sampled_history")


def test_import_does_not_load_native_extension():
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "python")
    code = """
import importlib.abc
import sys
class RejectNative(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if fullname == 'polyfempy.polyfempy':
            raise AssertionError('Import eagerly loaded native bindings')
sys.meta_path.insert(0, RejectNative())
import polyfempy
assert callable(polyfempy.model)
assert callable(polyfempy.solve)
assert 'Path' not in polyfempy.__all__
assert 'builtins' not in polyfempy.__all__
"""
    subprocess.run([sys.executable, "-c", code], cwd=root, env=env, check=True)
