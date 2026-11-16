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

    with pytest.raises(NotImplementedError, match="sidesets_func"):
        solve(sidesets_func=lambda *_: 1)


def test_public_api_imports():
    from polyfempy.runtime import Result, solve

    assert solve is not None
    assert Result is not None


def test_runtime_imports_are_the_preferred_solve_surface():
    import polyfempy.runtime as runtime
    from polyfempy.runtime import Result, solve

    assert runtime.CORE_RUNTIME == ["solve", "Result"]
    assert runtime.__all__ == runtime.CORE_RUNTIME
    assert solve is not None
    assert Result is not None


def test_generated_api_package_is_the_preferred_generated_surface():
    from polyfempy.generated_api import generated_api as polyfem
    from polyfempy.generated_api import generated_class

    assert polyfem.Root is generated_class.Root


def test_runtime_recommended_surface_is_small():
    import polyfempy.runtime as runtime

    assert runtime.CORE_RUNTIME == ["solve", "Result"]
    assert runtime.__all__ == runtime.CORE_RUNTIME

    for name in (
        "SimulationConfig",
        "Output",
        "ParaviewOutput",
        "result_output",
        "solve_with_timing",
        "configure_windows_runtime",
        "batch_solve",
        "HistoryView",
    ):
        assert not hasattr(runtime, name)


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

    solve_module = importlib.import_module("polyfempy.runtime.solve")

    assert solve_module.__all__ == ["solve"]
    assert "sampled_vtu_fallback" not in inspect.signature(solve_module.solve).parameters
    assert not hasattr(solve_module, "populate_sampled_history")
