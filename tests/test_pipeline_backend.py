from __future__ import annotations

import json
from typing import Any

import numpy as np
import pytest


@pytest.mark.parametrize("array_mesh", [False, True])
@pytest.mark.parametrize("log_level,expected", [(4, 4), ("off", 6)])
def test_solve_respects_log_level_for_both_mesh_inputs(monkeypatch, array_mesh, log_level, expected):
    from polyfempy import solve
    from polyfempy.runtime import _solve_backend as backend

    calls = {}

    class FakeSolver:
        def set_settings(self, settings, strict_validation):
            calls["settings"] = json.loads(settings)

        def load_mesh_from_settings(self):
            calls["mesh"] = "file"

        def set_mesh(self, vertices, cells):
            calls["mesh"] = "array"

        def solve(self, log_level):
            calls["log_level"] = log_level
            return {"_result_bundle": True, "sol": np.zeros((6, 1))}

    monkeypatch.setattr(backend, "build_solver", FakeSolver)
    config: dict[str, Any] = {"output": {"log": {"level": log_level}}}
    if array_mesh:
        solve([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]], cfg=config)
    else:
        config["geometry"] = [{"mesh": "triangle.obj"}]
        solve(cfg=config)

    assert calls["mesh"] == ("array" if array_mesh else "file")
    assert calls["settings"]["output"]["log"]["level"] == log_level
    assert calls["log_level"] == expected


def test_backend_adapter_exports_pipeline_runtime_helpers():
    from polyfempy.runtime import _solve_backend as backend

    expected = [
        "build_solver",
        "configure_solver",
        "apply_sidesets",
        "run_solver_stage",
    ]

    for name in expected:
        assert hasattr(backend, name)


def test_solve_module_does_not_reexport_backend_helpers():
    import importlib

    pipeline = importlib.import_module("polyfempy.runtime.solve")

    for name in (
        "build_solver",
        "configure_solver",
        "apply_sidesets",
        "run_solver_stage",
    ):
        assert not hasattr(pipeline, name)


def test_backend_adapter_no_longer_exports_legacy_context_or_retouch_helpers():
    from polyfempy.runtime import _solve_backend as backend

    assert not hasattr(backend, "SolverConfigContext")
    assert not hasattr(backend, "_retouch_bc_after_mesh")


def test_run_solver_stage_leaves_basis_and_assembly_to_backend_solve():
    from polyfempy.runtime import _solve_backend as backend

    class FakeSolver:
        def build_basis(self):
            raise AssertionError("Python runtime should not call build_basis")

        def assemble(self):
            raise AssertionError("Python runtime should not call assemble")

        def solve(self, log_level=None):
            return {"_result_bundle": True, "log_level": log_level}

    result = backend.run_solver_stage(
        FakeSolver(),
        {"output": {"log": {"level": "debug"}}},
    )

    assert result == {"_result_bundle": True, "log_level": 1}


def test_run_solver_stage_accepts_integer_log_level_from_config():
    from polyfempy.runtime import _solve_backend as backend

    class FakeSolver:
        def solve(self, log_level=None):
            return {"_result_bundle": True, "log_level": log_level}

    result = backend.run_solver_stage(
        FakeSolver(),
        {"output": {"log": {"level": 4}}},
    )

    assert result == {"_result_bundle": True, "log_level": 4}


def test_run_solver_stage_requires_varform_solve_entry_point():
    from polyfempy.runtime import _solve_backend as backend

    class LegacyRunOnlySolver:
        def run(self):
            return {"_result_bundle": True}

    try:
        backend.run_solver_stage(LegacyRunOnlySolver(), {})
    except RuntimeError as exc:
        assert "requires solve" in str(exc)
    else:
        raise AssertionError("run-only legacy backends must not be accepted")


def test_apply_sidesets_rejects_legacy_mesh_mutation_path():
    from polyfempy.runtime import _solve_backend as backend

    try:
        backend.apply_sidesets(object(), lambda *_: 1)
    except NotImplementedError as exc:
        assert "VarForm" in str(exc)
        assert "sidesets_func" in str(exc)
    else:
        raise AssertionError("sidesets_func should be rejected for the VarForm runtime")


def test_backend_adapter_source_has_no_legacy_solver_fallbacks():
    from pathlib import Path

    source = Path("python/polyfempy/runtime/_solve_backend.py").read_text(encoding="utf-8")

    assert '"State"' not in source
    assert '"run"' not in source
    assert "solver.settings(settings_json" not in source
    assert "solver.mesh()" not in source
    assert "set_mesh_data" not in source
    assert "load_mesh_from_points" not in source
