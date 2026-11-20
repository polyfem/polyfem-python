"""Public solves configure the native solver once and retain its outputs."""

import sys
from types import SimpleNamespace
import json

import numpy as np
import pytest

from polyfempy import solve


@pytest.mark.parametrize("log_level", [4, "off", "debug"])
def test_solve_preserves_logging_and_runs_native_pipeline(monkeypatch, log_level):
    calls = []
    sol = np.arange(6.)

    class FakeSolver:
        def set_settings(self, settings, strict_validation):
            assert strict_validation is True
            assert json.loads(settings)["output"]["log"]["level"] == log_level
            calls.append("settings")

        def load_mesh_from_settings(self):
            calls.append("file")

        def solve(self):
            calls.append("solve")
            return sol

    monkeypatch.setitem(sys.modules, "polyfempy.polyfempy", SimpleNamespace(Solver=FakeSolver))
    cfg = {"geometry": [{"mesh": "triangle.obj"}], "output": {"log": {"level": log_level}}}
    result = solve(cfg=cfg)
    assert calls == ["settings", "file", "solve"]
    np.testing.assert_array_equal(result.sol, sol)
    assert result.meta == {}
