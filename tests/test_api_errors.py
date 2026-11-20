"""Invalid user inputs should produce useful errors without a native solve."""

import json

import numpy as np
import pytest

from polyfempy import solve
from polyfempy.runtime import _solve_backend as backend
from polyfempy.runtime._solve_contract import choose_mesh_source


VERTICES = [[0, 0], [1, 0], [0, 1]]


@pytest.mark.parametrize("vertices,cells,message", [
    (VERTICES, [[0, 1, 1.9]], "fractional"),
    (VERTICES, [[0, 1, 3]], "zero-based"),
    (VERTICES, [[0, 1, -1]], "zero-based"),
    (VERTICES, [[0, 1, 2**40]], "zero-based"),
    ([[0, 0], [float("nan"), 0]], [[0, 1, 1]], "finite real"),
    ([[0, 0], [float("inf"), 0]], [[0, 1, 1]], "finite real"),
    ([0, 1, 2], [[0, 1, 2]], "shape"),
    ([], [[0, 1, 2]], "shape"),
    (VERTICES, [0, 1, 2], "2D"),
    (VERTICES, [], "non-empty"),
])
def test_bad_arrays_fail_before_backend_creation(monkeypatch, vertices, cells, message):
    def unexpected_backend():
        pytest.fail("Invalid mesh reached the native backend")
    monkeypatch.setattr(backend, "build_solver", unexpected_backend)
    with pytest.raises(ValueError, match=message):
        solve(vertices=vertices, cells=cells, cfg={})


def test_integer_valued_connectivity_is_preserved():
    source = choose_mesh_source(VERTICES, [[0.0, 1.0, 2.0]], None)
    np.testing.assert_array_equal(source.cells, [[0, 1, 2]])
    assert source.cells is not None
    assert source.cells.dtype == np.int32


def test_config_path_errors_name_the_file(tmp_path):
    missing = tmp_path / "missing.json"
    with pytest.raises(FileNotFoundError, match="Configuration file not found") as error:
        solve(cfg=missing)
    assert str(missing) in str(error.value)
    bad = tmp_path / "bad.json"
    bad.write_text('{"geometry": }', encoding="utf-8")
    with pytest.raises(ValueError, match="line 1, column") as error:
        solve(cfg=bad)
    assert str(bad) in str(error.value)


def test_rejected_elastic_material_has_parameter_hint():
    class RejectingSolver:
        def set_settings(self, *args, **kwargs):
            raise RuntimeError("Invalid input json file")
    payload = {"materials": [{"type": "NeoHookean", "E": 1000}]}
    with pytest.raises(RuntimeError, match=r"materials\[0\].*both E and nu"):
        backend._apply_settings_json(RejectingSolver(), json.dumps(payload))
    payload["materials"][0].update({"lambda": 1, "mu": 2})
    with pytest.raises(RuntimeError, match="^Invalid input json file$"):
        backend._apply_settings_json(RejectingSolver(), json.dumps(payload))


def test_successful_backend_settings_are_not_revalidated():
    class AcceptingSolver:
        def set_settings(self, value, strict_validation):
            assert strict_validation is True
            assert json.loads(value) == {"materials": {"type": "CustomMaterial"}}
    backend._apply_settings_json(AcceptingSolver(), '{"materials": {"type": "CustomMaterial"}}')
