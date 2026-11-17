"""Keep model selections attached when solve receives an array mesh."""

import copy
import json

import numpy as np
import pytest

import polyfem
from polyfempy import solve
from polyfempy.runtime import _solve_backend as backend


VERTICES = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]], dtype=float)
CELLS = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])


class RecordingSolver:
    def __init__(self):
        self.calls = []

    def set_settings(self, settings, strict_validation):
        assert strict_validation is True
        self.settings = json.loads(settings)
        self.calls.append("settings")

    def set_mesh(self, vertices, cells):
        self.calls.append("raw arrays")

    def set_settings_with_mesh(self, settings, vertices, cells):
        self.set_settings(settings, strict_validation=True)
        self.calls.append("configured arrays")
        np.testing.assert_array_equal(vertices, VERTICES)
        np.testing.assert_array_equal(cells, CELLS)
        assert cells.dtype == np.int32

    def load_mesh_from_settings(self):
        self.calls.append("file")

    def solve(self, log_level):
        self.calls.append("solve")
        return {"_result_bundle": True, "sol": np.zeros((10, 1))}


@pytest.mark.parametrize("config_kind", ["model", "dict", "file"])
@pytest.mark.parametrize("whole_boundary", [False, True])
def test_solve_preserves_automatic_material_and_boundary_links(
    monkeypatch, tmp_path, config_kind, whole_boundary,
):
    model = polyfem.model()
    body = model.mesh(mesh="array.msh")
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    if whole_boundary:
        body.dirichlet(value=[0, 0])
    else:
        body.surface_axis(axis="-x", position=0, relative=True).dirichlet(value=[0, 0])
        body.surface_axis(axis="x", position=1, relative=True).neumann(value=[1, 0])
    config = model.config(rhs=[0, -1])
    original = config.as_dict()
    if config_kind == "dict":
        config = copy.deepcopy(original)
    elif config_kind == "file":
        config = tmp_path / "simulation.json"
        config.write_text(json.dumps(original), encoding="utf-8")
    solver = RecordingSolver()
    monkeypatch.setattr(backend, "build_solver", lambda: solver)

    solve(VERTICES, CELLS, cfg=config)

    assert solver.calls == ["settings", "configured arrays", "solve"]
    geometry = solver.settings["geometry"][0]
    assert geometry == original["geometry"][0]
    assert solver.settings["materials"][0]["id"] == geometry["volume_selection"]
    conditions = solver.settings["boundary_conditions"]
    selection = geometry["surface_selection"]
    if whole_boundary:
        assert conditions["dirichlet_boundary"][0]["id"] == selection
    else:
        assert conditions["dirichlet_boundary"][0]["id"] == selection[0]["id"]
        assert conditions["neumann_boundary"][0]["id"] == selection[1]["id"]
        assert selection[0]["id"] != selection[1]["id"]
    assert model.config(rhs=[0, -1]).as_dict() == original
    if config_kind == "dict":
        assert config == original


def test_array_mesh_keeps_obstacles_and_skips_disabled_bodies(monkeypatch):
    config = {
        "geometry": [
            {"mesh": "unused.obj", "enabled": False},
            {"mesh": "array.msh", "surface_selection": 4},
            {"type": "ground", "height": -1, "is_obstacle": True},
        ],
        "materials": {"type": "NeoHookean", "E": 1000, "nu": 0.3},
    }
    original = copy.deepcopy(config)
    solver = RecordingSolver()
    monkeypatch.setattr(backend, "build_solver", lambda: solver)

    solve(VERTICES, CELLS, cfg=config)

    assert solver.calls == ["settings", "configured arrays", "solve"]
    assert solver.settings["geometry"] == original["geometry"]
    assert config == original


@pytest.mark.parametrize("geometry,message", [
    ([{"mesh": "first.obj"}, {"mesh": "second.obj"}], "exactly one"),
    ([{"mesh": "unused.obj", "enabled": False}], "exactly one"),
    ([{"type": "ground", "height": 0, "is_obstacle": True}], "exactly one"),
    ([{"type": "mesh_array", "mesh": "part.obj", "array": {}}], "type 'mesh'"),
])
def test_ambiguous_array_geometry_fails_before_loading_settings(monkeypatch, geometry, message):
    solver = RecordingSolver()
    monkeypatch.setattr(backend, "build_solver", lambda: solver)

    with pytest.raises(ValueError, match=message):
        solve(VERTICES, CELLS, cfg={"geometry": geometry})

    assert solver.calls == []


def test_old_backend_reports_missing_array_selection_support(monkeypatch):
    class OldSolver:
        def set_settings(self, *args, **kwargs):
            pytest.fail("An incompatible backend must fail before settings are applied")

    monkeypatch.setattr(backend, "build_solver", OldSolver)
    with pytest.raises(RuntimeError, match="cannot apply geometry selections.*array meshes"):
        solve(VERTICES, CELLS, cfg={"geometry": {"mesh": "array.msh"}})


def test_arrays_without_geometry_keep_the_existing_backend_path(monkeypatch):
    class OldSolver(RecordingSolver):
        def __getattribute__(self, name):
            if name == "set_settings_with_mesh":
                raise AttributeError(name)
            return super().__getattribute__(name)

    solver = OldSolver()
    monkeypatch.setattr(backend, "build_solver", lambda: solver)
    solve(VERTICES, CELLS, cfg={"materials": {"type": "NeoHookean", "E": 1000, "nu": 0.3}})
    assert solver.calls == ["settings", "raw arrays", "solve"]


def test_file_mesh_still_uses_the_file_loader(monkeypatch):
    solver = RecordingSolver()
    monkeypatch.setattr(backend, "build_solver", lambda: solver)
    solve(cfg={"geometry": {"mesh": "part.msh"}})
    assert solver.calls == ["settings", "file", "solve"]


@pytest.mark.parametrize("whole_boundary", [False, True])
def test_native_array_selections_match_a_file_solve(tmp_path, monkeypatch, whole_boundary):
    import polyfempy

    if not polyfempy.cpp_backend_available():
        pytest.skip("C++ backend is unavailable; array selection integration requires a native build")

    monkeypatch.chdir(tmp_path)
    mesh = tmp_path / "square.obj"
    mesh.write_text(
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv .5 .5 0\n"
        "f 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n",
        encoding="utf-8",
    )
    model = polyfem.model()
    body = model.mesh(mesh=str(mesh))
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    if whole_boundary:
        body.dirichlet(value=[0, 0])
    else:
        body.surface_axis(axis="-x", position=0, relative=True).dirichlet(value=[0, 0])
    model.solver(linear=polyfem.linear(solver="Eigen::SimplicialLDLT"))
    config = model.config(rhs=[0, -1])

    file_result = solve(cfg=config)
    # The array path must not open the configured FEM file.
    mesh.unlink()
    array_result = solve(VERTICES, CELLS, cfg=config)

    assert np.isfinite(array_result.sol).all()
    assert abs(array_result.sol).max() > 1e-6
    # File decoding can change the numbering of the same finite-element nodes.
    np.testing.assert_allclose(
        np.sort(array_result.sol.ravel()), np.sort(file_result.sol.ravel()), rtol=1e-6, atol=1e-9,
    )
