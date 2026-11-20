"""PolyFEM resolves mesh and selection resources from the JSON file's directory."""

import json

import numpy as np
import pytest

from polyfempy import solve


@pytest.mark.parametrize("geometry_list", [False, True])
def test_backend_resolves_config_resources(tmp_path, monkeypatch, geometry_list):
    pytest.importorskip("polyfempy.polyfempy")
    config_dir = tmp_path / "inputs"
    config_dir.mkdir()
    (config_dir / "square.obj").write_text(
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv .5 .5 0\n"
        "f 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n")
    (config_dir / "volume.txt").write_text("1\n1\n1\n1\n")
    geometry = {"mesh": "square.obj", "volume_selection": "volume.txt"}
    cfg = {
        "geometry": [geometry] if geometry_list else geometry,
        "materials": [{"id": 1, "type": "NeoHookean", "E": 1000, "nu": 0.3}],
        "boundary_conditions": {"dirichlet_boundary": [{"id": "all", "value": [0, 0]}], "rhs": [0, -1]},
        "solver": {"linear": {"solver": "Eigen::SimplicialLDLT"}},
        "output": {"log": {"level": "off", "path": ""}},
    }
    path = config_dir / "simulation.json"
    path.write_text(json.dumps(cfg))
    monkeypatch.chdir(tmp_path)
    result = solve(cfg=path)
    assert np.isfinite(result.sol).all()
    assert abs(result.sol).max() > 1e-6
