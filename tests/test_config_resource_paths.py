"""Generated and JSON configurations resolve the same resource paths."""

import copy

import pytest

from polyfempy.runtime._solve_contract import (
    prepare_generated_backend_payload,
    process_json_config,
)


@pytest.mark.parametrize("prepare", [process_json_config, prepare_generated_backend_payload])
@pytest.mark.parametrize("geometry_list", [False, True])
def test_resource_paths_are_relative_to_config(tmp_path, prepare, geometry_list):
    resources = ["beam.msh", "volume.txt", "surface.txt", "point.txt", "proxy.ply", "map.hdf5"]
    for name in resources:
        (tmp_path / name).touch()
    geometry = {
        "mesh": "beam.msh",
        "volume_selection": "volume.txt",
        "surface_selection": {"file": "surface.txt"},
        "point_selection": [{"file": "point.txt"}],
    }
    payload = {
        "root_path": str(tmp_path / "config.json"),
        "geometry": [geometry] if geometry_list else geometry,
        "contact": {"collision_mesh": {"mesh": "proxy.ply", "linear_map": "map.hdf5"}},
    }
    original = copy.deepcopy(payload)

    prepared = prepare(payload)

    assert payload == original
    resolved = prepared["geometry"][0] if geometry_list else prepared["geometry"]
    assert resolved["mesh"] == str(tmp_path / "beam.msh")
    assert resolved["volume_selection"] == str(tmp_path / "volume.txt")
    assert resolved["surface_selection"]["file"] == str(tmp_path / "surface.txt")
    assert resolved["point_selection"][0]["file"] == str(tmp_path / "point.txt")
    assert prepared["contact"]["collision_mesh"] == {
        "mesh": str(tmp_path / "proxy.ply"),
        "linear_map": str(tmp_path / "map.hdf5"),
    }
