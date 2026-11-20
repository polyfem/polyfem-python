"""Configuration loading passes dictionaries and generated payloads through."""

import json

import pytest

from polyfempy.solve import _config

CONFIG = {"geometry": [{"mesh": "beam.msh"}],
          "materials": {"type": "LinearElasticity", "E": 20, "nu": 0.3}}


def test_generated_configuration_is_not_rewritten():
    class Config:
        def as_dict(self):
            return CONFIG

    assert _config(Config()) == CONFIG


def test_dictionary_is_not_rewritten():
    assert _config(CONFIG) == CONFIG


def test_json_path_sets_resource_root(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(CONFIG))
    assert _config(path) == dict(CONFIG, root_path=str(path.resolve()))


def test_explicit_resource_root_is_preserved(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(dict(CONFIG, root_path="other-inputs")))
    assert _config(path)["root_path"] == "other-inputs"


@pytest.mark.parametrize("cfg", [None, 123, 1.5, ["a"], object()])
def test_unsupported_configurations_fail(cfg):
    with pytest.raises(TypeError, match="cfg"):
        _config(cfg)


def test_configuration_errors_propagate():
    class Config:
        def as_dict(self):
            raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        _config(Config())
