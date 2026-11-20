"""Unit tests for JSON normalization helpers in `_solve_contract`.

Covers:
    - clean_json_for_cpp:    recursively drop None leaves from dicts and lists
    - process_json_config:   pop `common`, preserve backend output keys,
                             resolve relative mesh paths via `root_path`

These tests never touch the C++ backend. They only exercise pure Python stages.
"""

from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO / "python") not in sys.path:
    sys.path.insert(0, str(_REPO / "python"))

from polyfempy.runtime._solve_contract import (  # noqa: E402
    clean_json_for_cpp,
    process_json_config,
)


def _canonical_path_text(path: Path) -> str:
    return str(path.resolve())


class CleanJsonForCppTests(unittest.TestCase):
    def test_drops_top_level_none_values(self):
        self.assertEqual(
            clean_json_for_cpp({"a": 1, "b": None, "c": 2}),
            {"a": 1, "c": 2},
        )

    def test_drops_nested_dict_none_leaves(self):
        self.assertEqual(
            clean_json_for_cpp({"outer": {"x": None, "y": 2}}),
            {"outer": {"y": 2}},
        )

    def test_drops_none_items_from_lists(self):
        self.assertEqual(
            clean_json_for_cpp({"items": [None, 1, None, 2]}),
            {"items": [1, 2]},
        )

    def test_keeps_empty_dicts_and_lists(self):
        # clean_json_for_cpp only drops *None* leaves. An empty dict / list
        # is a legitimate value and must be preserved, otherwise solver blocks
        # like "boundary_conditions": {} could vanish.
        self.assertEqual(
            clean_json_for_cpp({"empty_dict": {}, "empty_list": [], "real": None}),
            {"empty_dict": {}, "empty_list": []},
        )

    def test_preserves_zero_and_false_and_empty_string(self):
        self.assertEqual(
            clean_json_for_cpp({"zero": 0, "flag": False, "text": "", "skip": None}),
            {"zero": 0, "flag": False, "text": ""},
        )

    def test_does_not_mutate_input(self):
        original = {"a": 1, "b": None, "nested": {"x": None}}
        snapshot = copy.deepcopy(original)
        clean_json_for_cpp(original)
        self.assertEqual(original, snapshot)

    def test_cleans_list_items_once(self):
        class OneShotDict(dict):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.items_calls = 0

            def items(self):
                self.items_calls += 1
                if self.items_calls > 1:
                    raise AssertionError("list item was cleaned more than once")
                return super().items()

        item = OneShotDict({"keep": 1})

        self.assertEqual(clean_json_for_cpp([item]), [{"keep": 1}])
        self.assertEqual(item.items_calls, 1)

    def test_passes_through_scalars(self):
        for v in (1, 1.5, "text", True, False, 0):
            with self.subTest(value=v):
                self.assertEqual(clean_json_for_cpp(v), v)


class ProcessJsonConfigTests(unittest.TestCase):
    def _base(self) -> dict:
        return {
            "pde": "LinearElasticity",
            "common": "legacy-should-be-dropped",
            "output": {
                "directory": "out",
                "json": True,
            },
            "materials": [{"type": "LinearElasticity", "E": 20, "nu": 0.3}],
        }

    def test_pops_common_key(self):
        processed = process_json_config(self._base())
        self.assertNotIn("common", processed)

    def test_keeps_backend_output_keys(self):
        processed = process_json_config(self._base())
        self.assertIn("output", processed)
        self.assertEqual(processed["output"].get("directory"), "out")
        self.assertTrue(processed["output"].get("json"))

    def test_leaves_output_schema_validation_to_backend(self):
        payload = self._base()
        payload["output"].update(
            {
                "future_backend_option": {"enabled": True},
                "nested_backend_option": {"mode": "raw"},
            }
        )

        processed = process_json_config(payload)

        self.assertEqual(processed["output"], payload["output"])

    def test_does_not_mutate_input(self):
        original = self._base()
        snapshot = copy.deepcopy(original)
        process_json_config(original)
        self.assertEqual(original, snapshot)

    def test_leaves_absolute_mesh_paths_untouched(self):
        d = self._base()
        d["geometry"] = [{"mesh": "/abs/path/to/mesh.msh"}]
        d["root_path"] = "/tmp/cfg.json"
        processed = process_json_config(d)
        self.assertEqual(processed["geometry"][0]["mesh"], "/abs/path/to/mesh.msh")

    def test_resolves_relative_mesh_against_root_path_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg_dir = Path(tmp)
            mesh_file = cfg_dir / "beam.msh"
            mesh_file.write_text("placeholder", encoding="utf-8")

            d = self._base()
            d["geometry"] = [{"mesh": "beam.msh"}]
            d["root_path"] = str(cfg_dir / "config.json")

            processed = process_json_config(d)
            self.assertEqual(
                processed["geometry"][0]["mesh"],
                _canonical_path_text(mesh_file),
            )

    def test_resolves_relative_mesh_against_sibling_meshes_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case_dir = root / "case"
            case_dir.mkdir()
            meshes_dir = root / "meshes"
            meshes_dir.mkdir()
            mesh_file = meshes_dir / "beam.msh"
            mesh_file.write_text("placeholder", encoding="utf-8")

            d = self._base()
            d["geometry"] = [{"mesh": "beam.msh"}]
            d["root_path"] = str(case_dir / "config.json")

            processed = process_json_config(d)
            # Not under case_dir directly, so falls back to ../meshes/beam.msh.
            self.assertEqual(
                processed["geometry"][0]["mesh"],
                _canonical_path_text(mesh_file),
            )

    def test_leaves_unresolvable_relative_mesh_alone(self):
        d = self._base()
        d["geometry"] = [{"mesh": "does_not_exist.msh"}]
        d["root_path"] = "/tmp/nonexistent_dir/config.json"
        processed = process_json_config(d)
        # When the file cannot be located, the original mesh string is kept
        # (the C++ solver will raise its own error later if that's wrong).
        self.assertEqual(processed["geometry"][0]["mesh"], "does_not_exist.msh")

    def test_handles_geometry_entries_that_are_not_dicts(self):
        d = self._base()
        # Non-dict entries should be ignored gracefully rather than crashing.
        d["geometry"] = ["not-a-dict", {"mesh": "beam.msh"}]
        d["root_path"] = "/tmp/cfg.json"
        processed = process_json_config(d)
        self.assertEqual(processed["geometry"][0], "not-a-dict")


if __name__ == "__main__":
    unittest.main()
