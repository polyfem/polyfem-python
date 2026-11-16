"""Unit tests for small solve-contract helpers.

Covers edges that the other pipeline tests don't exercise directly:
    - _promote_materials_to_list: dict->list promotion, optional type-inference

These tests never touch the C++ backend. They only exercise pure Python helpers.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO / "python") not in sys.path:
    sys.path.insert(0, str(_REPO / "python"))

from polyfempy.runtime._solve_contract import (  # noqa: E402
    _promote_materials_to_list,
)


class PromoteMaterialsToListTests(unittest.TestCase):
    def test_promotes_dict_to_singleton_list(self):
        payload = {"materials": {"E": 20, "nu": 0.3}}
        _promote_materials_to_list(payload)
        self.assertEqual(payload["materials"], [{"E": 20, "nu": 0.3}])

    def test_leaves_list_materials_untouched(self):
        payload = {"materials": [{"E": 1}, {"E": 2}]}
        original = payload["materials"]
        _promote_materials_to_list(payload)
        self.assertIs(payload["materials"], original)

    def test_noop_when_materials_missing(self):
        payload = {"pde": "LinearElasticity"}
        _promote_materials_to_list(payload)
        self.assertNotIn("materials", payload)

    def test_infer_type_fills_missing_type_for_linear_elasticity(self):
        payload = {"pde": "LinearElasticity", "materials": {"E": 20, "nu": 0.3}}
        _promote_materials_to_list(payload, infer_type_from_pde=True)
        self.assertEqual(payload["materials"][0]["type"], "LinearElasticity")

    def test_infer_type_uses_laplacian_for_poisson(self):
        payload = {"pde": "Poisson", "materials": {"k": 1.0}}
        _promote_materials_to_list(payload, infer_type_from_pde=True)
        self.assertEqual(payload["materials"][0]["type"], "Laplacian")

    def test_infer_type_does_not_override_existing_type(self):
        payload = {
            "pde": "Poisson",
            "materials": {"type": "NeoHookean", "E": 20, "nu": 0.3},
        }
        _promote_materials_to_list(payload, infer_type_from_pde=True)
        self.assertEqual(payload["materials"][0]["type"], "NeoHookean")

    def test_without_infer_type_flag_no_type_is_added(self):
        payload = {"pde": "LinearElasticity", "materials": {"E": 20, "nu": 0.3}}
        _promote_materials_to_list(payload, infer_type_from_pde=False)
        self.assertNotIn("type", payload["materials"][0])


if __name__ == "__main__":
    unittest.main()
