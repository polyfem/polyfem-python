"""Tests for the minimal runtime ``Result`` container."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO / "python") not in sys.path:
    sys.path.insert(0, str(_REPO / "python"))

from polyfempy.runtime.result import Result  # noqa: E402


def test_result_stores_raw_backend_solution_and_meta():
    sol = np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]])

    result = Result(sol, meta={"from_bundle": True})

    np.testing.assert_array_equal(result.sol, sol)
    assert result.meta == {"from_bundle": True}
    assert not hasattr(result, "vertices")
    assert not hasattr(result, "cells")
    assert not hasattr(result, "p")
    assert not hasattr(result, "u")


def test_result_keeps_raw_backend_solution_shape():
    result = Result(np.arange(6))

    np.testing.assert_array_equal(result.sol, np.arange(6))


def test_removed_runtime_output_and_convenience_helpers_are_absent():
    result = Result(np.zeros((3, 2)))

    for name in (
        "backend",
        "vertices",
        "cells",
        "V",
        "p",
        "u",
        "fields",
        "point_data",
        "cell_data",
        "history",
        "sampled_data",
        "body_ids",
        "von_mises",
        "stress",
        "strain",
        "has_field",
        "require_field",
        "point_field",
        "cell_field",
        "set_field",
        "remove_field",
        "field_names",
        "field",
        "n_vertices",
        "n_cells",
        "available_fields",
        "summary",
        "as_numpy",
        "to_backend",
        "to_torch",
        "set_sampled_field",
        "sampled_field",
        "field_by_body",
        "to_meshio",
        "from_meshio",
        "read",
        "write",
        "to_vtk",
        "get_percentile_from_von_mises",
    ):
        assert not hasattr(result, name)
        assert not hasattr(Result, name)
