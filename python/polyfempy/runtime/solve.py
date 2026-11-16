"""Public solve entry point and internal orchestration."""

from __future__ import annotations

from typing import Any, Callable, Optional

import numpy as np
from numpy.typing import ArrayLike, DTypeLike

from . import _solve_backend as _backend
from . import _solve_contract as _contract
from .result import Result

__all__ = ["solve"]


def _result_from_backend_return(ret) -> Result:
    """Convert the standard C++ backend result bundle into a Result."""
    if not isinstance(ret, dict) or not ret.get("_result_bundle"):
        raise RuntimeError(
            "PolyFEM backend did not return the standard result bundle. "
            "Expected a dict with _result_bundle and sol."
        )

    if "p" in ret:
        raise RuntimeError(
            "PolyFEM backend returned legacy pressure field 'p'. "
            "Rebuild the VarForm backend so Solver.solve() returns only sol."
        )
    if "u" in ret:
        raise RuntimeError(
            "PolyFEM backend returned legacy u field. "
            "Rebuild the VarForm backend so Solver.solve() returns raw sol."
        )

    required = ("sol",)
    missing = [key for key in required if key not in ret]
    if missing:
        raise RuntimeError(
            "PolyFEM backend result bundle is missing required keys: "
            + ", ".join(missing)
        )

    meta = {}
    if isinstance(ret.get("meta"), dict):
        meta.update(ret["meta"])

    return Result(
        np.asarray(ret["sol"]),
        meta=meta,
    )


def run_pipeline(
    vertices=None,
    cells=None,
    cfg=None,
    sidesets_func: Optional[Callable] = None,
    dtype=None,
) -> Result:
    """Drive the full solve pipeline."""
    _backend.apply_sidesets(None, sidesets_func)
    canonical = _contract.prepare_canonical_solve_input(
        vertices=vertices,
        cells=cells,
        cfg=cfg,
        dtype=dtype,
    )
    solver = _backend.build_solver()
    _backend.configure_solver(
        solver,
        canonical.config,
        canonical.full_json,
        canonical.mesh_source,
        backend_settings=canonical.backend_settings,
    )
    ret = _backend.run_solver_stage(solver, canonical.backend_settings)
    return _result_from_backend_return(ret)


def solve(
    vertices: ArrayLike | None = None,
    cells: ArrayLike | None = None,
    cfg: Any = None,
    sidesets_func: Optional[Callable] = None,
    dtype: DTypeLike | None = None,
) -> Result:
    """Run a forward simulation with the compiled PolyFEM backend.

    Args:
        vertices: Optional vertex coordinates. Supply both vertices and cells
            to use an array mesh instead of mesh files in the configuration.
        cells: Optional element connectivity with zero-based vertex indices.
            Connectivity is converted to int32 before calling the backend.
        cfg: Required simulation configuration. Use model.config(), a dict,
            or a JSON file path. Objects with as_dict(), to_dict(), or
            to_json_dict() are also accepted. A JSON file resolves relative
            mesh paths from its location.
        sidesets_func: Reserved for compatibility. Only None is supported.
            Set boundary conditions on the model or in the mesh/configuration.
        dtype: Optional NumPy dtype for input vertex coordinates. This does
            not select the precision of the native solver.

    Returns:
        A Result with the raw finite-element solution in result.sol and
        backend metadata in result.meta. Solution rows are not necessarily
        mesh vertices. result.sol is a NumPy array, including for tensor input.

    Raises:
        ValueError: The configuration or a required mesh input is missing.
        TypeError: The configuration has an unsupported type.
        NotImplementedError: A legacy sideset callback or unsupported array
            mesh tags are supplied.
        RuntimeError: The native backend cannot load or the solve fails.

    Example:
        result = solve(cfg=model.config())
    """
    return run_pipeline(
        vertices=vertices,
        cells=cells,
        cfg=cfg,
        sidesets_func=sidesets_func,
        dtype=dtype,
    )
