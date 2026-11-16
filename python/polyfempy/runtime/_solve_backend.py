"""Backend adapter helpers for the VarForm forward solve pipeline.

This module owns only the mechanics of talking to the compiled PolyFEM
VarForm binding: constructing the backend solver, applying settings, attaching
array meshes, and calling ``solve``. User-facing config and mesh-source
semantics live in ``_solve_contract.py``; the final backend bundle is converted
in ``solve.py``.
"""

from __future__ import annotations

import copy
import importlib
import json
from typing import Any, Callable, Dict, Optional

from ._solve_contract import MeshSource, build_canonical_solver_settings


def build_solver():
    """Instantiate the compiled VarForm Solver binding."""
    try:
        import polyfempy as pf
    except Exception as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "polyfempy bindings not found. Please install/compile them first."
        ) from exc

    if not getattr(pf, "cpp_backend_available", lambda: False)():
        err = getattr(pf, "cpp_backend_error", lambda: None)()
        raise RuntimeError(
            "C++ backend not loaded. Install a wheel for your Python version "
            "and platform, or follow BUILD.md to build from source. "
            f"Import error: {err}"
        )

    core = importlib.import_module("polyfempy.polyfempy")
    if not hasattr(core, "Solver"):
        raise RuntimeError("Compiled VarForm backend must expose Solver.")

    solver = core.Solver()
    missing = [
        name
        for name in ("set_settings", "load_mesh_from_settings", "set_mesh", "solve")
        if not hasattr(solver, name)
    ]
    if missing:
        raise RuntimeError(
            "Compiled VarForm Solver is missing required methods: "
            + ", ".join(missing)
        )
    return solver


def _apply_settings_json(solver, settings_json: str) -> None:
    if not hasattr(solver, "set_settings"):
        raise RuntimeError("Compiled VarForm Solver is missing set_settings(...).")
    try:
        solver.set_settings(settings_json, strict_validation=True)
    except RuntimeError as exc:
        # Add a focused hint only after the backend rejects the settings.
        materials = json.loads(settings_json).get("materials", [])
        materials = [materials] if isinstance(materials, dict) else materials
        if isinstance(materials, list):
            for index, material in enumerate(materials):
                if not isinstance(material, dict):
                    continue
                if material.get("type") not in ("LinearElasticity", "NeoHookean"):
                    continue
                if not ({"E", "nu"} <= material.keys() or {"lambda", "mu"} <= material.keys()):
                    raise RuntimeError(
                        f"{exc}. materials[{index}] ({material['type']}) requires "
                        "both E and nu, or both lambda and mu. "
                        "For example, set E=1000 and nu=0.3."
                    ) from exc
        raise


def _configure_json_mode(
    solver,
    full_json: Optional[Dict[str, Any]],
    cfg,
    *,
    backend_settings: Optional[Dict[str, Any]] = None,
) -> None:
    settings_dict = (
        copy.deepcopy(backend_settings)
        if backend_settings is not None
        else build_canonical_solver_settings(
            cfg,
            full_json=full_json,
            mesh_source=MeshSource(mode="json"),
        )
    )
    _apply_settings_json(solver, json.dumps(settings_dict))
    solver.load_mesh_from_settings()


def _configure_array_mode(
    solver,
    cfg,
    mesh_source: MeshSource,
    *,
    backend_settings: Optional[Dict[str, Any]] = None,
) -> None:
    settings_dict = (
        copy.deepcopy(backend_settings)
        if backend_settings is not None
        else build_canonical_solver_settings(
            cfg,
            full_json=None,
            mesh_source=mesh_source,
        )
    )
    _apply_settings_json(solver, json.dumps(settings_dict))
    _apply_mesh_array_mode(solver, mesh_source)


def _apply_mesh_array_mode(solver, mesh_source: MeshSource) -> None:
    if mesh_source.body_ids is not None:
        raise NotImplementedError(
            "VarForm array mode does not support Python-side body_ids. "
            "Put body tags in the mesh/config before calling solve()."
        )
    if mesh_source.boundary_ids is not None:
        raise NotImplementedError(
            "VarForm array mode does not support Python-side boundary_ids. "
            "Put boundary tags in the mesh/config before calling solve()."
        )
    if not hasattr(solver, "set_mesh"):
        raise RuntimeError("Compiled VarForm Solver is missing set_mesh(...).")

    solver.set_mesh(mesh_source.vertices, mesh_source.cells)


def configure_solver(
    solver,
    cfg,
    full_json: Optional[dict],
    mesh_source: MeshSource,
    *,
    backend_settings: Optional[Dict[str, Any]] = None,
) -> None:
    """Set settings and load or attach the mesh."""
    if mesh_source.mode == "json":
        _configure_json_mode(
            solver,
            full_json,
            cfg,
            backend_settings=backend_settings,
        )
        return

    _configure_array_mode(
        solver,
        cfg,
        mesh_source,
        backend_settings=backend_settings,
    )


def apply_sidesets(solver, sidesets_func: Optional[Callable]) -> None:
    """Reject the old Python-side sideset mutation hook.

    New VarForm solves should receive boundary information through settings and
    mesh data before ``load_mesh``/``set_mesh``. The old hook depended on
    mutating the backend mesh object after loading.
    """
    if sidesets_func is None:
        return
    raise NotImplementedError(
        "sidesets_func is not supported by the VarForm runtime. "
        "Encode boundary ids in the mesh/config before calling solve()."
    )


def _resolve_log_level(full_json: Optional[dict]) -> int:
    if not (full_json and "output" in full_json and "log" in full_json["output"]):
        return 2
    log_cfg = full_json["output"]["log"]
    # Keep terminal logging aligned with the generated OutputLog default.
    # ``OutputLog.to_dict()`` omits ``level`` when it is the default "debug",
    # so falling back to "info" here would silently make console output less
    # verbose than the file log.
    log_level = log_cfg.get("level", "debug")
    if isinstance(log_level, int) and not isinstance(log_level, bool):
        return max(0, min(6, log_level))

    if isinstance(log_level, str) and log_level.strip().isdigit():
        return max(0, min(6, int(log_level.strip())))

    log_level_map = {
        "trace": 0,
        "debug": 1,
        "info": 2,
        "warn": 3,
        "warning": 3,
        "error": 4,
        "critical": 5,
        "off": 6,
    }
    if isinstance(log_level, str):
        return log_level_map.get(log_level.lower(), 2)
    return 2


def run_solver_stage(solver, full_json: Optional[dict]):
    """Run the VarForm backend solve entry point and return its raw result."""
    if not hasattr(solver, "solve"):
        raise RuntimeError("Compiled VarForm backend requires solve entry point.")

    return solver.solve(log_level=_resolve_log_level(full_json))


__all__ = [
    "apply_sidesets",
    "build_solver",
    "configure_solver",
    "run_solver_stage",
]
