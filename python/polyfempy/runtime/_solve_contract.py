"""Shared solve contract helpers.

This internal module owns the parts of ``solve`` that define user-facing
semantics before the C++ backend is touched: config normalization, mesh-source
selection, and backend settings construction.
"""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, overload

import numpy as np

from . import tensor as T


class NoMeshSourceError(ValueError):
    """Raised when a solve call has no usable JSON or array mesh source."""


@dataclass
class MeshSource:
    """Resolved mesh source for a solve call."""

    mode: str
    vertices: Optional[np.ndarray] = None
    cells: Optional[np.ndarray] = None
    body_ids: Optional[np.ndarray] = None
    boundary_ids: Optional[np.ndarray] = None
    v_backend: str = "numpy"


@dataclass
class CanonicalSolveInput:
    """Normalized config, mesh, and backend settings for a solve pipeline."""

    config: Any
    full_json: Optional[Dict[str, Any]]
    mesh_source: MeshSource
    backend_settings: Dict[str, Any]
    metadata: Dict[str, Any]


def is_generated_config(obj: Any) -> bool:
    """Return True for objects that look like generated config objects."""
    return callable(getattr(obj, "as_dict", None))


def generated_payload_from_config(cfg: Any) -> Dict[str, Any]:
    """Return the plain dict emitted by a generated config object."""
    as_dict = getattr(cfg, "as_dict", None)
    if not callable(as_dict):
        raise TypeError(f"cfg must expose as_dict(), got {type(cfg).__name__}")

    payload = as_dict()
    if not isinstance(payload, dict):
        raise TypeError(
            f"cfg.as_dict() must return dict, got {type(payload).__name__}"
        )
    return payload


def prepare_generated_backend_payload(
    payload: Dict[str, Any],
    *,
    root_path: str | None = None,
) -> Dict[str, Any]:
    """Prepare a generated ``Root.as_dict()`` payload for the backend."""
    backend = copy.deepcopy(payload)
    resolved_root_path = root_path or backend.get("root_path")
    if resolved_root_path:
        _resolve_geometry_mesh_paths_in_place(backend, str(resolved_root_path))
        _resolve_collision_mesh_paths_in_place(backend, str(resolved_root_path))
    _restore_backend_payload_semantics(backend)
    return _drop_none(backend)


_SOLVER_ADVANCED_GENERATED_DEFAULTS = {
    "cache_size": 900000,
    "lump_mass_matrix": False,
    "lagged_regularization_weight": 0.0,
    "lagged_regularization_iterations": 1,
    "check_inversion": "Discrete",
    "jacobian_threshold": 0.0,
    "characteristic_length": -1.0,
    "characteristic_force_density": 10000.0,
}

_NONLINEAR_GENERATED_DEFAULTS = {
    "rel_grad_norm_tol": 1e-10,
    "newton_decrement_tol": 0.0,
    "rel_x_delta_tol": 0.0,
    "first_grad_norm_tol": 1e-12,
    "norm_type": "L2",
    "max_iterations": 500,
    "allow_out_of_iterations": False,
}

_LINE_SEARCH_GENERATED_DEFAULTS = {
    "use_grad_norm_tol": 1e-06,
    "min_step_size": 1e-10,
    "max_step_size_iter": 30,
    "min_step_size_final": 1e-20,
    "max_step_size_iter_final": 100,
    "default_init_step_size": 1.0,
    "step_ratio": 0.5,
}

_NONLINEAR_ADVANCED_GENERATED_DEFAULTS = {
    "f_delta": 0.0,
    "f_delta_step_tol": 100,
    "derivative_along_delta_x_tol": 0,
    "apply_gradient_fd": "None",
    "gradient_fd_eps": 1e-07,
}


def _restore_backend_payload_semantics(payload: Dict[str, Any]) -> None:
    _restore_materials_semantics(payload)
    _restore_transformation_semantics(payload)
    _restore_nonlinear_solver_semantics(payload)
    _drop_solver_advanced_generated_defaults(payload)


def _restore_materials_semantics(payload: Dict[str, Any]) -> None:
    materials = payload.get("materials")
    if not isinstance(materials, list) or len(materials) != 1:
        return

    material = materials[0]
    if isinstance(material, dict) and "id" not in material:
        payload["materials"] = material


def _restore_transformation_semantics(payload: Dict[str, Any]) -> None:
    geometry = payload.get("geometry")
    if not isinstance(geometry, list):
        return

    for entry in geometry:
        if not isinstance(entry, dict):
            continue

        transformation = entry.get("transformation")
        if not isinstance(transformation, dict):
            continue

        for scalar_or_vector_key in ("scale", "rotation"):
            value = transformation.get(scalar_or_vector_key)
            if isinstance(value, list) and len(value) == 1:
                transformation[scalar_or_vector_key] = value[0]

        for empty_list_key in ("translation", "rotation", "scale"):
            if transformation.get(empty_list_key) == []:
                transformation.pop(empty_list_key)

        if transformation.get("rotation_mode") == "xyz":
            transformation.pop("rotation_mode")

        if not transformation:
            entry.pop("transformation")


def _drop_solver_advanced_generated_defaults(payload: Dict[str, Any]) -> None:
    solver = payload.get("solver")
    if not isinstance(solver, dict):
        return

    advanced = solver.get("advanced")
    if not isinstance(advanced, dict):
        return

    for key, default_value in _SOLVER_ADVANCED_GENERATED_DEFAULTS.items():
        if advanced.get(key) == default_value:
            advanced.pop(key)

    if not advanced:
        solver.pop("advanced")


def _restore_nonlinear_solver_semantics(payload: Dict[str, Any]) -> None:
    solver = payload.get("solver")
    if not isinstance(solver, dict):
        return

    nonlinear = solver.get("nonlinear")
    if not isinstance(nonlinear, dict):
        return

    for key, default_value in _NONLINEAR_GENERATED_DEFAULTS.items():
        if nonlinear.get(key) == default_value:
            nonlinear.pop(key)

    line_search = nonlinear.get("line_search")
    if isinstance(line_search, dict):
        for key, default_value in _LINE_SEARCH_GENERATED_DEFAULTS.items():
            if line_search.get(key) == default_value:
                line_search.pop(key)
        if not line_search:
            nonlinear.pop("line_search")

    advanced = nonlinear.get("advanced")
    if isinstance(advanced, dict):
        _rename_key(advanced, "f_delta_tol", "f_delta")
        for key, default_value in _NONLINEAR_ADVANCED_GENERATED_DEFAULTS.items():
            if advanced.get(key) == default_value:
                advanced.pop(key)
        if not advanced:
            nonlinear.pop("advanced")


def _rename_key(payload: Dict[str, Any], old: str, new: str) -> None:
    if old in payload and new not in payload:
        payload[new] = payload.pop(old)
    else:
        payload.pop(old, None)


def _drop_none(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {
            key: cleaned
            for key, value in obj.items()
            if (cleaned := _drop_none(value)) is not None
        }
    if isinstance(obj, list):
        return [cleaned for value in obj if (cleaned := _drop_none(value)) is not None]
    return obj


def _resolve_geometry_mesh_paths_in_place(payload: Dict[str, Any], root_path: str) -> None:
    geometry = payload.get("geometry")
    if not isinstance(geometry, list):
        return

    root = Path(root_path).resolve()
    root_dir = root if root.is_dir() else root.parent
    for entry in geometry:
        if not isinstance(entry, dict):
            continue

        _resolve_geometry_selection_files_in_place(entry, root_path)

        mesh_path = entry.get("mesh")
        if not isinstance(mesh_path, str) or not mesh_path.strip():
            continue

        mesh_file = Path(mesh_path)
        if mesh_file.is_absolute():
            continue

        direct = root_dir / mesh_file
        if direct.exists():
            entry["mesh"] = str(direct)
            continue

        sibling = root_dir.parent / "meshes" / mesh_file.name
        if sibling.exists():
            entry["mesh"] = str(sibling)


def _resolve_geometry_selection_files_in_place(entry: Dict[str, Any], root_path: str) -> None:
    for key in ("volume_selection", "surface_selection", "point_selection"):
        selection = entry.get(key)
        if isinstance(selection, str) and selection.strip():
            entry[key] = _resolve_existing_path_from_root(root_path, selection)
        else:
            _resolve_selection_file_paths_in_place(selection, root_path)


def _resolve_selection_file_paths_in_place(selection: Any, root_path: str) -> None:
    if isinstance(selection, list):
        for item in selection:
            _resolve_selection_file_paths_in_place(item, root_path)
        return

    if not isinstance(selection, dict):
        return

    file_path = selection.get("file")
    if isinstance(file_path, str) and file_path.strip():
        selection["file"] = _resolve_existing_path_from_root(root_path, file_path)


def _resolve_existing_path_from_root(root_path: str, path_value: str) -> str:
    path = Path(path_value)
    if path.is_absolute():
        return path_value

    root = Path(root_path).resolve()
    root_dir = root if root.is_dir() else root.parent
    resolved = root_dir / path
    if resolved.exists():
        return str(resolved)
    return path_value


def _resolve_collision_mesh_paths_in_place(payload: Dict[str, Any], root_path: str) -> None:
    contact = payload.get("contact")
    if not isinstance(contact, dict):
        return

    collision_mesh = contact.get("collision_mesh")
    if not isinstance(collision_mesh, dict):
        return

    for key in ("mesh", "linear_map"):
        value = collision_mesh.get(key)
        if isinstance(value, str) and value.strip():
            collision_mesh[key] = _resolve_existing_path_from_root(root_path, value)


def cfg_array_mesh_payload(cfg: Any) -> Optional[Dict[str, Any]]:
    if isinstance(cfg, dict):
        payload = cfg.get("_mesh_array_mode")
        if isinstance(payload, dict):
            return payload
        extras = cfg.get("extras")
        if isinstance(extras, dict):
            payload = extras.get("_mesh_array_mode")
            if isinstance(payload, dict):
                return payload
        return None

    extras = getattr(cfg, "extras", None)
    if isinstance(extras, dict):
        payload = extras.get("_mesh_array_mode")
        if isinstance(payload, dict):
            return payload
    return None


def _load_json_payload(path: str | Path) -> Dict[str, Any]:
    json_path = Path(path).resolve()
    try:
        with json_path.open("r", encoding="utf-8") as f:
            payload = json.load(f)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Configuration file not found: {json_path}. "
            "Check the path or pass cfg=model.config()."
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in {json_path}, line {exc.lineno}, column {exc.colno}: {exc.msg}"
        ) from exc
    if not isinstance(payload, dict):
        raise TypeError(f"JSON config must contain an object, got {type(payload).__name__}")
    payload.setdefault("root_path", str(json_path))
    return payload


def _payload_dict_or_raise(cfg: Any, *, context: str) -> Dict[str, Any]:
    if isinstance(cfg, dict):
        return copy.deepcopy(cfg)

    for method_name in ("to_json_dict", "to_dict"):
        method = getattr(cfg, method_name, None)
        if not callable(method):
            continue
        try:
            payload = method()
        except Exception as exc:
            raise RuntimeError(f"cfg.{method_name}() failed while {context}") from exc
        if not isinstance(payload, dict):
            raise TypeError(
                f"cfg.{method_name}() must return dict while {context}, "
                f"got {type(payload).__name__}"
            )
        return copy.deepcopy(payload)

    raise TypeError(
        "cfg must be dict, str/Path JSON file path, generated object with as_dict(), "
        f"or expose to_dict()/to_json_dict(); got {type(cfg).__name__}"
    )


def normalize_config(cfg: Any) -> Dict[str, Any]:
    """Accept a backend-shaped dict, JSON path, or generated-style config object."""
    if cfg is None:
        raise ValueError("cfg (configuration) is required")
    if isinstance(cfg, dict):
        return copy.deepcopy(cfg)
    if isinstance(cfg, (str, Path)):
        return _load_json_payload(cfg)
    if is_generated_config(cfg):
        return copy.deepcopy(generated_payload_from_config(cfg))
    return _payload_dict_or_raise(cfg, context="normalizing config")


def validate_config(cfg: Any) -> None:
    if hasattr(cfg, "validate"):
        cfg.validate()


def _ensure_i32(cells: np.ndarray) -> np.ndarray:
    return cells.astype(np.int32, copy=False) if cells.dtype != np.int32 else cells


def _validate_mesh_arrays(vertices: np.ndarray, cells: np.ndarray) -> None:
    if vertices.ndim != 2 or vertices.shape[1] not in (2, 3) or not len(vertices):
        raise ValueError("vertices must be a non-empty array with shape (n, 2) or (n, 3).")
    if vertices.dtype.kind not in "iuf" or not np.isfinite(vertices).all():
        raise ValueError("vertices must contain finite real coordinates; remove NaN and infinity.")
    if cells.ndim != 2 or not cells.size:
        raise ValueError("cells must be a non-empty 2D array, with one element per row.")
    if cells.dtype.kind not in "iuf" or not np.isfinite(cells).all():
        raise ValueError("cells must contain finite integer vertex indices.")
    if cells.dtype.kind == "f" and not np.equal(cells, np.floor(cells)).all():
        raise ValueError("cells must contain integer vertex indices; fractional indices are invalid.")
    if cells.min() < 0 or cells.max() >= len(vertices):
        raise ValueError(
            f"cells contains a vertex index outside 0..{len(vertices) - 1}. "
            "Use zero-based indices into vertices."
        )
    if cells.max() > np.iinfo(np.int32).max:
        raise ValueError("cells indices exceed the int32 range supported by the backend.")


def _promote_materials_to_list(payload: Dict[str, Any], *, infer_type_from_pde: bool = False) -> None:
    materials = payload.get("materials")
    if not isinstance(materials, dict):
        return
    if infer_type_from_pde and "type" not in materials:
        pde = payload.get("pde", "LinearElasticity")
        materials["type"] = "Laplacian" if pde == "Poisson" else "LinearElasticity"
    payload["materials"] = [materials]


def merge_user_cfg_over_full_json(cfg: Any, full_json: Dict[str, Any]) -> Dict[str, Any]:
    cfg_dict = _payload_dict_or_raise(cfg, context="merging Python-side config over full JSON")

    merged = copy.deepcopy(full_json) if isinstance(full_json, dict) else {}
    extras = cfg_dict.pop("extras", None)
    for key, value in cfg_dict.items():
        merged[key] = value

    if isinstance(extras, dict) and "_root_path" in extras:
        merged["root_path"] = extras["_root_path"]
    elif isinstance(full_json, dict) and "root_path" in full_json:
        merged["root_path"] = full_json["root_path"]

    return merged


def build_full_json(cfg: Any) -> Optional[Dict[str, Any]]:
    validate_config(cfg)

    if cfg_array_mesh_payload(cfg) is not None:
        return None

    cfg_dict = _payload_dict_or_raise(cfg, context="building full JSON settings")
    extras = cfg_dict.get("extras")
    if isinstance(extras, dict) and "_full_json_config" in extras:
        return merge_user_cfg_over_full_json(cfg_dict, extras["_full_json_config"])

    if not (isinstance(cfg_dict, dict) and "geometry" in cfg_dict):
        return None

    full_json = cfg_dict
    full_json.pop("extras", None)
    full_json.pop("_mesh_array_mode", None)
    return full_json


def choose_mesh_source(
    vertices: Any,
    cells: Any,
    full_json: Optional[Dict[str, Any]],
    *,
    dtype: Any = None,
    cfg: Any = None,
) -> MeshSource:
    """Resolve the real mesh source for this solve.

    Explicit partial array input is an error even when JSON geometry exists,
    because silently falling back would change the user's requested data binding.
    """
    if (vertices is None) ^ (cells is None):
        raise ValueError("array mode requires both vertices and cells")

    array_payload = cfg_array_mesh_payload(cfg)
    used_cfg_array_payload = False
    if vertices is None and cells is None and array_payload is not None:
        vertices = array_payload.get("vertices")
        cells = array_payload.get("cells")
        used_cfg_array_payload = True
        if vertices is None or cells is None:
            raise ValueError("guided array mesh payload requires both vertices and cells")

    if vertices is not None and cells is not None:
        V_np, v_backend = T.as_numpy(vertices, dtype=dtype)
        C_np, _ = T.as_numpy(cells)
        _validate_mesh_arrays(V_np, C_np)
        body_ids_np = None
        boundary_ids_np = None
        if used_cfg_array_payload and array_payload is not None:
            if array_payload.get("body_ids") is not None:
                body_ids_np = np.asarray(array_payload["body_ids"], dtype=np.int32).reshape(-1)
            if array_payload.get("boundary_ids") is not None:
                boundary_ids_np = np.asarray(array_payload["boundary_ids"], dtype=np.int32).reshape(-1)
        return MeshSource(
            mode="guided_array" if used_cfg_array_payload else "array",
            vertices=V_np,
            cells=_ensure_i32(C_np),
            body_ids=body_ids_np,
            boundary_ids=boundary_ids_np,
            v_backend=v_backend,
        )

    if full_json is not None and "geometry" in full_json:
        return MeshSource(mode="json", v_backend="numpy")

    raise NoMeshSourceError(
        "Either provide vertices/cells arrays, or use JSON config with geometry (mesh files)"
    )


@overload
def clean_json_for_cpp(obj: Dict[str, Any]) -> Dict[str, Any]: ...


@overload
def clean_json_for_cpp(obj: Any) -> Any: ...


def clean_json_for_cpp(obj: Any) -> Any:
    if isinstance(obj, dict):
        cleaned = {}
        for key, value in obj.items():
            cleaned_value = clean_json_for_cpp(value)
            if cleaned_value is not None:
                cleaned[key] = cleaned_value
        return cleaned
    if isinstance(obj, list):
        cleaned_items = []
        for item in obj:
            cleaned_item = clean_json_for_cpp(item)
            if cleaned_item is not None:
                cleaned_items.append(cleaned_item)
        return cleaned_items
    return obj


def process_json_config(full_json: Dict[str, Any]) -> Dict[str, Any]:
    processed = copy.deepcopy(full_json)
    processed.pop("common", None)
    processed.pop("extras", None)
    processed.pop("_mesh_array_mode", None)

    root_path = processed.get("root_path")
    if root_path:
        _resolve_geometry_mesh_paths_in_place(processed, str(root_path))
    return processed


def _strip_guided_array_placeholders(settings: Dict[str, Any]) -> None:
    geometry = settings.get("geometry")
    if not isinstance(geometry, list):
        return
    kept = []
    for entry in geometry:
        if isinstance(entry, dict) and str(entry.get("mesh", "")).startswith("__array_body__:"):
            continue
        kept.append(entry)
    if kept:
        settings["geometry"] = kept
    else:
        settings.pop("geometry", None)


def build_canonical_solver_settings(
    cfg: Any,
    *,
    full_json: Optional[Dict[str, Any]],
    mesh_source: MeshSource,
) -> Dict[str, Any]:
    """Build backend JSON/settings after config and mesh source are known."""
    if mesh_source.mode == "json":
        if full_json is None:
            raise ValueError("JSON mesh mode requires full_json settings")
        return clean_json_for_cpp(process_json_config(full_json))

    settings = _payload_dict_or_raise(cfg, context="building array-mode backend settings")

    if mesh_source.mode == "guided_array":
        _strip_guided_array_placeholders(settings)

    settings.pop("extras", None)
    settings.pop("_mesh_array_mode", None)

    if "geometry" not in settings:
        settings["geometry"] = [
            {"type": "ground", "height": 0.0, "enabled": True, "is_obstacle": False}
        ]

    _promote_materials_to_list(settings, infer_type_from_pde=True)
    return clean_json_for_cpp(settings)


def _prepare_generated_canonical_solve_input(
    *,
    vertices: Any,
    cells: Any,
    cfg: Any,
    dtype: Any = None,
) -> CanonicalSolveInput:
    payload = generated_payload_from_config(cfg)
    mesh_source = choose_mesh_source(
        vertices,
        cells,
        payload,
        dtype=dtype,
        cfg=cfg,
    )
    backend_settings = prepare_generated_backend_payload(payload)
    return CanonicalSolveInput(
        config=cfg,
        full_json=payload,
        mesh_source=mesh_source,
        backend_settings=backend_settings,
        metadata={
            "mesh_source": mesh_source.mode,
            "config_source": "generated",
        },
    )


def prepare_canonical_solve_input(
    *,
    vertices: Any,
    cells: Any,
    cfg: Any,
    dtype: Any = None,
) -> CanonicalSolveInput:
    """Normalize config, resolve mesh source, and build backend settings once."""
    if is_generated_config(cfg):
        return _prepare_generated_canonical_solve_input(
            vertices=vertices,
            cells=cells,
            cfg=cfg,
            dtype=dtype,
        )

    config = normalize_config(cfg)
    full_json = build_full_json(config)
    mesh_source = choose_mesh_source(
        vertices,
        cells,
        full_json,
        dtype=dtype,
        cfg=config,
    )
    backend_settings = build_canonical_solver_settings(
        config,
        full_json=full_json,
        mesh_source=mesh_source,
    )
    return CanonicalSolveInput(
        config=config,
        full_json=full_json,
        mesh_source=mesh_source,
        backend_settings=backend_settings,
        metadata={"mesh_source": mesh_source.mode},
    )
