"""Smoke and source-parity tests for classic generated API examples."""

from __future__ import annotations

import copy
import importlib.util
import json
import py_compile
import sys
from pathlib import Path
from typing import Any

import pytest

from polyfempy.runtime._solve_contract import prepare_canonical_solve_input
from polyfempy.generated_api import generated_api as polyfem


ROOT = Path(__file__).resolve().parents[1]
CLASSIC_EXAMPLES = ROOT / "examples" / "classic_example"
POLYFEM_DATA_CONTACT = ROOT / "polyfem-data" / "contact"
POLYFEM_DATA_EXAMPLES = POLYFEM_DATA_CONTACT / "examples"

EXAMPLE_CASES = [
    (
        CLASSIC_EXAMPLES / "2D" / "contact_2d_golf_ball_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "2D" / "golf-ball.json",
    ),
    (
        CLASSIC_EXAMPLES
        / "2D"
        / "contact_2d_golf_ball_deformable_wall_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "2D" / "golf-ball-doformable-wall.json",
    ),
    (
        CLASSIC_EXAMPLES / "2D" / "contact_2d_friction_circle_rollers_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "2D" / "friction" / "circle-rollers.json",
    ),
    (
        CLASSIC_EXAMPLES / "3D" / "contact_3d_golf_ball_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "3D" / "golf-ball.json",
    ),
    (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_friction_high_school_slopetest_generated_api.py",
        POLYFEM_DATA_EXAMPLES
        / "3D"
        / "friction"
        / "high-school-physics-slopetest-mu=0.50.json",
    ),
    (
        CLASSIC_EXAMPLES / "3D" / "contact_3d_large_ratios_sphere_mat_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "3D" / "large-ratios" / "sphere-mat.json",
    ),
]

DOLPHIN_MATERIAL_SHAPE_CASES = [
    (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_stress_tests_dolphin_funnel_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "3D" / "stress-tests" / "dolphin-funnel.json",
    ),
    (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_stress_tests_dolphin_funnel_linf_generated_api.py",
        POLYFEM_DATA_EXAMPLES / "3D" / "stress-tests" / "dolphin-funnel-Linf.json",
    ),
]

EXPECTED_GENERATED_API_FALLBACK_SOURCES: set[str] = set()


def _ordinary_contact_source_paths() -> set[Path]:
    sources: set[Path] = set()
    for dim in ("2D", "3D"):
        for path in (POLYFEM_DATA_EXAMPLES / dim).rglob("*.json"):
            if _is_auxiliary_contact_json(path):
                continue
            sources.add(path.resolve())
    return sources


def _is_auxiliary_contact_json(path: Path) -> bool:
    return (
        path.name == "common.json"
        or path.name.endswith("-common.json")
    )


def _generated_example_files() -> list[Path]:
    return sorted(CLASSIC_EXAMPLES.glob("[23]D/contact_*_generated_api.py"))


def _declared_source_paths() -> set[Path]:
    return {
        Path(_import_example(path).SOURCE_JSON).resolve()
        for path in _generated_example_files()
    }


def _all_example_cases() -> list[tuple[Path, Path]]:
    cases = []
    for example_path in _generated_example_files():
        example = _import_example(example_path)
        cases.append((example_path, Path(example.SOURCE_JSON).resolve()))
    return cases


def _builder_safe_source_paths() -> set[Path]:
    safe = set()
    for source_path in _ordinary_contact_source_paths():
        config = _load_source_config(source_path)
        try:
            polyfem.config(**copy.deepcopy(config))
        except TypeError:
            continue
        if _is_builder_safe_config(config):
            safe.add(source_path)
    return safe


def _is_builder_safe_config(config: dict[str, Any]) -> bool:
    geometry = config.get("geometry")
    if not isinstance(geometry, list):
        return False
    if any(
        not isinstance(item, dict)
        or item.get("type") not in ("mesh", "mesh_array", "mesh_sequence")
        for item in geometry
    ):
        return False

    normal_geometry = [item for item in geometry if not item.get("is_obstacle")]
    if any(
        "volume_selection" in item
        and (
            not isinstance(item.get("volume_selection"), int)
            or item.get("volume_selection") <= 0
        )
        for item in normal_geometry
    ):
        return False
    volume_geometry = [
        item for item in normal_geometry
        if isinstance(item.get("volume_selection"), int)
        and item.get("volume_selection") > 0
    ]
    obstacle_surface_ids = {
        item.get("surface_selection")
        for item in geometry
        if item.get("is_obstacle")
        and isinstance(item.get("surface_selection"), int)
        and item.get("surface_selection") > 0
    }
    if any(
        item.get("is_obstacle")
        and item.get("point_selection") not in (None, 0)
        for item in geometry
    ):
        return False

    volume_ids = {
        item.get("volume_selection")
        for item in volume_geometry
        if isinstance(item.get("volume_selection"), int)
        and item.get("volume_selection") > 0
    }
    if not all(item.get("volume_selection") in volume_ids for item in volume_geometry):
        return False

    materials = _optional_items(config.get("materials"))
    if any(not isinstance(item, dict) for item in materials):
        return False
    id_materials = [item for item in materials if "id" in item]
    if any(item.get("id") not in volume_ids for item in id_materials):
        return False

    selection_ids = _builder_selection_ids(geometry)
    if selection_ids is None:
        return False
    if not _builder_safe_boundary_conditions(
        config.get("boundary_conditions"),
        selection_ids,
        obstacle_surface_ids,
    ):
        return False
    return _builder_safe_initial_conditions(
        config.get("initial_conditions"),
        volume_ids,
    )


def _builder_selection_ids(geometry: list[dict[str, Any]]) -> set[int] | None:
    selection_ids = set()
    for item in geometry:
        for field in ("surface_selection", "point_selection"):
            value = item.get(field)
            if value is None:
                continue
            if isinstance(value, int):
                if value > 0:
                    selection_ids.add(value)
                continue
            if isinstance(value, dict):
                items = [value]
            elif isinstance(value, list):
                items = value
            else:
                return None
            for selection in items:
                if not isinstance(selection, dict):
                    return None
                selection_id = selection.get("id")
                if not isinstance(selection_id, int) or selection_id <= 0:
                    return None
                selection_ids.add(selection_id)
    return selection_ids


def _builder_safe_boundary_conditions(
    boundary_conditions: Any,
    selection_ids: set[int],
    obstacle_surface_ids: set[int],
) -> bool:
    if boundary_conditions is None:
        return True
    if not isinstance(boundary_conditions, dict):
        return False

    bindable = {
        "dirichlet_boundary",
        "neumann_boundary",
        "normal_aligned_neumann_boundary",
        "pressure_boundary",
        "pressure_cavity",
        "obstacle_displacements",
    }
    if any(key not in {"rhs", *bindable} for key in boundary_conditions):
        return False

    for section in bindable:
        for item in _optional_items(boundary_conditions.get(section)):
            if not isinstance(item, dict):
                return False
            expected_ids = (
                obstacle_surface_ids
                if section == "obstacle_displacements"
                else selection_ids
            )
            if item.get("id") not in expected_ids:
                return False
    return True


def _builder_safe_initial_conditions(
    initial_conditions: Any,
    volume_ids: set[int],
) -> bool:
    if initial_conditions is None:
        return True
    if not isinstance(initial_conditions, dict):
        return False

    bindable = {"velocity", "solution", "acceleration"}
    if any(key not in bindable for key in initial_conditions):
        return False

    for section in bindable:
        for item in _optional_items(initial_conditions.get(section)):
            if not isinstance(item, dict) or item.get("id") not in volume_ids:
                return False
    return True


def _optional_items(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _case_id(case: tuple[Path, Path]) -> str:
    example_path, _source_path = case
    return example_path.stem


def _require_examples_submodule() -> None:
    if not CLASSIC_EXAMPLES.exists():
        pytest.skip(
            "examples submodule is not initialized; run "
            "`git submodule update --init examples`"
        )


def _require_polyfem_data_submodule() -> None:
    if not POLYFEM_DATA_EXAMPLES.exists():
        pytest.skip(
            "polyfem-data submodule is not initialized; run "
            "`git submodule update --init polyfem-data`"
        )


@pytest.fixture(autouse=True)
def _skip_without_required_submodules():
    _require_examples_submodule()
    _require_polyfem_data_submodule()


def test_classic_generated_examples_cover_all_contact_source_jsons():
    assert _declared_source_paths() == _ordinary_contact_source_paths()


def test_classic_generated_examples_show_api_configuration():
    offenders = []
    fallback_sources = {
        (POLYFEM_DATA_EXAMPLES / source).resolve()
        for source in EXPECTED_GENERATED_API_FALLBACK_SOURCES
    }
    for example_path, source_path in _all_example_cases():
        text = example_path.read_text(encoding="utf-8")
        has_constructor_config = (
            "polyfem.config(" in text or "model.config(" in text
        )
        has_known_fallback = (
            source_path in fallback_sources and "SourcePayloadConfig(" in text
        )
        if "load_polyfem_config(SOURCE_JSON)" in text:
            offenders.append(example_path.relative_to(ROOT).as_posix())
        elif not has_constructor_config and not has_known_fallback:
            offenders.append(example_path.relative_to(ROOT).as_posix())

    assert offenders == []


def test_builder_safe_classic_examples_use_model_builder():
    offenders = []
    builder_safe_sources = _builder_safe_source_paths()
    for example_path, source_path in _all_example_cases():
        if source_path not in builder_safe_sources:
            continue
        text = example_path.read_text(encoding="utf-8")
        if "model = polyfem.model()" not in text or "polyfem_config = model.config(" not in text:
            offenders.append(example_path.relative_to(ROOT).as_posix())

    assert offenders == []


def test_obstacle_mesh_classic_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "2D"
        / "contact_2d_unit_tests_erleben_spikes_generated_api.py"
    )
    source_path = (
        POLYFEM_DATA_EXAMPLES
        / "2D"
        / "unit-tests"
        / "erleben"
        / "spikes.json"
    )

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert "polyfem_config = model.config(" in text


def test_model_builder_obstacle_mesh_does_not_create_volume_selection():
    model = polyfem.model()

    obstacle = model.mesh(is_obstacle=True,
        mesh="obstacle.obj",
        surface_selection=0,
    )
    config = model.config()

    assert obstacle is None
    assert config.as_dict()["geometry"] == [
        {
            "mesh": "obstacle.obj",
            "surface_selection": 0,
            "is_obstacle": True,
            "type": "mesh",
        }
    ]


@pytest.mark.parametrize("factory", ["mesh", "mesh_array"])
def test_unified_geometry_can_bind_material_and_obstacle_motion(factory):
    model = polyfem.model()
    kwargs = {"mesh": "part.msh"}
    if factory == "mesh_array":
        kwargs["array"] = polyfem.array(size=[2, 1], offset=1)
    body = getattr(model, factory)(**kwargs)
    obstacle = getattr(model, factory)(**kwargs, is_obstacle=True)
    assert "volume_selection" not in model.config().as_dict()["geometry"][0]
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    obstacle.dirichlet(value=[1, 0])
    cfg = model.config().as_dict()
    assert cfg["materials"][0]["id"] == cfg["geometry"][0]["volume_selection"]
    assert "volume_selection" not in cfg["geometry"][1]
    assert cfg["boundary_conditions"]["dirichlet_boundary"][0]["id"] == cfg["geometry"][1]["surface_selection"]


def test_unified_geometry_can_share_a_volume_without_explicit_ids():
    model = polyfem.model()
    first = model.mesh(mesh="first.msh")
    second = model.mesh_array(
        mesh="second.msh", array=polyfem.array(size=[2, 1], offset=1),
        same_volume_as=first,
    )
    second.material(polyfem.neo_hookean(E=1000, nu=0.3))
    cfg = model.config().as_dict()
    assert cfg["geometry"][0]["volume_selection"] == cfg["geometry"][1]["volume_selection"]
    assert cfg["materials"][0]["id"] == cfg["geometry"][0]["volume_selection"]


@pytest.mark.parametrize(
    "method,value",
    [("dirichlet", [0, 0]), ("neumann", [1, 0]), ("pressure", 2.0)],
)
def test_body_boundary_shortcut_matches_explicit_whole_surface(method, value):
    direct = polyfem.model()
    body = direct.mesh(mesh="part.msh")
    getattr(body, method)(value=value)
    getattr(body, method)(value=value)
    explicit = polyfem.model()
    surface = explicit.mesh(mesh="part.msh").surface_all()
    getattr(surface, method)(value=value)
    assert direct.config().as_dict() == explicit.config().as_dict()


def test_body_boundary_shortcut_keeps_bodies_and_volume_loads_separate():
    model = polyfem.model()
    first = model.mesh(mesh="first.msh")
    second = model.mesh(mesh="second.msh", surface_selection=1)
    first.rhs(value=[0, 1])
    first.dirichlet(value=[0, 0])
    second.dirichlet(value=[1, 0])
    cfg = model.config().as_dict()
    first_id, second_id = [g["surface_selection"] for g in cfg["geometry"]]
    assert first_id != second_id
    assert second_id == 1
    assert cfg["boundary_conditions"]["dirichlet_boundary"] == [
        {"id": first_id, "value": [0, 0]},
        {"id": second_id, "value": [1, 0]},
    ]
    assert cfg["boundary_conditions"]["rhs"] == [
        {"id": cfg["geometry"][0]["volume_selection"], "value": [0, 1]}
    ]


@pytest.mark.parametrize("whole_first", [True, False])
def test_body_boundary_shortcut_rejects_overwriting_partial_selections(whole_first):
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    if whole_first:
        body.dirichlet(value=[0, 0])
    else:
        bottom = body.surface_axis(axis=-1, position=-1)
        bottom.dirichlet(value=[0, 0])
    before = model.config().as_dict()
    with pytest.raises(ValueError, match="surface|boundary"):
        if whole_first:
            body.surface_axis(axis=1, position=1)
        else:
            body.dirichlet(value=[1, 0])
    assert model.config().as_dict() == before


def test_body_boundary_shortcut_requires_mesh_and_preserves_file_selection():
    model = polyfem.model()
    with pytest.raises(ValueError, match="body.mesh"):
        model.body().dirichlet(value=[0, 0])
    body = model.mesh(
        mesh="part.msh",
        surface_selection=polyfem.surface_selection_file(file="sidesets.txt"),
    )
    before = model.config().as_dict()
    with pytest.raises(ValueError, match="file-based"):
        body.dirichlet(value=[0, 0])
    assert model.config().as_dict() == before


def test_top_level_config_with_model_delegates_to_model_builder():
    model = polyfem.model()

    material = polyfem.neo_hookean(E=100000.0, nu=0.3)
    body = model.mesh(mesh="beam.msh")
    body.material(material)
    body.surface_all(id=1).dirichlet(value=[0.0, 0.0, 0.0])

    solver = polyfem.solver(max_threads=1)
    output = polyfem.output(directory="out")
    config_kwargs = {
        "rhs": [10, 100, 0],
        "space": polyfem.space(
            discr_order=1,
            advanced=polyfem.space_advanced(quadrature_order=4),
        ),
        "solver": solver,
        "output": output,
    }

    builder_config = model.config(**config_kwargs)
    top_level_config = polyfem.config(model=model, **config_kwargs)

    assert top_level_config.as_dict() == builder_config.as_dict()

    builder_canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=builder_config,
        dtype=None,
    )
    top_level_canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=top_level_config,
        dtype=None,
    )
    assert top_level_canonical.backend_settings == builder_canonical.backend_settings


def test_top_level_body_mesh_and_model_rhs_match_legacy_builder_payload():
    legacy_model = polyfem.model()
    legacy_body = legacy_model.mesh(mesh="beam.msh")
    legacy_body.material(polyfem.neo_hookean(E=100000.0, nu=0.3))
    legacy_body.surface_all(id=1).dirichlet(value=[0.0, 0.0, 0.0])
    legacy_config = polyfem.config(model=legacy_model, rhs=[10, 100, 0])

    model = polyfem.model()
    mesh_descriptor = polyfem.mesh(mesh="beam.msh")
    mesh_payload_before_attach = mesh_descriptor.as_dict()
    body = polyfem.body(model=model)
    body.mesh(mesh_descriptor)
    body.material(polyfem.neo_hookean(E=100000.0, nu=0.3))
    body.surface_all(id=1).dirichlet(value=[0.0, 0.0, 0.0])
    model.rhs([10, 100, 0])
    config = polyfem.config(model=model)

    assert config.as_dict() == legacy_config.as_dict()
    assert mesh_descriptor.as_dict() == mesh_payload_before_attach

    legacy_canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=legacy_config,
        dtype=None,
    )
    canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=config,
        dtype=None,
    )
    assert canonical.backend_settings == legacy_canonical.backend_settings


def test_model_rhs_rejects_duplicate_rhs_sources():
    model = polyfem.model()
    model.rhs([10, 100, 0])

    with pytest.raises(TypeError, match="rhs"):
        polyfem.config(model=model, rhs=[0, 0, 0])
    with pytest.raises(TypeError, match="rhs"):
        polyfem.config(
            model=model,
            boundary_conditions=polyfem.boundary_conditions(rhs=[0, 0, 0]),
        )


def test_body_rhs_preserves_each_body_selection_and_rejects_conflicts():
    model = polyfem.model()
    first = model.mesh(mesh="first.obj")
    second = model.mesh(mesh="second.obj")
    first.material(polyfem.neo_hookean(E=1000, nu=0.3))
    second.material(polyfem.neo_hookean(E=2000, nu=0.3))
    first.rhs(value=[1, 0])
    second.rhs(value=[0, 2])
    payload = polyfem.config(model=model).as_dict()
    ids = [item["volume_selection"] for item in payload["geometry"]]
    assert ids[0] != ids[1]
    assert payload["boundary_conditions"]["rhs"] == [
        {"id": ids[0], "value": [1, 0]},
        {"id": ids[1], "value": [0, 2]},
    ]
    assert [item["id"] for item in payload["materials"]] == ids
    with pytest.raises(TypeError, match="rhs"):
        first.rhs(value=[3, 4])
    with pytest.raises(TypeError, match="rhs"):
        model.rhs([1, 1])
    with pytest.raises(TypeError, match="rhs"):
        polyfem.config(model=model, rhs=[1, 1])
    with pytest.raises(TypeError, match="rhs"):
        polyfem.config(model=model, boundary_conditions={"rhs": [1, 1]})

    global_model = polyfem.model()
    global_model.rhs([1, 1])
    global_body = global_model.mesh(mesh="first.obj")
    with pytest.raises(TypeError, match="rhs"):
        global_body.rhs(value=[1, 0])
    unbound = polyfem.body(model=polyfem.model())
    with pytest.raises(ValueError, match="body.mesh"):
        unbound.rhs(value=[1, 0])


def test_global_material_geometry_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "2D"
        / "contact_2d_static_friction_slope_generated_api.py"
    )
    source_path = POLYFEM_DATA_EXAMPLES / "2D" / "static" / "friction-slope.json"

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert "model.mesh(" in text
    assert "polyfem_config = model.config(" in text


def test_model_builder_geometry_mesh_does_not_create_volume_selection():
    model = polyfem.model()

    geometry = model.mesh(
        mesh="static-wall.obj",
    )
    surface = geometry.surface_all(id=1)
    surface.dirichlet(value=[0, 0])
    config = model.config()

    assert config.as_dict()["geometry"] == [
        {
            "mesh": "static-wall.obj",
            "surface_selection": 1,
            "type": "mesh",
        }
    ]
    assert config.as_dict()["boundary_conditions"]["dirichlet_boundary"] == [
        {
            "id": 1,
            "value": [0, 0],
        }
    ]


def test_single_dict_selection_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_unit_tests_2_cubes_generated_api.py"
    )
    source_path = POLYFEM_DATA_EXAMPLES / "3D" / "unit-tests" / "2-cubes.json"

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert ".surface_axis(" in text
    assert "polyfem_config = model.config(" in text


def test_mesh_sequence_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_mesh_sequence_kick_generated_api.py"
    )
    source_path = POLYFEM_DATA_EXAMPLES / "3D" / "mesh-sequence" / "kick.json"

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert "model.mesh_sequence(" in text
    assert "polyfem_config = model.config(" in text


def test_model_builder_mesh_sequence_preserves_payload():
    model = polyfem.model()

    obstacle = model.mesh_sequence(is_obstacle=True,
        mesh_sequence="kick-sequence/",
        fps=24,
    )
    config = model.config()

    assert obstacle is None
    assert config.as_dict()["geometry"] == [
        {
            "mesh_sequence": "kick-sequence/",
            "fps": 24,
            "is_obstacle": True,
            "type": "mesh_sequence",
        }
    ]



def test_implicit_obstacle_motion_targets_only_the_bound_obstacle():
    model = polyfem.model()
    model.mesh(is_obstacle=True, mesh="lower.obj")
    obstacle = model.mesh(is_obstacle=True, mesh="upper.obj")
    obstacle.dirichlet(value=[1, 0, 0])
    moved = model.config().as_dict()
    assert "surface_selection" not in moved["geometry"][0]
    surface_id = moved["geometry"][1]["surface_selection"]
    assert moved["boundary_conditions"]["dirichlet_boundary"] == [
        {"id": surface_id, "value": [1, 0, 0]}
    ]


def test_mesh_array_pile_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_pile_cubes_generated_api.py"
    )
    source_path = POLYFEM_DATA_EXAMPLES / "3D" / "pile" / "cubes.json"

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert "model.mesh_array(" in text
    assert "polyfem_config = model.config(" in text


def test_model_builder_mesh_array_without_volume_selection_preserves_payload():
    model = polyfem.model()

    mesh_array = model.mesh_array(
        mesh="cube.msh",
        array=polyfem.array(
            size=[2, 3, 4],
            offset=1.25,
        ),
    )
    config = model.config()

    assert callable(mesh_array.material)
    assert config.as_dict()["geometry"] == [
        {
            "mesh": "cube.msh",
            "array": {
                "size": [2, 3, 4],
                "offset": 1.25,
            },
            "type": "mesh_array",
        }
    ]


def test_obstacle_displacement_example_uses_model_builder():
    example_path = (
        CLASSIC_EXAMPLES
        / "3D"
        / "contact_3d_friction_ball_rollers_generated_api.py"
    )
    source_path = (
        POLYFEM_DATA_EXAMPLES
        / "3D"
        / "friction"
        / "ball-rollers.json"
    )

    assert source_path.resolve() in _builder_safe_source_paths()

    text = example_path.read_text(encoding="utf-8")
    assert "model = polyfem.model()" in text
    assert ".obstacle_displacement(" in text


def test_model_builder_can_bind_obstacle_surface_displacement():
    model = polyfem.model()

    surface = model.mesh(is_obstacle=True,
        mesh="roller.obj",
        surface_selection=1000,
    )
    surface.obstacle_displacement(value=[1, 0, 0])
    config = model.config()

    assert config.as_dict()["geometry"] == [
        {
            "mesh": "roller.obj",
            "surface_selection": 1000,
            "is_obstacle": True,
            "type": "mesh",
        }
    ]
    assert config.as_dict()["boundary_conditions"]["obstacle_displacements"] == [
        {
            "id": 1000,
            "value": [1, 0, 0],
        }
    ]


def _import_example(example_path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(example_path.stem, example_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _import_module(module_path: Path, module_name: str) -> Any:
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)
    return module


def _deep_merge(base: Any, override: Any) -> Any:
    if isinstance(base, dict) and isinstance(override, dict):
        merged = copy.deepcopy(base)
        for key, value in override.items():
            merged[key] = (
                _deep_merge(merged[key], value)
                if key in merged
                else copy.deepcopy(value)
            )
        return merged
    return copy.deepcopy(override)


def _load_source_config(source_path: Path) -> dict[str, Any]:
    source_path = source_path.resolve()
    with source_path.open(encoding="utf-8") as f:
        source = json.load(f)

    common_ref = source.pop("common", None)
    source.pop("tests", None)
    source.pop("default_params", None)
    patch = source.pop("patch", None)

    if common_ref is None:
        merged = source
    else:
        common_path = (source_path.parent / common_ref).resolve()
        merged = _deep_merge(_load_source_config(common_path), source)

    if patch:
        merged = _apply_json_patch(merged, patch)

    _resolve_mesh_paths(merged, source_path.parent)
    _normalize_solver_fallbacks(merged)
    _normalize_generated_mesh_payloads(merged, source_path.parent)
    return merged


def _expected_config_from_source(source_path: Path) -> dict[str, Any]:
    config = _load_source_config(source_path)
    output = config.setdefault("output", {})
    case = source_path.relative_to(POLYFEM_DATA_EXAMPLES).with_suffix("")
    output.setdefault("directory", f"results/{case.as_posix()}")
    # The generated configuration represents these singleton entries as lists.
    if isinstance(config.get("materials"), dict):
        config["materials"] = [config["materials"]]
    boundary = config.get("boundary_conditions", {})
    if isinstance(boundary.get("dirichlet_boundary"), dict):
        boundary["dirichlet_boundary"] = [boundary["dirichlet_boundary"]]
    return config


def _apply_json_patch(config: dict[str, Any], patch: list[dict[str, Any]]) -> dict[str, Any]:
    patched = copy.deepcopy(config)
    for item in patch:
        assert item.get("op") == "replace"
        _replace_json_pointer(patched, item["path"], item.get("value"))
    return patched


def _replace_json_pointer(config: Any, pointer: str, value: Any) -> None:
    target = config
    tokens = _json_pointer_tokens(pointer)
    for token in tokens[:-1]:
        target = target[int(token)] if isinstance(target, list) else target[token]

    last = tokens[-1]
    if isinstance(target, list):
        target[int(last)] = copy.deepcopy(value)
    else:
        target[last] = copy.deepcopy(value)


def _json_pointer_tokens(pointer: str) -> list[str]:
    assert pointer.startswith("/")
    return [
        token.replace("~1", "/").replace("~0", "~")
        for token in pointer[1:].split("/")
    ]


def _resolve_mesh_paths(value: Any, source_dir: Path) -> None:
    if isinstance(value, list):
        for item in value:
            _resolve_mesh_paths(item, source_dir)
    elif isinstance(value, dict):
        for key in ("mesh", "linear_map"):
            path_value = value.get(key)
            if isinstance(path_value, str):
                value[key] = str((source_dir / path_value).resolve())
        for child in value.values():
            _resolve_mesh_paths(child, source_dir)


def _normalize_solver_fallbacks(config: dict[str, Any]) -> None:
    solver = config.get("solver")
    if not isinstance(solver, dict):
        return
    linear = solver.get("linear")
    if not isinstance(linear, dict):
        return
    solver_name = linear.get("solver")
    if isinstance(solver_name, list):
        linear["solver"] = solver_name[0]


def _normalize_generated_mesh_payloads(config: dict[str, Any], source_dir: Path) -> None:
    geometry = config.get("geometry")
    if not isinstance(geometry, list):
        return

    for item in geometry:
        if not isinstance(item, dict) or "mesh" not in item:
            continue
        item.setdefault("type", "mesh")
        _normalize_file_selection(item, "volume_selection", source_dir)
        _normalize_file_selection(item, "surface_selection", source_dir)

        transformation = item.get("transformation")
        if not isinstance(transformation, dict):
            continue
        for field in ("rotation", "scale"):
            if isinstance(transformation.get(field), (int, float)):
                transformation[field] = [float(transformation[field])]


def _normalize_file_selection(item: dict[str, Any], field: str, source_dir: Path) -> None:
    selection = item.get(field)
    if isinstance(selection, str):
        item[field] = {"file": _resolve_source_file_path(selection, source_dir)}
    elif isinstance(selection, dict):
        file_path = selection.get("file")
        if isinstance(file_path, str):
            selection["file"] = _resolve_source_file_path(file_path, source_dir)


def _resolve_source_file_path(path_value: str, source_dir: Path) -> str:
    path = Path(path_value)
    if path.is_absolute():
        return path_value
    return str((source_dir / path).resolve())


@pytest.mark.parametrize("case", EXAMPLE_CASES, ids=_case_id)
def test_classic_generated_example_compiles(case: tuple[Path, Path]):
    example_path, _source_path = case
    py_compile.compile(str(example_path), doraise=True)


@pytest.mark.parametrize("case", EXAMPLE_CASES, ids=_case_id)
def test_classic_generated_example_declares_polyfem_data_source(case: tuple[Path, Path]):
    example_path, source_path = case
    example = _import_example(example_path)

    assert Path(example.SOURCE_JSON).resolve() == source_path.resolve()
    assert POLYFEM_DATA_CONTACT in Path(example.SOURCE_JSON).resolve().parents


def test_classic_common_helpers_support_polyfem_data_root_env(monkeypatch, tmp_path):
    data_root = tmp_path / "polyfem-data-alt"
    monkeypatch.setenv("POLYFEM_DATA_ROOT", str(data_root))

    common_2d = _import_module(
        CLASSIC_EXAMPLES / "2D" / "_contact_2d_common.py",
        "contact_2d_common_env_test",
    )
    common_3d = _import_module(
        CLASSIC_EXAMPLES / "3D" / "_contact_3d_common.py",
        "contact_3d_common_env_test",
    )

    assert common_2d.POLYFEM_DATA_ROOT == data_root.resolve()
    assert common_2d.POLYFEM_DATA_CONTACT == data_root.resolve() / "contact"
    assert common_2d.CONTACT_EXAMPLES_DIR == data_root.resolve() / "contact" / "examples"
    assert common_2d.MESHES_DIR == data_root.resolve() / "contact" / "meshes"
    assert common_3d.POLYFEM_DATA_ROOT == data_root.resolve()
    assert common_3d.POLYFEM_DATA_CONTACT == data_root.resolve() / "contact"


def test_contact_source_loader_resolves_file_selection_paths():
    source_path = POLYFEM_DATA_EXAMPLES / "3D" / "hex_tet" / "params.json"
    loader = _import_module(
        CLASSIC_EXAMPLES / "_contact_source_loader.py",
        "contact_source_loader_file_selection_test",
    )

    payload = loader.load_contact_source_payload(source_path)

    assert payload["geometry"][0]["surface_selection"] == {
        "file": str((source_path.parent / "surface_sidesets11.txt").resolve())
    }


def _assert_source_parity(actual, expected):
    """Compare geometry labels and their consumers under a one-to-one renumbering."""
    actual = copy.deepcopy(actual)
    maps = {"volume": {}, "boundary": {}}

    def match(label, source_label, group):
        if not isinstance(label, int) or not isinstance(source_label, int) or min(label, source_label) <= 0:
            assert label == source_label
            return label
        mapping = maps[group]
        assert mapping.setdefault(label, source_label) == source_label
        return source_label

    def align(value, source, group):
        if isinstance(value, list):
            assert len(value) == len(source)
            return [align(v, s, group) for v, s in zip(value, source)]
        if isinstance(value, dict):
            if "id" in value:
                value["id"] = match(value["id"], source["id"], group)
            return value
        return match(value, source, group)

    expected = copy.deepcopy(expected)
    conditions = expected.get("boundary_conditions", {})
    for geometry, source in zip(actual["geometry"], expected["geometry"]):
        if (source.get("is_obstacle") and source.get("type", "mesh") in ("mesh", "mesh_array")
                and source.get("surface_selection") == 0 and "surface_selection" not in geometry):
            # GeometryLoader only consults an obstacle tag when choosing its displacement.
            for key in ("dirichlet_boundary", "obstacle_displacements"):
                entries = conditions.get(key, [])
                for entry in entries if isinstance(entries, list) else [entries]:
                    target = entry["id"]
                    assert target not in (0, "all")
                    assert not isinstance(target, list) or 0 not in target
            del source["surface_selection"]
        for field, group in (("volume_selection", "volume"), ("surface_selection", "boundary"), ("point_selection", "boundary")):
            if field in geometry:
                if (field == "surface_selection" and not source.get("is_obstacle")
                        and source[field] == 0 and isinstance(geometry[field], int) and geometry[field] > 0):
                    label = geometry[field]
                    assert maps[group].setdefault(label, 0) == 0
                    geometry[field] = 0
                else:
                    geometry[field] = align(geometry[field], source[field], group)
    for mapping in maps.values():
        assert len(set(mapping.values())) == len(mapping)

    def translate(items, group):
        for item in items if isinstance(items, list) else [items]:
            if isinstance(item, dict) and isinstance(item.get("id"), int):
                label = item["id"]
                mapping = maps[group]
                if label in mapping:
                    item["id"] = mapping[label]
                else:
                    # An external mesh tag can stay literal. An old, reassigned ID cannot.
                    assert label not in mapping.values(), f"Stale {group} ID {label}"

    translate(actual.get("materials"), "volume")
    translate(actual.get("space", {}).get("discr_order"), "volume")
    for items in actual.get("initial_conditions", {}).values():
        translate(items, "volume")
    for items in actual.get("boundary_conditions", {}).values():
        translate(items, "boundary")
    assert actual == expected


def test_source_parity_rejects_stale_boundary_id():
    source = {
        "geometry": [{"surface_selection": 1000}],
        "boundary_conditions": {"dirichlet_boundary": [{"id": 1000, "value": [1, 0]}]},
    }
    actual = copy.deepcopy(source)
    actual["geometry"][0]["surface_selection"] = 1
    with pytest.raises(AssertionError, match="Stale boundary ID 1000"):
        _assert_source_parity(actual, source)


@pytest.mark.parametrize("case", EXAMPLE_CASES, ids=_case_id)
def test_classic_generated_example_matches_polyfem_data_source(case: tuple[Path, Path]):
    example_path, source_path = case
    example = _import_example(example_path)

    _assert_source_parity(example.polyfem_config.as_dict(), _expected_config_from_source(source_path))


def test_all_classic_generated_examples_match_polyfem_data_sources():
    for example_path, source_path in _all_example_cases():
        example = _import_example(example_path)
        try:
            _assert_source_parity(example.polyfem_config.as_dict(), _expected_config_from_source(source_path))
        except AssertionError as error:
            raise AssertionError(str(source_path)) from error



@pytest.mark.parametrize("case", DOLPHIN_MATERIAL_SHAPE_CASES, ids=_case_id)
def test_dolphin_generated_backend_payload_keeps_source_materials_shape(
    case: tuple[Path, Path],
):
    example_path, source_path = case
    source_config = _load_source_config(source_path)
    example = _import_example(example_path)

    canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=example.polyfem_config,
        dtype=None,
    )

    assert isinstance(source_config["materials"], dict)
    assert canonical.backend_settings["materials"] == source_config["materials"]


def test_classic_generated_api_fallbacks_are_limited_to_known_schema_gaps():
    actual = set()
    for _example_path, source_path in _all_example_cases():
        example = _import_example(_example_path)
        fallback_reason = getattr(example.polyfem_config, "fallback_reason", None)
        if fallback_reason:
            actual.add(source_path.relative_to(POLYFEM_DATA_EXAMPLES).as_posix())

    assert actual == EXPECTED_GENERATED_API_FALLBACK_SOURCES


@pytest.mark.parametrize("case", EXAMPLE_CASES, ids=_case_id)
def test_classic_generated_example_uses_generated_solve_path(case: tuple[Path, Path]):
    example_path, _source_path = case
    example = _import_example(example_path)

    canonical = prepare_canonical_solve_input(
        vertices=None,
        cells=None,
        cfg=example.polyfem_config,
        dtype=None,
    )

    assert canonical.metadata == {"mesh_source": "json", "config_source": "generated"}
    assert canonical.mesh_source.mode == "json"
    assert "geometry" in canonical.backend_settings


def test_all_classic_generated_examples_use_generated_solve_path():
    failures = []
    for example_path, _source_path in _all_example_cases():
        example = _import_example(example_path)
        canonical = prepare_canonical_solve_input(
            vertices=None,
            cells=None,
            cfg=example.polyfem_config,
            dtype=None,
        )
        if canonical.metadata != {"mesh_source": "json", "config_source": "generated"}:
            failures.append(example_path.relative_to(ROOT).as_posix())

    assert failures == []
