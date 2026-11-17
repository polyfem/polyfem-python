"""Check model associations without loading meshes or running the native solver."""

import pytest

import polyfem


def select_left(body):
    return body.surface_axis(axis="-x", position=0.0, relative=True)


def test_repeated_partial_selection_keeps_one_geometry_rule():
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    first = select_left(body)
    first.dirichlet(value=[0, 0])
    second = select_left(body)
    second.neumann(value=[1, 0])
    cfg = model.config().as_dict()
    selections = cfg["geometry"][0]["surface_selection"]
    assert len(selections) == 1
    for key in ("dirichlet_boundary", "neumann_boundary"):
        assert cfg["boundary_conditions"][key][0]["id"] == selections[0]["id"]


@pytest.mark.parametrize("direct_first", [False, True])
def test_repeated_whole_surface_reuses_its_id(direct_first):
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    if direct_first:
        first = body.dirichlet(value=[0, 0])
    else:
        first = body.surface_all().dirichlet(value=[0, 0])
    second = body.surface_all().neumann(value=[1, 0])
    assert first.selection_ref.backend_id == second.selection_ref.backend_id
    cfg = model.config().as_dict()
    assert cfg["geometry"][0]["surface_selection"] == first.selection_ref.backend_id
    assert cfg["boundary_conditions"]["neumann_boundary"][0]["id"] == first.selection_ref.backend_id


@pytest.mark.parametrize("field", ["point_selection", "surface_selection", "volume_selection"])
def test_existing_selection_rules_reserve_ids_before_auto_allocation(field):
    model = polyfem.model()
    model.mesh(mesh="tagged.msh", **{field: [{"id": 1, "axis": "-x", "position": 0.0}]})
    body = model.mesh(mesh="automatic.msh")
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    body.dirichlet(value=[0, 0])
    cfg = model.config().as_dict()
    produced_field = "volume_selection" if field == "volume_selection" else "surface_selection"
    assert cfg["geometry"][1][produced_field] != 1
    assert cfg["geometry"][0][field][0]["id"] == 1


def test_bind_rejects_a_foreign_ref_even_when_backend_ids_match():
    first_model = polyfem.model()
    foreign = first_model.mesh(mesh="first.msh").surface_all()
    second_model = polyfem.model()
    local = second_model.mesh(mesh="second.msh").surface_all()
    assert foreign.selection_ref.backend_id == local.selection_ref.backend_id
    before = second_model.config().as_dict()
    with pytest.raises(ValueError, match="this model"):
        second_model.bind("surface.dirichlet", foreign.selection_ref, {"value": [9, 9]})
    assert second_model.config().as_dict() == before
    second_model.bind("surface.dirichlet", local.selection_ref, {"value": [0, 0]})
    assert second_model.config().as_dict()["boundary_conditions"]["dirichlet_boundary"] == [
        {"id": local.selection_ref.backend_id, "value": [0, 0]},
    ]


@pytest.mark.parametrize("reverse", [False, True])
def test_multiple_bodies_keep_their_own_material_and_boundary_ids(reverse):
    model = polyfem.model()
    bodies = [model.mesh(mesh="first.msh"), model.mesh(mesh="second.msh")]
    indices = [1, 0] if reverse else [0, 1]
    for index in indices:
        body = bodies[index]
        body.material(polyfem.neo_hookean(E=1000 * (index + 1), nu=0.3))
        body.velocity(value=[index, 0])
        select_left(body).dirichlet(value=[index, 0])
        body.point_axis(axis="-y", position=0.0, relative=True).neumann(value=[index + 1, 0])
    cfg = model.config().as_dict()
    boundary_ids = []
    for index, geometry in enumerate(cfg["geometry"]):
        material = next(item for item in cfg["materials"] if item["id"] == geometry["volume_selection"])
        assert material["E"] == 1000 * (index + 1)
        for field, consumer, value in (
            ("surface_selection", "dirichlet_boundary", [index, 0]),
            ("point_selection", "neumann_boundary", [index + 1, 0]),
        ):
            tag = geometry[field][0]["id"]
            boundary_ids.append(tag)
            assert next(item for item in cfg["boundary_conditions"][consumer] if item["id"] == tag)["value"] == value
    assert len(set(boundary_ids)) == 4
    assert len({g["volume_selection"] for g in cfg["geometry"]}) == 2


def test_shared_surfaces_and_volumes_work_without_explicit_ids():
    model = polyfem.model()
    first = model.mesh(mesh="first.msh")
    second = model.mesh(mesh="second.msh", same_volume_as=first)
    first.material(polyfem.neo_hookean(E=1000, nu=0.3))
    second.material(polyfem.neo_hookean(E=1000, nu=0.3))
    surface = first.surface_all()
    second.surface_all(same_as=surface).dirichlet(value=[0, 0])
    surface.dirichlet(value=[0, 0])
    cfg = model.config().as_dict()
    assert cfg["geometry"][0]["volume_selection"] == cfg["geometry"][1]["volume_selection"]
    assert cfg["geometry"][0]["surface_selection"] == cfg["geometry"][1]["surface_selection"]
    assert len(cfg["materials"]) == 1
    assert len(cfg["boundary_conditions"]["dirichlet_boundary"]) == 1


def test_cross_model_sharing_is_rejected_without_changing_the_model():
    first = polyfem.model().mesh(mesh="first.msh")
    surface = first.surface_all()
    model = polyfem.model()
    local = model.mesh(mesh="local.msh")
    before = model.config().as_dict()
    with pytest.raises(ValueError, match="this model"):
        model.mesh(mesh="second.msh", same_volume_as=first)
    with pytest.raises(ValueError, match="this model"):
        local.surface_all(same_as=surface)
    assert model.config().as_dict() == before


def test_point_reference_cannot_be_used_as_a_surface_reference():
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    point = body.point_axis(axis="-x", position=0.0, relative=True)
    before = model.config().as_dict()
    with pytest.raises(ValueError):
        body.surface_all(same_as=point)
    with pytest.raises(ValueError, match="expects surface"):
        model.bind("surface.dirichlet", point.selection_ref, {"value": [0, 0]})
    assert model.config().as_dict() == before


def test_explicit_shared_ids_remain_supported():
    model = polyfem.model()
    for mesh in ("first.msh", "second.msh"):
        body = model.mesh(mesh=mesh, volume_selection=7)
        body.material(polyfem.neo_hookean(E=1000, nu=0.3))
        body.surface_all(id=9).dirichlet(value=[0, 0])
    cfg = model.config().as_dict()
    assert [g["volume_selection"] for g in cfg["geometry"]] == [7, 7]
    assert [g["surface_selection"] for g in cfg["geometry"]] == [9, 9]
    assert len(cfg["materials"]) == 1
    assert len(cfg["boundary_conditions"]["dirichlet_boundary"]) == 1


def test_file_surface_refs_belong_to_the_model_and_reserve_all_tags(tmp_path):
    sideset = tmp_path / "boundary.txt"
    sideset.write_text("1 0\n5 1\n", encoding="utf-8")
    model = polyfem.model()
    tagged = model.mesh(mesh="tagged.msh")
    tagged.surface_from_file(sideset, tag=1).dirichlet(value=[0, 0])
    tagged.surface_from_file(sideset, tag=5).neumann(value=[1, 0])
    other = model.mesh(mesh="other.msh").surface_all()
    other.dirichlet(value=[1, 0])
    assert other.selection_ref.backend_id > 5
    boundary = model.config().as_dict()["boundary_conditions"]
    assert boundary["dirichlet_boundary"][0]["id"] == 1
    assert boundary["neumann_boundary"][0]["id"] == 5


def test_obstacle_and_body_surface_refs_have_distinct_ids():
    model = polyfem.model()
    obstacle = model.mesh(mesh="wall.obj", is_obstacle=True)
    first = obstacle.dirichlet(value=[0, 0])
    obstacle.surface_all().dirichlet(value=[0, 0])
    second = model.mesh(mesh="part.msh").surface_all().dirichlet(value=[1, 0])
    assert first.selection_ref.backend_id != second.selection_ref.backend_id
    cfg = model.config().as_dict()
    assert len(cfg["boundary_conditions"]["dirichlet_boundary"]) == 2
