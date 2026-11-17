"""Model settings must survive export without changing existing config calls."""

import pytest

import polyfem


@pytest.mark.parametrize(
    "section, fields",
    [
        ("space", {"discr_order": 2}),
        ("initial_conditions", {"velocity": [{"id": 7, "value": [1, 0]}]}),
        ("boundary_conditions", {"rhs": [0, -1]}),
        ("time", {"dt": 0.01, "tend": 0.1}),
        ("output", {"directory": "results"}),
    ],
)
def test_model_settings_match_explicit_config(section, fields):
    model = polyfem.model()
    returned = getattr(model, section)(**fields)
    if section == "boundary_conditions":
        assert isinstance(returned, dict)
    else:
        assert returned is model
    expected = polyfem.model().config(**{section: getattr(polyfem, section)(**fields)})
    assert model.config().as_dict() == expected.as_dict()
    assert model.config().as_dict() == expected.as_dict()


@pytest.mark.parametrize("as_list", [False, True])
def test_material_setting_accepts_one_material_or_a_list(as_list):
    material = polyfem.neo_hookean(E=1000, nu=0.3)
    value = [material] if as_list else material
    model = polyfem.model()
    assert model.materials(value) is model
    assert model.config().as_dict() == polyfem.model().config(materials=value).as_dict()


@pytest.mark.parametrize("settings_first", [False, True])
def test_global_fields_and_body_lists_survive_in_either_order(settings_first):
    model = polyfem.model()
    body = model.mesh(mesh="part.msh", volume_selection=7)

    def settings():
        model.space(discr_orderq=4)
        model.initial_conditions(velocity=[{"id": 8, "value": [1, 0]}])
        model.materials(polyfem.neo_hookean(id=8, E=2000, nu=0.3))
        model.boundary_conditions(rhs=[0, -1])

    if settings_first:
        settings()
    body.discr_order(order=2)
    body.velocity(value=[0, 0])
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    body.dirichlet(value=[0, 0])
    if not settings_first:
        settings()

    config = model.config().as_dict()
    assert config["space"] == {"discr_orderq": 4, "discr_order": [{"id": 7, "order": 2}]}
    assert [entry["id"] for entry in config["materials"]] == [8, 7]
    assert [entry["id"] for entry in config["initial_conditions"]["velocity"]] == [8, 7]
    assert config["boundary_conditions"]["rhs"] == [0, -1]
    assert len(config["boundary_conditions"]["dirichlet_boundary"]) == 1
    for section in ("space", "materials", "initial_conditions", "boundary_conditions"):
        assert getattr(model, section)() == config[section]


def test_replacing_settings_preserves_body_values_without_duplicates():
    model = polyfem.model()
    body = model.mesh(mesh="part.msh", volume_selection=7)
    body.velocity(value=[0, 0])
    body.dirichlet(value=[0, 0])
    model.initial_conditions(velocity=[{"id": 8, "value": [1, 0]}])
    model.initial_conditions(velocity=[{"id": 9, "value": [2, 0]}])
    model.space(discr_order=1)
    model.space(discr_order=2)
    model.boundary_conditions(rhs=[0, -1])
    model.boundary_conditions(rhs=[0, -2])
    first = model.config().as_dict()
    assert first == model.config().as_dict()
    assert first["space"] == {"discr_order": 2}
    assert [entry["id"] for entry in first["initial_conditions"]["velocity"]] == [9, 7]
    assert first["boundary_conditions"]["rhs"] == [0, -2]
    assert len(first["boundary_conditions"]["dirichlet_boundary"]) == 1


def test_boundary_dictionary_can_still_be_passed_to_config():
    model = polyfem.model()
    model.mesh(mesh="part.msh").dirichlet(value=[0, 0])
    boundary = model.boundary_conditions(rhs=[0, -1])
    assert isinstance(boundary, dict)
    expected = model.config().as_dict()
    assert model.config(boundary_conditions=boundary).as_dict() == expected
    assert len(expected["boundary_conditions"]["dirichlet_boundary"]) == 1


def test_explicit_boundary_list_and_body_conditions_are_merged_once():
    model = polyfem.model()
    model.mesh(mesh="part.msh").dirichlet(value=[0, 0])
    boundary = model.boundary_conditions(
        dirichlet_boundary=[{"id": 9, "value": [1, 0]}],
    )
    assert len(boundary["dirichlet_boundary"]) == 2
    assert model.config(boundary_conditions=boundary).as_dict()["boundary_conditions"] == boundary
    assert model.config().as_dict()["boundary_conditions"] == boundary


def test_body_and_selection_ids_stay_automatic_with_model_settings():
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    body.material(polyfem.neo_hookean(E=1000, nu=0.3))
    body.velocity(value=[0, 0])
    body.surface_axis(axis="-x", position=0.0, relative=True).dirichlet(value=[0, 0])
    model.space(discr_order=2)
    model.boundary_conditions(rhs=[0, -1])
    config = model.config().as_dict()
    geometry = config["geometry"][0]
    assert config["materials"][0]["id"] == geometry["volume_selection"]
    assert config["initial_conditions"]["velocity"][0]["id"] == geometry["volume_selection"]
    assert config["boundary_conditions"]["dirichlet_boundary"][0]["id"] == geometry["surface_selection"][0]["id"]


@pytest.mark.parametrize("section", ["space", "materials", "initial_conditions", "boundary_conditions"])
def test_conflicting_config_argument_is_not_silently_ignored(section):
    model = polyfem.model()
    if section == "materials":
        model.materials(polyfem.neo_hookean(E=1000, nu=0.3))
    elif section == "space":
        model.space(discr_order=2)
    elif section == "initial_conditions":
        model.initial_conditions(velocity=[{"id": 7, "value": [0, 0]}])
    else:
        model.boundary_conditions(rhs=[0, -1])
    with pytest.raises(TypeError, match="already set"):
        model.config(**{section: {}})


@pytest.mark.parametrize("settings_first", [False, True])
@pytest.mark.parametrize("section", ["space", "boundary_conditions"])
def test_global_and_body_values_cannot_replace_each_other(section, settings_first):
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    fields = {"discr_order": 2} if section == "space" else {"rhs": [0, -1]}
    if settings_first:
        getattr(model, section)(**fields)
    if section == "space":
        body.discr_order(order=1)
    else:
        body.rhs(value=[0, -2])
    if settings_first:
        with pytest.raises((TypeError, ValueError)):
            model.config()
    else:
        before = model.config().as_dict()
        with pytest.raises((TypeError, ValueError)):
            getattr(model, section)(**fields)
        assert model.config().as_dict() == before


def test_reading_or_exporting_does_not_change_settings():
    model = polyfem.model()
    fields = {"velocity": [{"id": 7, "value": [1, 0]}]}
    model.initial_conditions(**fields)
    expected = model.config().as_dict()
    fields["velocity"][0]["value"][0] = 99
    model.initial_conditions()["velocity"][0]["value"][0] = 88
    assert model.config().as_dict() == expected


def test_explicit_config_is_still_a_one_time_override():
    model = polyfem.model()
    body = model.mesh(mesh="part.msh")
    body.discr_order(order=1)
    before = model.config().as_dict()
    assert model.config(space=polyfem.space(discr_order=2)).as_dict()["space"] == {"discr_order": 2}
    assert model.config().as_dict() == before
    model.config(rhs=[0, -1])
    assert model.config().as_dict() == before
