"""Schema-generated time factories and builder payload compatibility."""

import pytest

from polyfempy.generated_api import generated_api as pf


@pytest.mark.parametrize("fields", [
    {"dt": 0.01, "tend": 0.5},
    {"dt": 0.01, "time_steps": 50},
    {"tend": 0.5, "time_steps": 50},
])
def test_time_factory_and_builder_match_schema(fields):
    fields = dict(fields, t0=0, integrator="BDF")
    expected = pf.config(time=fields).as_dict()["time"]
    assert pf.time(**fields).as_dict() == expected
    model = pf.model()
    assert model.time(**fields) is model
    assert model.config().as_dict()["time"] == expected


@pytest.mark.parametrize("fields", [
    {"dt": 0.01},
    {"dt": 0.01, "tend": 0.5, "typo": True},
])
def test_invalid_time_parameters_fail(fields):
    with pytest.raises(TypeError):
        pf.time(**fields)
    with pytest.raises(TypeError):
        pf.model().time(**fields)


def test_existing_time_wrapper_and_model_configuration():
    old = pf.time(pf.object1(t0=0, dt=0.01, tend=0.5, integrator="BDF"))
    model = pf.model()
    model.time(t0=0, dt=0.01, tend=0.5, integrator="BDF")
    model.units(length="cm", mass="g", time="s")
    model.output(directory="runs/example")
    config = model.config().as_dict()
    assert config["time"] == old.as_dict()
    assert config["units"] == {"length": "cm", "mass": "g", "time": "s"}
    assert config["output"]["directory"] == "runs/example"
    assert "time" not in pf.model().config().as_dict()
    with pytest.raises(TypeError, match="already set"):
        model.config(time=old)
