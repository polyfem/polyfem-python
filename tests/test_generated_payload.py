"""Generated configurations reach PolyFEM with their chosen values intact."""

import copy

import polyfempy as pf

from polyfempy.solve import _config


def test_generated_material_union_is_not_rewritten():
    cfg = pf.config(materials=[pf.neo_hookean(E=1000, nu=0.3)])
    assert _config(cfg) == cfg.as_dict()


def test_tagged_materials_remain_a_list():
    cfg = pf.config(materials=[pf.neo_hookean(E=1000, nu=0.3, id=7)])
    assert _config(cfg)["materials"] == [{"id": 7, "type": "NeoHookean", "E": 1000, "nu": 0.3}]


def test_explicit_solver_values_pass_through_unchanged():
    cfg = pf.config(solver=pf.solver(
        advanced=pf.solver_advanced(cache_size=900000, lump_mass_matrix=False),
        nonlinear=pf.nonlinear(max_iterations=500, advanced=pf.nonlinear_advanced(f_delta_tol=0.0)),
    ))
    before = cfg.as_dict()
    assert _config(cfg) == before


def test_unset_solver_options_stay_unset():
    cfg = pf.config(solver=pf.solver(advanced=pf.solver_advanced(), nonlinear=pf.nonlinear()))
    assert _config(cfg)["solver"] == {"advanced": {}, "nonlinear": {}}


def test_backend_options_and_paths_are_not_rewritten():
    cfg = {"geometry": [{"mesh": "relative.msh"}], "root_path": "/tmp/config.json",
           "output": {"future_option": {"enabled": False}}, "common": "parameters.json"}
    original = copy.deepcopy(cfg)
    assert _config(cfg) == original
    assert cfg == original
