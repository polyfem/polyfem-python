from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_forward_solver_binding_targets_varform_state_only():
    source = (ROOT / "src" / "state" / "state.cpp").read_text(
        encoding="utf-8",
        errors="ignore",
    )

    assert "#include <polyfem/State.hpp>" in source
    assert "solve_problem" not in source
    assert ".build_basis(" not in source
    assert ".assemble_rhs(" not in source
    assert ".assemble_mass_mat(" not in source
    assert "s.assembler" not in source
    assert "s.bases" not in source
    assert "set_per_element_material" not in source
    assert 'bundle["sol"] = sol' in source
    assert 'bundle["u"]' not in source
    assert 'bundle["p"]' not in source
    assert 'bundle["vertices"]' not in source
    assert 'bundle["cells"]' not in source


def test_default_cpp_module_does_not_register_differentiable_bindings():
    src_cmake = (ROOT / "src" / "CMakeLists.txt").read_text()
    binding_cpp = (ROOT / "src" / "binding.cpp").read_text()

    assert "add_subdirectory(differentiable)" not in src_cmake
    assert "add_subdirectory(solver)" not in src_cmake
    assert "differentiable/binding.hpp" not in binding_cpp
    assert "solver/binding.hpp" not in binding_cpp
    assert "define_differentiable_cache" not in binding_cpp
    assert "define_adjoint" not in binding_cpp
    assert "define_objective" not in binding_cpp
    assert "define_opt_utils" not in binding_cpp
    assert "define_nonlinear_problem" not in binding_cpp
    assert not (ROOT / "src" / "solver" / "binding.hpp").exists()
    assert not (ROOT / "src" / "solver" / "nl_problem.cpp").exists()
    assert not (ROOT / "src" / "solver" / "CMakeLists.txt").exists()
