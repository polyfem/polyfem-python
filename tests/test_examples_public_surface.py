from __future__ import annotations

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
CLASSIC_EXAMPLES = EXAMPLES / "classic_example"


def _require_examples_submodule() -> None:
    if not CLASSIC_EXAMPLES.exists():
        pytest.skip(
            "examples submodule is not initialized; run "
            "`git submodule update --init examples`"
        )


def _example_source_files() -> list[Path]:
    _require_examples_submodule()
    return sorted(path for path in EXAMPLES.rglob("*.py") if path.is_file())


def test_examples_do_not_reintroduce_removed_impact_template_name():
    offenders = []
    for path in _example_source_files():
        text = path.read_text(encoding="utf-8")
        if "make_impact_template" in text:
            offenders.append(path.relative_to(ROOT).as_posix())

    assert offenders == []


def test_classic_examples_use_generated_api_helpers_only():
    expected = {
        "examples/classic_example/2D/_contact_2d_common.py",
        "examples/classic_example/2D/contact_2d_friction_circle_rollers_generated_api.py",
        "examples/classic_example/2D/contact_2d_golf_ball_deformable_wall_generated_api.py",
        "examples/classic_example/2D/contact_2d_golf_ball_generated_api.py",
        "examples/classic_example/3D/_contact_3d_common.py",
        "examples/classic_example/3D/contact_3d_friction_high_school_slopetest_generated_api.py",
        "examples/classic_example/3D/contact_3d_golf_ball_generated_api.py",
        "examples/classic_example/3D/contact_3d_large_ratios_sphere_mat_generated_api.py",
    }
    actual = {path.relative_to(ROOT).as_posix() for path in _example_source_files()}
    assert expected <= actual


def test_classic_examples_do_not_import_legacy_api_package():
    offenders = []
    for path in _example_source_files():
        text = path.read_text(encoding="utf-8")
        if "polyfempy.api" in text:
            offenders.append(path.relative_to(ROOT).as_posix())

    assert offenders == []


def test_classic_examples_import_public_polyfem():
    _require_examples_submodule()
    for example in CLASSIC_EXAMPLES.glob("[23]D/*_generated_api.py"):
        lines = example.read_text(encoding="utf-8").splitlines()
        assert "import polyfem" in lines
        assert "from polyfempy.runtime import solve" in lines


def test_classic_examples_keep_solver_and_output_explicit():
    offenders = []
    for path in _example_source_files():
        text = path.read_text(encoding="utf-8")
        if "common_solver" in text or "common_output" in text:
            offenders.append(path.relative_to(ROOT).as_posix())

    assert offenders == []
