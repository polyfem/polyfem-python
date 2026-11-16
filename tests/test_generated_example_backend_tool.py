from __future__ import annotations

import json
import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOL_PATH = ROOT / "tools" / "check_generated_example_backend.py"
BATCH_TOOL_PATH = ROOT / "tools" / "run_generated_contact_backend_checks.py"


def _load_tool_module():
    spec = importlib.util.spec_from_file_location(
        "check_generated_example_backend_for_test",
        TOOL_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        return module
    finally:
        sys.modules.pop(spec.name, None)


def _load_batch_tool_module():
    spec = importlib.util.spec_from_file_location(
        "run_generated_contact_backend_checks_for_test",
        BATCH_TOOL_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        return module
    finally:
        sys.modules.pop(spec.name, None)


def test_generated_example_backend_tool_exists_and_documents_manual_use():
    text = TOOL_PATH.read_text(encoding="utf-8")

    assert "generated example" in text
    assert "--example" in text
    assert "--source-json" in text
    assert "--require-tests-match" in text


def test_compare_metric_maps_reports_pass_and_fail_rows():
    tool = _load_tool_module()

    rows = tool.compare_metric_maps(
        left={"err_l2": 1.0, "err_h1": 2.0},
        right={"err_l2": 1.0 + 1e-12, "err_h1": 2.1},
        tolerance=1e-9,
    )

    by_key = {row.key: row for row in rows}
    assert by_key["err_l2"].passed
    assert not by_key["err_h1"].passed


def test_backend_tool_applies_polyfem_test_time_steps_to_tend_dt_config(tmp_path):
    tool = _load_tool_module()

    source_json = tmp_path / "example.json"
    source_json.write_text(
        json.dumps(
            {
                "time": {"tend": 2.0, "dt": 0.25},
                "tests": {"time_steps": 3},
            }
        ),
        encoding="utf-8",
    )

    payload, test_directive = tool.load_reduced_backend_payload(source_json)

    assert test_directive == 3
    assert payload["time"] == {"dt": 0.25, "time_steps": 3}


def test_backend_tool_applies_polyfem_test_time_steps_to_tend_steps_config(tmp_path):
    tool = _load_tool_module()

    source_json = tmp_path / "example.json"
    source_json.write_text(
        json.dumps(
            {
                "time": {"tend": 2.0, "time_steps": 8},
                "tests": {"time_steps": 2},
            }
        ),
        encoding="utf-8",
    )

    payload, test_directive = tool.load_reduced_backend_payload(source_json)

    assert test_directive == 2
    assert payload["time"] == {"dt": 0.25, "time_steps": 2}


def test_backend_tool_defaults_dynamic_examples_to_one_step(tmp_path):
    tool = _load_tool_module()

    source_json = tmp_path / "example.json"
    source_json.write_text(
        json.dumps({"time": {"dt": 0.1, "time_steps": 10}}),
        encoding="utf-8",
    )

    payload, test_directive = tool.load_reduced_backend_payload(source_json)

    assert test_directive == 1
    assert payload["time"] == {"dt": 0.1, "time_steps": 1}


def test_backend_tool_adds_json_output_to_reduced_payload(tmp_path):
    tool = _load_tool_module()

    source_json = tmp_path / "example.json"
    source_json.write_text(
        json.dumps({"time": {"dt": 0.1, "time_steps": 10}}),
        encoding="utf-8",
    )

    payload, _test_directive = tool.load_reduced_backend_payload(source_json)

    assert payload["output"]["json"] == "sim.json"


def test_backend_tool_applies_polyfem_verify_run_linear_solver_patch(tmp_path):
    tool = _load_tool_module()

    default_source = tmp_path / "example.json"
    default_source.write_text(
        json.dumps(
            {
                "solver": {"linear": {"solver": ["A", "B"]}},
                "time": {"dt": 0.1, "time_steps": 10},
            }
        ),
        encoding="utf-8",
    )
    navier_source = tmp_path / "navier-example.json"
    navier_source.write_text(
        json.dumps(
            {
                "solver": {"linear": {"solver": ["A", "B"]}},
                "time": {"dt": 0.1, "time_steps": 10},
            }
        ),
        encoding="utf-8",
    )

    default_payload, _ = tool.load_reduced_backend_payload(default_source)
    navier_payload, _ = tool.load_reduced_backend_payload(navier_source)

    assert (
        default_payload["solver"]["linear"]["solver"]
        == "Eigen::SimplicialLDLT"
    )
    assert navier_payload["solver"]["linear"]["solver"] == "Eigen::SparseLU"


def test_generated_run_isolates_output_from_example(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from polyfempy import runtime

    tool = _load_tool_module()
    cfg = {"output": {"directory": "results/sample", "log": {"path": "sample.log"}}}
    monkeypatch.setattr(tool, "_import_example", lambda path: SimpleNamespace(polyfem_config=cfg))
    captured = {}
    monkeypatch.setattr(runtime, "solve", lambda *, cfg: captured.update(cfg))
    workspace = tool.run_generated_example(
        Path("sample.py"), Path("sample.json"), tmp_path, test_time_steps="static", keep_visual_output=True
    )
    assert captured["output"] == {
        "directory": str(workspace), "json": "sim.json", "log": {"path": str(workspace / "polyfem.log")}
    }
    assert cfg == {"output": {"directory": "results/sample", "log": {"path": "sample.log"}}}


def test_backend_tool_restores_generated_payload_semantics_before_backend_run(tmp_path):
    tool = _load_tool_module()

    class GeneratedConfig:
        def as_dict(self):
            return {
                "geometry": [
                    {
                        "mesh": str(tmp_path / "mesh.obj"),
                        "transformation": {
                            "rotation": [90.0],
                            "scale": [0.5],
                            "translation": [],
                        },
                    }
                ],
            }

    payload = tool._config_payload_for_backend_run(GeneratedConfig())

    assert payload["geometry"][0]["transformation"] == {
        "rotation": 90.0,
        "scale": 0.5,
    }


def test_backend_tool_keeps_all_and_static_test_directives_unpatched(tmp_path):
    tool = _load_tool_module()

    all_source = tmp_path / "all.json"
    all_source.write_text(
        json.dumps(
            {
                "time": {"dt": 0.1, "time_steps": 10},
                "tests": {"time_steps": "all"},
            }
        ),
        encoding="utf-8",
    )
    static_source = tmp_path / "static.json"
    static_source.write_text(json.dumps({"tests": {}}), encoding="utf-8")

    all_payload, all_directive = tool.load_reduced_backend_payload(all_source)
    static_payload, static_directive = tool.load_reduced_backend_payload(static_source)

    assert all_directive == "all"
    assert all_payload["time"] == {"dt": 0.1, "time_steps": 10}
    assert static_directive == "static"
    assert "time" not in static_payload


def test_backend_tool_prunes_visual_output_files(tmp_path):
    tool = _load_tool_module()

    keep = tmp_path / "sim.json"
    keep.write_text("{}", encoding="utf-8")
    log = tmp_path / "polyfem.log"
    log.write_text("log", encoding="utf-8")
    for suffix in (".vtu", ".vtm", ".pvd"):
        (tmp_path / f"sim{suffix}").write_text("large", encoding="utf-8")

    removed = tool.prune_visual_outputs(tmp_path)

    assert sorted(path.suffix for path in removed) == [".pvd", ".vtm", ".vtu"]
    assert keep.exists()
    assert log.exists()
    assert not any(tmp_path.glob("*.vtu"))
    assert not any(tmp_path.glob("*.vtm"))
    assert not any(tmp_path.glob("*.pvd"))


def test_backend_tool_uses_short_backend_workspace_names(tmp_path):
    tool = _load_tool_module()

    example_path = (
        tmp_path
        / "contact_2d_friction_high_school_physics_slopetest_mu_0_49_generated_api.py"
    )
    source_json = tmp_path / "high-school-physics-slopetest-mu=0.49.json"
    output_root = tmp_path / "generated-contact-backend-check-20260813-model-full"

    generated_workspace = tool._generated_workspace(output_root, example_path)
    source_workspace = tool._source_json_workspace(output_root, source_json)

    assert generated_workspace == output_root / "generated" / "run"
    assert source_workspace == output_root / "source-json" / "run"
    assert example_path.stem not in str(generated_workspace)
    assert source_json.stem not in str(source_workspace)


def test_generated_contact_workflow_runs_batch_backend_checker():
    workflow = (
        ROOT / ".github" / "workflows" / "generated-contact-backend.yml"
    ).read_text(
        encoding="utf-8"
    )

    assert "tools/run_generated_contact_backend_checks.py" in workflow
    assert "--require-tests-match" in workflow
    assert "tools/generated_contact_expected_failures.json" in workflow


def test_batch_backend_tool_reads_polyfem_active_contact_lists():
    tool = _load_batch_tool_module()

    cases = tool.iter_active_contact_cases(tool.CONTACT_LISTS)
    by_source = {case.source_rel: case for case in cases}

    assert len(cases) == 75
    assert all(Path(case.generated_example).exists() for case in cases)
    assert "contact/examples/3D/pile/cubes.json" in by_source
    assert "contact/examples/3D/hex_tet/params.json" in by_source
    assert "contact/examples/3D/rigid/proxy/screw.json" in by_source
    assert "contact/examples/3D/static/two-cubes.json" not in by_source
    assert by_source[
        "contact/examples/3D/friction/high-school-physics-slopetest-mu=0.50.json"
    ].generated_example.endswith(
        "contact_3d_friction_high_school_slopetest_generated_api.py"
    )


def test_batch_backend_tool_writes_logs_outside_checkout(tmp_path, monkeypatch):
    tool = _load_batch_tool_module()
    case = tool.ContactCase(
        list_file="contact_2d.txt",
        line_number=1,
        source_rel="contact/examples/2D/example.json",
        source_json=str(tmp_path / "source.json"),
        generated_example=str(ROOT / "examples" / "example.py"),
    )
    completed = tool.subprocess.CompletedProcess(
        args=[], returncode=0, stdout="Reference metrics passed.", stderr=""
    )
    monkeypatch.setattr(tool.subprocess, "run", lambda *args, **kwargs: completed)

    result = tool.run_case(
        case,
        tmp_path,
        generated_source_tolerance=1e-5,
        require_tests_match=True,
        keep_visual_output=False,
        log_level=2,
        max_threads=1,
    )
    log_file = Path(result.log_file)
    assert log_file.is_absolute()
    assert "Reference metrics passed." in log_file.read_text(encoding="utf-8")
    tool.write_summaries(tmp_path, [case], [result])
    summary = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    assert summary["passed"] == 1
    assert summary["failed"] == 0


def test_batch_backend_tool_loads_expected_failure_config(tmp_path):
    tool = _load_batch_tool_module()

    config_path = tmp_path / "expected_failures.json"
    config_path.write_text(
        json.dumps(
            {
                "ignored": [
                    {
                        "source": "contact/examples/3D/rigid/proxy/screw.json",
                        "reason": "teacher-approved HDF5 issue",
                        "approved": True,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    expected_failures = tool.load_expected_failures(config_path)

    assert set(expected_failures) == {
        "contact/examples/3D/rigid/proxy/screw.json"
    }
    ignored = expected_failures["contact/examples/3D/rigid/proxy/screw.json"]
    assert ignored.source == "contact/examples/3D/rigid/proxy/screw.json"
    assert ignored.reason == "teacher-approved HDF5 issue"
    assert ignored.approved is True


def test_default_expected_failure_config_does_not_ignore_passing_screw_case():
    tool = _load_batch_tool_module()

    expected_failures = tool.load_expected_failures(tool.DEFAULT_EXPECTED_FAILURES)

    assert "contact/examples/3D/rigid/proxy/screw.json" not in expected_failures


def test_batch_backend_tool_classifies_expected_failures_in_summary(tmp_path):
    tool = _load_batch_tool_module()
    expected_failures = {
        "contact/examples/3D/rigid/proxy/screw.json": tool.ExpectedFailure(
            source="contact/examples/3D/rigid/proxy/screw.json",
            reason="teacher-approved HDF5 issue",
            approved=True,
        ),
        "contact/examples/3D/expected/pass-now.json": tool.ExpectedFailure(
            source="contact/examples/3D/expected/pass-now.json",
            reason="old expected failure",
            approved=True,
        ),
    }
    raw_results = [
        tool.CaseResult(
            source_rel="contact/examples/2D/friction/card-house.json",
            generated_example="examples/card-house.py",
            status="PASS",
            returncode=0,
            log_file="logs/card-house.log",
        ),
        tool.CaseResult(
            source_rel="contact/examples/3D/rigid/proxy/screw.json",
            generated_example="examples/screw.py",
            status="FAIL",
            returncode=1,
            log_file="logs/screw.log",
        ),
        tool.CaseResult(
            source_rel="contact/examples/3D/expected/pass-now.json",
            generated_example="examples/pass-now.py",
            status="PASS",
            returncode=0,
            log_file="logs/pass-now.log",
        ),
        tool.CaseResult(
            source_rel="contact/examples/3D/pile/cubes.json",
            generated_example="examples/cubes.py",
            status="FAIL",
            returncode=1,
            log_file="logs/cubes.log",
        ),
    ]

    results = tool.apply_expected_failures(raw_results, expected_failures)

    assert [result.status for result in results] == [
        "PASS",
        "IGNORED",
        "UNEXPECTED_PASS",
        "FAIL",
    ]
    assert results[1].raw_status == "FAIL"
    assert results[1].reason == "teacher-approved HDF5 issue"
    assert results[2].raw_status == "PASS"
    assert results[2].reason == "old expected failure"

    tool.write_summaries(tmp_path, [], results)

    summary = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    assert summary["passed"] == 1
    assert summary["ignored"] == 1
    assert summary["failed"] == 1
    assert summary["unexpected_pass"] == 1
    assert summary["unexpected_fail"] == 1

    summary_text = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "PASS:            1" in summary_text
    assert "IGNORED:         1" in summary_text
    assert "FAIL:            1" in summary_text
    assert "Unexpected pass: 1" in summary_text
    assert "Unexpected fail: 1" in summary_text
