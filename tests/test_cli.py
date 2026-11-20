from __future__ import annotations

import subprocess
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_backend_info_reports_backend_status(capsys):
    from polyfempy.cli import main

    exit_code = main(["backend-info"])

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "backend_available=" in captured.out


def test_backend_info_can_require_backend(monkeypatch, capsys):
    from polyfempy import cli

    fake_polyfempy = types.SimpleNamespace(
        cpp_backend_available=lambda: False,
        cpp_backend_error=lambda: RuntimeError("not built"),
    )
    monkeypatch.setattr(cli, "_load_polyfempy", lambda: fake_polyfempy)

    exit_code = cli.main(["backend-info", "--require"])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "C++ backend not loaded" in captured.err


def test_solve_reports_missing_backend_cleanly(monkeypatch, capsys, tmp_path):
    from polyfempy import cli

    fake_polyfempy = types.SimpleNamespace(
        cpp_backend_available=lambda: False,
        cpp_backend_error=lambda: RuntimeError("not built"),
    )
    monkeypatch.setattr(cli, "_load_polyfempy", lambda: fake_polyfempy)
    config = tmp_path / "input.json"
    config.write_text("{}", encoding="utf-8")

    exit_code = cli.main(["solve", str(config)])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "C++ backend not loaded" in captured.err


def test_solve_dispatches_json_config_to_backend(monkeypatch, tmp_path):
    from polyfempy import cli

    calls = []

    def polyfem_command(**kwargs):
        calls.append(kwargs)

    fake_polyfempy = types.SimpleNamespace(
        cpp_backend_available=lambda: True,
        cpp_backend_error=lambda: None,
        polyfem_command=polyfem_command,
    )
    monkeypatch.setattr(cli, "_load_polyfempy", lambda: fake_polyfempy)

    config = tmp_path / "input.json"
    config.write_text("{}", encoding="utf-8")
    output_dir = tmp_path / "out"

    exit_code = cli.main(
        [
            "solve",
            str(config),
            "--output-dir",
            str(output_dir),
            "--log-level",
            "4",
            "--max-threads",
            "2",
            "--no-strict-validation",
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "json": str(config),
            "yaml": "",
            "log_level": 4,
            "strict_validation": False,
            "max_threads": 2,
            "output_dir": str(output_dir),
        }
    ]


def test_solve_dispatches_yaml_config_to_backend(monkeypatch, tmp_path):
    from polyfempy import cli

    calls = []

    def polyfem_command(**kwargs):
        calls.append(kwargs)

    fake_polyfempy = types.SimpleNamespace(
        cpp_backend_available=lambda: True,
        cpp_backend_error=lambda: None,
        polyfem_command=polyfem_command,
    )
    monkeypatch.setattr(cli, "_load_polyfempy", lambda: fake_polyfempy)

    config = tmp_path / "input.yaml"
    config.write_text("problem: test\n", encoding="utf-8")

    exit_code = cli.main(["solve", str(config)])

    assert exit_code == 0
    assert calls == [
        {
            "json": "",
            "yaml": str(config),
            "log_level": 2,
            "strict_validation": True,
            "max_threads": 1,
            "output_dir": "",
        }
    ]


def test_python_module_help_works_without_backend(monkeypatch):
    monkeypatch.setenv("PYTHONPATH", str(ROOT / "python"))
    result = subprocess.run(
        [sys.executable, "-m", "polyfempy", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "backend-info" in result.stdout
    assert "solve" in result.stdout


def test_pyproject_declares_console_script():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert "[project.scripts]" in pyproject
    assert 'polyfempy = "polyfempy.cli:main"' in pyproject
    assert "cmake>=3.25" in pyproject


def test_setup_py_always_passes_conda_prefix_to_cmake():
    setup_py = (ROOT / "setup.py").read_text(encoding="utf-8")

    conda_prefix_line = next(
        line for line in setup_py.splitlines() if "CONDA_PREFIX" in line
    )
    assert conda_prefix_line.startswith("        conda_prefix = ")


def test_setup_py_disables_heavy_optional_polysolve_dependencies():
    setup_py = (ROOT / "setup.py").read_text(encoding="utf-8")

    assert "-DPOLYSOLVE_WITH_SPECTRA=OFF" in setup_py
    assert "-DPOLYSOLVE_WITH_AMGCL=OFF" in setup_py
    assert "-DPOLYSOLVE_WITH_UMFPACK=OFF" in setup_py
    assert "-DPOLYSOLVE_WITH_MPI=OFF" in setup_py
    assert "-DPOLYSOLVE_WITH_HYPRE=OFF" in setup_py


def test_github_actions_generates_ignored_api_before_tests():
    workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8"
    )

    assert "matrix:" in workflow
    assert '"3.10"' in workflow
    assert '"3.11"' in workflow
    assert '"3.12"' in workflow
    assert "submodules: recursive" in workflow
    assert "python tools/generate_polyfem_api.py" in workflow
    assert "python -m pytest tests -q" in workflow


def test_github_actions_include_python_backend_and_release_workflows():
    workflow_dir = ROOT / ".github" / "workflows"
    workflow_names = {path.name for path in workflow_dir.glob("*.yml")}

    assert {
        "generated-contact-backend.yml",
        "release-packages.yml",
        "test.yml",
    } <= workflow_names


def test_backend_workflow_checks_main_and_ready_main_pull_requests():
    workflow = (
        ROOT / ".github" / "workflows" / "generated-contact-backend.yml"
    ).read_text(encoding="utf-8")

    assert "push:" in workflow
    assert "branches:" in workflow
    assert "- main" in workflow
    assert "- jingyao" not in workflow
    assert "pull_request:\n    branches: [main]" in workflow
    assert "types: [opened, synchronize, reopened, ready_for_review, edited]" in workflow
    assert "workflow_dispatch:" in workflow
    assert "GENERATED_SOURCE_TOLERANCE" in workflow
    assert "github.event.inputs.generated_source_tolerance || '1e-5'" in workflow
    assert "MAX_THREADS" in workflow
    assert "github.event.inputs.max_threads || '1'" in workflow
    assert "submodules: recursive" in workflow
    assert "python tools/generate_polyfem_api.py" in workflow
    assert "python -m build --no-isolation" in workflow
    assert "python -m pip install --force-reinstall --no-deps dist/*.whl" in workflow
    assert "python -I -m polyfempy backend-info --require" in workflow
    smoke, sweep = workflow.split("- name: Run active generated contact backend sweep")
    assert "tools/check_installed_package.py" in smoke
    assert "--readme README.md --reference polyfem-data/units/neohookean.json" in smoke
    assert "github.event.pull_request.draft == false" in workflow
    assert "github.event.action != 'edited' || github.event.changes.base != null" in workflow
    assert "if: always()" in sweep
    assert "tools/run_generated_contact_backend_checks.py" in sweep
    assert "--require-tests-match" in workflow
    assert "tools/generated_contact_expected_failures.json" in workflow
    assert "actions/upload-artifact@v7" in workflow


def test_github_actions_use_node24_compatible_action_versions():
    workflow_dir = ROOT / ".github" / "workflows"
    workflows = {
        path.name: path.read_text(encoding="utf-8")
        for path in workflow_dir.glob("*.yml")
    }

    assert workflows
    for workflow in workflows.values():
        assert "actions/checkout@v4" not in workflow
        assert "actions/setup-python@v5" not in workflow
        assert "actions/checkout@v7" in workflow
        assert "actions/setup-python@v7" in workflow
