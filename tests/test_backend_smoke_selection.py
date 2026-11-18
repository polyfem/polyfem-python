"""Check smoke-test backend selection without compiling native code."""

import hashlib
import importlib.util
from importlib.machinery import EXTENSION_SUFFIXES
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def smoke():
    path = Path(__file__).with_name("test_backend_smoke.py")
    spec = importlib.util.spec_from_file_location("backend_smoke_for_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_package(root, marker="selected", *, native_path=None, available=True):
    package = root / "polyfempy"
    package.mkdir(parents=True)
    native = native_path or package / ("polyfempy" + EXTENSION_SUFFIXES[0])
    native.write_bytes(b"native test fixture")
    (package / "__init__.py").write_text(
        "from types import SimpleNamespace\n"
        f"_core = SimpleNamespace(__file__={str(native)!r})\n"
        f"def cpp_backend_available(): return {available!r}\n"
        "def cpp_backend_error(): return 'test backend unavailable'\n",
        encoding="utf-8",
    )
    (package / "cli.py").write_text(
        "def main(args):\n"
        f"    print({marker!r})\n"
        "    return 0\n",
        encoding="utf-8",
    )
    return native


def test_selected_build_wins_over_pythonpath_and_working_directory(smoke, tmp_path, monkeypatch):
    selected = tmp_path / "current build"
    native = _fake_package(selected)
    old_env = tmp_path / "old-env"
    old_cwd = tmp_path / "old-cwd"
    _fake_package(old_env, "OLD_ENV_USED")
    _fake_package(old_cwd, "OLD_CWD_USED")
    monkeypatch.setenv("PYTHONPATH", str(old_env))

    probe = smoke._run_backend(selected, [], cwd=old_cwd)
    assert probe.returncode == 0, probe.stderr
    info = json.loads(probe.stdout)
    assert Path(info["native"]) == native.resolve()
    assert info["sha256"] == hashlib.sha256(native.read_bytes()).hexdigest()

    solve = smoke._run_backend(selected, ["solve"], cwd=old_cwd, expected_backend=info)
    assert solve.returncode == 0, solve.stderr
    assert "selected" in solve.stdout
    assert "OLD_ENV_USED" not in solve.stdout
    assert "OLD_CWD_USED" not in solve.stdout


def test_missing_selection_skips_without_probing_a_backend(smoke, tmp_path, monkeypatch):
    def unexpected_probe(*args, **kwargs):
        raise AssertionError("A backend must not be selected automatically.")

    monkeypatch.setattr(smoke, "_run_backend", unexpected_probe)
    config = SimpleNamespace(getoption=lambda name: None)
    with pytest.raises(pytest.skip.Exception, match="--backend-package-dir"):
        smoke.test_backend_forward_solve_smoke(tmp_path, config)


def test_missing_selected_package_fails_instead_of_falling_back(smoke, tmp_path, monkeypatch):
    _fake_package(tmp_path / "old")
    monkeypatch.setenv("PYTHONPATH", str(tmp_path / "old"))
    config = SimpleNamespace(getoption=lambda name: str(tmp_path / "missing"))

    with pytest.raises(AssertionError, match="Selected package is missing"):
        smoke.test_backend_forward_solve_smoke(tmp_path, config)


def test_unavailable_selected_backend_fails_instead_of_skipping(smoke, tmp_path):
    selected = tmp_path / "selected"
    _fake_package(selected, available=False)
    config = SimpleNamespace(getoption=lambda name: str(selected))

    with pytest.raises(AssertionError, match="Native backend is unavailable"):
        smoke.test_backend_forward_solve_smoke(tmp_path, config)


def test_backend_outside_selected_package_is_rejected(smoke, tmp_path):
    selected = tmp_path / "selected"
    old_native = tmp_path / ("old" + EXTENSION_SUFFIXES[0])
    _fake_package(selected, native_path=old_native)

    result = smoke._run_backend(selected, [], cwd=tmp_path)

    assert result.returncode != 0
    assert "Loaded backend outside selected package" in result.stderr


def test_python_module_cannot_stand_in_for_native_backend(smoke, tmp_path):
    selected = tmp_path / "selected"
    _fake_package(selected, native_path=selected / "polyfempy" / "polyfempy.py")

    result = smoke._run_backend(selected, [], cwd=tmp_path)

    assert result.returncode != 0
    assert "Expected a compiled native extension" in result.stderr


def test_changed_binary_is_rejected_before_solve(smoke, tmp_path):
    selected = tmp_path / "selected"
    native = _fake_package(selected, "SOLVE_MUST_NOT_RUN")
    probe = smoke._run_backend(selected, [], cwd=tmp_path)
    assert probe.returncode == 0, probe.stderr
    info = json.loads(probe.stdout)
    native.write_bytes(b"a different native build")

    result = smoke._run_backend(selected, ["solve"], cwd=tmp_path, expected_backend=info)

    assert result.returncode != 0
    assert "Native backend changed between probe and solve" in result.stderr
    assert "SOLVE_MUST_NOT_RUN" not in result.stdout
