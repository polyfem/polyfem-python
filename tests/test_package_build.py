"""Check build failures without compiling the native solver."""

import runpy
import subprocess
import sys
import sysconfig
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from setuptools import Distribution


@pytest.fixture
def build_script():
    with patch("setuptools.setup"):
        return runpy.run_path(str(Path(__file__).resolve().parents[1] / "setup.py"))


def test_missing_cmake_has_installation_hint(build_script):
    command = build_script["CMakeBuild"](Distribution())
    with patch("subprocess.check_output", side_effect=FileNotFoundError):
        with pytest.raises(RuntimeError, match="CMake 3.25 or later"):
            command.run()


@pytest.mark.parametrize("runtime", ["nanobind.dll", "libnanobind.so", "libnanobind.so.1", "libnanobind.dylib"])
def test_editable_mapping_includes_shared_runtime(build_script, tmp_path, runtime):
    from setuptools.command.build_ext import build_ext

    source = tmp_path / "python" / "polyfempy"
    output = tmp_path / "build" / "polyfempy"
    source.mkdir(parents=True)
    (source / runtime).write_bytes(b"runtime")
    (source / "debug.pdb").write_bytes(b"debug symbols")
    extension = "polyfempy.pyd"
    mapping = {str(output / extension): str(source / extension)}
    command = build_script["CMakeBuild"](Distribution())
    with patch.object(build_ext, "get_output_mapping", return_value=mapping.copy()):
        assert command.get_output_mapping() == {
            **mapping, str(output / runtime): str(source / runtime),
        }


@pytest.mark.parametrize(
    "system,generator",
    [("Windows", "Ninja"), ("Linux", "Ninja"), ("Windows", "Visual Studio 17 2022")],
)
def test_build_uses_cmake_parallel_options(
    build_script, monkeypatch, tmp_path, system, generator
):
    monkeypatch.setenv("CMAKE_GENERATOR", generator)
    monkeypatch.setenv("N_THREADS", "2")
    monkeypatch.setitem(
        sys.modules,
        "nanobind",
        SimpleNamespace(__file__=str(tmp_path / "nanobind/__init__.py")),
    )
    command = build_script["CMakeBuild"](Distribution())
    command.build_temp = str(tmp_path / "build")
    command.build_lib = str(tmp_path / "lib")
    with (
        patch("platform.system", return_value=system),
        patch("subprocess.check_call") as call,
    ):
        command.build_extension(build_script["CMakeExtension"]("polyfempy.polyfempy"))
    configure, build = [entry.args[0] for entry in call.call_args_list]
    assert ("-A" in configure) == (
        generator.startswith("Visual Studio") and sys.maxsize > 2**32
    )
    visual_studio = system == "Windows" and generator.startswith("Visual Studio")
    assert build[-3:] == ["--parallel", "1" if visual_studio else "2", "--verbose"]
    assert "--" not in build
    assert "-DPOLYFEM_PORTABLE_BUILD=ON" in configure
    assert f"-DPython_EXTENSION_MODULE_SUFFIX={sysconfig.get_config_var('EXT_SUFFIX')}" in configure
    if visual_studio:
        assert call.call_args.kwargs["env"]["_CL_"].endswith("/MP2")


def test_failed_build_is_not_repeated(build_script, monkeypatch, tmp_path):
    monkeypatch.setitem(
        sys.modules,
        "nanobind",
        SimpleNamespace(__file__=str(tmp_path / "nanobind/__init__.py")),
    )
    command = build_script["CMakeBuild"](Distribution())
    command.build_temp = str(tmp_path / "build")
    command.build_lib = str(tmp_path / "lib")
    with patch(
        "subprocess.check_call",
        side_effect=[0, subprocess.CalledProcessError(1, "cmake")],
    ) as call:
        with pytest.raises(RuntimeError, match="CMake build failed"):
            command.build_extension(
                build_script["CMakeExtension"]("polyfempy.polyfempy")
            )
    assert call.call_count == 2


def test_existing_ninja_build_does_not_receive_visual_studio_options(
    build_script, monkeypatch, tmp_path
):
    monkeypatch.delenv("CMAKE_GENERATOR", raising=False)
    monkeypatch.setenv("N_THREADS", "2")
    monkeypatch.setitem(
        sys.modules,
        "nanobind",
        SimpleNamespace(__file__=str(tmp_path / "nanobind/__init__.py")),
    )
    (tmp_path / "CMakeCache.txt").write_text("CMAKE_GENERATOR:INTERNAL=Ninja\n")
    command = build_script["CMakeBuild"](Distribution())
    command.build_temp = str(tmp_path)
    command.build_lib = str(tmp_path / "lib")
    with (
        patch("platform.system", return_value="Windows"),
        patch("subprocess.check_call") as call,
    ):
        command.build_extension(build_script["CMakeExtension"]("polyfempy.polyfempy"))
    configure, build = [entry.args[0] for entry in call.call_args_list]
    assert "-A" not in configure
    assert build[-3:] == ["--parallel", "2", "--verbose"]


def test_missing_generated_sources_fail_before_packaging(build_script, tmp_path):
    check = build_script["check_package_sources"]
    with patch.dict(check.__globals__, ROOT=tmp_path):
        with pytest.raises(RuntimeError, match="tools/generate_polyfem_api.py"):
            check()


@pytest.mark.parametrize("workers", ["0", "-1", "many"])
def test_invalid_worker_count_has_clear_error(
    build_script, monkeypatch, tmp_path, workers
):
    monkeypatch.setenv("N_THREADS", workers)
    command = build_script["CMakeBuild"](Distribution())
    command.build_lib = str(tmp_path)
    with pytest.raises(RuntimeError, match="N_THREADS must be a positive integer"):
        command.build_extension(build_script["CMakeExtension"]("polyfempy.polyfempy"))


def test_package_discovery_preserves_installed_names():
    root = Path(__file__).resolve().parents[1]
    with patch("setuptools.setup") as setup:
        runpy.run_path(str(root / "setup.py"))
    options = setup.call_args.kwargs
    packages = options["packages"]
    assert {"polyfem", "polyfempy", "polyfempy.runtime", "polyfempy.generated_api"} <= set(packages)
    assert not any(name.startswith("python.") for name in packages)
    distribution = Distribution({"packages": packages, "package_dir": options["package_dir"]})
    command = distribution.get_command_obj("build_py")
    command.ensure_finalized()
    for package in ("polyfem", "polyfempy", "polyfempy.runtime"):
        assert Path(command.get_package_dir(package)) == Path("python", *package.split("."))
    assert Path(command.get_package_dir("polyfempy.generated_api.generator")) == Path("python-from-jse/generator")
    assert options["ext_modules"][0].name == "polyfempy.polyfempy"
