"""Configure the parent adapters without native dependencies or compilation."""

from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def write(root, name, content):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


@pytest.mark.parametrize("python_enabled", ["ON", "OFF"])
def test_parent_adapters_configure_without_backend_patches(tmp_path, python_enabled):
    cmake = shutil.which("cmake")
    if cmake is None:
        pytest.skip("CMake is required for adapter configuration tests")

    # Spaces exercise patch paths without relying on shell quoting.
    source = tmp_path / "parent with spaces"
    for name in ("polyfem.cmake", "libigl.cmake"):
        write(source, f"cmake/recipes/{name}",
              (ROOT / "cmake/recipes" / name).read_text(encoding="utf-8"))
    write(source, "CMakeLists.txt", """
cmake_minimum_required(VERSION 3.25)
project(AdapterTest NONE)
# These rules are only needed by generators that require a language rule.
# No build command is run.
set(CMAKE_CXX_COMPILE_OBJECT "unused")
set(CMAKE_CXX_ARCHIVE_CREATE "unused")
set(CMAKE_CXX_ARCHIVE_FINISH "unused")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes")
include(polyfem)
get_target_property(links polyfem INTERFACE_LINK_LIBRARIES)
get_target_property(definitions polyfem COMPILE_DEFINITIONS)
file(WRITE "${CMAKE_BINARY_DIR}/result.txt"
    "option=${POLYFEM_WITH_PYTHON}\nlinks=${links}\ndefinitions=${definitions}\n")
""")
    write(source, "polyfem/CMakeLists.txt", """
option(POLYFEM_WITH_PYTHON "Upstream interpreter support" ON)
if(POLYFEM_WITH_PYTHON)
    message(FATAL_ERROR "The backend must not configure interpreter embedding")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes")
include(libigl)
include(libigl)
add_library(polyfem STATIC dummy.cpp)
set_target_properties(polyfem PROPERTIES LINKER_LANGUAGE CXX)
add_library(polyfem::polyfem ALIAS polyfem)
""")
    write(source, "polyfem/dummy.cpp", "")
    write(source, "polyfem/cmake/recipes/libigl.cmake",
          'message(FATAL_ERROR "The parent libigl recipe was not selected")')
    write(source, "polyfem/cmake/patches/igl_first.patch", "test patch")
    write(source, "polyfem/cmake/patches/igl_second.patch", "test patch")
    write(source, "cmake/recipes/FindPython3.cmake", """
if(NOT Python3_FIND_COMPONENTS STREQUAL "Interpreter;Development.Module")
    message(FATAL_ERROR "Only interpreter and extension development are needed")
endif()
set(Python3_EXECUTABLE "test python/python.exe")
add_library(Python3::Module INTERFACE IMPORTED)
""")
    write(source, "polyfem/cmake/recipes/pybind11.cmake", """
if(NOT TARGET Python3::Module)
    message(FATAL_ERROR "Python extension support was not configured")
endif()
add_library(pybind11::module INTERFACE IMPORTED)
set_property(TARGET pybind11::module PROPERTY
    INTERFACE_LINK_LIBRARIES Python3::Module)
""")
    write(source, "cmake/recipes/CPM.cmake", """
function(CPMAddPackage)
    cmake_parse_arguments(PACKAGE "" "NAME;GITHUB_REPOSITORY;GIT_TAG" "PATCH_COMMAND" ${ARGN})
    file(APPEND "${CMAKE_BINARY_DIR}/dependency.txt"
        "${PACKAGE_NAME}|${PACKAGE_GITHUB_REPOSITORY}|${PACKAGE_GIT_TAG}\n")
    foreach(arg IN LISTS PACKAGE_PATCH_COMMAND)
        file(APPEND "${CMAKE_BINARY_DIR}/patch-args.txt" "${arg}\n")
    endforeach()
    add_library(igl::core INTERFACE IMPORTED GLOBAL)
endfunction()
""")
    build = tmp_path / "build"
    result = subprocess.run(
        [cmake, "-S", str(source), "-B", str(build),
         f"-DPOLYFEM_WITH_PYTHON={python_enabled}"],
        capture_output=True, text=True, timeout=45,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    config = (build / "result.txt").read_text()
    assert f"option={python_enabled}" in config
    assert "pybind11::embed" not in config
    assert "Python3::Python" not in config
    if python_enabled == "ON":
        assert "pybind11::module" in config
        assert "POLYFEM_WITH_PYTHON" in config
        assert 'POLYFEM_PYTHON_EXECUTABLE="test python/python.exe"' in config
    else:
        assert "pybind11::module" not in config
        assert "POLYFEM_WITH_PYTHON" not in config
    assert (build / "dependency.txt").read_text().splitlines() == [
        "libigl|libigl/libigl|aeeea9b416a5b5474cb9d12c07ff1238a4c83ba1"
    ]
    patch_args = (build / "patch-args.txt").read_text().splitlines()
    assert Path(patch_args[0]).name.lower() in ("git", "git.exe")
    assert patch_args[1:4] == ["apply", "--ignore-whitespace", "--whitespace=nowarn"]
    assert patch_args[4:] == [
        (source / f"polyfem/cmake/patches/igl_{name}.patch").as_posix()
        for name in ("first", "second")
    ]
