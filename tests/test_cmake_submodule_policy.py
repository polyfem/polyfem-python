from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def test_cmake_defaults_to_polyfem_data_submodule():
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")

    assert '${CMAKE_CURRENT_SOURCE_DIR}/polyfem-data' in cmake
    assert '${CMAKE_CURRENT_SOURCE_DIR}/data/' not in cmake
    assert "file(MAKE_DIRECTORY ${POLYFEMPY_DATA_ROOT})" not in cmake


def test_polyfem_data_recipe_does_not_fetch_data_repo():
    recipe = (ROOT / "cmake" / "recipes" / "polyfem_data.cmake").read_text(
        encoding="utf-8"
    )

    assert "FetchContent" not in recipe
    assert "GIT_REPOSITORY https://github.com/polyfem/polyfem-data" not in recipe
    assert "git submodule update --init --recursive" in recipe
    assert "FATAL_ERROR" in recipe


def test_polyfem_recipe_uses_polyfem_submodule_instead_of_fetching_repo():
    recipe = (ROOT / "cmake" / "recipes" / "polyfem.cmake").read_text(
        encoding="utf-8"
    )

    assert "CPMAddPackage" not in recipe
    assert "FetchContent" not in recipe
    assert "GIT_REPOSITORY https://github.com/polyfem/polyfem" not in recipe
    assert "/polyfem" in recipe
    assert "add_subdirectory" in recipe
    assert "git submodule update --init --recursive" in recipe
    assert "FATAL_ERROR" in recipe


def test_root_cpm_version_matches_polyfem_submodule_cpm_version():
    root_cpm = (ROOT / "cmake" / "recipes" / "CPM.cmake").read_text(
        encoding="utf-8"
    )
    polyfem_cpm = (ROOT / "polyfem" / "cmake" / "recipes" / "CPM.cmake").read_text(
        encoding="utf-8"
    )

    root_version = re.search(r"CPM_DOWNLOAD_VERSION ([0-9.]+)", root_cpm)
    polyfem_version = re.search(r"CPM_DOWNLOAD_VERSION ([0-9.]+)", polyfem_cpm)

    assert root_version is not None
    assert polyfem_version is not None
    assert root_version.group(1) == polyfem_version.group(1)
