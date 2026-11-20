"""Protect against broken GitHub links and incomplete release archives."""

import pytest

from tools.check_distribution import BUILDER_FILES, PACKAGE_FILES, SOURCE_FILES, check_files


def source_files():
    return {f"python/{name}" for name in PACKAGE_FILES} | SOURCE_FILES | {f"python-from-jse/generator/{n}" for n in BUILDER_FILES}


def test_source_archive_requires_generator_config_and_builder():
    names = source_files()
    check_files(names, wheel=False)
    for name in ["generator-config/api_aliases.json", "python-from-jse/generator/model_builder.py"]:
        with pytest.raises(ValueError, match="Missing distribution files"):
            check_files(names - {name}, wheel=False)


@pytest.mark.parametrize("path", ["tests/test_cli.py", "build/native.obj", "polyfempy/__pycache__/a.pyc"])
def test_source_archive_rejects_development_files(path):
    with pytest.raises(ValueError, match="distribution"):
        check_files(source_files() | {path}, wheel=False)


def test_wheel_must_include_native_extension():
    names = PACKAGE_FILES | {"polyfempy-0.8.dist-info/licenses/LICENSE"}
    names |= {f"polyfempy/generated_api/generator/{n}" for n in BUILDER_FILES}
    with pytest.raises(ValueError, match="native"):
        check_files(names, wheel=True)
    check_files(names | {"polyfempy/polyfempy.cp311-win_amd64.pyd"}, wheel=True)


def test_source_archive_requires_python_source_layout():
    names = source_files()
    names.difference_update(f"python/{name}" for name in PACKAGE_FILES)
    names.update(PACKAGE_FILES)
    with pytest.raises(ValueError, match="Missing distribution files"):
        check_files(names, wheel=False)
