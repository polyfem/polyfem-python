"""Protect against broken GitHub links and incomplete release archives."""

from pathlib import Path
import subprocess

import pytest

from tools.check_docs import check_markdown, headings
from tools.check_distribution import BUILDER_FILES, PACKAGE_FILES, SOURCE_FILES, check_files


def test_submodule_file_exists_but_relative_github_link_is_rejected(tmp_path):
    examples = tmp_path / "examples"
    examples.mkdir()
    (examples / "beam.py").write_text("pass\n")
    readme = tmp_path / "README.md"
    readme.write_text("[Beam](examples/beam.py)\n")
    with pytest.raises(ValueError, match="pinned GitHub URL"):
        check_markdown(readme, tmp_path, [examples])


def test_local_anchor_and_python_syntax_are_checked(tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text("# Run a model\n[Run](#run-a-model)\n```python\nprint(1)\n```\n")
    assert check_markdown(readme, tmp_path, []) == set()
    readme.write_text("[Run](#missing)\n")
    with pytest.raises(ValueError, match="heading is missing"):
        check_markdown(readme, tmp_path, [])
    readme.write_text("```python\nif:\n```\n")
    with pytest.raises(SyntaxError):
        check_markdown(readme, tmp_path, [])


def test_repeated_headings_and_code_examples():
    assert headings("# Setup\n# Setup\n```sh\n# Not a heading\n```\n") == {"setup", "setup-1"}


def source_files():
    return {f"python/{name}" for name in PACKAGE_FILES} | SOURCE_FILES | {f"python-from-jse/generator/{n}" for n in BUILDER_FILES}


def test_source_archive_requires_quickstart_mesh_and_builder():
    names = source_files()
    check_files(names, wheel=False)
    for name in ["examples/quickstart/beam.obj", "python-from-jse/generator/model_builder.py"]:
        with pytest.raises(ValueError, match="Missing distribution files"):
            check_files(names - {name}, wheel=False)


@pytest.mark.parametrize("path", ["tests/test_cli.py", "build/native.obj", "polyfempy/__pycache__/a.pyc"])
def test_source_archive_rejects_development_files(path):
    with pytest.raises(ValueError, match="distribution"):
        check_files(source_files() | {path}, wheel=False)


def test_wheel_must_include_native_extension():
    names = PACKAGE_FILES | {"polyfem/__init__.py", "polyfempy-0.8.dist-info/licenses/LICENSE"}
    names |= {f"polyfempy/generated_api/generator/{n}" for n in BUILDER_FILES}
    with pytest.raises(ValueError, match="native"):
        check_files(names, wheel=True)
    check_files(names | {"polyfempy/polyfempy.cp311-win_amd64.pyd"}, wheel=True)


@pytest.mark.parametrize("path, ignored", [
    ("docs/guide.md", False), ("tests/fixtures/mesh.obj", False),
    ("tests/fixtures/reference.csv", False), ("tests/fixtures/reference.npz", False),
    ("build/native.obj", True), ("dist/package.whl", True),
    ("results/solution.vtu", True), (".venv/pyvenv.cfg", True),
])
def test_repository_ignore_rules_preserve_source_data(path, ignored):
    result = subprocess.run(
        ["git", "check-ignore", "--no-index", "-q", path],
        cwd=Path(__file__).resolve().parents[1], check=False,
    )
    assert result.returncode == (0 if ignored else 1)


def test_source_archive_requires_python_source_layout():
    names = source_files()
    names.difference_update(f"python/{name}" for name in PACKAGE_FILES)
    names.update(PACKAGE_FILES)
    with pytest.raises(ValueError, match="Missing distribution files"):
        check_files(names, wheel=False)
