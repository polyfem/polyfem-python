"""Check the contents of a source archive or wheel without extracting it."""

from __future__ import annotations

import argparse
import glob
from pathlib import Path, PurePosixPath
import tarfile
import zipfile


PACKAGE_FILES = {
    "polyfempy/__init__.py", "polyfempy/__main__.py", "polyfempy/cli.py", "polyfempy/runtime/solve.py",
    "polyfempy/runtime/result.py", "polyfempy/runtime/_solve_contract.py",
    "polyfempy/runtime/__init__.py", "polyfempy/runtime/_solve_backend.py",
    "polyfempy/runtime/_runtime.py", "polyfempy/runtime/tensor.py",
    "polyfempy/generated_api/__init__.py",
    "polyfempy/generated_api/generated_api.py",
    "polyfempy/generated_api/generated_class.py",
    "polyfempy/generated_api/builder_api_manifest.json",
    "polyfempy/generated_api/class_tree_manifest.json",
    "polyfempy/generated_api/generated_api_manifest.json",
}
BUILDER_FILES = {"__init__.py", "model_builder.py", "id_relationships.py", "selection_refs.py"}
SOURCE_FILES = {
    "pyproject.toml", "setup.py", "CMakeLists.txt", "LICENSE", "README.md", "INSTALL.md", "API.md",
    "src/binding.cpp", "polyfem/CMakeLists.txt", "polyfem/LICENSE",
    "polyfem/json-specs/input-spec.json", "python-from-jse/LICENSE",
    "python/polyfem/__init__.py", "tools/check_installed_package.py",
    "examples/quickstart/beam.py", "examples/quickstart/beam.obj",
    "examples/quickstart/static_elasticity.py", "examples/quickstart/transient_elasticity.py",
}


def check_files(names: set[str], *, wheel: bool) -> None:
    for name in names:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe archive path: {name}")
        if set(path.parts) & {"__pycache__", ".git", ".pytest_cache", ".ruff_cache"}:
            raise ValueError(f"Cache or repository file in distribution: {name}")
        if path.parts[0] in {"tests", "build", "dist", "results"}:
            raise ValueError(f"Development or output file in distribution: {name}")
    builder = "polyfempy/generated_api/generator" if wheel else "python-from-jse/generator"
    package_files = PACKAGE_FILES if wheel else {f"python/{name}" for name in PACKAGE_FILES}
    required = package_files | {f"{builder}/{name}" for name in BUILDER_FILES}
    required |= {"polyfem/__init__.py"} if wheel else SOURCE_FILES
    missing = required - names
    if missing:
        raise ValueError("Missing distribution files: " + ", ".join(sorted(missing)))
    if wheel:
        if not any(n.startswith("polyfempy/polyfempy") and n.endswith((".so", ".pyd")) for n in names):
            raise ValueError("Wheel has no native polyfempy extension.")
        if not any(".dist-info/licenses/" in n and n.endswith("LICENSE") for n in names):
            raise ValueError("Wheel has no license files.")


def check_archive(path: Path) -> None:
    wheel = path.suffix == ".whl"
    if wheel:
        with zipfile.ZipFile(path) as archive:
            names = {m.filename for m in archive.infolist() if not m.is_dir()}
    else:
        with tarfile.open(path) as archive:
            members = archive.getmembers()
            roots = {PurePosixPath(m.name).parts[0] for m in members}
            if len(roots) != 1:
                raise ValueError("Source archive must have one top-level directory.")
            names = {m.name.split("/", 1)[1] for m in members if m.isfile()}
    check_files(names, wheel=wheel)
    print(f"Distribution files passed: {path.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", nargs="+", help="Archive paths or glob patterns")
    for pattern in parser.parse_args().archives:
        matches = glob.glob(pattern)
        if not matches:
            parser.error(f"No archives match {pattern}")
        for path in matches:
            check_archive(Path(path))


if __name__ == "__main__":
    main()
