"""Generate the packaged configuration API; CMake runs this during installation."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re
import subprocess
import sys
import urllib.request

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR_ROOT = REPO_ROOT / "python-from-jse"
GENERATOR_CONFIG_DIR = REPO_ROOT / "generator-config"
GENERATED_DIR = REPO_ROOT / "python/polyfempy/generated_api"
POLYSOLVE_SPEC_CACHE_DIR = REPO_ROOT / "build/polysolve-json-specs"
LINKED_SOLVER_SPEC_FILES = ("linear-solver-spec.json", "nonlinear-solver-spec.json")


def _polysolve_revision_from_recipe(polyfem_source_dir: Path) -> str | None:
    recipe = polyfem_source_dir / "cmake" / "recipes" / "polysolve.cmake"
    if not recipe.exists():
        return None

    match = re.search(
        r"polyfem/polysolve#([A-Za-z0-9_.-]+)",
        recipe.read_text(encoding="utf-8"),
    )
    if match is None:
        return None
    return match.group(1)


def _cache_polysolve_solver_specs(polyfem_source_dir: Path) -> Path | None:
    revision = _polysolve_revision_from_recipe(polyfem_source_dir)
    if revision in (None, '.', '..'):
        return None

    cache_dir = POLYSOLVE_SPEC_CACHE_DIR / revision
    cache_dir.mkdir(parents=True, exist_ok=True)
    for spec_file in LINKED_SOLVER_SPEC_FILES:
        target = cache_dir / spec_file
        if target.exists():
            continue

        url = (
            "https://raw.githubusercontent.com/polyfem/polysolve/"
            f"{revision}/{spec_file}"
        )
        try:
            urllib.request.urlretrieve(url, target)
        except Exception as exc:
            target.unlink(missing_ok=True)
            print(
                f"Warning: could not download {spec_file} from polysolve "
                f"revision {revision}: {exc}",
                file=sys.stderr,
            )
            return None

    return cache_dir


def write_public_exports() -> None:
    api_file = GENERATED_DIR / "generated_api.py"
    source = api_file.read_text(encoding="utf-8").split("\n__all__ = ", 1)[0]
    exports = ["Root", *(
        node.name for node in ast.parse(source).body
        if isinstance(node, (ast.FunctionDef, ast.ClassDef))
        and not node.name.startswith("_")
    )]
    api_file.write_text(source + "\n__all__ = " + repr(exports) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polyfem-source-dir", type=Path)
    parser.add_argument("--resolved-schema-file", type=Path)
    parser.add_argument("--include-spec-dir", type=Path, action="append", default=[])
    parser.add_argument("--check", action="store_true", help="Run example parity tests after generation")
    args = parser.parse_args(argv)
    if args.resolved_schema_file is None and args.polyfem_source_dir is None:
        parser.error("Install with pip to generate automatically, or pass --polyfem-source-dir")
    schema = args.resolved_schema_file or args.polyfem_source_dir / "json-specs/input-spec.json"
    generator = GENERATOR_ROOT / "tools/generate_with_overrides.py"
    for path in (schema, generator):
        if not path.is_file():
            parser.error(f"Missing {path}; initialize the pinned submodules with git submodule update --init --recursive")
    includes = list(args.include_spec_dir)
    search = [schema.parent, *includes]
    if args.resolved_schema_file is None and not all(
        any((folder / name).is_file() for folder in search)
        for name in LINKED_SOLVER_SPEC_FILES
    ):
        cache = _cache_polysolve_solver_specs(args.polyfem_source_dir)
        if cache is None:
            parser.error("Linked solver schemas are missing; pass --include-spec-dir with their directory")
        includes.append(cache)
    command = [sys.executable, str(generator), "--schema-file", str(schema),
               "--output-file", str(GENERATED_DIR / "generated_class.py"),
               "--api-output-file", str(GENERATED_DIR / "generated_api.py"),
               "--manifest-dir", str(GENERATED_DIR),
               "--relationships", str(GENERATOR_CONFIG_DIR / "id_relationships.json"),
               "--api-aliases", str(GENERATOR_CONFIG_DIR / "api_aliases.json"),
               "--model-entry", "polyfempy.model"]
    for include in includes:
        command.extend(["--include-spec-dir", str(include)])
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    write_public_exports()
    if args.check:
        subprocess.run([sys.executable, "-m", "pytest", "tests/test_generated_api_example.py"],
                       cwd=REPO_ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
