"""Generate PolyFEM's packaged Python API from the repo root."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import urllib.request


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR_ROOT = REPO_ROOT / "python-from-jse"
GENERATOR_CONFIG_DIR = REPO_ROOT / "generator-config"
DEFAULT_POLYFEM_SOURCE_DIR = REPO_ROOT / "polyfem"
POLYFEMPY_PACKAGE_DIR = REPO_ROOT / "python" / "polyfempy"
GENERATED_DIR = POLYFEMPY_PACKAGE_DIR / "generated_api"
POLYSOLVE_SPEC_CACHE_DIR = REPO_ROOT / "build" / "polysolve-json-specs"
LINKED_SOLVER_SPEC_FILES = (
    "linear-solver-spec.json",
    "nonlinear-solver-spec.json",
)


def polyfem_schema_file(polyfem_source_dir: Path) -> Path:
    return polyfem_source_dir / "json-specs" / "input-spec.json"


def workflow_steps(
    run_checks: bool,
    *,
    schema_file: Path,
    include_spec_dirs: list[Path] | None = None,
) -> list[tuple[list[str], Path, dict[str, str] | None]]:
    include_spec_dirs = include_spec_dirs or []
    generator_command = [
        sys.executable,
        str(GENERATOR_ROOT / "tools" / "generate_with_overrides.py"),
        "--schema-file",
        str(schema_file),
    ]
    for include_spec_dir in include_spec_dirs:
        generator_command.extend(["--include-spec-dir", str(include_spec_dir)])
    generator_command.extend([
        "--output-file",
        str(GENERATED_DIR / "generated_class.py"),
        "--api-output-file",
        str(GENERATED_DIR / "generated_api.py"),
        "--manifest-dir",
        str(GENERATED_DIR),
        "--relationships",
        str(GENERATOR_CONFIG_DIR / "id_relationships.json"),
        "--api-aliases",
        str(GENERATOR_CONFIG_DIR / "api_aliases.json"),
        "--model-entry",
        "polyfem.model",
    ])
    steps: list[tuple[list[str], Path, dict[str, str] | None]] = [
        (
            generator_command,
            REPO_ROOT,
            None,
        ),
    ]
    if run_checks:
        steps.extend([
            (
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    "tests/test_examples_public_surface.py",
                    "tests/test_generated_api_example.py",
                ],
                REPO_ROOT,
                None,
            ),
        ])
    return steps


def missing_required_paths(
    schema_file: Path,
    include_spec_dirs: list[Path],
    *,
    require_linked_specs: bool = True,
) -> list[Path]:
    required = [
        GENERATOR_ROOT,
        GENERATOR_ROOT / "tools" / "generate_with_overrides.py",
        schema_file,
        GENERATOR_CONFIG_DIR,
        GENERATOR_CONFIG_DIR / "api_aliases.json",
        GENERATOR_CONFIG_DIR / "id_relationships.json",
        POLYFEMPY_PACKAGE_DIR,
    ]
    for include_spec_dir in include_spec_dirs:
        required.extend([
            include_spec_dir,
        ])
    if require_linked_specs:
        search_dirs = [schema_file.parent, *include_spec_dirs]
        for spec_file in LINKED_SOLVER_SPEC_FILES:
            if not any((search_dir / spec_file).exists() for search_dir in search_dirs):
                required.append(schema_file.parent / spec_file)
    return [path for path in required if not path.exists()]


def _linked_solver_specs_are_available(search_dirs: list[Path]) -> bool:
    return all(
        any((search_dir / spec_file).exists() for search_dir in search_dirs)
        for spec_file in LINKED_SOLVER_SPEC_FILES
    )


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


def resolve_include_spec_dirs(
    schema_file: Path,
    include_spec_dirs: list[Path],
    polyfem_source_dir: Path,
) -> list[Path]:
    """Return include spec dirs, auto-caching PolySolve linked specs if needed."""
    resolved_include_dirs = list(include_spec_dirs)
    search_dirs = [schema_file.parent, *resolved_include_dirs]
    if _linked_solver_specs_are_available(search_dirs):
        return resolved_include_dirs

    cache_dir = _cache_polysolve_solver_specs(polyfem_source_dir)
    if cache_dir is not None:
        resolved_include_dirs.append(cache_dir)

    return resolved_include_dirs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate python/polyfempy/generated_api from python-from-jse and "
            "generator-config."
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run backend-free generated API and example parity checks.",
    )
    parser.add_argument(
        "--polyfem-source-dir",
        default=DEFAULT_POLYFEM_SOURCE_DIR,
        type=Path,
        help=(
            "Path to a PolyFEM source checkout. Defaults to the "
            "polyfem submodule."
        ),
    )
    parser.add_argument(
        "--resolved-schema-file",
        default=None,
        type=Path,
        help=(
            "Path to a complete PolyFEM schema JSON dumped from the PolyFEM "
            "build/embedded spec. When set, automatic PolySolve linked-spec "
            "cache resolution is skipped."
        ),
    )
    parser.add_argument(
        "--include-spec-dir",
        action="append",
        default=[],
        type=Path,
        help=(
            "Additional directory to search for included spec files. "
            "Can be passed more than once."
        ),
    )
    args = parser.parse_args(argv)

    if args.resolved_schema_file is not None:
        schema_file = args.resolved_schema_file
        include_spec_dirs = list(args.include_spec_dir)
        require_linked_specs = False
    else:
        schema_file = polyfem_schema_file(args.polyfem_source_dir)
        include_spec_dirs = resolve_include_spec_dirs(
            schema_file,
            args.include_spec_dir,
            args.polyfem_source_dir,
        )
        require_linked_specs = True

    missing = missing_required_paths(
        schema_file,
        include_spec_dirs,
        require_linked_specs=require_linked_specs,
    )
    if missing:
        for path in missing:
            print(f"Missing required path: {path}", file=sys.stderr)
        if args.resolved_schema_file is not None:
            print(
                "If using --resolved-schema-file, build PolyFEM and dump the "
                "embedded/resolved spec first, then pass that JSON path.",
                file=sys.stderr,
            )
        else:
            print(
                "If this is a fresh checkout, run "
                "`git submodule update --init --recursive` so polyfem "
                "and python-from-jse are present. "
                "If the missing file is a linked solver spec, pass "
                "`--include-spec-dir <dir-containing-linked-specs>`.",
                file=sys.stderr,
            )
        return 1

    for command, cwd, extra_env in workflow_steps(
        args.check,
        schema_file=schema_file,
        include_spec_dirs=include_spec_dirs,
    ):
        print("+ " + " ".join(command), flush=True)
        env = None
        if extra_env is not None:
            env = os.environ.copy()
            env.update(extra_env)
        result = subprocess.run(command, cwd=cwd, env=env)
        if result.returncode != 0:
            return result.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
