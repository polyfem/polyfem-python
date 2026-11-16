"""Compare a generated example backend run against its source JSON run.

This is a manual diagnostic tool for generated example backend coverage. It
does two comparisons:

1. generated example output vs. source JSON output
2. generated example output vs. the source JSON ``tests`` block

By default only the first comparison is required to pass. Use
``--require-tests-match`` after the PolyFEM backend and polyfem-data tests are
known to be synchronized.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXAMPLE = (
    ROOT
    / "examples"
    / "classic_example"
    / "2D"
    / "contact_2d_golf_ball_deformable_wall_generated_api.py"
)
DEFAULT_SOURCE_JSON = (
    ROOT
    / "polyfem-data"
    / "contact"
    / "examples"
    / "2D"
    / "golf-ball-doformable-wall.json"
)
DEFAULT_OUTPUT_ROOT = ROOT / "build" / "generated-example-backend-check"
VISUAL_OUTPUT_SUFFIXES = {".pvd", ".vtm", ".vtu"}


@dataclass(frozen=True)
class MetricDiff:
    key: str
    left: float
    right: float
    diff: float
    tolerance: float

    @property
    def passed(self) -> bool:
        return self.diff <= self.tolerance


def _error_metrics(payload: dict[str, Any]) -> dict[str, float]:
    return {
        key: float(value)
        for key, value in payload.items()
        if key.startswith("err_") and isinstance(value, (int, float))
    }


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise TypeError(f"expected JSON object in {path}")
    return payload


def _deep_merge(base: Any, override: Any) -> Any:
    if isinstance(base, dict) and isinstance(override, dict):
        merged = copy.deepcopy(base)
        for key, value in override.items():
            merged[key] = (
                _deep_merge(merged[key], value)
                if key in merged
                else copy.deepcopy(value)
            )
        return merged
    return copy.deepcopy(override)


def _json_pointer_tokens(pointer: str) -> list[str]:
    if not pointer.startswith("/"):
        raise ValueError(f"JSON patch path must start with '/': {pointer!r}")
    return [
        token.replace("~1", "/").replace("~0", "~")
        for token in pointer[1:].split("/")
    ]


def _replace_json_pointer(payload: Any, pointer: str, value: Any) -> None:
    target = payload
    tokens = _json_pointer_tokens(pointer)
    for token in tokens[:-1]:
        target = target[int(token)] if isinstance(target, list) else target[token]

    last = tokens[-1]
    if isinstance(target, list):
        target[int(last)] = copy.deepcopy(value)
    else:
        target[last] = copy.deepcopy(value)


def _apply_json_patch(
    payload: dict[str, Any],
    patch: list[dict[str, Any]],
) -> dict[str, Any]:
    patched = copy.deepcopy(payload)
    for item in patch:
        if item.get("op") != "replace":
            raise ValueError(f"unsupported JSON patch op: {item!r}")
        _replace_json_pointer(patched, item["path"], item.get("value"))
    return patched


def _load_effective_source_payload(source_json: Path) -> dict[str, Any]:
    source_json = source_json.resolve()
    source = _load_json(source_json)
    common_ref = source.pop("common", None)
    patch = source.pop("patch", None)

    if common_ref is None:
        payload = source
    else:
        common_path = (source_json.parent / common_ref).resolve()
        payload = _deep_merge(_load_effective_source_payload(common_path), source)

    if patch:
        payload = _apply_json_patch(payload, patch)
    return payload


def _resolve_mesh_paths(value: Any, source_dir: Path) -> None:
    if isinstance(value, list):
        for item in value:
            _resolve_mesh_paths(item, source_dir)
        return

    if not isinstance(value, dict):
        return

    mesh = value.get("mesh")
    if isinstance(mesh, str):
        value["mesh"] = str((source_dir / mesh).resolve())

    for child in value.values():
        _resolve_mesh_paths(child, source_dir)


def verify_run_linear_solver(source_json: Path) -> str:
    source_text = str(source_json).lower()
    if any(
        marker in source_text
        for marker in ("navier", "bilaplace", "thermoelastic")
    ):
        return "Eigen::SparseLU"
    return "Eigen::SimplicialLDLT"


def apply_verify_run_linear_solver(payload: dict[str, Any], source_json: Path) -> None:
    solver = payload.setdefault("solver", {})
    if not isinstance(solver, dict):
        solver = {}
        payload["solver"] = solver

    linear = solver.setdefault("linear", {})
    if not isinstance(linear, dict):
        linear = {}
        solver["linear"] = linear

    linear["solver"] = verify_run_linear_solver(source_json)


def _test_time_steps(payload: dict[str, Any]) -> int | str:
    if "time" not in payload:
        return "static"

    tests = payload.get("tests")
    if not isinstance(tests, dict) or "time_steps" not in tests:
        return 1

    time_steps = tests["time_steps"]
    if isinstance(time_steps, bool):
        raise TypeError("tests.time_steps must be an integer, 'all', or 'static'")
    if isinstance(time_steps, (int, float)):
        return int(time_steps)
    if time_steps in {"all", "static"}:
        return time_steps
    raise TypeError(f"unsupported tests.time_steps value: {time_steps!r}")


def apply_test_time_steps(
    payload: dict[str, Any],
    time_steps: int | str,
) -> dict[str, Any]:
    reduced = copy.deepcopy(payload)
    if not isinstance(time_steps, int):
        return reduced

    time_args = reduced.get("time")
    if not isinstance(time_args, dict):
        raise ValueError("numeric tests.time_steps requires a time block")

    time_args = copy.deepcopy(time_args)
    if "tend" in time_args and "dt" in time_args:
        time_args.pop("tend", None)
        time_args["time_steps"] = time_steps
    elif "tend" in time_args and "time_steps" in time_args:
        original_steps = int(time_args["time_steps"])
        time_args["dt"] = float(time_args["tend"]) / original_steps
        time_args["time_steps"] = time_steps
        time_args.pop("tend", None)
    elif "dt" in time_args and "time_steps" in time_args:
        time_args["time_steps"] = time_steps
    else:
        raise ValueError("time block must contain two of tend, dt, and time_steps")

    reduced["time"] = time_args
    return reduced


def load_reduced_backend_payload(source_json: Path) -> tuple[dict[str, Any], int | str]:
    source_json = source_json.resolve()
    payload = _load_effective_source_payload(source_json)
    time_steps = _test_time_steps(payload)

    backend_payload = copy.deepcopy(payload)
    backend_payload.pop("tests", None)
    backend_payload.pop("default_params", None)
    backend_payload["root_path"] = str(source_json)
    _resolve_mesh_paths(backend_payload, source_json.parent)
    apply_verify_run_linear_solver(backend_payload, source_json)
    ensure_json_output(backend_payload)
    return apply_test_time_steps(backend_payload, time_steps), time_steps


def ensure_json_output(payload: dict[str, Any]) -> None:
    output = payload.setdefault("output", {})
    if not isinstance(output, dict):
        raise TypeError("output block must be an object")
    output.setdefault("json", "sim.json")


def _config_payload_for_backend_run(cfg: Any) -> dict[str, Any]:
    if isinstance(cfg, dict):
        return copy.deepcopy(cfg)

    as_dict = getattr(cfg, "as_dict", None)
    if callable(as_dict):
        payload = as_dict()
        if not isinstance(payload, dict):
            raise TypeError(
                f"cfg.as_dict() must return dict, got {type(payload).__name__}"
            )
        from polyfempy.runtime._solve_contract import prepare_generated_backend_payload

        return prepare_generated_backend_payload(payload)

    raise TypeError(
        "polyfem_config must be a dict or object "
        f"with as_dict(), got {type(cfg).__name__}"
    )


def write_backend_payload(
    payload: dict[str, Any],
    output_root: Path,
    source_json: Path,
) -> Path:
    input_dir = output_root / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    payload_path = input_dir / f"{source_json.stem}_reduced.json"
    payload_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return payload_path


def prune_visual_outputs(output_dir: Path) -> list[Path]:
    removed: list[Path] = []
    if not output_dir.exists():
        return removed

    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.suffix.lower() in VISUAL_OUTPUT_SUFFIXES:
            path.unlink()
            removed.append(path)
    return removed


def load_output_metrics(output_dir: Path) -> dict[str, float]:
    sim_json = output_dir / "sim.json"
    if not sim_json.exists():
        raise FileNotFoundError(f"backend output is missing {sim_json}")
    return _error_metrics(_load_json(sim_json))


def load_test_metrics(source_json: Path) -> tuple[dict[str, float], float]:
    payload = _load_effective_source_payload(source_json)
    tests = payload.get("tests")
    if not isinstance(tests, dict):
        return {}, 1e-8
    tolerance = float(tests.get("margin", 1e-8))
    return _error_metrics(tests), tolerance


def compare_metric_maps(
    *,
    left: dict[str, float],
    right: dict[str, float],
    tolerance: float,
) -> list[MetricDiff]:
    keys = sorted(set(left) | set(right))
    rows: list[MetricDiff] = []
    for key in keys:
        if key not in left or key not in right:
            rows.append(MetricDiff(key, left.get(key, float("nan")), right.get(key, float("nan")), float("inf"), tolerance))
            continue
        diff = abs(left[key] - right[key])
        rows.append(MetricDiff(key, left[key], right[key], diff, tolerance))
    return rows


def _print_rows(title: str, rows: Sequence[MetricDiff]) -> bool:
    print(title)
    all_passed = True
    for row in rows:
        status = "PASS" if row.passed else "FAIL"
        all_passed = all_passed and row.passed
        print(
            f"  {status} {row.key}: "
            f"left={row.left:.17g} right={row.right:.17g} "
            f"diff={row.diff:.3g} tolerance={row.tolerance:.3g}"
        )
    return all_passed


def _import_example(example_path: Path):
    if str(example_path.parent) not in sys.path:
        sys.path.insert(0, str(example_path.parent))
    spec = importlib.util.spec_from_file_location(example_path.stem, example_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load example module from {example_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _generated_workspace(output_root: Path, example_path: Path) -> Path:
    return output_root / "generated" / "run"


def _source_json_workspace(output_root: Path, source_json: Path) -> Path:
    return output_root / "source-json" / "run"


def run_generated_example(
    example_path: Path,
    source_json: Path,
    output_root: Path,
    *,
    test_time_steps: int | str,
    keep_visual_output: bool,
) -> Path:
    from polyfempy.runtime import solve

    module = _import_example(example_path)
    workspace = _generated_workspace(output_root, example_path)
    workspace.mkdir(parents=True, exist_ok=True)
    payload = _config_payload_for_backend_run(module.polyfem_config)
    output = payload.setdefault("output", {})
    output["directory"] = str(workspace)
    output.setdefault("json", "sim.json")
    if isinstance(output.get("log"), dict):
        output["log"]["path"] = str(workspace / "polyfem.log")
    payload = apply_test_time_steps(payload, test_time_steps)
    apply_verify_run_linear_solver(payload, source_json)
    solve(cfg=payload)
    if not keep_visual_output:
        removed = prune_visual_outputs(workspace)
        if removed:
            print(f"pruned generated visual outputs: {len(removed)} files")
    return workspace


def run_source_json(
    reduced_source_json: Path,
    source_json: Path,
    output_root: Path,
    *,
    log_level: int,
    max_threads: int,
    keep_visual_output: bool,
) -> Path:
    workspace = _source_json_workspace(output_root, source_json)
    workspace.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        "-I",
        "-m",
        "polyfempy",
        "solve",
        str(reduced_source_json),
        "--output-dir",
        str(workspace),
        "--log-level",
        str(log_level),
        "--max-threads",
        str(max_threads),
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(
            "source JSON backend run failed with exit code "
            f"{result.returncode}: {' '.join(command)}"
        )
    if not keep_visual_output:
        removed = prune_visual_outputs(workspace)
        if removed:
            print(f"pruned source JSON visual outputs: {len(removed)} files")
    return workspace


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run a generated example and its source JSON through the backend, "
            "then compare their err_* metrics."
        )
    )
    parser.add_argument(
        "--example",
        type=Path,
        default=DEFAULT_EXAMPLE,
        help="generated example Python file to run",
    )
    parser.add_argument(
        "--source-json",
        type=Path,
        default=DEFAULT_SOURCE_JSON,
        help="polyfem-data source JSON mirrored by the generated example",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="directory for generated and source JSON backend outputs",
    )
    parser.add_argument(
        "--generated-source-tolerance",
        type=float,
        default=1e-8,
        help="required tolerance for generated output vs source JSON output",
    )
    parser.add_argument(
        "--require-tests-match",
        action="store_true",
        help=(
            "also fail if generated output does not match the source JSON "
            "tests block"
        ),
    )
    parser.add_argument(
        "--keep-visual-output",
        action="store_true",
        help="keep .vtu/.vtm/.pvd files instead of pruning them after comparison",
    )
    parser.add_argument("--log-level", type=int, default=2)
    parser.add_argument("--max-threads", type=int, default=1)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    example_path = args.example.resolve()
    source_json = args.source_json.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"generated example: {example_path}")
    print(f"source JSON:       {source_json}")
    print(f"output root:       {output_root}")

    source_payload, test_time_steps = load_reduced_backend_payload(source_json)
    reduced_source_json = write_backend_payload(source_payload, output_root, source_json)
    print(f"tests.time_steps:  {test_time_steps}")
    print(f"reduced JSON:      {reduced_source_json}")

    generated_output = run_generated_example(
        example_path,
        source_json,
        output_root,
        test_time_steps=test_time_steps,
        keep_visual_output=args.keep_visual_output,
    )
    source_output = run_source_json(
        reduced_source_json,
        source_json,
        output_root,
        log_level=args.log_level,
        max_threads=args.max_threads,
        keep_visual_output=args.keep_visual_output,
    )

    generated_metrics = load_output_metrics(generated_output)
    source_metrics = load_output_metrics(source_output)
    test_metrics, test_tolerance = load_test_metrics(source_json)

    print(f"generated output:  {generated_output}")
    print(f"source output:     {source_output}")

    generated_source_rows = compare_metric_maps(
        left=generated_metrics,
        right=source_metrics,
        tolerance=args.generated_source_tolerance,
    )
    generated_matches_source = _print_rows(
        "generated output vs source JSON output:",
        generated_source_rows,
    )

    generated_matches_tests = True
    if test_metrics:
        tests_rows = compare_metric_maps(
            left=generated_metrics,
            right=test_metrics,
            tolerance=test_tolerance,
        )
        generated_matches_tests = _print_rows(
            "generated output vs source JSON tests:",
            tests_rows,
        )
    else:
        print("source JSON has no tests block; skipping tests comparison")

    if not generated_matches_source:
        return 1
    if args.require_tests_match and not generated_matches_tests:
        return 1
    if not generated_matches_tests:
        print(
            "WARNING: generated output matches the source JSON run, but does "
            "not match the source JSON tests block."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
