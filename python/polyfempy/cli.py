"""Small command-line entry points for local PolyFEM-Python builds."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence


def _load_polyfempy():
    import polyfempy as pf

    return pf


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="polyfempy",
        description="PolyFEM-Python command-line helpers.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    backend_info = subparsers.add_parser(
        "backend-info",
        help="print whether the compiled C++ backend is available",
    )
    backend_info.add_argument(
        "--require",
        action="store_true",
        help="return a non-zero exit code when the C++ backend is unavailable",
    )

    solve = subparsers.add_parser(
        "solve",
        help="run PolyFEM with a JSON or YAML config file through the C++ backend",
    )
    solve.add_argument("config", help="path to a .json, .yaml, or .yml config file")
    solve.add_argument(
        "--output-dir",
        default="",
        help="optional output directory passed to the backend",
    )
    solve.add_argument(
        "--log-level",
        type=int,
        default=2,
        help="PolyFEM log level, from 0 (all logs) to 6 (off)",
    )
    solve.add_argument(
        "--max-threads",
        type=int,
        default=1,
        help="maximum number of backend threads",
    )
    solve.add_argument(
        "--strict-validation",
        dest="strict_validation",
        action="store_true",
        default=True,
        help="enable backend settings validation",
    )
    solve.add_argument(
        "--no-strict-validation",
        dest="strict_validation",
        action="store_false",
        help="disable backend settings validation",
    )

    return parser


def _backend_status(pf) -> tuple[bool, object]:
    available = bool(getattr(pf, "cpp_backend_available", lambda: False)())
    error = getattr(pf, "cpp_backend_error", lambda: None)()
    return available, error


def _print_backend_info(pf, *, require: bool) -> int:
    available, error = _backend_status(pf)
    print(f"backend_available={available}")
    print(f"backend_error={error}")
    if require and not available:
        print(
            "C++ backend not loaded. Install a wheel for your Python version "
            "and platform, or follow BUILD.md to build from source.",
            file=sys.stderr,
        )
        return 1
    return 0


def _solve_config(pf, args: argparse.Namespace) -> int:
    available, error = _backend_status(pf)
    if not available:
        print(
            "C++ backend not loaded. Real solves require the compiled backend. "
            f"Error: {error}",
            file=sys.stderr,
        )
        return 1

    if not hasattr(pf, "polyfem_command"):
        print(
            "C++ backend is loaded but does not expose polyfem_command.",
            file=sys.stderr,
        )
        return 1

    config = Path(args.config)
    suffix = config.suffix.lower()
    yaml_path = str(config) if suffix in {".yaml", ".yml"} else ""
    json_path = "" if yaml_path else str(config)

    try:
        pf.polyfem_command(
            json=json_path,
            yaml=yaml_path,
            log_level=args.log_level,
            strict_validation=args.strict_validation,
            max_threads=args.max_threads,
            output_dir=args.output_dir,
        )
    except Exception as exc:
        print(f"polyfempy solve failed: {exc}", file=sys.stderr)
        return 1

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    pf = _load_polyfempy()

    if args.command == "backend-info":
        return _print_backend_info(pf, require=args.require)

    if args.command == "solve":
        return _solve_config(pf, args)

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
