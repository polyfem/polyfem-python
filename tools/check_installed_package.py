"""Check an installed package from a directory outside the source checkout."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


def first_python_example(path: Path) -> str:
    blocks = re.findall(r"```python\n(.*?)```", path.read_text(encoding="utf-8"), re.DOTALL)
    if not blocks:
        raise RuntimeError(f"{path} has no Python example.")
    return blocks[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readme", type=Path, required=True)
    parser.add_argument("--reference", type=Path, help="PolyFEM cube reference JSON")
    parser.add_argument("--example", type=Path, action="append", default=[])
    args = parser.parse_args()
    readme = args.readme.resolve()
    examples = [
        ("README", first_python_example(readme), 0.003648, 5e-7),
        (
            "Installation",
            first_python_example(readme.with_name("INSTALL.md")),
            0.00011234567901234568,
            1e-10,
        ),
        (
            "Python expression",
            '''
from pathlib import Path
from polyfempy import solve

load = Path("body_load.py").resolve()
load.write_text("def value(x, y, z, t, index):\\n    return -1.0\\n", encoding="utf-8")
vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]])
cells = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
config = {
    "materials": {"type": "NeoHookean", "E": 1000, "nu": 0.3},
    "boundary_conditions": {
        "dirichlet_boundary": [{"id": "all", "value": [0, 0]}],
        "rhs": [0, {"file_name": str(load), "function_name": "value"}],
    },
    "solver": {"linear": {"solver": "Eigen::SimplicialLDLT"}},
}
result = solve(vertices, cells, cfg=config)
''',
            0.00011234567901234568,
            1e-10,
        ),
    ]

    with TemporaryDirectory(prefix="polyfem-installed-") as folder:
        shutil.copyfile(
            readme.parent / "examples/quickstart/beam.obj", Path(folder) / "beam.obj"
        )
        # -I excludes the checkout, PYTHONPATH, and user site-packages.
        check = """
import os
import pathlib
import sys
streams = (sys.stdout, sys.stderr)
openmp_setting = os.environ.get("KMP_DUPLICATE_LIB_OK")
import polyfem
import polyfempy
if (sys.stdout, sys.stderr) != streams or os.environ.get("KMP_DUPLICATE_LIB_OK") != openmp_setting:
    raise RuntimeError("Import changed the host streams or OpenMP settings.")
for module in (polyfem, polyfempy):
    path = pathlib.Path(module.__file__).resolve()
    if not path.is_relative_to(pathlib.Path(sys.prefix).resolve()):
        raise RuntimeError(f"Expected an installed package, imported {path}")
    print(f"{module.__name__}: {path}")
if not polyfempy.cpp_backend_available():
    raise RuntimeError(f"Native backend is unavailable: {polyfempy.cpp_backend_error()}")
from importlib.metadata import distribution
dist = distribution("polyfempy")
if not any(entry.name == "polyfempy" and entry.value == "polyfempy.cli:main" for entry in dist.entry_points):
    raise RuntimeError("The package is missing the polyfempy console entry point.")
print(f"Installed version: {dist.version}")
import numpy as np
"""
        for label, code, expected, atol in examples:
            validation = (
                f"\nnp.testing.assert_allclose(abs(result.sol).max(), {expected!r}, "
                f"rtol=1e-5, atol={atol!r})\n"
                f"print({label!r} + ' example passed.')\n"
            )
            subprocess.run(
                [sys.executable, "-I", "-c", check + code + validation],
                cwd=folder,
                check=True,
                timeout=300,
            )
        subprocess.run(
            [sys.executable, "-I", "-m", "polyfempy", "backend-info", "--require"],
            cwd=folder,
            check=True,
            timeout=300,
        )
        if args.reference is not None:
            reference = args.reference.resolve()
            expected = json.loads(reference.read_text(encoding="utf-8"))["tests"]
            output = Path(folder) / "cube"
            subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-m",
                    "polyfempy",
                    "solve",
                    str(reference),
                    "--output-dir",
                    str(output),
                    "--log-level",
                    "2",
                    "--max-threads",
                    "1",
                ],
                cwd=folder,
                check=True,
                timeout=300,
            )
            actual = json.loads((output / "sim.json").read_text(encoding="utf-8"))
            metrics = {
                key: value for key, value in expected.items() if key.startswith("err_")
            }
            if not metrics:
                raise RuntimeError("The reference JSON has no error metrics.")
            for key, value in metrics.items():
                if not abs(actual[key] - value) <= expected.get("margin", 1e-8):
                    raise AssertionError(f"{key}: got {actual[key]}, expected {value}")
            print(f"Cube reference passed ({len(metrics)} metrics).")
        for example in args.example:
            subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-c",
                    "import runpy, sys, numpy as np; "
                    "result = runpy.run_path(sys.argv[1])['main'](); "
                    "assert np.isfinite(result.sol).all(); "
                    "assert abs(result.sol).max() > 0",
                    str(example.resolve()),
                ],
                cwd=folder,
                check=True,
                timeout=300,
            )
            print(f"Example passed: {example.name}")


if __name__ == "__main__":
    main()
