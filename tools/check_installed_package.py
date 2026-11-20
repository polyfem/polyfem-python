"""Smoke-test an installed wheel outside the source checkout."""

import subprocess
import sys
from tempfile import TemporaryDirectory


CHECK = r"""
from pathlib import Path
import sys
import numpy as np
import polyfempy as pf
from polyfempy.polyfempy import Solver

assert Path(pf.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [.5, .5]], dtype=float)
cells = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
mesh = Path("square.obj")
mesh.write_text("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv .5 .5 0\nf 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n")
model = pf.model()
body = model.mesh(mesh=str(mesh.resolve()))
body.material(pf.neo_hookean(E=1000, nu=.3))
body.surface_axis(axis="-x", position=0, relative=True).dirichlet(value=[0, 0])
model.solver(linear=pf.linear(solver="Eigen::SimplicialLDLT"))
config = model.config(rhs=[0, -1])
file_result = pf.solve(cfg=config)
mesh.unlink()
array_result = pf.solve(vertices, cells, cfg=config)
assert np.isfinite(array_result.sol).all()
assert abs(array_result.sol).max() > 1e-6
np.testing.assert_allclose(np.sort(array_result.sol.ravel()), np.sort(file_result.sol.ravel()), rtol=1e-6, atol=1e-9)
load = Path("body_load.py").resolve()
load.write_text("def value(x, y, z, t, index):\n    return -1.0\n")
function_result = pf.solve(vertices, cells, cfg=model.config(rhs=[0, {"file_name": str(load), "function_name": "value"}]))
np.testing.assert_allclose(function_result.sol, array_result.sol, rtol=1e-6, atol=1e-9)
print("Installed package: file, NumPy array, and Python load-function solves passed.")
"""


def main():
    with TemporaryDirectory(prefix="polyfem-installed-") as folder:
        subprocess.run([sys.executable, "-I", "-c", CHECK], cwd=folder, check=True, timeout=300)


if __name__ == "__main__":
    main()
