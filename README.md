# PolyFEM for Python

`polyfempy` provides model builders and Python bindings to the PolyFEM solver.

## Install

From a checkout with its submodules initialized:

```sh
python -m pip install .
```

Pip handles dependencies, API generation, and compilation. For development,
use `python -m pip install -e ".[dev]"`. See [BUILD.md](BUILD.md) for source-build
requirements and [INSTALL.md](INSTALL.md) for wheel installation.

## Run a simulation

```python
import numpy as np
import polyfempy as pf

vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [.5, .5]], dtype=float)
cells = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
model = pf.model()
body = model.mesh(mesh="array")
body.material(pf.neo_hookean(E=1000, nu=.3))
body.surface_axis(axis="-x", position=0, relative=True).dirichlet(value=[0, 0])
model.solver(linear=pf.linear(solver="Eigen::SimplicialLDLT"))
result = pf.solve(vertices, cells, cfg=model.config(rhs=[0, -1]))
print(result.sol)
```

For a mesh file, pass its path to `model.mesh()` and call
`pf.solve(cfg=model.config(...))` without vertex and connectivity arrays.
`result.sol` contains the raw finite-element solution; rows are not necessarily
mesh vertices.

See [API.md](API.md) for materials, selections, loads, and time settings, and
[TESTING.md](TESTING.md) for validation commands.

## Source layout

- `python/polyfempy/`: generated model builders and simulation runtime.
- `src/`: nanobind bindings to the native solver.
- `cmake/recipes/`: CPM, nanobind, and PolyFEM integration.
- `python-from-jse/`: pinned generator used during builds.
- `polyfem-data/`: reference data; CPM fetches the pinned backend and nanobind sources.
- `tools/`: API generation and release checks.
