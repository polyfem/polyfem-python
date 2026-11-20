# PolyFEM for Python

`polyfempy` provides model builders and Python bindings to the PolyFEM solver.
`import polyfem` is an alias with the same public API.

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
import polyfempy as pf

model = pf.model()
body = model.mesh(mesh="part.msh")
body.material(pf.neo_hookean(E=1000, nu=.3))
body.surface_axis(axis="-x", position=0, relative=True).dirichlet(value=[0, 0])
model.solver(linear=pf.linear(solver="Eigen::SimplicialLDLT"))
result = pf.solve(model.config(rhs=[0, -1]))
print(result.sol)
```

Supply a 2D mesh file at `part.msh`, or pass another path to `model.mesh()`.
`result.sol` contains the raw finite-element solution; rows are not necessarily
mesh vertices.

See [API.md](API.md) for materials, selections, loads, and time settings, and
[TESTING.md](TESTING.md) for validation commands.

## Source layout

- `python/polyfempy/solve.py`: configuration preparation, solving, and results.
- `python/polyfempy/generated_api/`: generated model builders.
- `python/polyfem/`: automatic alias for `polyfempy`.
- `src/`: nanobind bindings to the native solver.
- `cmake/recipes/`: CPM, nanobind, and PolyFEM integration.
- `python-from-jse/`: pinned generator used during builds.
- `polyfem-data/`: reference data; CPM fetches a pinned PolyFEM commit and nanobind release.
- `tools/`: API generation and release checks.
