# PolyFEM for Python

Build a finite element model with `import polyfem` and run it with `from polyfempy import solve`. The Python API configures the simulation and the compiled PolyFEM backend computes the solution.

## Installation

This forward API is not yet published on PyPI or conda-forge.
[Build this branch from source](BUILD.md), or install a tested wheel supplied
by a maintainer. See the [installation guide](INSTALL.md) for commands and
platform status.

## Run a forward simulation

Simulate a 1 m x 0.2 m elastic beam. Fix its left end, apply a downward body
force, and compute its displacement. This example uses meters and pascals.

Download [beam.obj](https://raw.githubusercontent.com/polyfem/python_data/6b4751b9c3f40fb9653601efe32ce6b6d0ead96f/quickstart/beam.obj)
(12 KB) and save the following code as `beam.py` in the same folder.

```python
import polyfem
from polyfempy import solve

model = polyfem.model()
beam = model.mesh(mesh="beam.obj")
beam.material(polyfem.linear_elasticity(E=1e7, nu=0.3))

# Fix the left end at x = 0.
left = beam.surface_axis(axis="-x", position=0.0)
left.dirichlet(value=[0, 0])

# Apply a downward body force (N/m^3).
model.rhs([0, -1000])
model.solver(linear=polyfem.linear(solver="Eigen::SimplicialLDLT"))

result = solve(cfg=model.config())
print(f"Largest displacement component: {abs(result.sol).max():.6f} m")
```

Run it from that folder:

```sh
python beam.py
```

The last line should be approximately:

```text
Largest displacement component: 0.003648 m
```

Try changing `E` to adjust the material stiffness, or change the load in
`model.rhs(...)`. The [complete script](https://github.com/polyfem/python_data/blob/6b4751b9c3f40fb9653601efe32ce6b6d0ead96f/quickstart/beam.py)
also works when run from another directory, using the mesh beside it.

`result.sol` is the raw finite-element solution; do not assume that each row
corresponds to a mesh vertex. For a small check that creates its own mesh,
run the [installation self-check](INSTALL.md#check-the-installation).

For larger models, see the [classic example guide](https://github.com/polyfem/python_data/blob/ffd670690bf20021e1c380d68b0afe845fb275fd/README.md) and the
[pinned example sources](https://github.com/polyfem/python_data/tree/ffd670690bf20021e1c380d68b0afe845fb275fd/classic_example).
These reference examples require their mesh data and may select linear solvers
that are not enabled in every backend build.

Download either of these scripts and run it with Python in the same environment:
[static elasticity](https://github.com/polyfem/python_data/blob/ffd670690bf20021e1c380d68b0afe845fb275fd/quickstart/static_elasticity.py) and
[transient elasticity](https://github.com/polyfem/python_data/blob/ffd670690bf20021e1c380d68b0afe845fb275fd/quickstart/transient_elasticity.py).
Both create their own mesh and use `Eigen::SimplicialLDLT`. The transient example
also shows time settings and an output directory.

## Source layout

- `python/polyfem/` exposes the model builder and generated configuration helpers.
- `python/polyfempy/runtime/` implements `solve` and returns the solution in `result.sol`.
- `src/` binds the Python runtime to the C++ solver in the `polyfem/` submodule.
- `python-from-jse/` generates the configuration API from PolyFEM schemas and the rules in `generator-config/`.
- `polyfem-data/` supplies the official meshes, input configurations, and reference values.

## Forward solver validation

Push checks run Python and API tests on Windows, Linux, and macOS. Native builds
and simulation checks run in separate workflows. On a ready-for-review pull
request targeting `main`, **Generated Contact Backend** builds and installs a
wheel, runs the short examples and cube reference check, then runs the 75 active
contact cases. See [TESTING.md](TESTING.md) for triggers, commands, and artifacts.

See the [API guide](API.md) for geometry, materials, boundary conditions and loads.
