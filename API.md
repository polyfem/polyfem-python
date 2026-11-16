# Forward API

For a complete example with a supplied mesh, use the [README beam](README.md#run-a-forward-simulation).
This guide shows API patterns. Supply your own mesh paths; `part.msh` below
must be a 2D mesh. Later snippets are alternatives, not one script to run in
sequence. They assume `import polyfem` and a fresh `model` with the named bodies.

```python
import polyfem
from polyfempy import solve

model = polyfem.model()
body = model.mesh(mesh="part.msh")
body.material(polyfem.neo_hookean(E=1000, nu=0.3))
body.dirichlet(value=[0, 0])

model.solver(linear=polyfem.linear(solver="Eigen::SimplicialLDLT"))
config = model.config(rhs=[0, -1])
result = solve(cfg=config)
```

Use vectors with two components for 2D and three for 3D. Choose values consistent with the model's units.

## Run an existing configuration

`solve(cfg=...)` accepts `model.config()`, a configuration dictionary, or a JSON
file path. For an existing JSON file:

```python
from polyfempy import solve

result = solve(cfg="simulation.json")
```

By default, relative mesh paths in that file are resolved from the JSON file's
directory. Relative paths in a model or dictionary use the working directory
unless you set `root_path`.

## Use an array mesh

Pass vertex coordinates and zero-based cell indices directly to `solve` when
the mesh is already in Python. This complete example needs no mesh file:

```python
import numpy as np
from polyfempy import solve

vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]])
cells = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
config = {
    "materials": {"type": "NeoHookean", "E": 1000, "nu": 0.3},
    "boundary_conditions": {
        "dirichlet_boundary": [{"id": "all", "value": [0, 0]}],
        "rhs": [0, -1],
    },
    "solver": {"linear": {"solver": "Eigen::SimplicialLDLT"}},
}

result = solve(vertices, cells, cfg=config)
print(f"Largest displacement component: {abs(result.sol).max():.8g}")
```

The result is about `0.00011234568`. Each vertex row contains two coordinates
for 2D or three for 3D; each cell row contains its vertex indices. Supply both
arrays together. They take precedence over mesh files in `cfg`.

## Geometry

| Method | Purpose |
| --- | --- |
| `model.mesh(mesh=...)` | Add a mesh |
| `model.mesh_array(mesh=..., array=...)` | Repeat a mesh in an array |
| `model.mesh_sequence(mesh_sequence=..., fps=...)` | Add a sequence of meshes |

Set `is_obstacle=True` for an obstacle. Geometry methods use the names and parameters generated from the spec.

```python
obstacles = model.mesh_array(
    mesh="edge.obj",
    array=polyfem.array(offset=1, size=[21, 1]),
    is_obstacle=True,
)
obstacles.dirichlet(value=[0, 0])
```

Use these geometry methods for both deformable meshes and obstacles.

## Materials and body settings

```python
body.material(polyfem.neo_hookean(E=1000, nu=0.3, rho=1000))
body.velocity(value=[0, 0])
body.discr_order(order=1)
```

The returned body connects these settings to its mesh. The builder writes the volume ID when a body-level setting needs it. Creating a mesh alone preserves its existing volume tags.

## All or part of a boundary

Apply a condition directly to the body to select its entire outer boundary.

```python
body.dirichlet(value=[0, 0])
```

Alternatively, select part of the boundary before applying the condition:

```python
# Select the left side of the mesh bounding box.
boundary = body.surface_axis(axis="-x", position=0.0, relative=True)
boundary.dirichlet(value=[0, 0])
```

With `relative=True`, position 0 is the minimum and position 1 is the maximum
of the bounding box along that axis. `"-x"` selects the side at or below the
chosen x position. With `relative=False`, position uses mesh coordinates.
The selection parameters follow the PolyFEM spec. The API assigns the boundary ID and connects the condition to it. The same binding mechanism supports `neumann`, `normal_aligned_neumann` and `pressure`.

A direct whole-boundary condition and a partial or file-based surface selection on the same body are rejected when they would overwrite each other's selection. Apply conditions to explicit selections when a body has multiple boundary regions.

For an existing sideset file:

```python
boundary = body.surface_from_file("surface_sidesets.txt")
boundary.dirichlet(value=[0, 0])
```

If the file contains several groups, supply `tag=...` to identify the intended group.

Use `surface_all()` when you need an explicit whole-boundary object. It also supports shared labels:

```python
first_surface = first.surface_all()
second_surface = second.surface_all(same_as=first_surface)
first_surface.dirichlet(value=[0, 0])
```

Here `first` and `second` are bodies in the same model. Both surfaces share
one boundary ID, so the condition applies to both.

## Body loads

For a uniform load across the model:

```python
config = model.config(rhs=[0, -1])
```

For two bodies named `first` and `second` in the same model:

```python
first.rhs(value=[0, -1])
second.rhs(value=[0, 0])
config = model.config()
```

Choose one of these load scopes per model. Supplying both a global RHS and body-specific RHS raises an error.

## Time, contact and output

Add these settings to a model with geometry, materials, loads, and constraints
already configured. A transient elastic model also needs material density.
For a complete transient example, use the
[transient square](https://github.com/polyfem/python_data/blob/ffd670690bf20021e1c380d68b0afe845fb275fd/quickstart/transient_elasticity.py).

```python
config = model.config(
    time=polyfem.time(dt=0.01, tend=0.5, integrator="BDF"),
    contact=polyfem.contact(enabled=True, dhat=0.001),
    output=polyfem.output(
        directory="results/example",
        json="sim.json",
        paraview=polyfem.output_paraview(file_name="sim.pvd"),
    ),
)
result = solve(cfg=config)
```

A relative output directory is resolved from the working directory. `result.sol`
is the raw finite-element solution array; its rows are not necessarily mesh
vertices. Open the saved `.pvd` file in ParaView to view the result on the mesh.

`result.sol` is a NumPy array. For a transient simulation, it contains the final
state; use the saved output files for earlier time steps. `result.meta` contains
backend metadata, whose keys may depend on the solve.

Configuration sections can also be set with `model.time(...)`, `model.output(...)` and other generated section methods. Use one configuration style for each section. Passing a section to `model.config(...)` after setting it on the model raises an error.

To show only errors in the terminal, set
`model.output(log=polyfem.output_log(level="error"))`. For a dictionary, use
`config["output"] = {"log": {"level": "error"}}`. Log levels apply to both
file meshes and array meshes. Add the log settings to the existing output
section if you also want to save results.

## Input errors

`solve` checks array meshes before calling the native solver. Vertex coordinates
must be finite, and cell indices must be integers starting at zero and smaller
than the number of vertices. JSON configuration errors include the file path
and line number. If the backend rejects an incomplete `LinearElasticity` or
`NeoHookean` material, the exception identifies the required parameter pairs.

## Generation and maintenance

Configuration factories and geometry entry points are generated from the pinned specs. Shared builder code manages scene objects and IDs. `generator-config/id_relationships.json` defines how bodies and selections connect to materials, loads and boundary conditions.

PolySolve specs are cached by the revision selected by the pinned PolyFEM source. After updating the source or generator, regenerate and test using [BUILD.md](BUILD.md) and [TESTING.md](TESTING.md).
