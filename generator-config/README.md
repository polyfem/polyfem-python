# PolyFEM generator configuration

The API generator reads two files in this directory:

- `api_aliases.json` defines public factory names and aliases.
- `id_relationships.json` connects geometry IDs to materials, boundary conditions and contact settings.

The generic generator lives in the `python-from-jse/` submodule.

## Generated global configuration helpers

The generic generator accepts keyword arguments directly on polymorphic
factories, for example `polyfem.time(dt=0.01, tend=0.5, integrator="BDF")`.
It selects a schema branch using its required fields; explicit variant objects
remain supported for compatibility.

Root configuration sections also generate methods on the model builder:

```python
import polyfem

model = polyfem.model()
model.time(dt=0.01, tend=0.5, integrator="BDF")
model.units(length="cm", mass="g", time="s")
model.output(directory="runs/example")
config = model.config()
```

These setters return the model and replace the same section when called again.
Passing that section again to `model.config(...)` raises an error instead of
silently overriding it. Shared builder methods also support setting `space`,
`materials`, and `initial_conditions`, while preserving their no-argument reads.
`boundary_conditions(...)` saves keyword settings but keeps its dictionary return
value for existing callers. Passing that unchanged dictionary to `model.config(...)`
is supported without duplicating boundary entries.

Explicit section settings are combined with body settings when exporting the
model. Lists are appended; conflicting non-list fields raise an error. Repeating
a setter replaces only its explicit settings. ID assignment remains in the
shared builder and relationship rules, so users can work with body and selection
objects. See [Model settings](../API.md#model-settings) for examples and merge rules.

Maintain this behavior in `python-from-jse/generator/JsonToTreeClass.py` and
`python-from-jse/generator/model_builder.py`, not the generated output. From the
repository root, regenerate with:

```powershell
python tools/generate_polyfem_api.py
```

## Geometry and ID binding

Geometry factories also generate model methods such as `model.mesh(...)`,
`model.mesh_array(...)` and `model.mesh_sequence(...)`. Use
`is_obstacle=True` to mark an obstacle.

Shared builder logic registers geometry and binds IDs. Body-level consumers
such as `body.material(...)` use volume relationships. Direct boundary
conditions such as `body.dirichlet(...)` use the configured surface
relationships and select the whole boundary. Partial selections retain their
own handles.

See [API.md](../API.md) for user-facing examples. New geometry types use the
same generated dispatch. New relationships must be described in
`id_relationships.json` and tested against the backend configuration.
