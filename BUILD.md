# Build the Python bindings

For an installed binary package, start with [INSTALL.md](INSTALL.md).
This guide covers building and changing the library.

Use a Git checkout of the intended release branch. GitHub's **Download ZIP** does
not include submodule contents. Keep the recorded submodule revisions together;
updating each dependency to its latest branch is not the same configuration.

Clone the forward branch (use a short parent directory on Windows):

```sh
git clone --branch jingyaoyu_forward --recurse-submodules https://github.com/peteryu0131/polyfem-python.git
cd polyfem-python
```

For an existing checkout, select the intended revision first, then run
`git submodule update --init --recursive` from its root. This checks out the
commits recorded by the parent repository. Do not add `--remote` when reproducing
a release: it follows remote branches instead of the recorded commits.
Maintainers should review dependency updates and commit the new submodule
pointers explicitly.

Use an environment with Python 3.10 or later, CMake 3.25 or later, and a C++ compiler. On Windows, use an x64 Visual Studio development environment. To use an existing conda environment named `polyfem`, activate it with `conda activate polyfem` first.

On Windows, use a short checkout path such as `C:\src\polyfem-python`. Some
third-party MSBuild projects still fail when generated paths exceed 260 characters.

## Development install

```sh
python -m pip install "setuptools>=77" wheel build "cmake>=3.25" "nanobind>=2" numpy pytest
python tools/generate_polyfem_api.py
python -m pip install -e . --no-build-isolation
python -m polyfempy backend-info --require
```

Both Python packages live under `python/`. After moving from an older checkout
layout, repeat the development install above to refresh editable-install paths.

Regenerate the API after updating the schemas, generator or generator configuration. Generation writes to `python/polyfempy/generated_api/` and may download the linked solver schemas. PolySolve schemas are cached under `build/polysolve-json-specs/<revision>/`, using the revision in the pinned PolyFEM source. Switching that revision selects a separate cache. The first native build also needs access to the pinned third-party dependencies.

Generation is currently an explicit prerequisite, not an automatic installation
step. Run it before building or installing from source. Installed wheel users
receive the generated API and its small builder runtime inside the wheel.
They do not need a separate generator checkout. A build fails early if these
files are missing.

The distribution is named `polyfempy`. Check its public model-building import
with the same interpreter used for installation:

```sh
python -c "import polyfem; print(type(polyfem.model()).__name__)"
```

## Build and check a wheel

From a fresh checkout with the development dependencies and generated API in place, run:

```sh
python -m build --no-isolation
```

This creates a source archive, then builds the wheel from the extracted archive.
Set `N_THREADS` to a positive integer to limit parallel compilation. The default
is at most two workers. Visual Studio builds one project at a time and assigns
the worker limit to its compiler, so parallel settings do not multiply. Release builds enable
PolyFEM's portable CPU settings.

macOS builds use Apple Accelerate and Eigen solvers on both Intel and Apple
Silicon. MKL is disabled on macOS. The examples select
`Eigen::SimplicialLDLT`, which is available on all supported platforms.

In a separate environment with a compatible Python version and architecture,
install the resulting wheel by its full path. Run the installed-package check
from the checkout:

```sh
python tools/check_installed_package.py --readme README.md --reference polyfem-data/units/neohookean.json
```

The check starts isolated Python processes in a temporary directory. It verifies
the import location, native backend, README beam, installation self-check, Python
load function, and six cube error metrics. Run from the checkout root so the
check can find `INSTALL.md`, `examples/quickstart/beam.obj`, and the cube reference
in `polyfem-data`.
The cube data is not included in the source archive; omit `--reference` there
to check the beam, installation self-check, and Python load function.
For the full example sweep, see [TESTING.md](TESTING.md).

Check both `import polyfem` and `python -I -m polyfempy backend-info --require` there.
An editable install or a wheel assembled with an existing local binary does not
verify a fresh native build of the release revision.

## Build adapters

The parent repository owns the Python extension build adapters in
`cmake/recipes/polyfem.cmake` and `cmake/recipes/libigl.cmake`. They do not edit
the PolyFEM submodule.

The Python adapter disables the backend's interpreter-embedding setup, then
links its target with `pybind11::module` and enables Python expressions.
`POLYFEM_WITH_PYTHON=OFF` still disables Python expressions. Nanobind remains
the binding library for the public extension.

The libigl adapter uses the pinned dependency version and the patches shipped
with PolyFEM. It applies them with `git apply`, without shell redirection or a
separate `patch` program. Review this adapter when upgrading the PolyFEM pin.

The adapter tests configure small CMake fixtures without downloading dependencies
or compiling C++. They do not replace native build and installed-package checks.

## Source archives

The source archive includes the binding sources, pinned PolyFEM sources,
generated API, and builder runtime. It does not need Git submodules after
extraction. A compiler and network access for pinned native dependencies are
still required. Mesh data and the full example collection are not build inputs.

Use a Git checkout when changing schemas or regenerating the API. The archive
contains the generated package, not the complete generator development tools
or the repository test suite. Run development tests from a Git checkout.
The archive retains `tools/check_installed_package.py` and the quick-start
scripts for checking an installed build.

## Windows compatibility settings

Importing `polyfempy` leaves Python's output streams and OpenMP environment
settings unchanged by default. If an existing application needs the previous
console and OpenMP workarounds, set `POLYFEMPY_WINDOWS_AUTOCONFIG=1` before
importing. `POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG=1` always disables automatic changes.
Use this option only for an application that needs those settings.

## What the build tools do

| Component | Role |
| --- | --- |
| `python-from-jse` | Generates Python configuration classes and model-building helpers from schemas. |
| `python/polyfempy/runtime/` | Prepares solve inputs and wraps the native result as `Result`. |
| nanobind | Exposes the C++ solver and its methods to Python and converts supported argument and result types. |
| CMake and a C++ compiler | Build PolyFEM and the Python extension. |
| pip and conda | Install packages into the selected Python environment. |

The call path is `model.config()` -> Python `solve()` -> the bound C++ `Solver`
-> PolyFEM -> `Result`. The bindings are declared in `src/binding.cpp` and
`src/state/state.cpp`. They are compiled, not generated by `python-from-jse`.

The build uses nanobind's shared runtime. On Windows, the reviewed wheel contains
`polyfempy.pyd` and `nanobind.dll`. Release packaging must include and check the
corresponding native libraries on each platform. End users should not have to
install the binding toolchain to use a binary package.
