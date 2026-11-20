# Build the Python bindings

Clone this repository with its pinned submodules:

```sh
git clone --recurse-submodules https://github.com/polyfem/polyfem-python.git
cd polyfem-python
```

Install with the Python interpreter you want to use:

```sh
python -m pip install .
```

For development:

```sh
python -m pip install -e ".[dev]"
```

Pip installs the build requirements declared in `pyproject.toml`. The build
backend runs CMake, generates the Python API, builds the native extension, and
installs it alongside the Python package. No separate generation or dependency
installation command is needed. A source build requires a C++17 compiler and
network access for dependencies that have not been cached.

An existing conda environment works too; activate it before running pip.
Python 3.10 or later is required. On Windows, use an x64 Visual Studio developer
shell and a short checkout path.

## Wheels and source archives

```sh
python -m build
python tools/check_distribution.py "dist/*.tar.gz" "dist/*.whl"
```

The `build` frontend is included in the development extra. Source archives
contain the generator and generated API; CPM fetches pinned backend sources. CMake installs
only the extension and the builder runtime into wheels; build-time schema tools
stay in the build directory. Nanobind supplies NumPy/Eigen conversion and the
platform-specific extension suffix.

The three recipes load CPM and fetch pinned nanobind and PolyFEM sources.
Build defaults are declared in `pyproject.toml`. A small patch to the pinned
backend fixes the libigl patch command and links Python as an extension.
Reference data is used by tests and needs no CMake setup. The backend submodule
is no longer required.

The native build uses `build/<wheel-tag>/`. Override CMake options through pip,
for example `python -m pip install . -Ccmake.define.POLYFEM_WITH_TRIANGLE=ON`.
The default build uses two compiler workers. Override that with
`-Cbuild.tool-args=-j1` when memory is limited.
