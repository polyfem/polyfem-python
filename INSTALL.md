# Install PolyFEM for Python

The distribution is named `polyfempy`. Use `import polyfem` to build a model and
`from polyfempy import solve` to run it.

## Release status

The forward release is being prepared. Its reviewed Windows x64 wheel works
with CPython 3.11, including in a new conda environment. Linux and macOS native
packages still need successful release builds and installation checks.

The existing public packages contain older APIs. As checked on 2026-11-15,
[conda-forge](https://anaconda.org/conda-forge/polyfempy) provides version 0.5.2
for Windows x64, Linux x64, and Intel macOS. [PyPI](https://pypi.org/project/polyfempy/)
provides a 0.2.3 source archive. Installing either does not install this forward
branch or its new model builder.

The new release keeps conda, pip, and source installation as release targets.
The initial installation checks use CPython 3.11 on Windows x64, Linux x64,
Intel macOS, and Apple Silicon macOS. A configured job is not proof that a
platform works.

This branch requires Python 3.10 or later. Python-only tests cover CPython
3.10-3.12. Native package checks for Python 3.10 and Python 3.12 are deferred
until release preparation. The legacy project's Python 3.7/3.8
support does not apply to the current API.

## Install a reviewed wheel with pip

This section assumes a maintainer has supplied a tested `.whl` file.
There is no published forward release to download from the package indexes.
If you do not have a wheel, follow [BUILD.md](BUILD.md).
Choose a wheel for your operating system, CPU architecture, and Python version.
A Windows wheel cannot be used on Linux or macOS. Intel and Apple Silicon Macs
also need matching binaries.

Create a new environment. For Windows with CPython 3.11 installed:

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
```

For Linux or macOS with CPython 3.11 installed:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
```

From the directory containing the wheel, install it with the active environment's
Python. This filename is the reviewed Windows x64 / CPython 3.11 build:

```sh
python -m pip install --only-binary=:all: ./polyfempy-0.8-cp311-cp311-win_amd64.whl
python -m pip check
```

Use the corresponding filename for another platform when its checked wheel is
available. `--only-binary=:all:` prevents an unexpected source build. Wheel users
do not need CMake, a compiler, a separate nanobind installation, or Git submodules.

## Use pip inside a conda environment

Create the environment first, then install the matching wheel as shown above:

```sh
conda create -n polyfem-user --override-channels -c conda-forge python=3.11 pip numpy
conda activate polyfem-user
```

Use `python -m pip` after activation so installation uses this environment.
See conda's [environment guide](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html#using-pip-in-an-environment)
for guidance on mixing conda and pip packages.

## Install a native conda package

The original library's supported command is:

```sh
conda install -c conda-forge polyfempy
```

Until the feedstock is updated, this installs the legacy release. The new native
conda package must pass the recipe's installation and solve checks before it is
submitted to the existing [feedstock](https://github.com/conda-forge/polyfempy-feedstock).
Maintainers can build and test it using [conda/README.md](conda/README.md).

A native conda package and a pip wheel installed inside conda are separate
distribution routes. Both are part of the release checks.

## Check the installation

```sh
python -I -c "import polyfem; from polyfempy import solve; print(type(polyfem.model()).__name__)"
python -I -m polyfempy backend-info --require
```

The second command must report `backend_available=True`.

Save the following code as `check_install.py` and run `python check_install.py`.
It creates a small square mesh, fixes its boundary, and applies a downward body
force. Only the center vertex can move. No example data download is needed.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

import polyfem
from polyfempy import solve

with TemporaryDirectory() as folder:
    mesh = Path(folder) / "square.obj"
    mesh.write_text(
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv .5 .5 0\n"
        "f 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n",
        encoding="utf-8",
    )

    model = polyfem.model()
    square = model.mesh(mesh=str(mesh))
    square.material(polyfem.neo_hookean(E=1000, nu=0.3))
    square.surface_all().dirichlet(value=[0, 0])
    model.rhs([0, -1])
    model.solver(linear=polyfem.linear(solver="Eigen::SimplicialLDLT"))

    result = solve(cfg=model.config())
    print(f"Largest displacement component: {abs(result.sol).max():.8f}")
```

The last line should report a value of about `0.00011235`. This checks that the
installed package can run a complete simulation. Then try the
[README cantilever beam](README.md#run-a-forward-simulation) with its supplied mesh.

If pip rejects a wheel, check the Python version and architecture before trying
another file. If an import fails, check the active interpreter with
`python -c "import sys; print(sys.executable)"`. A successful Python import alone
does not establish that the native solver can load.

## Build from source

Use [BUILD.md](BUILD.md) when changing the library or building for a platform
without a compatible binary. Source builds need a C++ compiler and native build
dependencies. The current equivalents of the legacy setup.py commands are
`python -m pip install .` and `python -m pip install -e .`, after the preparation
steps in that guide.
