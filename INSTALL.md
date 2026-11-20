# Installation

For this branch, install from a checkout with its submodules initialized:

```sh
python -m pip install .
```

Pip obtains Python and build dependencies from `pyproject.toml` and runs the
build automatically. See [BUILD.md](BUILD.md) for cloning and compiler details.

To install a wheel supplied by a maintainer:

```sh
python -m pip install ./polyfempy-<version>-<python>-<abi>-<platform>.whl
```

Choose a wheel matching your Python version, operating system, and CPU
architecture. Wheels include the generated API and native extension. No
compiler or generator is needed to install one.

A conda environment can be used by activating it before running pip. This
repository does not maintain a separate native conda recipe.

## Check the installation

```sh
python -c "import polyfempy; from polyfempy.polyfempy import Solver; print(polyfempy.model())"
python tools/check_installed_package.py
```

The second command runs small file and NumPy-array simulations using the
installed package in an isolated process outside the checkout.
