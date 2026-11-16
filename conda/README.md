# Build the conda package

This recipe builds a native conda package from the same source archive used for
the wheels. It does not wrap a prebuilt wheel. The existing conda-forge release
supports `win-64`, `linux-64`, and `osx-64`; those remain required targets.
The initial installation checks use Python 3.11 on these three platforms
and `osx-arm64` for M-series Macs. Other Python versions need validation
before publishing.

The recipe is prepared for local and CI validation. It has not been published
to conda-forge. Keep its version aligned with `project.version` in
`pyproject.toml` when preparing a release.

## Local build

Activate a separate conda environment containing `conda-build`. Prepare a source
archive using [BUILD.md](../BUILD.md), then set its URL and SHA-256 checksum.
For the current version, in PowerShell:

```powershell
$archive = Get-Item dist/polyfempy-0.8.tar.gz
$env:POLYFEMPY_SDIST_URL = ([Uri]$archive.FullName).AbsoluteUri
$env:POLYFEMPY_SDIST_SHA256 = (Get-FileHash $archive -Algorithm SHA256).Hash.ToLower()
$env:CONDA_BLD_PATH = "C:/pf-cb"
$env:CPM_SOURCE_CACHE = "C:/pf-cpm"
conda-build conda/recipe --python 3.11 --no-anaconda-upload --override-channels -c conda-forge --output-folder conda-out
```

On Linux or macOS:

```sh
export CPM_SOURCE_CACHE="$PWD/.cache/CPM"
export POLYFEMPY_SDIST_URL="$(python -c 'from pathlib import Path; print(Path("dist/polyfempy-0.8.tar.gz").resolve().as_uri())')"
export POLYFEMPY_SDIST_SHA256="$(python -c 'import hashlib; from pathlib import Path; print(hashlib.sha256(Path("dist/polyfempy-0.8.tar.gz").read_bytes()).hexdigest())')"
conda-build conda/recipe --python 3.11 --no-anaconda-upload --override-channels -c conda-forge --output-folder conda-out
```

The source archive includes the generated API, builder runtime, and native
sources. The native build still downloads pinned third-party sources. The
recipe installs build tools separately from runtime dependencies and tests the
installed package in conda-build's test environment. It runs the README beam,
installation self-check, and all three quick-start scripts using the mesh
included in the source archive; an import-only check is not sufficient.

Windows also needs the Visual Studio 2022 C++ tools. Use short build and cache
paths. The recipe selects Ninja so the same build commands work across systems.

## Release workflow

Run the package matrix explicitly as described in
[TESTING.md](../TESTING.md#release-package-checks). Ordinary pushes do not build
the complete release matrix.

For manual dispatch, first make the workflow available on the repository's
default branch. Then use
**Actions -> Release Packages -> Run workflow** and select the release branch.
See [TESTING.md](../TESTING.md#release-package-checks) for the prerequisite.
Each native conda job builds for one platform and Python version, runs the
installed-package checks, and uploads the resulting package as an artifact.
The workflow does not upload to Anaconda.org or PyPI.

After all required jobs pass, update the existing
[polyfempy feedstock](https://github.com/conda-forge/polyfempy-feedstock).
Replace the local source parameters with the published archive URL and checksum,
preserve the reviewed build and test requirements, and rerender the feedstock
with its normal conda-forge tooling. Review native dependency and license
requirements there before publishing. A successful local build alone does not
complete the conda-forge release.
