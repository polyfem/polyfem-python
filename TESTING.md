# Testing the forward API

Run development tests from a Git checkout after initializing its submodules.
The source archive omits the repository test suite and its development-only
inputs. It retains the installed-package check and quick-start examples.
See [BUILD.md](BUILD.md) for the native build.

## Python and example tests

Generate the Python API before running the tests:

```sh
python -m pip install numpy pytest "setuptools>=77"
python tools/generate_polyfem_api.py
python -m pytest -q -rs
```

The pytest configuration selects `tests/` and adds the public package directory
to the import path. Submodule tests run separately.
The suite includes comparisons between classic example configurations and
official JSON inputs. This forward branch does not require PyTorch.

On a backend-free machine, the native smoke test is expected to skip. Windows
also skips the tests for the non-Windows runtime branch. Check the reasons with
`-rs`; a skipped native test does not establish that the solver works.

The generator has its own standalone suite. From `python-from-jse/`, run
`python tools/regenerate_and_test.py`. CI uses this same command to generate
the generic test files before compiling and testing them. Generating the parent
package API does not create these files in `python-from-jse/generated/`.
The parent wrapper's `--check` option runs the parent API/example checks,
not that standalone suite.

## Native solver

With a wheel installed in a clean environment, run the README beam,
installation self-check, and Neo-Hookean cube reference check:

```sh
python tools/check_installed_package.py --readme README.md --reference polyfem-data/units/neohookean.json --example examples/quickstart/beam.py --example examples/quickstart/static_elasticity.py --example examples/quickstart/transient_elasticity.py
```

This check runs outside the checkout and fails if the backend cannot load.
It checks the README beam and the small square in `INSTALL.md` against their
expected results, repeats the square with a Python load function, checks all six
cube metrics, then runs the three quick-start
scripts. The beam mesh is copied from `examples/quickstart/beam.obj`; the check
does not download it. An editable install is not accepted as a wheel validation.

For a development install, `python -m pytest tests/test_backend_smoke.py -q`
remains available. It may skip when the backend is absent; that is not release
evidence.

## Code checks

```sh
python -m pip install "ruff>=0.17.5,<0.18"
python -m ruff check python setup.py tools tests
```

These checks cover syntax, undefined names, and unused imports in handwritten
Python code. Generated output and external submodules are excluded.

Check documentation links, Python snippet syntax, and built archive contents:

```sh
python tools/check_docs.py --online
python tools/check_distribution.py "dist/*.tar.gz" "dist/*.whl"
```

Omit `--online` to check local links without network access. The documentation
check rejects relative links that cross a Git submodule boundary, even when the
file exists locally. Archive checks require runtime files and reject repository
tests, build output, and caches. Use only the archive patterns you have built.
The Test workflow checks documentation and a fresh source archive on every push.

## GitHub Actions

| When | What runs |
| --- | --- |
| Push or pull request | Python, generated API, and generator unit tests on Windows, macOS and Linux with Python 3.10-3.12 |
| Manual release-package check before publishing | Source archive, pip wheels, native conda packages, and installed-package checks across the release matrix |
| Draft PR, or PR targeting another branch | Fast Python and documentation checks only |
| Ready-for-review PR targeting `main` | Clean wheel installation, beam, installation self-check, cube and quick-start checks, then the full contact sweep on open, update, reopen, marking ready or retargeting to `main` |
| Manual run of Generated Contact Backend, or push to main | The same wheel and simulation checks, then the full contact sweep |

Before merging, check **Forward solver checks (Ubuntu, Python 3.11)** on the latest PR revision. The full sweep runs the 75 active cases and compares the generated configurations with the official inputs and their reference metrics. For an on-demand check, open **Actions -> Generated Contact Backend -> Run workflow** and select the branch. To run the same check locally:

```sh
python tools/run_generated_contact_backend_checks.py --output-root build/contact-check --generated-source-tolerance 1e-5 --max-threads 1 --require-tests-match
```

Read `build/contact-check/summary.txt` for the outcome and the per-case logs for failures. In GitHub Actions, download the **generated-contact-backend-results** artifact from the workflow run to inspect those files.

The sweep follows each reference input's `tests.time_steps`; a transient case
without that field is reduced to one step. It also selects a portable linear
solver for both runs, generally `Eigen::SimplicialLDLT`. Passing the sweep checks
these regression configurations, not every example's original full-duration
command or its original solver choice.

## Release-package checks

The **Release Packages** workflow builds one source archive, then uses
it for pip wheels and native conda packages. The initial installation checks
use CPython 3.11 on Windows x64, Linux x64, Intel macOS, and Apple Silicon
macOS. This preserves the legacy conda package's platforms and adds native
support for M-series Macs. There are nine jobs: one source archive,
four wheel builds, and four conda builds.
Expand and validate the matrix for Python 3.10 and Python 3.12
before publishing packages for those targets.
Each job installs its output and runs the README beam, installation self-check,
Python load function, and all three quick-start scripts. Each simulation has a
five-minute timeout so a stalled native call fails the check.
Wheel builds use cibuildwheel's dependency repair step before installation.

Pushes and pull requests do not run the package matrix. Use **Actions -> Release Packages
-> Run workflow** when checking installation packages before a release. Manual
runs become available once the workflow exists on the repository's default
branch; see
[GitHub's manual-run instructions](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow).

Before that workflow reaches the default branch, use the local build and
installation commands above. No separate validation branch is needed.
Review every job for the exact revision being released.
Successful jobs produce downloadable artifacts without publishing to package indexes. See
[INSTALL.md](INSTALL.md) for the verified user installation routes and
[conda/README.md](conda/README.md) for the native conda recipe.

The separate contact-backend workflow builds a source archive and wheel on Linux. It
installs the wheel in a fresh environment and checks it outside the checkout.
This does not establish native wheel support on Windows or macOS. Before
publishing, validate each intended platform and Python version:

1. Confirm `import polyfem`, `from polyfempy import solve`, and backend loading.
2. Run the README beam and installation self-check against their expected values.
3. Run the cube reference check with the pinned data and inspect all six metrics.
4. Run the required contact sweep on the same release revision.
5. If publishing a source archive, build a wheel from that extracted archive.
   See [BUILD.md](BUILD.md#source-archives) for the archive contents.

A configured GitHub workflow is a test plan; review a successful run for the
actual release commit before treating it as release evidence.
