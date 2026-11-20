# Testing

Install development dependencies through the project:

```sh
python -m pip install -e ".[dev]"
python -m pytest tests -q
python -m ruff check python tools tests
```

The normal installation generates the configuration API automatically. After
changing schemas or generator configuration, regenerate it by
reinstalling the project.

Python tests cover configuration normalization, model selection IDs, array
validation, and solver dispatch. Native array-selection tests run when the
extension is available. Reference examples require the pinned `examples` and
`polyfem-data` submodules. Their historical import name is adapted only inside
tests; the installed public package is `polyfempy`.

To test an installed wheel outside the source directory:

```sh
python tools/check_installed_package.py
```

This checks imports and compares a file-mesh solve with the equivalent NumPy
array solve, including material and boundary selections.

For a selected existing native build and its reference data:

```sh
python -m pytest tests/test_backend_smoke.py --backend-package-dir /path/to/site-packages -s
```

Full contact regressions are an explicit test run:

```sh
python -m pytest tests/test_contact_regressions.py -m contact -q
```

CI has one workflow that installs the project and runs Python tests, lint, and
installed native checks on Ubuntu, macOS, and Windows.
The separate Release Packages workflow builds a source archive and platform
wheels with cibuildwheel. It checks archive contents and each installed wheel
without publishing them.
