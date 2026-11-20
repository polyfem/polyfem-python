"""Options for tests that use an existing native build."""


def pytest_addoption(parser):
    parser.addoption(
        "--backend-package-dir",
        help="Package root containing polyfempy and the native backend to test.",
    )
