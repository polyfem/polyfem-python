from pathlib import Path
import json
import subprocess
import sys

import pytest


_BACKEND_RUNNER = """
import hashlib
import json
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
import sys

package_dir = Path(sys.argv[1]).resolve()
expected_package = package_dir / "polyfempy" / "__init__.py"
if not expected_package.is_file():
    raise RuntimeError(f"Selected package is missing: {expected_package}")
sys.path.insert(0, str(package_dir))
import polyfempy as pf

package_path = Path(pf.__file__).resolve()
if package_path != expected_package.resolve():
    raise RuntimeError(f"Loaded package {package_path}; expected {expected_package}")
if not pf.cpp_backend_available():
    raise RuntimeError(f"Native backend is unavailable: {pf.cpp_backend_error()}")
native_path = Path(pf._core.__file__).resolve()
if native_path.parent != expected_package.parent.resolve():
    raise RuntimeError(f"Loaded backend outside selected package: {native_path}")
if not any(native_path.name.endswith(suffix) for suffix in EXTENSION_SUFFIXES):
    raise RuntimeError(f"Expected a compiled native extension: {native_path}")

digest = hashlib.sha256()
with native_path.open("rb") as stream:
    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(chunk)
info = {"package": str(package_path), "native": str(native_path), "sha256": digest.hexdigest()}
expected = json.loads(sys.argv[2])
if expected is not None and info != expected:
    raise RuntimeError(f"Native backend changed between probe and solve: expected {expected}, got {info}")
print(json.dumps(info), flush=True)

if len(sys.argv) > 3:
    from polyfempy.cli import main
    raise SystemExit(main(sys.argv[3:]))
"""


def _run_backend(package_dir, args, *, cwd, expected_backend=None):
    # -I ignores PYTHONPATH, the working directory, and user site-packages.
    return subprocess.run(
        [
            sys.executable, "-I", "-c", _BACKEND_RUNNER,
            str(Path(package_dir).resolve()), json.dumps(expected_backend), *args,
        ],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def _load_expected_error_metrics(config_path: Path) -> tuple[dict[str, float], float]:
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    tests = payload.get("tests")
    if not isinstance(tests, dict):
        raise AssertionError(f"backend smoke config has no tests block: {config_path}")

    margin = float(tests.get("margin", 1e-8))
    expected = {
        key: float(value)
        for key, value in tests.items()
        if key.startswith("err_") and isinstance(value, (int, float))
    }
    if not expected:
        raise AssertionError(f"backend smoke config has no err_* test metrics: {config_path}")
    return expected, margin


def _load_actual_error_metrics(output_dir: Path) -> dict[str, float]:
    stats_path = output_dir / "sim.json"
    if not stats_path.exists():
        raise AssertionError(f"backend did not write expected stats file: {stats_path}")

    payload = json.loads(stats_path.read_text(encoding="utf-8"))
    return {
        key: float(value)
        for key, value in payload.items()
        if key.startswith("err_") and isinstance(value, (int, float))
    }


def test_backend_forward_solve_smoke(tmp_path, pytestconfig):
    package_dir = pytestconfig.getoption("backend_package_dir")
    if package_dir is None:
        pytest.skip("Select the build to test with --backend-package-dir; no backend is chosen automatically.")
    root = Path(__file__).resolve().parents[1]
    backend_info = _run_backend(package_dir, [], cwd=tmp_path)
    assert backend_info.returncode == 0, backend_info.stdout + backend_info.stderr
    selected_backend = json.loads(backend_info.stdout)
    print(f"Native backend: {selected_backend['native']}")
    print(f"SHA-256: {selected_backend['sha256']}")

    config_path = root / "polyfem-data" / "units" / "neohookean.json"
    assert config_path.is_file(), f"Backend smoke config is missing: {config_path}. Initialize polyfem-data."
    expected_metrics, margin = _load_expected_error_metrics(config_path)

    output_dir = tmp_path / "backend_smoke"
    result = _run_backend(
        package_dir,
        [
            "solve",
            str(config_path),
            "--output-dir",
            str(output_dir),
            "--log-level",
            "2",
            "--max-threads",
            "1",
        ],
        cwd=tmp_path,
        expected_backend=selected_backend,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert output_dir.exists()
    assert any(output_dir.iterdir())

    actual_metrics = _load_actual_error_metrics(output_dir)
    for key, expected_value in expected_metrics.items():
        assert key in actual_metrics
        assert abs(actual_metrics[key] - expected_value) <= margin
