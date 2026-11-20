"""Opt-in native parity checks for the backend's active contact cases."""

from functools import lru_cache
import json
from pathlib import Path
import subprocess
import sys

import pytest

from polyfempy.solve import _config
from test_generated_api_example import _all_example_cases, _import_example, _load_source_config

ROOT = Path(__file__).resolve().parents[1]
BACKEND = None
for cache in sorted((ROOT / "build").glob("*/CMakeCache.txt"), key=lambda p: p.stat().st_mtime, reverse=True):
    for line in cache.read_text().splitlines():
        if line.startswith("CPM_PACKAGE_polyfem_SOURCE_DIR:INTERNAL="):
            BACKEND = Path(line.partition("=")[2])
            break
    if BACKEND is not None:
        break
CASES = []
for name in ("contact_2d.txt", "contact_3d.txt"):
    listing = BACKEND / "tests" / name if BACKEND else ROOT / "missing-backend"
    if listing.exists():
        CASES.extend(line.strip() for line in listing.read_text().splitlines()
                     if line.strip().startswith("contact/examples/"))


@lru_cache(maxsize=1)
def example_map():
    return {source: example for example, source in _all_example_cases()}


RUN = """
import json, sys
from pathlib import Path
from polyfempy import solve
from polyfempy.polyfempy import Solver
workspace = Path(sys.argv[1])
for label in ('generated', 'source'):
    payload = json.loads((workspace / (label + '.json')).read_text())
    if label == 'generated':
        solve(cfg=payload)
    else:
        solver = Solver()
        solver.set_settings(json.dumps(payload), strict_validation=True)
        solver.load_mesh_from_settings()
        solver.solve()
"""


@pytest.mark.contact
@pytest.mark.parametrize("case", CASES, ids=lambda case: case.removeprefix("contact/examples/"))
def test_contact_api_matches_source_json(case, tmp_path):
    source = (ROOT / "polyfem-data" / case).resolve()
    if not source.exists():
        pytest.skip("Initialize polyfem-data for native contact regressions")
    probe = subprocess.run([sys.executable, "-I", "-c", "from polyfempy.polyfempy import Solver"],
                           capture_output=True, text=True)
    assert probe.returncode == 0, "Install the native package before running contact regressions: " + probe.stderr
    generated = _import_example(example_map()[source]).polyfem_config.as_dict()
    reference = _load_source_config(source)
    steps = json.loads(source.read_text()).get("tests", {}).get("time_steps", 1)
    for label, payload in (("generated", generated), ("source", reference)):
        payload = _config(dict(payload, root_path=str(source)))
        time = payload.get("time")
        if time:
            dt = time.get("dt")
            if dt is None:
                dt = (time["tend"] - time.get("t0", 0)) / time["time_steps"]
            time.update(dt=dt, time_steps=steps, tend=time.get("t0", 0) + steps * dt)
        payload.setdefault("solver", {}).update(linear={"solver": "Eigen::SimplicialLDLT"}, max_threads=1)
        payload["output"] = {"directory": str(tmp_path / label), "json": "sim.json"}
        (tmp_path / (label + ".json")).write_text(json.dumps(payload))
    result = subprocess.run([sys.executable, "-I", "-c", RUN, str(tmp_path)],
                            capture_output=True, text=True, timeout=360, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    generated_metrics = json.loads((tmp_path / "generated/sim.json").read_text())
    source_metrics = json.loads((tmp_path / "source/sim.json").read_text())
    keys = [key for key in source_metrics if key.startswith("err_")]
    assert keys, "Native solver did not report error metrics"
    for key in keys:
        assert generated_metrics[key] == pytest.approx(source_metrics[key], rel=1e-5, abs=1e-12)
