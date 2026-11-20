from __future__ import annotations

from pathlib import Path


def test_removed_output_adapter_module_is_absent():
    package_dir = Path(__file__).resolve().parents[1] / "python" / "polyfempy"

    assert not (package_dir / "runtime" / "_solve_outputs.py").exists()


def test_solve_module_does_not_reexport_output_helpers():
    import importlib

    pipeline = importlib.import_module("polyfempy.runtime.solve")

    assert not hasattr(pipeline, "NativeOutputs")
    assert not hasattr(pipeline, "extract_native_outputs")
    assert not hasattr(pipeline, "populate_sampled_history")
    assert not hasattr(pipeline, "apply_sampled_vtu_fallback")
    assert not hasattr(pipeline, "finalize_result")
