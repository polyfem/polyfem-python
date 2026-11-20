"""Keep generated exports deterministic and linked schemas pinned."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def workflow():
    spec = importlib.util.spec_from_file_location("generate_api", ROOT / "tools/generate_polyfem_api.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_schema_cache_tracks_backend_revision(tmp_path, monkeypatch):
    module = workflow()
    monkeypatch.setattr(module, "POLYSOLVE_SPEC_CACHE_DIR", tmp_path / "cache")
    recipe = tmp_path / "polyfem/cmake/recipes/polysolve.cmake"
    recipe.parent.mkdir(parents=True)
    downloads = []
    def download(url, target):
        downloads.append(url)
        Path(target).write_text("{}")
    monkeypatch.setattr(module.urllib.request, "urlretrieve", download)
    for revision in ("abc123", "def456", "abc123"):
        recipe.write_text(f'CPMAddPackage("gh:polyfem/polysolve#{revision}")')
        assert module._cache_polysolve_solver_specs(tmp_path / "polyfem") == tmp_path / "cache" / revision
    assert len(downloads) == 4
    assert all("/abc123/" in url or "/def456/" in url for url in downloads)


def test_generated_exports_hide_implementation_imports_and_are_idempotent(tmp_path, monkeypatch):
    module = workflow()
    monkeypatch.setattr(module, "GENERATED_DIR", tmp_path)
    api = tmp_path / "generated_api.py"
    api.write_text("from pathlib import Path\nfrom .generated_class import Root\ndef model(): pass\ndef _private(): pass\n")
    module.write_public_exports()
    first = api.read_text()
    module.write_public_exports()
    assert api.read_text() == first
    assert "__all__ = ['Root', 'model']" in first
