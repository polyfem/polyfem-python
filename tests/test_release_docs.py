from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_local_doc_folder_is_not_tracked_release_content():
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")

    assert "/doc/" in gitignore
