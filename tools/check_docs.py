"""Check project Markdown links and Python syntax without building the solver."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import unquote, urldefrag, urlsplit
from urllib.request import Request, urlopen


def headings(text: str) -> set[str]:
    """Return GitHub-style anchors for the plain headings used in these guides."""
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    counts: Counter[str] = Counter()
    anchors = set()
    for title in re.findall(r"^#{1,6}\s+(.+)$", text, re.MULTILINE):
        slug = re.sub(r"[^\w\- ]", "", title.strip().lower()).replace(" ", "-")
        anchors.add(f"{slug}-{counts[slug]}" if counts[slug] else slug)
        counts[slug] += 1
    return anchors


def check_markdown(path: Path, repo: Path, submodules: list[Path]) -> set[str]:
    text = path.read_text(encoding="utf-8")
    for code in re.findall(r"```python\s*\n(.*?)```", text, re.DOTALL):
        ast.parse(code, filename=str(path))
    prose = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    urls = set()
    for target in re.findall(r"\[[^\]]*\]\(([^)\s]+)\)", prose):
        url, fragment = urldefrag(target)
        scheme = urlsplit(url).scheme
        if scheme in ("http", "https"):
            urls.add(url)
            continue
        if scheme == "mailto":
            continue
        resolved = (path.parent / unquote(url)).resolve() if url else path.resolve()
        if not resolved.is_relative_to(repo.resolve()):
            raise ValueError(f"{path}: link leaves its repository: {target}")
        if any(resolved.is_relative_to(sub.resolve()) for sub in submodules):
            raise ValueError(f"{path}: use a pinned GitHub URL for submodule link: {target}")
        if not resolved.exists():
            raise ValueError(f"{path}: link target is missing: {target}")
        if fragment and resolved.suffix == ".md":
            if unquote(fragment) not in headings(resolved.read_text(encoding="utf-8")):
                raise ValueError(f"{path}: heading is missing: {target}")
    return urls


def check_url(url: str) -> None:
    for attempt in range(3):
        try:
            request = Request(url, headers={"User-Agent": "PolyFEM-doc-check"})
            with urlopen(request, timeout=30) as response:
                if response.status != 200:
                    raise ValueError(f"HTTP {response.status}")
            return
        except (OSError, ValueError) as exc:
            if attempt == 2:
                raise RuntimeError(f"Could not read documentation link {url}: {exc}") from exc
            time.sleep(1 + attempt)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--online", action="store_true", help="Also request public HTTP links")
    args = parser.parse_args()
    root = args.root.resolve()
    submodules = [root / name for name in ("examples", "python-from-jse", "polyfem", "polyfem-data")]
    urls = set()
    count = 0
    for repo in [root, root / "examples", root / "python-from-jse"]:
        names = subprocess.check_output(
            ["git", "ls-files", "-z", "*.md"], cwd=repo
        ).decode("utf-8").split("\0")
        for name in filter(None, names):
            urls.update(check_markdown(repo / name, repo, submodules if repo == root else []))
            count += 1
    if args.online:
        for url in sorted(urls):
            check_url(url)
        print(f"Checked {len(urls)} public URLs.")
    print(f"Checked links and Python syntax in {count} Markdown files.")


if __name__ == "__main__":
    main()
