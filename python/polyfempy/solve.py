"""Forward solves from model configurations, dictionaries, or JSON files."""

import json
from pathlib import Path

import numpy as np

__all__ = ["Result", "solve"]


class Result:
    """Raw finite-element solution coefficients and optional metadata."""

    def __init__(self, sol, *, meta=None):
        self.sol = np.ascontiguousarray(sol)
        self.meta = {} if meta is None else dict(meta)


def _config(cfg):
    if isinstance(cfg, (str, Path)):
        path = Path(cfg).resolve()
        cfg = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(cfg, dict):
            cfg = dict(cfg, root_path=cfg.get("root_path", str(path)))
    elif callable(getattr(cfg, "as_dict", None)):
        cfg = cfg.as_dict()
    if not isinstance(cfg, dict):
        raise TypeError("cfg must be a model configuration, dictionary, or JSON path")
    return cfg


def solve(cfg):
    """Solve a model configuration, dictionary, or JSON file using PolyFEM."""
    settings = _config(cfg)
    from .polyfempy import Solver

    solver = Solver()
    solver.set_settings(json.dumps(settings), strict_validation=True)
    solver.load_mesh_from_settings()
    return Result(solver.solve())
