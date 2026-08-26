from __future__ import annotations

from pathlib import Path

import yaml


def load_config(path: str | Path) -> dict:
    """Load a YAML experiment config (see configs/default.yaml)."""
    with open(path) as fh:
        return yaml.safe_load(fh)
