"""Standard figure export: every figure is saved as a self-contained bundle.

results/figures/<name>/
    <name>.png  <name>.svg  <name>.pdf     the figure in all three formats
    <key>.csv                              source data behind the figure
    <script name>.py                       copy of the script that created it
"""
from __future__ import annotations

import shutil
from pathlib import Path


def save_figure(fig, name: str, out_root: str | Path,
                source_data=None, script: str | Path | None = None,
                dpi: int = 150) -> Path:
    """Save `fig` under out_root/<name>/ as png+svg+pdf with source data + script.

    source_data: a DataFrame, or a dict {stem: DataFrame} for multiple tables.
    script: pass __file__ from the calling script.
    """
    d = Path(out_root) / name
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(d / f"{name}.{ext}", dpi=dpi)
    if source_data is not None:
        tables = source_data if isinstance(source_data, dict) else {"source_data": source_data}
        for key, df in tables.items():
            df.to_csv(d / f"{key}.csv", index=False)
    if script is not None:
        shutil.copy2(script, d / Path(script).name)
    return d
