"""Shared pseudobulk-delta construction, additive priors, and metrics.

Single source of truth for Phase 2/3 model comparisons: identical deltas
(log1p CPM-1e6, plate-matched, cell-weighted, plate14 excluded), identical
gene sets, identical metrics. Mirrors scripts/eval_baselines.py conventions.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

MIN_CELLS = 50
KEY = ["cell_line_id", "drug", "conc"]


def load_pseudobulk(pb_dir: str | Path):
    """Return (lognorm matrix [conds x genes], condition DataFrame with 'row')."""
    pb_dir = Path(pb_dir)
    counts = np.load(pb_dir / "pseudobulk_counts.npz")["counts"]
    cond = pd.read_csv(pb_dir / "conditions.csv").reset_index().rename(
        columns={"index": "row"})
    X = np.log1p(counts / (counts.sum(axis=1, keepdims=True) + 1e-9) * 1e6)
    return X, cond


def build_deltas(X: np.ndarray, cond: pd.DataFrame,
                 min_cells: int = MIN_CELLS,
                 exclude_plates: tuple = ("plate14",),
                 keep_plate: bool = False):
    """Plate-matched, cell-weighted deltas per (line, drug, conc) triple.

    Returns (G, DELTA): G is a DataFrame with KEY columns (reset index aligns
    with DELTA rows), DELTA is float32 [n_triples x n_genes].

    With keep_plate=True the per-plate rows are returned instead of the
    cell-weighted average over plates, and G carries a 'plate' column. Needed
    to test whether within-(line, drug) agreement is inflated by shared plate
    effects rather than real biology.
    """
    cond = cond[~cond.plate.isin(exclude_plates)]
    ctrl = cond[(cond.drug == "DMSO_TF") & (cond.n_cells >= min_cells)]
    ctrl_idx = {(r.cell_line_id, r.plate): r.row for r in ctrl.itertuples()}

    rows, keys, weights, plates = [], [], [], []
    treated = cond[(cond.drug != "DMSO_TF") & (cond.n_cells >= min_cells)]
    for r in treated.itertuples():
        ci = ctrl_idx.get((r.cell_line_id, r.plate))
        if ci is None:
            continue
        keys.append((r.cell_line_id, r.drug, float(r.conc)))
        weights.append(r.n_cells)
        plates.append(r.plate)
        rows.append(X[r.row] - X[ci])
    key_df = pd.DataFrame(keys, columns=KEY)
    key_df["w"] = weights
    deltas_all = np.stack(rows)

    if keep_plate:
        out = key_df[KEY].copy()
        out["plate"] = plates
        return out.reset_index(drop=True), deltas_all.astype(np.float32)

    key_df["gid"] = key_df.groupby(KEY, observed=True).ngroup()
    n = key_df.gid.nunique()
    W = np.zeros((n, deltas_all.shape[1]))
    wsum = np.zeros(n)
    for i, (gid, w) in enumerate(zip(key_df.gid, key_df.w)):
        W[gid] += w * deltas_all[i]
        wsum[gid] += w
    DELTA = (W / wsum[:, None]).astype(np.float32)
    G = (key_df.drop_duplicates("gid").sort_values("gid")[KEY]
         .reset_index(drop=True))
    return G, DELTA


def additive_prior(G: pd.DataFrame, DELTA: np.ndarray,
                   train_mask: np.ndarray, loo: bool = False) -> np.ndarray:
    """Additive-shift prior per row of G: mean TRAIN delta of the same
    (drug, conc). With loo=True a train row's own delta is excluded from its
    prior (prevents target leakage when the prior is a model input).
    Rows with no available train delta get a zero prior."""
    sums: dict = {}
    counts: dict = {}
    for i in np.where(train_mask)[0]:
        k = (G.drug[i], G.conc[i])
        sums[k] = sums.get(k, 0.0) + DELTA[i]
        counts[k] = counts.get(k, 0) + 1
    P = np.zeros_like(DELTA)
    for i in range(len(G)):
        k = (G.drug[i], G.conc[i])
        if k not in sums:
            continue
        s, c = sums[k], counts[k]
        if loo and train_mask[i]:
            if c > 1:
                P[i] = (s - DELTA[i]) / (c - 1)
        else:
            P[i] = s / c
    return P


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a - a.mean(), b - b.mean()
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / den) if den > 0 else np.nan


def responsive_genes(DELTA: np.ndarray, train_mask: np.ndarray,
                     n: int = 2000) -> np.ndarray:
    return np.sort(np.argsort(DELTA[train_mask].var(axis=0))[-n:])


def delta_metrics(pred: np.ndarray, true: np.ndarray,
                  resp: np.ndarray) -> dict:
    de100 = np.argsort(np.abs(true))[-100:]
    return {
        "r_hvg": pearson(pred[resp], true[resp]),
        "r_de100": pearson(pred[de100], true[de100]),
        "rmse_hvg": float(np.sqrt(np.mean((pred[resp] - true[resp]) ** 2))),
    }
