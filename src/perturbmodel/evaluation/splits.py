"""Train/test splits shared by every model and baseline.

The split is defined at the (cell_line_id, drug, conc) triple level and is
deterministic given (seed, test_frac) regardless of row order — triples are
canonically sorted before drawing, so pseudobulk baselines and single-cell
models hold out the SAME conditions.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEY = ["cell_line_id", "drug", "conc"]


def held_out_condition_triples(df: pd.DataFrame, seed: int = 0,
                               test_frac: float = 0.2) -> set[tuple]:
    """Return the set of test (cell_line_id, drug, conc) triples.

    Controls (drug == DMSO_TF) are never held out. A drawn triple is kept in
    the test set only if its (drug, conc) is still observed in at least one
    training line — the additive-shift baseline and context-transfer models
    need the perturbation observed somewhere.
    """
    t = (df.loc[df.drug != "DMSO_TF", KEY].drop_duplicates()
         .sort_values(KEY).reset_index(drop=True))
    rng = np.random.default_rng(seed)
    test = rng.random(len(t)) < test_frac
    train_pairs = set(zip(t.drug[~test], t.conc[~test]))
    keep = test & np.array([(d, c) in train_pairs
                            for d, c in zip(t.drug, t.conc)])
    return {tuple(r) for r in t[keep].itertuples(index=False)}


def split_column(df: pd.DataFrame, test_triples: set[tuple],
                 seed: int = 0, val_frac: float = 0.1) -> pd.Series:
    """Assign 'train' / 'val' / 'test' per row of df (any granularity).

    Rows whose triple is in test_triples -> 'test'. Remaining rows are split
    row-wise into train/val (controls always follow train/val, never test).
    """
    keys = list(zip(df.cell_line_id, df.drug, df.conc))
    is_test = np.array([k in test_triples for k in keys])
    rng = np.random.default_rng(seed + 1)
    is_val = (rng.random(len(df)) < val_frac) & ~is_test
    out = np.where(is_test, "test", np.where(is_val, "val", "train"))
    return pd.Series(out, index=df.index, name="split")
