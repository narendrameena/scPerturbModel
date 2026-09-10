"""Read a perturbation atlas's DESIGN from its metadata, without its counts.

The design calculation needs three integers -- contexts, perturbations,
replicates per condition -- and nothing else. Those live in an ``.h5ad``'s ``obs``
table, so they can be read from a 500 MB file in under a second without touching
the expression matrix. That is what makes it practical to score a whole
collection: RESULTS.md sec.48 applies this to all 38 scPerturb datasets.

Two decisions in here carry the weight, and both are about what does NOT count.

**A replicate must be independently treated.** ``REP`` excludes ``nperts`` (which
counts perturbations per cell, not replicates) and bare ``well``/``lane``, which
in 10x datasets label a capture rather than an independent treatment. Splitting
one well's cells in half yields two measurements sharing that well's dose, batch
and dropout; their covariance carries the shared noise straight into the
interaction term. That is the confound sec.38 found inflating published indices,
and it is why this module counts annotated replicates rather than cells.

**A missing label is not a level.** Unlabelled cells surface as ``''`` when a
categorical code is -1, and some datasets write the literal string ``'None'``.
Counting either inflates the design -- sci-Plex 3 reads as 4 cell lines and 3
replicates when it has 3 and 2, overstating its pair count by half.
"""
from __future__ import annotations

import warnings
from pathlib import Path

__all__ = ["CTX", "PERT", "REP", "MISSING", "read_obs", "pick", "nlev",
           "describe"]

# candidate obs columns, in the order scPerturb's harmonised schema prefers
CTX = ["cell_line", "celltype", "cell_type", "tissue_type", "tissue", "organ",
       "disease"]
PERT = ["perturbation", "perturbation_1", "guide_id", "target", "drug"]
REP = ["replicate", "batch", "sample", "donor", "hto", "library", "plate"]

MISSING = {"", "nan", "none", "na", "null", "n/a", "unknown", "-1"}


def read_obs(path):
    """The ``obs`` table alone -- no expression matrix is loaded."""
    import h5py
    import numpy as np
    import pandas as pd

    out = {}
    with h5py.File(path, "r") as f:
        if "obs" not in f:
            return None
        g = f["obs"]
        for k in g.keys():
            try:
                node = g[k]
                if isinstance(node, h5py.Dataset):
                    v = node[...]
                elif "categories" in node and "codes" in node:
                    cats = node["categories"][...]
                    codes = node["codes"][...]
                    cats = np.array([c.decode() if isinstance(c, bytes) else c
                                     for c in cats])
                    v = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
                else:
                    continue
                if v.ndim != 1:
                    continue
                if v.dtype.kind == "S":
                    v = np.array([x.decode(errors="ignore") for x in v])
                out[k] = v
            except Exception:
                continue
    if not out:
        return None
    n = min(len(v) for v in out.values())
    return pd.DataFrame({k: v[:n] for k, v in out.items()})


def pick(cols, prefer):
    """First of ``prefer`` present in ``cols``, matched case-insensitively."""
    low = {c.lower(): c for c in cols}
    for p in prefer:
        if p in low:
            return low[p]
    return None


def nlev(s):
    """Distinct real levels of a column, missing markers excluded."""
    v = s.astype(str).str.strip()
    return v[~v.str.lower().isin(MISSING)]


def describe(path):
    """One atlas's design, or ``{'dataset', 'error'}`` if it cannot be read."""
    path = Path(path)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            o = read_obs(path)
    except Exception as e:
        return {"dataset": path.stem, "error": type(e).__name__}
    if o is None or not len(o):
        return {"dataset": path.stem, "error": "no obs"}
    cc = pick(o.columns, CTX)
    pc = pick(o.columns, PERT)
    rc = pick(o.columns, REP)
    if pc is None:
        return {"dataset": path.stem, "error": "no perturbation column"}
    # An unlabelled context column means one unnamed context, not zero.
    n_ctx = max(int(nlev(o[cc]).nunique()), 1) if cc else 1
    n_pert = max(int(nlev(o[pc]).nunique()), 1)
    # Replicates PER CONDITION -- not the number of batches in the dataset. An
    # atlas with 8 batches spread across 8 disjoint conditions has 8 batches and
    # one replicate, and contributes no pair.
    if rc is not None:
        w = o.copy()
        w["_r"] = o[rc].astype(str).str.strip()
        w = w[~w["_r"].str.lower().isin(MISSING)]
        keys = [c for c in (cc, pc) if c]
        if not len(w):
            n_rep = 1
        elif keys:
            n_rep = int(w.groupby(keys, observed=True)["_r"].nunique().median())
        else:
            n_rep = int(w["_r"].nunique())
    else:
        n_rep = 1
    # A CRISPR screen in one cell line is single-context BY DESIGN and it would
    # be unfair to score it as a failed atlas; a drug screen in one cell line is
    # a different statement. The breakdown separates them.
    ptype = "unknown"
    for c in ("perturbation_type", "perturbation_type_1"):
        if c in o.columns:
            v = nlev(o[c])
            if len(v):
                ptype = str(v.mode().iloc[0])
            break
    # EXACT pair count, not the one implied by the median replicate count. An
    # unbalanced atlas -- most conditions measured once, some many times -- has a
    # median of 1 and therefore an implied zero pairs, while actually carrying
    # thousands. Ten of the 38 scPerturb datasets are in exactly that state, with
    # up to 18,007 real pairs, and the tool told each of them "UNRESOLVABLE AT ANY
    # EFFECT SIZE". The median is kept for reporting; `n_pairs_exact` is what the
    # calculation should use when it is available.
    n_pairs_exact = None
    if rc is not None and len(w) and keys:
        per = w.groupby(keys, observed=True)["_r"].nunique()
        n_pairs_exact = int((per * (per - 1) // 2).sum())
    return {"dataset": path.stem, "n_cells": int(len(o)),
            "context_col": cc, "pert_col": pc, "rep_col": rc,
            "pert_type": ptype,
            "n_contexts": n_ctx, "n_perturbations": n_pert,
            "n_replicates": max(n_rep, 1),
            "n_pairs_exact": n_pairs_exact}
