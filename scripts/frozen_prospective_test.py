#!/usr/bin/env python3
"""FROZEN at the pre-registration of 2026-09-07. Do not modify.

`docs/PREREGISTRATION.md` Part B commits to running this file, unmodified, on the
first three qualifying perturbation atlases released on or after 2026-09-07. Its
SHA-256 is recorded in that document, so any edit is detectable. If a future data
format makes an edit unavoidable, the diff goes in
`docs/PREREGISTRATION_OUTCOMES.md` and the run is flagged as modified.

The point of freezing it is that the analysis cannot be adjusted once the answer
is visible. Every choice that could be tuned after the fact -- the noise constant,
the permutation count, the definition of a replicate, the rule for calling an
interaction resolved -- is fixed here as a constant rather than a flag.

    python scripts/frozen_prospective_test.py ATLAS.h5ad --out outcome.json

It prints, and writes as JSON:

  * the design read from the object's own metadata;
  * the detection floor the calculation predicts from that design ALONE, computed
    and printed BEFORE the interaction is estimated, so the prediction is on the
    record before its test;
  * the measured interaction share and a permutation test of it;
  * whether the prediction `share > floor <=> permutation rejects` held.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

# ---- frozen constants. Changing any of these invalidates the registration. ----
SNR = 0.20            # shared noise constant of RESULTS.md sec.46
N_FEAT = 2000         # features entering the averaged covariance
N_PERM = 200          # permutation draws, as in sec.35 and sec.39
ALPHA = 0.05          # rejection threshold for the permutation test
MIN_CELLS = 50        # fewest cells for a condition to be pseudobulked
CONTROL_HINTS = ("dmso", "control", "ctrl", "untreated", "vehicle", "none",
                 "non-targeting", "nontargeting", "scrambled", "safe-harbour")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    src = Path(__file__).resolve()
    print(f"frozen script sha256 = {sha256(src)}")
    print(f"  (must match docs/PREREGISTRATION.md; any mismatch means this run "
          f"is modified)\n")

    from perturbmodel.atlas_meta import describe
    from perturbmodel.design import min_detectable_share, n_pairs

    # ---------- 1. the prediction, from the design alone ----------
    d = describe(a.h5ad)
    if "error" in d:
        print(f"cannot read design: {d['error']}")
        return 2
    floor = min_detectable_share(d["n_contexts"], d["n_perturbations"],
                                 d["n_replicates"], SNR, N_FEAT)
    pairs = n_pairs(d["n_contexts"], d["n_perturbations"], d["n_replicates"])
    print("PREDICTION (from design only, before any expression is read)")
    print(f"  {d['n_contexts']} contexts x {d['n_perturbations']} perturbations "
          f"x {d['n_replicates']} replicates")
    print(f"  cross-replicate pairs : {pairs:,}")
    print(f"  detection floor       : {floor:.4f}")
    eligible = d["n_contexts"] >= 2 and d["n_perturbations"] >= 20
    print(f"  eligible under Part B : {eligible}")
    if pairs == 0:
        print("\n  PREDICTS: not separable from zero at any effect size "
              "(assertion 4).")
    print()

    out = {"dataset": d["dataset"], "design": d, "n_pairs": pairs,
           "predicted_floor": floor, "snr": SNR, "eligible": eligible,
           "script_sha256": sha256(src)}

    # ---------- 2. the measurement ----------
    try:
        import anndata as ad
        import pandas as pd
        A = ad.read_h5ad(a.h5ad)
    except Exception as e:
        print(f"could not load matrix ({type(e).__name__}); prediction stands "
              f"on the record, measurement deferred.")
        out["measured"] = None
        _write(out, a.out)
        return 0

    cc, pc, rc = d["context_col"], d["pert_col"], d["rep_col"]
    if rc is None or d["n_replicates"] < 2 or d["n_contexts"] < 2:
        print("MEASUREMENT: the design admits no cross-replicate pair inside a "
              "context,\n  so the interaction is not estimable. Assertion 4 is "
              "tested by whether any\n  reproducible interaction can be "
              "demonstrated here by ANY method; this\n  script reports none is "
              "available from the released annotation.")
        out["measured"] = None
        out["assertion4_applies"] = True
        _write(out, a.out)
        return 0

    obs = A.obs
    keys = [cc, pc, rc]
    X = A.X
    # pseudobulk: mean log1p-CPM per (context, perturbation, replicate)
    import scipy.sparse as sp
    lab, rows = {}, []
    for i, k in enumerate(zip(*[obs[c].astype(str) for c in keys])):
        lab.setdefault(k, []).append(i)
    groups = {k: v for k, v in lab.items() if len(v) >= MIN_CELLS}
    print(f"MEASUREMENT: {len(groups):,} conditions with >= {MIN_CELLS} cells")
    if len(groups) < 8:
        print("  too few to estimate; measurement deferred.")
        out["measured"] = None
        _write(out, a.out)
        return 0
    M, idx = [], []
    for k, v in groups.items():
        sub = X[v]
        s = np.asarray(sub.sum(0)).ravel() if sp.issparse(sub) else sub.sum(0)
        tot = s.sum()
        M.append(np.log1p(s / (tot if tot else 1.0) * 1e4))
        idx.append(k)
    M = np.asarray(M, dtype=np.float32)
    T = pd.DataFrame(idx, columns=["ctx", "pert", "rep"])
    # keep the N_FEAT most variable features, as the calculation assumes
    if M.shape[1] > N_FEAT:
        keep = np.argsort(M.var(0))[::-1][:N_FEAT]
        M = M[:, keep]

    def interaction_share(pert_labels):
        """Covariance between independent replicates of the same condition.

        Both main effects are removed out of fold: the perturbation effect from
        other contexts, the context effect from other perturbations within the
        same replicate. Noise contributes zero in expectation because the two
        members of each pair are independent.
        """
        R = M.copy()
        ctx, rep = T.ctx.to_numpy(), T.rep.to_numpy()
        pert = np.asarray(pert_labels)
        # perturbation main effect, from OTHER contexts
        for p_ in np.unique(pert):
            m = pert == p_
            for c_ in np.unique(ctx[m]):
                o = m & (ctx != c_)
                if o.sum():
                    R[m & (ctx == c_)] -= M[o].mean(0)
        # context main effect, from OTHER perturbations WITHIN a replicate
        for r_ in np.unique(rep):
            for c_ in np.unique(ctx):
                m = (rep == r_) & (ctx == c_)
                if not m.sum():
                    continue
                for p_ in np.unique(pert[m]):
                    o = m & (pert != p_)
                    if o.sum():
                        R[m & (pert == p_)] -= R[o].mean(0)
        num = den = 0.0
        for (c_, p_), g in T.groupby(["ctx", "pert"], observed=True):
            ii = g.index.to_numpy()
            if len(ii) < 2:
                continue
            for x in range(len(ii)):
                for y in range(x + 1, len(ii)):
                    num += float(R[ii[x]] @ R[ii[y]])
                    den += float(np.linalg.norm(R[ii[x]])
                                 * np.linalg.norm(R[ii[y]]))
        return num / den if den else np.nan

    share = interaction_share(T.pert.to_numpy())
    rng = np.random.default_rng(0)
    null = []
    for _ in range(N_PERM):
        q = T.pert.to_numpy().copy()
        for c_ in np.unique(T.ctx):
            m = T.ctx.to_numpy() == c_
            q[m] = rng.permutation(q[m])
        v = interaction_share(q)
        if np.isfinite(v):
            null.append(v)
    null = np.asarray(null)
    p_emp = float((1 + (null >= share).sum()) / (1 + len(null))) \
        if len(null) else float("nan")
    rejects = bool(p_emp < ALPHA)
    predicted = bool(share > floor)
    held = predicted == rejects

    print(f"  interaction share : {share:.4f}")
    print(f"  permutation null  : {null.mean():.4f} "
          f"[{np.percentile(null, 2.5):.4f}, {np.percentile(null, 97.5):.4f}]"
          if len(null) else "  permutation null  : n/a")
    print(f"  P                 : {p_emp:.4f}")
    print(f"\n  predicted resolvable (share > floor) : {predicted}")
    print(f"  permutation rejects                  : {rejects}")
    print(f"  PREDICTION HELD                      : {held}")
    out["measured"] = {"share": share, "p": p_emp, "rejects": rejects,
                       "predicted_resolvable": predicted, "held": held,
                       "n_perm": int(len(null))}
    _write(out, a.out)
    return 0


def _write(out, path):
    if path:
        Path(path).write_text(json.dumps(out, indent=2, default=str))
        print(f"\nwrote {path}")


if __name__ == "__main__":
    raise SystemExit(main())
