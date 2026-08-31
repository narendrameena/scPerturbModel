"""pertdecomp — replicate-validated decomposition of perturbation response.

Splits a perturbation atlas into the part of the response that is shared across
cellular contexts and the part that is context-specific, and — unlike a plain
variance partition — establishes that the context-specific part is *reproducible
signal* rather than noise, using the atlas's own replicate structure.

The tool encodes three lessons that are easy to get wrong:

1. **Replicates are mandatory.** Context-specific variance is estimated as the
   covariance of residuals between INDEPENDENT replicates of the same
   context x perturbation. Independent noise does not covary, so this separates
   real interaction from noise. Without a replicate axis the two are
   indistinguishable and the tool refuses to report an interaction share.

2. **Batch structure must be excluded, not just regressed.** Comparisons sharing
   a batch are dropped from every reproducibility estimate, and the tool reports
   the within- versus cross-batch inflation so the user can see what pooling
   would have cost. In Tahoe-100M that ratio is ~7x *after* batch-matched
   control normalisation.

3. **Stratified estimates need stratified power.** Per-perturbation indices are
   reported with the number of replicate pairs behind each, and flagged when
   that number is too small to support an estimate — the difference between
   "not context-dependent" and "not estimable here".

Typical use::

    from perturbmodel.decompose import decompose
    res = decompose(adata, context="cell_line", perturbation="drug",
                    replicate="plate", dose="dose", control="DMSO")
    print(res.report())
    res.per_perturbation.head()

Command line::

    python -m perturbmodel.decompose --h5ad atlas.h5ad \\
        --context cell_line --perturbation drug --replicate plate \\
        --dose dose --control DMSO
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

__all__ = ["decompose", "Decomposition"]

MIN_PAIRS_GLOBAL = 30
MIN_PAIRS_PER_PERT = 25


def _lognorm(X: np.ndarray, scale: float = 1e6) -> np.ndarray:
    tot = X.sum(1, keepdims=True)
    return np.log1p(X / np.where(tot > 0, tot, 1) * scale)


def _corr_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = A - A.mean(1, keepdims=True)
    B = B - B.mean(1, keepdims=True)
    d = np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)
    return np.where(d > 0, (A * B).sum(1) / np.where(d > 0, d, 1), np.nan)


@dataclass
class Decomposition:
    """Result of :func:`decompose`."""

    additive: float
    interaction: float
    n_contexts: int
    n_perturbations: int
    n_replicates: int
    n_pairs: int
    residual_correlations: pd.DataFrame
    batch_inflation: Optional[float]
    per_perturbation: pd.DataFrame
    n_genes: int
    warnings: list = field(default_factory=list)

    @property
    def interaction_share(self) -> float:
        tot = self.additive + self.interaction
        return float(self.interaction / tot) if tot > 0 else float("nan")

    def report(self) -> str:
        L = ["pertdecomp — replicate-validated response decomposition", ""]
        L.append(f"  contexts {self.n_contexts}   perturbations "
                 f"{self.n_perturbations}   replicates {self.n_replicates}"
                 f"   genes {self.n_genes}")
        L.append(f"  cross-replicate pairs used: {self.n_pairs}")
        L.append("")
        L.append(f"  shared / additive        {self.additive:.5f}  "
                 f"({1 - self.interaction_share:.0%} of reproducible)")
        L.append(f"  context x perturbation   {self.interaction:.5f}  "
                 f"({self.interaction_share:.0%} of reproducible)")
        if self.batch_inflation is not None:
            L.append("")
            L.append(f"  batch check: within-batch agreement is "
                     f"{self.batch_inflation:.1f}x cross-batch")
            if self.batch_inflation > 2:
                L.append("    -> batch structure is substantial; estimates here "
                         "exclude same-batch pairs, but any analysis that pools "
                         "them will overstate reproducibility by about this "
                         "factor")
        if len(self.residual_correlations):
            L.append("")
            L.append("  residual structure (median r, cross-replicate only):")
            for r in self.residual_correlations.itertuples():
                L.append(f"    {r.comparison:38s} {r.median_r:+.3f}   "
                         f"(null {r.null_r:+.3f}, n={r.n})")
            L.append("")
            L.append("    a residual that reproduces across replicates but not "
                     "across perturbations")
            L.append("    is an interaction, not a property of the context.")
        if len(self.per_perturbation):
            ok = self.per_perturbation.estimable.sum()
            L.append("")
            L.append(f"  per-perturbation index: {ok}/{len(self.per_perturbation)}"
                     f" estimable (>= {MIN_PAIRS_PER_PERT} replicate pairs)")
        for w in self.warnings:
            L.append(f"\n  WARNING: {w}")
        return "\n".join(L)

    def __repr__(self) -> str:
        return (f"<Decomposition interaction_share={self.interaction_share:.2f} "
                f"contexts={self.n_contexts} pairs={self.n_pairs}>")


def decompose(adata, context: str, perturbation: str, control,
              replicate: Optional[str] = None, dose: Optional[str] = None,
              batch: Optional[str] = None, n_genes: int = 2000,
              layer: Optional[str] = None, max_pairs: int = 4000,
              seed: int = 0) -> Decomposition:
    """Decompose perturbation response into shared and context-specific parts.

    Parameters
    ----------
    adata
        AnnData whose ``.X`` holds counts or log-normalised expression, one row
        per profile (pseudobulk recommended; single cells are aggregated by the
        caller).
    context, perturbation
        ``.obs`` columns naming the cellular context and the perturbation.
    control
        Value of ``perturbation`` marking untreated controls.
    replicate
        ``.obs`` column identifying INDEPENDENT repeats of the same
        context x perturbation — plates, donors, experimental replicates. Without
        it the interaction cannot be separated from noise and is not reported.
    dose, batch
        Optional. ``batch`` defaults to ``replicate``; comparisons sharing a
        batch are excluded from all reproducibility estimates.
    n_genes
        Number of most response-variable genes to analyse.
    """
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    warnings: list[str] = []
    obs = adata.obs.copy()
    X = adata.layers[layer] if layer else adata.X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    if X.min() >= 0 and X.max() > 50:
        X = _lognorm(X)

    for c in (context, perturbation):
        obs[c] = obs[c].astype(str)
    rep = replicate
    if rep is None:
        warnings.append(
            "no replicate column given: the interaction term cannot be "
            "separated from noise and is reported as NaN. Supply a plate, "
            "donor or experimental-replicate column.")
        obs["__rep"] = "single"; rep = "__rep"
    else:
        obs[rep] = obs[rep].astype(str)
    bat = batch or rep
    obs[bat] = obs[bat].astype(str)
    dcol = dose or "__dose"
    obs[dcol] = obs[dose].astype(str) if dose else "1"

    # replicate-matched controls
    ctrl = {k: g.index.to_numpy() for k, g in
            obs[obs[perturbation] == str(control)].groupby([context, rep],
                                                           observed=True)}
    if not ctrl:
        raise ValueError(f"no rows with {perturbation} == {control!r}")
    pos = {n: i for i, n in enumerate(obs.index)}
    rows, keys = [], []
    for name, r in obs[obs[perturbation] != str(control)].iterrows():
        c = ctrl.get((r[context], r[rep]))
        if c is None or len(c) == 0:
            continue
        rows.append(X[pos[name]] - X[[pos[x] for x in c]].mean(0))
        keys.append((r[context], r[perturbation], r[dcol], r[rep], r[bat]))
    if not rows:
        raise ValueError("no treated profiles had a replicate-matched control")
    D = np.stack(rows).astype(np.float32)
    K = pd.DataFrame(keys, columns=["ctx", "pert", "dose", "rep", "batch"])
    sel = np.sort(np.argsort(D.var(0))[-min(n_genes, D.shape[1]):])
    D = D[:, sel]

    # additive prior: leave-one-context-out mean for each (perturbation, dose)
    P = np.zeros_like(D)
    for _, g in K.groupby(["pert", "dose"], observed=True):
        ii = g.index.to_numpy()
        ctxs = K.ctx[ii].to_numpy()
        for i, ci in zip(ii, ctxs):
            other = ii[ctxs != ci]
            if len(other):
                P[i] = D[other].mean(0)
    R = D - P

    # reproducible interaction, cross-batch pairs only
    def pair_cov(same_batch: bool):
        cov, n = 0.0, 0
        for _, g in K.groupby(["ctx", "pert"], observed=True):
            v = g.index.to_numpy()
            rp, bt = K.rep[v].to_numpy(), K.batch[v].to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if rp[a] == rp[b]:
                        continue
                    if (bt[a] == bt[b]) != same_batch:
                        continue
                    cov += float(np.mean(R[v[a]] * R[v[b]])); n += 1
        return cov, n

    cov, npair = pair_cov(same_batch=False)
    interaction = max(cov / npair, 0.0) if npair else float("nan")
    additive = float(np.mean(P ** 2))
    if npair < MIN_PAIRS_GLOBAL:
        warnings.append(
            f"only {npair} cross-replicate pairs: the interaction estimate is "
            "unreliable. More replicates, not more cells, would fix this.")

    # batch inflation diagnostic
    inflation = None
    scov, sn = pair_cov(same_batch=True)
    if sn >= MIN_PAIRS_GLOBAL and npair >= MIN_PAIRS_GLOBAL and cov > 0:
        inflation = float((scov / sn) / (cov / npair))

    # residual structure, cross-replicate pairs only
    def sample(match, differ):
        out = []
        for _, g in K.groupby(match, observed=True):
            v = g.index.to_numpy()
            if len(v) < 2:
                continue
            for _ in range(min(len(v), 60)):
                i, j = rng.choice(v, 2, replace=False)
                if K.loc[i, differ] != K.loc[j, differ] and \
                        K.rep[i] != K.rep[j] and K.batch[i] != K.batch[j]:
                    out.append((i, j))
        return np.array(out[:max_pairs]) if out else np.empty((0, 2), int)

    crecs = []
    for label, (m, d) in {
            "across replicates (same ctx+pert)": (["ctx", "pert"], "rep"),
            "across perturbations (same ctx)": (["ctx"], "pert"),
            "across contexts (same pert)": (["pert"], "ctx")}.items():
        pr = sample(m, d)
        if len(pr) < MIN_PAIRS_GLOBAL:
            continue
        r = _corr_rows(R[pr[:, 0]], R[pr[:, 1]])
        rn = _corr_rows(R[rng.integers(0, len(R), len(pr))],
                        R[rng.integers(0, len(R), len(pr))])
        crecs.append({"comparison": label, "median_r": float(np.nanmedian(r)),
                      "null_r": float(np.nanmedian(rn)), "n": int(len(r))})

    # per-perturbation index, with power flags
    precs = []
    for pert, gd in K.groupby("pert", observed=True):
        cons, resid = 0.0, {}
        for _, gc in gd.groupby("dose", observed=True):
            ii = gc.index.to_numpy()
            mu = D[ii].mean(0)
            cons += float(np.mean(mu ** 2)) * len(ii)
            for i in ii:
                resid[i] = D[i] - mu
        cons /= max(len(gd), 1)
        c2, n2 = 0.0, 0
        for _, gc in gd.groupby("ctx", observed=True):
            v = gc.index.to_numpy()
            rp, bt = K.rep[v].to_numpy(), K.batch[v].to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if rp[a] == rp[b] or bt[a] == bt[b]:
                        continue
                    c2 += float(np.mean(resid[v[a]] * resid[v[b]])); n2 += 1
        ctx_var = max(c2 / n2, 0.0) if n2 else np.nan
        den = cons + (ctx_var if np.isfinite(ctx_var) else 0)
        precs.append({"perturbation": pert, "n_contexts": gd.ctx.nunique(),
                      "n_pairs": n2, "shared": cons, "context_specific": ctx_var,
                      "index": ctx_var / den if den > 0 else np.nan,
                      "estimable": n2 >= MIN_PAIRS_PER_PERT})
    per = pd.DataFrame(precs).sort_values("index", ascending=False)
    if len(per) and not per.estimable.any():
        warnings.append(
            "no perturbation has enough replicate pairs for its own index; "
            "report the pooled decomposition only.")

    return Decomposition(
        additive=additive, interaction=interaction,
        n_contexts=K.ctx.nunique(), n_perturbations=K.pert.nunique(),
        n_replicates=K.rep.nunique(), n_pairs=npair,
        residual_correlations=pd.DataFrame(crecs),
        batch_inflation=inflation, per_perturbation=per,
        n_genes=D.shape[1], warnings=warnings)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="pertdecomp",
        description="Replicate-validated decomposition of perturbation response")
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--context", required=True)
    ap.add_argument("--perturbation", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--replicate")
    ap.add_argument("--dose")
    ap.add_argument("--batch")
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--out", help="write the per-perturbation table here")
    a = ap.parse_args(argv)

    import anndata as ad
    res = decompose(ad.read_h5ad(a.h5ad), context=a.context,
                    perturbation=a.perturbation, control=a.control,
                    replicate=a.replicate, dose=a.dose, batch=a.batch,
                    n_genes=a.n_genes)
    print(res.report())
    if a.out:
        res.per_perturbation.to_csv(a.out, index=False)
        print(f"\nper-perturbation table -> {a.out}")


if __name__ == "__main__":
    main()
