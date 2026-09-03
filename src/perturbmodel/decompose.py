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
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

__all__ = ["decompose", "Decomposition", "find_replicate_batches"]

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


def find_replicate_batches(obs, context: str, perturbation: str,
                           batch: str, dose: Optional[str] = None,
                           min_jaccard: float = 0.5):
    """Detect batches that are near-duplicates of one another.

    A "replicate plate" is a batch whose set of (context, perturbation, dose)
    conditions substantially repeats another batch's. Atlases build these
    deliberately -- Tahoe-100M's plate 14 duplicates plate 6 across 50 lines and
    95 drugs -- and then commonly withhold them from model training, after which
    downstream pipelines inherit the exclusion and lose the only same-dose
    replicates in the dataset. That is not hypothetical: it is the mistake this
    package itself shipped with.

    Returns a DataFrame of batch pairs with their Jaccard overlap and the number
    of shared conditions, most overlapping first.
    """
    keys = [context, perturbation] + ([dose] if dose else [])
    sets = {b: set(map(tuple, g[keys].drop_duplicates().to_numpy()))
            for b, g in obs.groupby(batch, observed=True)}
    rows = []
    names = sorted(sets)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = sets[names[i]], sets[names[j]]
            if not a or not b:
                continue
            inter = len(a & b)
            jac = inter / len(a | b)
            if jac >= min_jaccard:
                rows.append({"batch_a": names[i], "batch_b": names[j],
                             "shared_conditions": inter,
                             "jaccard": round(jac, 3),
                             "frac_of_a": round(inter / len(a), 3),
                             "frac_of_b": round(inter / len(b), 3)})
    return (pd.DataFrame(rows).sort_values("jaccard", ascending=False)
            if rows else pd.DataFrame(
                columns=["batch_a", "batch_b", "shared_conditions", "jaccard",
                         "frac_of_a", "frac_of_b"]))


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

    # SPLIT PRIOR. Two independent leave-one-context-out priors per condition,
    # built from disjoint halves of the other contexts. This removes two biases
    # that a single shared prior creates and that were both live in earlier
    # versions of this code:
    #
    #   numerator   both members of a same-context pair subtracted the SAME
    #               prior, so its estimation noise appeared in their covariance
    #               as signal. With disjoint priors the noise cannot covary, and
    #               the estimator returns ~0 on data with no interaction at any
    #               noise level rather than only at the one it was tested at.
    #   denominator mean(P**2) estimates var(shared) PLUS noise/n_out and
    #               batch/n_rep, so it inflated with nuisance variance and
    #               shrank the reported share by a dataset-specific factor
    #               (recovery slope 0.47-1.01 depending on the nuisances).
    #               E[<P_A, P_B>] over disjoint halves is var(shared) alone.
    #
    # Conditions in a (pert, dose) group with too few other contexts to split
    # are dropped outright: previously their prior stayed zero and the raw
    # response entered the pair and null sums un-residualised.
    PA = np.zeros_like(D)
    PB = np.zeros_like(D)
    usable = np.zeros(len(D), dtype=bool)
    for _, g in K.groupby(["pert", "dose"], observed=True):
        ii = g.index.to_numpy()
        ctxs = K.ctx[ii].to_numpy()
        uc = np.unique(ctxs)
        if len(uc) < 3:                      # need self + two disjoint halves
            continue
        for i, ci in zip(ii, ctxs):
            other_c = uc[uc != ci]
            perm = rng.permutation(len(other_c))
            h = max(len(other_c) // 2, 1)
            ca = set(other_c[perm[:h]]); cb = set(other_c[perm[h:]])
            if not ca or not cb:
                continue
            ia = ii[np.array([c in ca for c in ctxs])]
            ib = ii[np.array([c in cb for c in ctxs])]
            if not len(ia) or not len(ib):
                continue
            PA[i] = D[ia].mean(0); PB[i] = D[ib].mean(0)
            usable[i] = True
    RA = D - PA
    RB = D - PB
    # keep only conditions that have both priors
    K = K[usable].copy()
    D, RA, RB, PA, PB = (D[usable], RA[usable], RB[usable],
                         PA[usable], PB[usable])
    K.index = np.arange(len(K))
    if not len(K):
        raise ValueError("no (perturbation, dose) group has >=3 contexts; the "
                         "split-prior estimator cannot be formed")

    # ---- remove each context's GENERAL response before calling anything an
    # interaction. Subtracting the perturbation prior leaves the context's own
    # response to being perturbed at all -- growth rate, stress tolerance, drug
    # metabolism -- inside the residual. That term is shared between replicate
    # batches exactly as a genuine interaction is, so it lands in the pair
    # covariance and is counted as context-dependence. In PRISM the same
    # quantity reproduces across disjoint compound halves at r = 0.989 and
    # carries a third of what had been reported as interaction.
    #
    # It is estimated TWICE, on disjoint halves of the OTHER perturbations, with
    # half A subtracted from RA and half B from RB. Subtracting one shared
    # estimate from both would inject its estimation noise into <RA, RB> as a
    # positive term -- replacing one bias with another. Disjoint halves make
    # those two errors independent, so the covariance stays clean.
    # Grouped by (context, BATCH), not context alone. The control wells are per
    # (context, batch), so every condition in a context on a batch carries the
    # same control noise. A correction averaged over batches carries a share of
    # each, and the cross-terms then subtract var(control)/2 from the pair
    # covariance -- large enough in simulation to drive the estimate to zero.
    # Within a batch the correction cancels exactly the nuisance its own
    # residual carries, and the two members of a cross-batch pair are corrected
    # using independent batches. Simulation: with a planted general response of
    # 0, 0.5 and 1.0 SD the corrected estimator returns 0.211, 0.212, 0.211
    # against a truth of 0.20, where the uncorrected one returns 0.170, 0.347,
    # 0.595.
    CA = np.zeros_like(RA)
    CB = np.zeros_like(RB)
    ctx_ok = np.zeros(len(K), dtype=bool)
    for _, g in K.groupby(["ctx", "batch"], observed=True):
        ii = g.index.to_numpy()
        perts = K.pert[ii].to_numpy()
        up = np.unique(perts)
        if len(up) < 3:          # need self plus two disjoint halves
            continue
        for i, pi in zip(ii, perts):
            other_p = up[up != pi]
            perm = rng.permutation(len(other_p))
            h = max(len(other_p) // 2, 1)
            pa, pb = set(other_p[perm[:h]]), set(other_p[perm[h:]])
            ja = ii[np.array([p in pa for p in perts])]
            jb = ii[np.array([p in pb for p in perts])]
            if not len(ja) or not len(jb):
                continue
            CA[i] = RA[ja].mean(0); CB[i] = RB[jb].mean(0)
            ctx_ok[i] = True
    n_drop = int((~ctx_ok).sum())
    if ctx_ok.sum() < MIN_PAIRS_GLOBAL:
        warnings.warn(
            f"only {int(ctx_ok.sum())} conditions have >=3 perturbations in "
            f"their context, so the general-context response cannot be removed; "
            f"the interaction estimate will include it and is an upper bound.")
    else:
        RA = RA - CA
        RB = RB - CB
        K = K[ctx_ok].copy()
        D, RA, RB, PA, PB = (D[ctx_ok], RA[ctx_ok], RB[ctx_ok],
                             PA[ctx_ok], PB[ctx_ok])
        K.index = np.arange(len(K))
        if n_drop:
            warnings.warn(f"dropped {n_drop} conditions whose context had <3 "
                          f"perturbations, so no general-response estimate")
    R = 0.5 * (RA + RB)          # for the descriptive residual correlations

    # reproducible interaction, cross-batch pairs only
    # Reservoir sample of cross-replicate pair covariances, for the bootstrap
    # below. Keeping the FIRST n pairs was tried and is wrong: pairs arrive in
    # groupby order, so the prefix is a handful of contexts and perturbations
    # rather than a sample of them, and the resulting interval did not even
    # bracket the point estimate. Reservoir sampling keeps a uniform draw.
    kept_pairs: list[float] = []
    kept_blocks: list = []          # perturbation id, for the block bootstrap
    _seen = [0]
    RESERVOIR = 200000

    def pair_cov(same_batch: bool):
        # Grouped by (ctx, pert, DOSE): different doses are different
        # conditions, not replicates -- pairing across them is the error this
        # package exists to warn about, and the previous groupby(["ctx","pert"])
        # committed it. One residual is taken against prior A and the other
        # against prior B so their prior-estimation noise is independent.
        cov, n = 0.0, 0
        for _, g in K.groupby(["ctx", "pert", "dose"], observed=True):
            v = g.index.to_numpy()
            rp, bt = K.rep[v].to_numpy(), K.batch[v].to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if rp[a] == rp[b]:
                        continue
                    if (bt[a] == bt[b]) != same_batch:
                        continue
                    c_ = float(np.mean(RA[v[a]] * RB[v[b]]))
                    cov += c_; n += 1
                    if not same_batch:
                        _seen[0] += 1
                        pk = K.pert.iat[v[a]]
                        if len(kept_pairs) < RESERVOIR:
                            kept_pairs.append(c_); kept_blocks.append(pk)
                        else:
                            j = int(rng.integers(0, _seen[0]))
                            if j < RESERVOIR:
                                kept_pairs[j] = c_; kept_blocks[j] = pk
        return cov, n

    def null_cov():
        """Matched null: pairs from DIFFERENT contexts, same perturbation.

        Residuals are taken against a mean estimated from a finite number of
        other contexts, so the mean subtracted from one residual still contains
        the other. That leaves a negative offset of order -2 sigma^2/(n_ctx-1)
        in EVERY pair, signal or not — with 47 contexts it is not negligible,
        and it is what made the raw covariance come out negative and the
        clamped per-perturbation index collapse to exact zeros.

        With the split-prior construction above the shared-noise term is
        already gone, so this offset should now be near zero; it is retained as
        a diagnostic and a safeguard, and the report prints it so a large value
        is visible rather than silently absorbed.
        """
        cov, n = 0.0, 0
        for _, g in K.groupby("pert", observed=True):
            v = g.index.to_numpy()
            if len(v) < 2:
                continue
            cx, bt = K.ctx[v].to_numpy(), K.batch[v].to_numpy()
            take = min(len(v) * 4, 4000)
            for _ in range(take):
                a, b = rng.integers(0, len(v), 2)
                if cx[a] == cx[b] or bt[a] == bt[b]:
                    continue
                cov += float(np.mean(RA[v[a]] * RB[v[b]])); n += 1
        return (cov / n) if n else 0.0, n

    # Flag deliberate replicate batches, since losing them is the single most
    # consequential processing mistake for this estimate.
    # K carries the internal column names ("ctx", "pert", "batch"); obs carries
    # the caller's. Passing obs here made decompose() raise KeyError on every
    # real dataset, and the tests missed it because they happen to name their
    # columns "ctx"/"pert".
    rep_pairs = find_replicate_batches(K, "ctx", "pert", "batch",
                                       "dose" if dose else None)
    if len(rep_pairs):
        top = rep_pairs.iloc[0]
        warnings.append(
            f"{len(rep_pairs)} batch pair(s) look like designed replicates, "
            f"strongest {top.batch_a} / {top.batch_b} (Jaccard "
            f"{top.jaccard}, {top.shared_conditions} shared conditions). "
            f"These carry the same-dose replication this estimate depends on; "
            f"do not exclude them when measuring, even if they are held out of "
            f"model training.")

    cov, npair = pair_cov(same_batch=False)
    offset, n_null = null_cov()
    raw = (cov / npair) if npair else float("nan")
    interaction = max(raw - offset, 0.0) if npair else float("nan")
    # var(shared) from two disjoint prior halves; mean(P**2) also contained
    # noise/n_out and batch/n_rep, which is what set the recovery slope to 0.74
    additive = float(np.mean(PA * PB))
    # 95% interval on the share, by resampling the cross-replicate pairs. The
    # pooled shares were previously reported as point estimates, which hid how
    # much of the difference between datasets is resolvable.
    share_ci = (float("nan"), float("nan"))
    if len(kept_pairs) >= 100:
        # BLOCK bootstrap over perturbations. Resampling pairs i.i.d. treats
        # every pair as independent, but each residual vector appears in many
        # pairs and every pair is a mean over the same genes; a coverage
        # experiment put the implied sd 3.9x below the actual spread of the
        # estimate across independent datasets, which is why intervals came out
        # as narrow as +/-0.25 points on a 57% quantity.
        arr = np.asarray(kept_pairs)
        blocks = np.asarray(kept_blocks)
        uniq = np.unique(blocks)
        idx_by_block = {b_: np.where(blocks == b_)[0] for b_ in uniq}
        b = []
        for _ in range(300):
            pick = rng.choice(uniq, len(uniq), replace=True)
            sel = np.concatenate([idx_by_block[b_] for b_ in pick])
            b.append(arr[sel].mean())
        b = np.array(b)
        lo, hi = np.percentile(b - offset, [2.5, 97.5])
        share_ci = (float(max(lo, 0) / (additive + max(lo, 0))),
                    float(max(hi, 0) / (additive + max(hi, 0))))
        warnings.append(
            f"interaction share 95% interval (bootstrap over "
            f"{len(arr):,} sampled cross-replicate pairs): "
            f"{share_ci[0]:.1%}-{share_ci[1]:.1%}")
    if npair and n_null:
        warnings.append(
            f"interaction is reported against a matched cross-context null: "
            f"raw same-context covariance {raw:+.5f}, cross-context offset "
            f"{offset:+.5f} (n={n_null}), difference {raw - offset:+.5f}. "
            f"Reporting the raw value would have given "
            f"{max(raw, 0.0) / (additive + max(raw, 0.0)):.0%} interaction "
            f"share instead of {interaction / (additive + interaction):.0%}.")
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
        # Truncating the prefix samples groupby ORDER, not pairs: the
        # "across replicates" row was drawn from six alphabetically-first
        # contexts while "across perturbations" covered all 47, so the defining
        # contrast compared statistics from different parts of the atlas. The
        # same bug was fixed for the bootstrap and missed here.
        if not out:
            return np.empty((0, 2), int)
        out = np.array(out)
        if len(out) > max_pairs:
            out = out[rng.choice(len(out), max_pairs, replace=False)]
        return out

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

    # per-perturbation index, with power flags.
    # The shared response must be the leave-one-CONTEXT-out mean, matching the
    # pooled estimate above. Two weaker choices both bias the residual
    # covariance downward, because the covariance is taken between two
    # replicates OF THE SAME CONTEXT:
    #   in-sample mean          -> residuals of n conditions sum to zero, so
    #                              E[r_a . r_b] = -sigma^2/n with no interaction
    #   leave-one-CONDITION-out -> worse, not better: the mean subtracted from
    #                              replicate a still contains replicate b (and
    #                              vice versa), giving roughly -2 sigma^2/(n-1)
    # Dropping the whole context removes both replicates from the mean, so the
    # two residuals share no data and the bias vanishes. This matters because
    # the index is clamped at zero: a downward bias silently turns weak
    # perturbations into exact zeros rather than small positive values, and it
    # scales as 1/n, so it penalises sparsely measured perturbations most.
    precs = []
    for pert, gd in K.groupby("pert", observed=True):
        cons, resid = 0.0, {}
        for _, gc in gd.groupby("dose", observed=True):
            ii = gc.index.to_numpy()
            ctx_i = K.ctx[ii].to_numpy()
            if len(np.unique(ctx_i)) < 2:
                continue
            tot = D[ii].sum(0)
            csum, ccnt = {}, {}
            for c in np.unique(ctx_i):
                m = ctx_i == c
                csum[c] = D[ii[m]].sum(0); ccnt[c] = int(m.sum())
            for i, c in zip(ii, ctx_i):
                n_out = len(ii) - ccnt[c]
                if n_out < 1:
                    continue
                loo = (tot - csum[c]) / n_out
                cons += float(np.mean(loo ** 2))
                resid[i] = D[i] - loo
        cons /= max(len(resid), 1)
        c2, n2 = 0.0, 0
        for _, gc in gd.groupby("ctx", observed=True):
            v = gc.index.to_numpy()
            rp, bt = K.rep[v].to_numpy(), K.batch[v].to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if rp[a] == rp[b] or bt[a] == bt[b]:
                        continue
                    if v[a] not in resid or v[b] not in resid:
                        continue          # dose had a single context; skipped
                    c2 += float(np.mean(resid[v[a]] * resid[v[b]])); n2 += 1
        # matched cross-context null for THIS perturbation, same construction
        # as the pooled estimate: removes the -2 sigma^2/(n_ctx-1) offset that
        # otherwise drives weak perturbations to exactly zero after clamping
        keys = np.array([k for k in gd.index.to_numpy() if k in resid])
        c3, n3 = 0.0, 0
        if len(keys) >= 2:
            kc, kb = K.ctx[keys].to_numpy(), K.batch[keys].to_numpy()
            for _ in range(min(len(keys) * 6, 3000)):
                a, b = rng.integers(0, len(keys), 2)
                if kc[a] == kc[b] or kb[a] == kb[b]:
                    continue
                c3 += float(np.mean(resid[keys[a]] * resid[keys[b]])); n3 += 1
        off = (c3 / n3) if n3 else 0.0
        ctx_var = max(c2 / n2 - off, 0.0) if n2 else np.nan
        den = cons + (ctx_var if np.isfinite(ctx_var) else 0)
        precs.append({"perturbation": pert, "n_contexts": gd.ctx.nunique(),
                      "n_pairs": n2, "n_null_pairs": n3,
                      "raw_cov": (c2 / n2) if n2 else np.nan,
                      "null_offset": off,
                      "shared": cons, "context_specific": ctx_var,
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
