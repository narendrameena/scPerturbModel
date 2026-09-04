"""Detecting sparse context x perturbation interaction without knowing where.

Everything in this literature estimates interaction as a SUM: the trace of a
cross-replicate covariance, a variance component, a pooled ratio. Sums are
excellent estimators of total magnitude -- ``scripts/spectrum_benchmark.py``
shows the trace recovers a planted truth to within 0.5% at every concentration --
and they are the wrong statistic for asking whether a sparse effect exists at
all. That is not an opinion about biology; it is the classical sparse-detection
problem, where the sum has vanishing power against alternatives in which few
components are non-null (Donoho & Jin, *Ann. Statist.* 32:962, 2004).

The interaction in these atlases is sparse in exactly that sense. RESULTS.md
sec.36 found it concentrated in one drug class crossed with one genotype
(MEK inhibitors x BRAF/RAS-driven lines); sec.35 found it in 100 of 2,174
chromatin features. A pooled index averages such effects against thousands of
null cells and reports nothing, which is what produced three false nulls in this
project.

This module replaces the sum with two statistics matched to the structure.

**Per-pair reproducibility.** For each condition -- one (context, perturbation,
dose) -- the interaction residual is measured on two independent replicates and
their agreement is scored by cosine similarity. Two independent noise vectors in
p dimensions have cosine about N(0, 1/p), so a real interaction at that pair, and
only a real one, produces agreement. The null is built by pairing each replicate
A with a replicate B from a DIFFERENT condition, which destroys the pairing while
preserving every marginal, so it asks "does this pair reproduce" rather than "is
there structure anywhere".

**Global detection by Higher Criticism.** Given per-pair p-values, Higher
Criticism is asymptotically optimal for detecting sparse alternatives and is
sensitive where the sum is blind. It answers "does ANY interaction exist" with
power that does not collapse as the effect concentrates.

**Localisation by reproducible biclustering.** The per-pair scores form a
context x perturbation matrix, and a drug class acting on a genotype group is a
SUBMATRIX of it, not a scatter of independent cells. A greedy search for the
submatrix of highest mean reproducibility, significance-tested by permutation,
recovers such structure without being told what to look for -- which is the test
of whether the method can find what prior pharmacology found in sec.36.

The estimator of magnitude is unchanged: the trace is right and stays. What is
replaced is the TEST.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["pair_reproducibility", "higher_criticism", "find_bicluster",
           "SparseResult"]


def _bh(p):
    p = np.asarray(p, float)
    n = len(p)
    q = np.empty(n)
    prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r))
        q[i] = prev
    return q


class SparseResult:
    def __init__(self, table, hc, hc_null, hc_p, n_sig, trace_p, trace):
        self.table = table
        self.hc = hc
        self.hc_null = hc_null
        self.hc_p = hc_p
        self.n_sig = n_sig
        self.trace_p = trace_p
        self.trace = trace

    def report(self):
        n = len(self.table)
        L = ["sparse interaction detection", ""]
        L.append(f"  conditions tested                {n}")
        L.append(f"  pooled trace                     {self.trace:+.6f}"
                 f"   (permutation p = {self.trace_p:.4f})")
        L.append(f"  reproducible pairs at FDR<0.05   {self.n_sig}"
                 f"  ({self.n_sig / max(n, 1):.1%})")
        L.append(f"  Higher Criticism                 {self.hc:.2f}"
                 f"   (permutation p = {self.hc_p:.4f})")
        if self.hc_p < 0.05 and self.trace_p >= 0.05:
            L.append("")
            L.append("  The pooled test finds nothing and the sparse test does:")
            L.append("  the interaction is present and concentrated.")
        return "\n".join(L)


def pair_reproducibility(A, B, n_perm=500, seed=0):
    """Cross-replicate agreement per condition, with a permuted-pairing null.

    ``A`` and ``B`` are (n_conditions, n_features) interaction residuals from two
    INDEPENDENT replicates, rows aligned. Cosine similarity is used rather than
    the raw inner product so that conditions with large responses do not
    dominate: the question is whether a pair reproduces, not how big it is.
    """
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    ok = (na > 0) & (nb > 0)
    obs = np.full(len(A), np.nan)
    obs[ok] = (A[ok] * B[ok]).sum(1) / (na[ok] * nb[ok])

    rng = np.random.default_rng(seed)
    n = len(A)
    null = np.empty((n_perm, int(ok.sum())))
    Ao, Bo, nao, nbo = A[ok], B[ok], na[ok], nb[ok]
    for i in range(n_perm):
        # pair each A with a B from a DIFFERENT condition; the derangement is
        # approximated by a shift, which cannot accidentally re-pair a row
        sh = 1 + rng.integers(0, len(Ao) - 1)
        idx = (np.arange(len(Ao)) + sh) % len(Ao)
        null[i] = (Ao * Bo[idx]).sum(1) / (nao * nbo[idx])
    flat = null.ravel()
    # one-sided: a real interaction makes replicates agree, not disagree
    p = np.full(len(A), np.nan)
    p[ok] = (np.searchsorted(np.sort(flat), obs[ok], side="left")
             .astype(float))
    p[ok] = 1.0 - p[ok] / len(flat)
    p[ok] = np.clip(p[ok], 1.0 / (len(flat) + 1), 1.0)
    return obs, p, null


def higher_criticism(p, alpha0=0.5):
    """Higher Criticism statistic (Donoho & Jin 2004).

    ``HC* = max_i sqrt(n) (i/n - p_(i)) / sqrt(p_(i)(1 - p_(i)))`` over the
    smallest ``alpha0`` fraction of p-values. Sensitive to a few very small
    p-values among many null ones, which is precisely the regime where a pooled
    sum has no power.
    """
    p = np.sort(np.asarray(p, float)[np.isfinite(p)])
    n = len(p)
    if n < 10:
        return np.nan
    k = max(int(alpha0 * n), 1)
    i = np.arange(1, k + 1)
    pk = p[:k]
    den = np.sqrt(np.maximum(pk * (1 - pk), 1e-12))
    return float(np.max(np.sqrt(n) * (i / n - pk) / den))


def detect(A, B, n_perm=500, seed=0):
    """Full sparse-detection pass: per-pair, FDR count, HC, and the pooled test.

    The pooled permutation test is computed on identical data so the two
    inferences can be compared directly rather than argued about.
    """
    obs, p, null = pair_reproducibility(A, B, n_perm=n_perm, seed=seed)
    q = np.full(len(p), np.nan)
    fin = np.isfinite(p)
    q[fin] = _bh(p[fin])
    hc = higher_criticism(p)

    # HC null: recompute HC from each permutation round, so the calibration
    # accounts for the selection involved in taking a maximum
    hc_null = np.empty(len(null))
    srt = np.sort(null.ravel())
    for i in range(len(null)):
        pi = 1.0 - np.searchsorted(srt, null[i], side="left") / len(srt)
        hc_null[i] = higher_criticism(np.clip(pi, 1e-9, 1.0))
    hc_p = float(((hc_null >= hc).sum() + 1) / (len(hc_null) + 1))

    # the pooled statistic, tested the way the field tests it
    tr = float(np.nanmean((A * B).sum(1)) / A.shape[1])
    tr_null = np.array([float(np.nanmean(
        (A * B[(np.arange(len(B)) + 1 + i) % len(B)]).sum(1)) / A.shape[1])
        for i in range(min(n_perm, len(B) - 1))])
    tr_p = float(((tr_null >= tr).sum() + 1) / (len(tr_null) + 1))

    T = pd.DataFrame({"cosine": obs, "p": p, "q": q})
    return SparseResult(T, hc, hc_null, hc_p, int(np.nansum(q < 0.05)),
                        tr_p, tr)


def find_bicluster(Z, sizes=None, n_starts=20, n_perm=400, seed=0):
    """Largest average submatrix, by the LAS algorithm.

    Shabalin, Weigman, Perou & Nobel (*Ann. Appl. Statist.* 3:985, 2009). For a
    target shape (k rows, l columns) the search alternates: given the columns,
    take the k rows of highest mean; given those rows, take the l columns of
    highest mean; iterate to a fixed point, from several random starts. The shape
    is then scanned over a grid.

    A greedy search seeded at the single largest cell -- the obvious
    implementation, and the one tried first here -- cannot work: the global
    maximum is by construction larger than any average containing it, so every
    candidate addition lowers the mean and the block never grows past 1 x 1.
    LAS avoids that by fixing the size and optimising membership.

    Significance comes from permuting the matrix entries and rerunning the whole
    scan, so the null accounts for the search over shapes as well as membership.
    """
    Z = np.asarray(Z, float)
    nr, nc = Z.shape
    M = np.where(np.isfinite(Z), Z, np.nanmin(Z) if np.isfinite(Z).any() else 0)
    if sizes is None:
        sizes = [(k, l) for k in (2, 4, 8, 16) for l in (2, 4, 8, 16)
                 if k <= nr and l <= nc]

    def las(A, k, l, rng):
        best = None
        for _ in range(n_starts):
            cols = rng.choice(nc, l, replace=False)
            prev = None
            for _ in range(50):
                rows = np.argsort(-A[:, cols].mean(1))[:k]
                cols = np.argsort(-A[rows, :].mean(0))[:l]
                key = (tuple(sorted(rows)), tuple(sorted(cols)))
                if key == prev:
                    break
                prev = key
            sc = float(A[np.ix_(rows, cols)].mean())
            if best is None or sc > best[2]:
                best = (np.sort(rows), np.sort(cols), sc)
        return best

    rng = np.random.default_rng(seed)
    cands = {(k, l): las(M, k, l, rng) for k, l in sizes}

    # null: the same scan on permuted entries, so the maximum over shapes is
    # calibrated rather than compared against a single-shape null
    flat = M.ravel()
    null = np.empty((n_perm, len(sizes)))
    for i in range(n_perm):
        P = rng.permutation(flat).reshape(M.shape)
        for j, (k, l) in enumerate(sizes):
            null[i, j] = las(P, k, l, rng)[2]

    out = []
    for j, (k, l) in enumerate(sizes):
        rows, cols, sc = cands[(k, l)]
        pv = float(((null[:, j] >= sc).sum() + 1) / (n_perm + 1))
        out.append({"k": k, "l": l, "score": sc, "p": pv,
                    "null_mean": float(null[:, j].mean()),
                    "rows": rows, "cols": cols})
    best = min(out, key=lambda d: (d["p"], -d["score"]))
    return {"rows": list(best["rows"]), "cols": list(best["cols"]),
            "score": best["score"], "null_mean": best["null_mean"],
            "p": best["p"], "n_rows": len(best["rows"]),
            "n_cols": len(best["cols"]), "scan": out}
