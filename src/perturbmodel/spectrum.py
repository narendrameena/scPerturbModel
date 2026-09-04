"""Reproducible-interaction spectrum: the index is a trace, so look at the rest.

Every context-dependence index in this literature, including this project's, has
the same shape. With two independent replicate measurements of the same
condition, ``gamma_A(i)`` and ``gamma_B(i)``, the interaction is estimated as

    I = (1/n) sum_i <gamma_A(i), gamma_B(i)>

which is exactly the **trace** of the cross-replicate gene-by-gene covariance

    M = (1/n) Gamma_A^T Gamma_B ,   I = tr(M).

Writing it that way makes the failure mode obvious. A trace is a sum of
eigenvalues, and it is small in two very different situations:

  * there is no reproducible interaction -- every eigenvalue is near zero;
  * the interaction is CONCENTRATED in a few directions and small as a FRACTION
    of everything measured.

**The trace is not biased by concentration.** ``scripts/spectrum_benchmark.py``
plants a constant amount of reproducible interaction and sweeps it from 2,000
directions down to 2; the trace recovers it correctly at every point. Anyone
expecting the published index to collapse under concentration -- as this module
originally claimed -- is wrong, and the simulation says so.

What the trace cannot do is distinguish those two cases, because both give a
small number, and a small VARIANCE FRACTION is not the same as no effect. That
conflation is what produced three false nulls in this project (RESULTS.md
sec.35-37). The spectrum separates them: it reports how many directions carry
reproducible interaction and which they are, so "small because absent" and
"small because concentrated" become different answers rather than the same
number.

It decomposes M and asks how much of the trace sits in directions that survive a
null in which the replicate pairing is destroyed:

    I_reproducible = sum over k of lambda_k, for lambda_k above the null edge

Noise cross-covariance has a spectrum roughly symmetric about zero, because two
independent measurements agree in no direction; real reproducible structure
produces large positive eigenvalues and no matching negative ones. The
eigenvectors are gene programmes, so the method also returns **where** the
interaction lives, which is what a trace can never do.

The estimator is unsupervised: no gene set, no pathway, no marker, no prior
pharmacology.
"""
from __future__ import annotations

import numpy as np

__all__ = ["interaction_spectrum", "SpectrumResult"]


class SpectrumResult:
    """Outcome of :func:`interaction_spectrum`."""

    def __init__(self, trace, eigvals, eigvecs, null_edge, null_max,
                 n_components, reproducible, total_positive, n_cond, n_feat):
        self.trace = trace
        self.eigvals = eigvals
        self.eigvecs = eigvecs
        self.null_edge = null_edge
        self.null_max = null_max
        self.n_components = n_components
        self.reproducible = reproducible
        self.total_positive = total_positive
        self.n_cond = n_cond
        self.n_feat = n_feat

    @property
    def concentration(self):
        """Fraction of the reproducible interaction held by the components.

        Near 1 means the interaction is concentrated in the few directions the
        null edge admits, so a small trace reflects a focused effect. Well below
        1 means it is spread across many weak directions, and the component
        count is then a floor rather than a description.
        """
        return (self.reproducible / self.trace) if self.trace > 0 else np.nan

    def report(self):
        L = ["reproducible-interaction spectrum", ""]
        L.append(f"  conditions {self.n_cond}   features {self.n_feat}")
        L.append(f"  trace (the published-style index)   {self.trace:+.6f}")
        L.append(f"  sum of positive eigenvalues         "
                 f"{self.total_positive:+.6f}")
        L.append(f"  null edge (permuted pairing)        {self.null_edge:+.6f}")
        L.append(f"  components above the null edge      {self.n_components}")
        L.append(f"  reproducible interaction            "
                 f"{self.reproducible:+.6f}")
        if self.trace > 0:
            L.append(f"  held by those components            "
                     f"{self.concentration:.0%} of the trace")
        return "\n".join(L)


def interaction_spectrum(A, B, n_perm=200, quantile=0.99, seed=0,
                         max_components=None):
    """Spectrum of the cross-replicate covariance, with a permutation null.

    Parameters
    ----------
    A, B : (n_conditions, n_features)
        Two INDEPENDENT measurements of the same conditions in the same order --
        different replicate plates, batches or samples. Independence is what
        makes the estimator noise-free: uncorrelated noise contributes zero to
        the cross-covariance in expectation, so no variance term has to be
        subtracted and no noise model is assumed.
    n_perm : int
        Permutations used to place the null edge. The permutation shuffles which
        row of ``B`` is paired with which row of ``A``, which destroys the
        replicate pairing while preserving both marginal distributions, every
        feature's variance, and all correlation structure within each matrix.
        It is therefore a null for "these two measure the same conditions", not
        for "these matrices have no structure".
    quantile : float
        Null quantile defining the edge above which an eigenvalue counts as
        reproducible.

    Returns
    -------
    SpectrumResult
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    if A.shape != B.shape:
        raise ValueError(f"A and B must match: {A.shape} vs {B.shape}")
    n, p = A.shape
    if n < 4:
        raise ValueError("need at least 4 conditions")

    def spec(Bm, k=None):
        # divided by n AND p so the trace is the per-element reproducible
        # variance -- the same scale the published-style index reports, which is
        # what makes the two directly comparable
        M = (A.T @ Bm) / (n * p)
        M = 0.5 * (M + M.T)          # the antisymmetric part carries no variance
        if k is None or k >= p:
            return np.linalg.eigvalsh(M)[::-1], M
        # only the top eigenvalues are needed for the null
        w = np.linalg.eigvalsh(M)[::-1]
        return w, M

    ev, M = spec(B)
    trace = float(np.trace(M))
    total_pos = float(ev[ev > 0].sum())

    rng = np.random.default_rng(seed)
    null_max = np.empty(n_perm)
    for i in range(n_perm):
        null_max[i] = spec(B[rng.permutation(n)])[0][0]
    edge = float(np.quantile(null_max, quantile))

    keep = ev > edge
    if max_components is not None:
        idx = np.where(keep)[0][:max_components]
        keep = np.zeros_like(keep)
        keep[idx] = True
    reproducible = float(ev[keep].sum())
    vecs = None
    if keep.any():
        w, V = np.linalg.eigh(M)
        order = np.argsort(w)[::-1]
        vecs = V[:, order[:int(keep.sum())]]
    return SpectrumResult(trace=trace, eigvals=ev, eigvecs=vecs,
                          null_edge=edge, null_max=null_max,
                          n_components=int(keep.sum()),
                          reproducible=reproducible, total_positive=total_pos,
                          n_cond=n, n_feat=p)


def condition_loadings(A, B, vecs):
    """Reproducible loading of each condition on each recovered component.

    ``<gamma_A(i), v> * <gamma_B(i), v>`` -- a product of two independent
    measurements, so a condition scores high only if it loads on the component
    in BOTH replicates. This is what identifies which (context, perturbation)
    pairs carry the interaction, without being told what to look for.
    """
    return (np.asarray(A) @ vecs) * (np.asarray(B) @ vecs)
