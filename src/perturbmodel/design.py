"""What an atlas must measure to resolve a context x perturbation interaction.

Every large perturbation atlas is designed by intuition: more cells, more
compounds, more contexts, and replicates last if at all. Tahoe-100M spends 95.6
million cells and replicates 13.5% of its conditions; Spear-ATAC replicates
almost everything and cannot resolve an interaction at all; LINCS replicates 97%
and can. Nobody publishes the calculation that would have said in advance which
of those designs could answer the question.

This module is that calculation, in two parts.

**Detectability.** The interaction is estimated as a covariance between
independent replicates, so its precision is set by the number of replicate PAIRS
and by the per-observation noise, not by the number of cells. For an atlas with
``n_ctx`` contexts, ``n_pert`` perturbations, ``n_rep`` replicates per condition
and noise-to-signal ``snr``, the standard error of the estimated interaction
share follows from the number of independent cross-replicate pairs
``n_ctx * n_pert * n_rep * (n_rep - 1) / 2``. Doubling the cells in a condition
reduces that condition's noise; adding a replicate creates a new pair. Only the
second buys detectability, and the module quantifies the exchange rate.

**Dimensionality.** RESULTS.md sec.41 and sec.43 measured the interaction to
occupy many directions -- about 13 of 48 contexts in Tahoe, 4 of 14 in LINCS --
and, within LINCS, to GROW with the number of contexts measured (slope 0.05
directions per context, r = +0.55, P = 1e-64). It does not saturate. That has a
consequence no amount of data repairs: a model with a fixed d-dimensional context
embedding can represent at most d directions, so as an atlas adds contexts the
fraction of interaction it can express falls. ``captured_fraction`` computes that
fraction, and it is a statement about model class rather than about training.

The two combine into a design question with a numeric answer: to resolve an
interaction of a given size and dimensionality, how many contexts, perturbations
and replicates are needed, and which of those is the binding constraint?
"""
from __future__ import annotations

import numpy as np

__all__ = ["n_pairs", "interaction_se", "min_detectable_share",
           "captured_fraction", "required_replicates", "audit_design"]


def n_pairs(n_ctx, n_pert, n_rep):
    """Independent cross-replicate pairs available to the estimator.

    Cells do not appear. A condition measured in one replicate contributes no
    pair however deeply it is sequenced, which is why replication and not cell
    count is the binding constraint on this measurement.
    """
    if n_rep < 2:
        return 0
    return int(n_ctx * n_pert * n_rep * (n_rep - 1) / 2)


# Fraction of the pair-average variance carried by the first-order (per-profile)
# term. The estimator is a U-statistic of order 2: averaging over pairs of k
# profiles, its variance is a/k + b/k^2, not b'/n_pairs. The two terms give
# different scaling laws -- SE ~ k^-1/2 if the first dominates, SE ~ n_pairs^-1/2
# if the second does -- and RESULTS.md sec.49 measured the real exponent at
# -0.370 in pair units, almost exactly midway, so BOTH terms matter at realistic
# design sizes. Assuming only the second (the original formula) understates the
# standard error whenever k is small, which is exactly the regime an atlas
# builder is in.
#
# U_FIRST_ORDER is calibrated once, out of sample, on LINCS phase 1 and is fixed
# thereafter. It is a property of how correlated two pairs sharing a profile are,
# which is a feature of the estimator rather than of any dataset.
U_FIRST_ORDER = 0.1013   # LINCS phase 1, fitted on half the compounds, validated on the other half (§49)


def n_profiles(n_ctx, n_pert, n_rep):
    """Independent replicate profiles, the U-statistic's actual sample size.

    A condition measured on ``n_rep`` plates contributes ``n_rep`` profiles and
    ``n_rep(n_rep-1)/2`` pairs. The pairs are not independent of one another --
    two pairs sharing a profile are correlated -- so the pair count overstates
    the information available.
    """
    return int(n_ctx * n_pert * max(n_rep, 0))


def interaction_se(n_ctx, n_pert, n_rep, snr, n_feat=2000,
                   u_first_order=None):
    """Standard error of the estimated interaction share.

    The estimator averages products of independent replicate measurements. Each
    product has variance dominated by ``(1 + 1/snr^2)`` when the signal is small,
    and averaging ``n_feat`` weakly dependent features contributes ``1/sqrt``.

    The sample-size term is a U-statistic variance, ``a/k + b/k^2`` for ``k``
    profiles, rather than ``1/n_pairs``. With ``u_first_order = 0`` this reduces
    exactly to the original ``1/sqrt(n_pairs)`` form, so the change is opt-in and
    the old behaviour is recoverable.
    """
    p = n_pairs(n_ctx, n_pert, n_rep)
    if p <= 0:
        return np.inf
    k = n_profiles(n_ctx, n_pert, n_rep)
    u = U_FIRST_ORDER if u_first_order is None else u_first_order
    # var ∝ u/k + (1-u)/p ; u = 0 recovers the pair-only assumption
    var = (u / max(k, 1)) + ((1.0 - u) / p)
    return float((1.0 + 1.0 / max(snr, 1e-6) ** 2)
                 * np.sqrt(var / max(n_feat, 1)))


def min_detectable_share(n_ctx, n_pert, n_rep, snr, n_feat=2000, z=1.96,
                         u_first_order=None):
    """Smallest interaction share this design can distinguish from zero.

    Expressed as a fraction of reproducible variance, so it is comparable with
    the indices this literature reports.
    """
    se = interaction_se(n_ctx, n_pert, n_rep, snr, n_feat,
                        u_first_order=u_first_order)
    return float(min(z * se, 1.0))


def captured_fraction(d_model, n_ctx, slope=0.05, intercept=0.7):
    """Fraction of the interaction a fixed d-dimensional context model can hold.

    The interaction's effective dimensionality grows with the number of contexts
    measured; ``slope`` and ``intercept`` default to the LINCS fit of sec.43
    (0.05 directions per context). A model with a d-dimensional context
    embedding can represent min(d, dim) directions, and -- taking the
    interaction's energy as roughly evenly spread over its effective dimensions,
    which the near-flat eigenspectra of sec.41 support -- captures
    ``min(d, dim) / dim``.

    The consequence is the point: for fixed ``d`` this falls as ``n_ctx`` grows.
    A model that fits an atlas of 20 contexts is under-specified for one of 200,
    and no amount of training data corrects a representational limit.
    """
    dim = max(intercept + slope * n_ctx, 1.0)
    return float(min(d_model, dim) / dim)


def required_replicates(target_share, n_ctx, n_pert, snr, n_feat=2000,
                        max_rep=12, u_first_order=None):
    """Replicates per condition needed to detect an interaction of that size."""
    for r in range(2, max_rep + 1):
        if min_detectable_share(n_ctx, n_pert, r, snr, n_feat,
                                u_first_order=u_first_order) <= target_share:
            return r
    return None


def audit_design(name, n_ctx, n_pert, n_rep, snr, n_feat=2000,
                 observed_share=None):
    """One atlas through the calculation, with what it could and could not do."""
    mds = min_detectable_share(n_ctx, n_pert, n_rep, snr, n_feat)
    out = {"atlas": name, "n_contexts": n_ctx, "n_perturbations": n_pert,
           "n_replicates": n_rep, "snr": snr,
           "n_pairs": n_pairs(n_ctx, n_pert, n_rep),
           "min_detectable_share": mds,
           "observed_share": observed_share,
           "resolvable": (None if observed_share is None
                          else bool(observed_share > mds)),
           "captured_by_d10_model": captured_fraction(10, n_ctx),
           "replicates_for_1pct": required_replicates(0.01, n_ctx, n_pert, snr,
                                                      n_feat)}
    return out
