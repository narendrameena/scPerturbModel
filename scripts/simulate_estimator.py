#!/usr/bin/env python3
"""Does the recommended estimator recover a KNOWN interaction share?

Everything in this project rests on one estimator, and until now it was justified
by argument plus the fact that the alternatives disagree with it. Disagreement is
not evidence of correctness. This script generates data whose interaction share
is set by construction and asks whether each estimator returns it.

The generative model mirrors the real design:

    d[c,p,q,r] = shared[p,q] + interaction[c,p] * a(q) + batch[r] + noise

for context c, perturbation p, dose q, replicate plate r, over G genes.
`shared` is the compound's average effect, `interaction` the context-specific
part, `batch` a plate-level offset shared by every condition on that plate, and
`a(q)` controls how much of the interaction persists across doses. The quantity
being estimated is

    true share = var(interaction) / (var(shared) + var(interaction))

with noise and batch excluded from the denominator, because both are nuisance
terms the estimator is supposed to remove rather than components of the
reproducible response.

Five estimators are compared on identical data:

  recommended    leave-one-context-out prior, covariance between replicate
                 plates at MATCHED dose, minus a matched cross-context null
  variance       mean squared residual instead of a covariance
  pooled batch   same as recommended but same-plate pairs are not excluded
  in-sample      the prior includes the condition being residualised
  cross-dose     pairs taken at different doses and treated as replicates

The sweep varies the true share, the noise level, the number of replicate
plates and the number of contexts, so the report is a bias curve rather than a
single number. A good estimator is one whose curve lies on the identity line
across the whole grid; the informative failure is an estimator that is unbiased
in one regime and badly biased in another, since that is invisible in any single
dataset.

Outputs: results/tables/estimator_simulation.csv
         figure bundle results/figures/00_manuscript/estimator_simulation/
"""
import argparse
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def simulate(true_share, n_ctx=30, n_pert=20, n_dose=3, n_rep=2, n_well=2,
             n_genes=300, sigma_total=1.0, sigma_noise=1.0, sigma_batch=0.0,
             dose_persistence=1.0, sigma_ctx=0.0, seed=0):
    """Return a long table of simulated responses with known components.

    ``sigma_ctx`` is each context's GENERAL response to being perturbed at all --
    growth rate, stress tolerance, drug metabolism. It is the nuisance the
    estimator must not count as interaction: like a real interaction it is
    identical across replicate plates, so it survives every noise-cancelling
    device and lands directly in the pair covariance. It does NOT cancel against
    the control well, because the control is untreated and this term is a
    response to treatment. Setting it to zero reproduces the earlier
    simulation exactly, so the two can be compared.
    """
    rng = np.random.default_rng(seed)
    s_i = np.sqrt(true_share) * sigma_total
    s_s = np.sqrt(1 - true_share) * sigma_total
    shared = rng.normal(0, s_s, (n_pert, n_dose, n_genes))
    ctx_eff = rng.normal(0, sigma_ctx, (n_ctx, n_genes)) if sigma_ctx > 0 \
        else np.zeros((n_ctx, n_genes))
    inter = rng.normal(0, s_i, (n_ctx, n_pert, n_genes))
    # a dose-specific component so that persistence < 1 is meaningful
    inter_d = rng.normal(0, s_i, (n_ctx, n_pert, n_dose, n_genes))
    batch = rng.normal(0, sigma_batch, (n_rep, n_genes)) if sigma_batch > 0 \
        else np.zeros((n_rep, n_genes))
    a = dose_persistence
    # Each condition is measured in n_well wells ON EACH plate. Without this
    # there are no same-plate pairs at all, and "pooling same-batch
    # comparisons" would be indistinguishable from excluding them -- the
    # failure mode would be untestable rather than absent.
    # A control well per (context, plate), subtracted from every treated well on
    # that plate -- exactly what build_deltas does to the real data. Without
    # this the simulation is unfaithful in a way that breaks the estimator
    # rather than testing it: an additive plate offset that the prior averages
    # over both plates contributes -(b1-b2)^2/4 to every cross-plate pair,
    # a negative term of size sigma_batch^2/2 that has no counterpart in
    # plate-matched data. `batch_residual` is the part that does NOT cancel,
    # which is what keeps the pooled-batch failure mode testable.
    ctrl = {}
    for c, r in product(range(n_ctx), range(n_rep)):
        ctrl[(c, r)] = batch[r] + rng.normal(0, sigma_noise, n_genes)
    bres = rng.normal(0, sigma_batch * 0.5, (n_rep, n_genes)) \
        if sigma_batch > 0 else np.zeros((n_rep, n_genes))
    rows, mats = [], []
    for c, p, q, r, w in product(range(n_ctx), range(n_pert), range(n_dose),
                                 range(n_rep), range(n_well)):
        v = (shared[p, q] + ctx_eff[c]
             + a * inter[c, p] + np.sqrt(max(1 - a * a, 0)) * inter_d[c, p, q]
             + batch[r] + bres[r]
             + rng.normal(0, sigma_noise, n_genes))
        rows.append((c, p, q, r, w)); mats.append(v - ctrl[(c, r)])
    K = pd.DataFrame(rows, columns=["ctx", "pert", "dose", "rep", "well"])
    return K, np.stack(mats).astype(np.float32)


def estimate(K, D, exclude_same_batch=True, in_sample=False, cross_dose=False,
             use_variance=False, subtract_null=True, split_prior=False,
             remove_ctx=True,
             rng=None):
    """Return the estimated interaction share under one set of choices."""
    rng = rng or np.random.default_rng(0)
    if split_prior:
        return _estimate_split(K, D, rng, cross_dose, remove_ctx)
    resid, shared_num, shared_den = {}, 0.0, 0
    for (p, q), g in K.groupby(["pert", "dose"], observed=True):
        ii = g.index.to_numpy(); cx = g.ctx.to_numpy()
        tot = D[ii].sum(0)
        if in_sample:
            mu = tot / len(ii)
            for i in ii:
                resid[i] = D[i] - mu
                shared_num += float(np.mean(mu ** 2)); shared_den += 1
        else:
            csum, ccnt = {}, {}
            for c in np.unique(cx):
                m = cx == c
                csum[c] = D[ii[m]].sum(0); ccnt[c] = int(m.sum())
            for i, c in zip(ii, cx):
                n_out = len(ii) - ccnt[c]
                if n_out < 1:
                    continue
                loo = (tot - csum[c]) / n_out
                resid[i] = D[i] - loo
                shared_num += float(np.mean(loo ** 2)); shared_den += 1
    shared = shared_num / max(shared_den, 1)

    if use_variance:
        inter = float(np.mean([np.mean(r ** 2) for r in resid.values()]))
        return inter / (shared + inter) if shared + inter > 0 else np.nan

    cov, n = 0.0, 0
    grp = ["ctx", "pert"] if cross_dose else ["ctx", "pert", "dose"]
    for _, g in K.groupby(grp, observed=True):
        v = [i for i in g.index.to_numpy() if i in resid]
        rp = K.rep.loc[v].to_numpy()
        dz = K.dose.loc[v].to_numpy()
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                if exclude_same_batch and rp[a] == rp[b]:
                    continue
                if cross_dose and dz[a] == dz[b]:
                    continue
                cov += float(np.mean(resid[v[a]] * resid[v[b]])); n += 1
    if n == 0:
        return np.nan
    raw = cov / n

    off = 0.0
    if subtract_null:
        keys = np.array(sorted(resid))
        kc = K.ctx.loc[keys].to_numpy(); kp = K.pert.loc[keys].to_numpy()
        kr = K.rep.loc[keys].to_numpy()
        c2, n2 = 0.0, 0
        for _ in range(4000):
            a, b = rng.integers(0, len(keys), 2)
            if kc[a] == kc[b] or kp[a] != kp[b] or kr[a] == kr[b]:
                continue
            c2 += float(np.mean(resid[keys[a]] * resid[keys[b]])); n2 += 1
        off = (c2 / n2) if n2 else 0.0
    inter = max(raw - off, 0.0)
    return inter / (shared + inter) if shared + inter > 0 else np.nan



def _estimate_split(K, D, rng, cross_dose=False, remove_ctx=True):
    """Split-prior estimator: two disjoint leave-one-context-out priors.

    Each condition gets priors from disjoint halves of the other contexts, so
    the two members of a replicate pair share no prior-estimation noise, and
    var(shared) can be estimated as <P_A, P_B> rather than mean(P**2) -- which
    also contains noise/n_out and batch/n_rep.
    """
    PA = np.zeros_like(D); PB = np.zeros_like(D)
    usable = np.zeros(len(D), bool)
    for _, g in K.groupby(["pert", "dose"], observed=True):
        ii = g.index.to_numpy(); cx = g.ctx.to_numpy()
        uc = np.unique(cx)
        if len(uc) < 3:
            continue
        for i, ci in zip(ii, cx):
            oc = uc[uc != ci]
            perm = rng.permutation(len(oc))
            h = max(len(oc) // 2, 1)
            ca, cb = set(oc[perm[:h]]), set(oc[perm[h:]])
            ia = ii[np.array([c in ca for c in cx])]
            ib = ii[np.array([c in cb for c in cx])]
            if len(ia) and len(ib):
                PA[i] = D[ia].mean(0); PB[i] = D[ib].mean(0); usable[i] = True
    if not usable.any():
        return np.nan
    RA, RB = D - PA, D - PB
    # remove each context's general response, estimated on disjoint halves of
    # the OTHER perturbations so the two subtractions carry independent error
    if remove_ctx:
        CA = np.zeros_like(RA); CB = np.zeros_like(RB)
        ok = np.zeros(len(D), bool)
        # Grouped by (ctx, REPLICATE PLATE), not ctx alone. The control well is
        # per (context, plate), so every condition in a context on a plate
        # carries the same control noise. A correction averaged over both plates
        # carries half of each, and the cross-terms E[RA_a CB_b] then subtract
        # var(ctrl)/2 from the pair covariance -- which drove the estimate to
        # zero. Estimated within a plate, the correction cancels exactly the
        # nuisance its own residual carries, and the two members of a
        # cross-plate pair get corrections from independent plates.
        for _, g in K.groupby(["ctx", "rep"], observed=True):
            ii = g.index.to_numpy(); px = g.pert.to_numpy()
            up = np.unique(px)
            if len(up) < 3:
                continue
            for i, pi in zip(ii, px):
                op = up[up != pi]
                perm = rng.permutation(len(op))
                h = max(len(op) // 2, 1)
                pa, pb = set(op[perm[:h]]), set(op[perm[h:]])
                ja = ii[np.array([p in pa for p in px])]
                jb = ii[np.array([p in pb for p in px])]
                if len(ja) and len(jb):
                    CA[i] = RA[ja].mean(0); CB[i] = RB[jb].mean(0); ok[i] = True
        if ok.any():
            RA, RB = RA - CA, RB - CB
            usable = usable & ok
    shared = float(np.mean(PA[usable] * PB[usable]))
    grp = ["ctx", "pert"] if cross_dose else ["ctx", "pert", "dose"]
    cov, n = 0.0, 0
    for _, g in K.groupby(grp, observed=True):
        v = [i for i in g.index.to_numpy() if usable[i]]
        rp = K.rep.loc[v].to_numpy(); dz = K.dose.loc[v].to_numpy()
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                if rp[a] == rp[b]:
                    continue
                if cross_dose and dz[a] == dz[b]:
                    continue
                cov += float(np.mean(RA[v[a]] * RB[v[b]])); n += 1
    if not n:
        return np.nan
    inter = max(cov / n, 0.0)
    return inter / (shared + inter) if shared + inter > 0 else np.nan

ESTIMATORS = {
    "recommended": dict(split_prior=True),
    "old (shared prior)": dict(),
    "residual variance": dict(use_variance=True),
    "pooled batch": dict(exclude_same_batch=False),
    "in-sample prior": dict(in_sample=True),
    "cross-dose as replicate": dict(cross_dose=True),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-seeds", type=int, default=3)
    ap.add_argument("--n-genes", type=int, default=300)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    shares = [0.0, 0.1, 0.2, 0.3, 0.5]
    rows = []
    # main sweep: recovery across the range, at a realistic noise level
    for share, seed in product(shares, range(args.n_seeds)):
        K, D = simulate(share, n_genes=args.n_genes, sigma_noise=1.0,
                        sigma_batch=0.8, dose_persistence=0.7, seed=seed)
        for name, kw in ESTIMATORS.items():
            est = estimate(K, D, rng=np.random.default_rng(seed), **kw)
            rows.append({"sweep": "share", "true_share": share,
                         "estimator": name, "estimate": est, "seed": seed,
                         "noise": 1.0, "n_rep": 2, "n_ctx": 40})
        print(f"  share={share:.1f} seed={seed} done", flush=True)

    # noise sweep at a fixed true share: does any estimator track noise?
    for noise, seed in product([0.5, 1.0, 2.0, 4.0], range(args.n_seeds)):
        K, D = simulate(0.2, n_genes=args.n_genes, sigma_noise=noise,
                        sigma_batch=0.5, dose_persistence=0.7, seed=100 + seed)
        for name, kw in ESTIMATORS.items():
            est = estimate(K, D, rng=np.random.default_rng(seed), **kw)
            rows.append({"sweep": "noise", "true_share": 0.2,
                         "estimator": name, "estimate": est, "seed": seed,
                         "noise": noise, "n_rep": 2, "n_ctx": 40})
        print(f"  noise={noise} seed={seed} done", flush=True)

    # general-sensitivity sweep: the failure mode found in PRISM, where a
    # context's response to everything reproduces across disjoint compound
    # halves at r = 0.989. It is shared between replicate plates exactly as a
    # real interaction is, so an estimator that does not remove it reports it as
    # context-dependence. True share is held at 0.2 throughout: any rise with
    # sigma_ctx is the estimator counting a context property as a relation.
    for sc, seed in product([0.0, 0.3, 0.6, 1.0], range(args.n_seeds)):
        K, D = simulate(0.2, n_genes=args.n_genes, sigma_noise=1.0,
                        sigma_batch=0.5, dose_persistence=0.7, sigma_ctx=sc,
                        seed=300 + seed)
        est = estimate(K, D, rng=np.random.default_rng(seed),
                       **ESTIMATORS["recommended"])
        rows.append({"sweep": "sigma_ctx", "true_share": 0.2,
                     "estimator": "recommended", "estimate": est, "seed": seed,
                     "noise": 1.0, "n_rep": 2, "n_ctx": 40, "sigma_ctx": sc})
        print(f"  sigma_ctx={sc} seed={seed} est={est:.3f}", flush=True)

    # context-count sweep: is the estimator stable at Tahoe's n?
    for n_ctx, seed in product([10, 20, 47, 100], range(args.n_seeds)):
        K, D = simulate(0.2, n_ctx=n_ctx, n_genes=args.n_genes,
                        sigma_noise=1.0, sigma_batch=0.5,
                        dose_persistence=0.7, seed=200 + seed)
        est = estimate(K, D, rng=np.random.default_rng(seed),
                       **ESTIMATORS["recommended"])
        rows.append({"sweep": "n_ctx", "true_share": 0.2,
                     "estimator": "recommended", "estimate": est, "seed": seed,
                     "noise": 1.0, "n_rep": 2, "n_ctx": n_ctx})
        print(f"  n_ctx={n_ctx} seed={seed} done", flush=True)

    R = pd.DataFrame(rows)
    R.to_csv(TAB / "estimator_simulation.csv", index=False)

    print("\n=== recovery of a known interaction share ===")
    piv = (R[R.sweep == "share"].groupby(["estimator", "true_share"]).estimate
           .mean().unstack())
    print(piv.round(3).to_string())
    print("\nbias (estimate − truth), averaged over the sweep:")
    b = (R[R.sweep == "share"].assign(bias=lambda d: d.estimate - d.true_share)
         .groupby("estimator").bias.agg(["mean", "std"]))
    print(b.round(4).to_string())
    rec = R[(R.sweep == "share") & (R.estimator == "recommended")]
    sl = stats.linregress(rec.true_share, rec.estimate)
    print(f"\nrecommended estimator: slope {sl.slope:.3f} "
          f"(1.0 = unbiased), intercept {sl.intercept:+.4f}, "
          f"R²={sl.rvalue ** 2:.4f}")

    print("\n=== does the estimate track NOISE (it should not) ===")
    pn = (R[R.sweep == "noise"].groupby(["estimator", "noise"]).estimate
          .mean().unstack())
    print(pn.round(3).to_string())

    print("\n=== stability with context count (true share 0.2) ===")
    pc = R[R.sweep == "n_ctx"].groupby("n_ctx").estimate.agg(["mean", "std"])
    print(pc.round(4).to_string())

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    cols = {"recommended": AQUA, "residual variance": ORANGE,
            "pooled batch": VIOLET, "in-sample prior": BLUE,
            "cross-dose as replicate": "#c25fb0"}
    for name, c in cols.items():
        s = (R[(R.sweep == "share") & (R.estimator == name)]
             .groupby("true_share").estimate.mean())
        ax[0].plot(s.index, s.to_numpy(), "o-", color=c, lw=1.9, ms=5,
                   label=name)
    ax[0].plot([0, 0.5], [0, 0.5], ls="--", color="#444", lw=1.3,
               label="truth")
    ax[0].set_xlabel("true interaction share")
    ax[0].set_ylabel("estimated share")
    ax[0].legend(frameon=False, fontsize=6.4)
    ax[0].set_title("a  Recovery of a known share", loc="left",
                    fontweight="bold", fontsize=9.5)

    for name, c in cols.items():
        s = (R[(R.sweep == "noise") & (R.estimator == name)]
             .groupby("noise").estimate.mean())
        ax[1].plot(s.index, s.to_numpy(), "o-", color=c, lw=1.9, ms=5)
    ax[1].axhline(0.2, ls="--", color="#444", lw=1.3)
    ax[1].set_xlabel("noise σ (true share fixed at 0.2)")
    ax[1].set_ylabel("estimated share")
    ax[1].set_title("b  Only the covariance estimator ignores noise",
                    loc="left", fontweight="bold", fontsize=9.5)

    g = R[R.sweep == "n_ctx"].groupby("n_ctx").estimate.agg(["mean", "std"])
    ax[2].errorbar(g.index, g["mean"], yerr=g["std"], fmt="o-", color=AQUA,
                   lw=1.9, ms=6, capsize=3)
    ax[2].axhline(0.2, ls="--", color="#444", lw=1.3, label="truth")
    ax[2].axvline(47, color=ORANGE, ls=":", lw=1.4, label="Tahoe (47)")
    ax[2].set_xscale("log")
    ax[2].set_xlabel("number of contexts")
    ax[2].set_ylabel("estimated share")
    ax[2].legend(frameon=False, fontsize=7)
    ax[2].set_title("c  Stable down to Tahoe's context count", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Estimator validation on data with a known interaction share",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "estimator_simulation", FIG,
                    source_data={"simulation": R,
                                 "recovery": piv.reset_index()},
                    script=__file__)
    print(f"\nfigure bundle -> {d}")


if __name__ == "__main__":
    main()
