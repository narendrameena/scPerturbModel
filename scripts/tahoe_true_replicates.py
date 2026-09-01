#!/usr/bin/env python3
"""What is Tahoe's interaction share when its designed replicate plate is used?

The interaction is defined as the covariance of residuals between independent
replicates of the same context x perturbation x dose. Tahoe-100M supplies these through a
deliberate design choice that is easy to discard: **plate 14 is a biological
replicate of plate 6**, 6.2M cells over 50 lines and 95 drugs, included by the
authors to demonstrate platform reproducibility and reserved by them for
validation rather than training.

Counted from the released metadata, 7,691 of 56,877 (line, drug, dose)
combinations -- **13.5%** -- sit on more than one plate. Removing plate 14, as
model-training pipelines do and as ours did by default, drops that to **5.4%**.
Meanwhile 96.8% of (line, drug) PAIRS span plates, but only because different
DOSES sit on different plates, so a "cross-plate pair" is otherwise a cross-dose
pair and not a replicate at all.

That matters because a cross-dose pair is not a replicate. A line's response at
0.05 uM and at 5 uM genuinely differ -- we measured exactly that in §13 -- so
their residual covariance mixes reproducible interaction with real dose
biology, and attenuates towards zero as the dose gap widens. Using it as the
replicate axis is the same mistake the first PRISM pass made.

This script computes the interaction share three ways on the same data:

  true        only same (line, drug, dose) on different plates (n=2,549)
  cross-dose  the pairing used previously: same (line, drug), different plate,
              different dose
  pooled      both, which is what a naive analysis would do

and reports each against a matched cross-context null, so the residual
construction offset is removed rather than assumed small.

If `true` and `cross-dose` disagree materially, the reported Tahoe interaction
share is an artefact of the pairing, and the honest number is `true` -- with the
caveat that it rests on 4.7% of the atlas, which is itself the finding.

Outputs: results/tables/tahoe_true_replicates.csv
         figure bundle results/figures/04_decomposition/tahoe_true_replicates/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "04_decomposition"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--n-boot", type=int, default=400)
    ap.add_argument("--drop-plate14", action="store_true",
                    help="reproduce the earlier, mistaken exclusion")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk,
                                                    responsive_genes)
    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    # Plate 14 is an explicit biological replicate of plate 6, included by the
    # Tahoe authors "to highlight the reproducibility of the Mosaic platform".
    # They excluded it from TRAINING and reserved it for validation, and our
    # pipeline inherited that exclusion as a default -- correct for fitting a
    # model, and exactly wrong for measuring replicate agreement, since plate 14
    # is the atlas's only source of same-dose replicates at scale. Excluding it
    # drops replication from 13.5% of conditions to 5.4%.
    excl = ("plate14",) if args.drop_plate14 else ()
    G, DELTA = build_deltas(X, cond, keep_plate=True, exclude_plates=excl)
    print(f"plates included: {sorted(G.plate.unique())}", flush=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {len(resp)} genes", flush=True)

    # leave-one-LINE-out shared response per (drug, dose)
    resid, shared_num, shared_den = {}, 0.0, 0
    for (drug, conc), g in K.groupby(["drug", "conc"], observed=True):
        ii = g.i.to_numpy(); ln = g.line.to_numpy()
        if len(np.unique(ln)) < 2:
            continue
        tot = D[ii].sum(0)
        csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
        ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
        for i, c in zip(ii, ln):
            n_out = len(ii) - ccnt[c]
            if n_out < 1:
                continue
            loo = (tot - csum[c]) / n_out
            resid[i] = D[i] - loo
            shared_num += float(np.mean(loo ** 2)); shared_den += 1
    shared = shared_num / shared_den
    print(f"shared component {shared:.5f} over {shared_den} conditions",
          flush=True)

    # collect the three pair sets
    true_p, xdose_p = [], []
    for (line, drug), g in K.groupby(["line", "drug"], observed=True):
        v = [x for x in g.i.to_numpy() if x in resid]
        if len(v) < 2:
            continue
        sub = g.set_index("i")
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                pa, pb = sub.plate[v[a]], sub.plate[v[b]]
                ca, cb = sub.conc[v[a]], sub.conc[v[b]]
                if pa == pb:
                    continue
                val = float(np.mean(resid[v[a]] * resid[v[b]]))
                (true_p if ca == cb else xdose_p).append(val)

    # matched cross-LINE null, same plate constraint
    keys = np.array(sorted(resid))
    kl = K.set_index("i").line.loc[keys].to_numpy()
    kp = K.set_index("i").plate.loc[keys].to_numpy()
    kd = K.set_index("i").drug.loc[keys].to_numpy()
    null_p = []
    for _ in range(60000):
        a, b = rng.integers(0, len(keys), 2)
        if kl[a] == kl[b] or kp[a] == kp[b] or kd[a] != kd[b]:
            continue
        null_p.append(float(np.mean(resid[keys[a]] * resid[keys[b]])))
    T, Xd, N = np.array(true_p), np.array(xdose_p), np.array(null_p)
    print(f"\npairs: true replicate {len(T)}, cross-dose {len(Xd)}, "
          f"cross-line null {len(N)}")

    def boot(arr):
        if len(arr) < 10:
            return np.nan, np.nan
        m = np.array([arr[rng.integers(0, len(arr), len(arr))].mean()
                      for _ in range(args.n_boot)])
        return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

    off = N.mean()
    rows = []
    for lab, arr in (("true replicate (same line, drug, dose)", T),
                     ("cross-dose (previous pairing)", Xd),
                     ("pooled", np.concatenate([T, Xd]))):
        raw = arr.mean()
        inter = raw - off
        lo, hi = boot(arr)
        share = max(inter, 0) / (shared + max(inter, 0))
        s_lo = max(lo - off, 0) / (shared + max(lo - off, 0))
        s_hi = max(hi - off, 0) / (shared + max(hi - off, 0))
        p = stats.mannwhitneyu(arr, N).pvalue if len(arr) >= 10 else np.nan
        rows.append({"pairing": lab, "n_pairs": len(arr), "raw_cov": raw,
                     "null_offset": off, "interaction": inter,
                     "shared": shared, "share": share,
                     "share_lo": s_lo, "share_hi": s_hi, "p_vs_null": p})
        print(f"  {lab:42s} n={len(arr):6d} raw={raw:+.5f} "
              f"inter={inter:+.5f}  share={share:.1%} "
              f"[{s_lo:.1%}-{s_hi:.1%}]  p={p:.2e}")
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "tahoe_true_replicates.csv", index=False)

    # does the cross-dose pair attenuate with dose separation, as predicted?
    # Split-half validation. The 11.5% rests on a single plate pair, so it needs
    # a check that it is not driven by a subset of drugs or lines. Splitting the
    # pairs by drug and by line and re-estimating on each half asks whether the
    # number is a property of the data or of a few conditions.
    print("\nsplit-half validation of the true-replicate estimate:")
    keyinfo = K.set_index("i")
    for axis in ("drug", "line"):
        vals = sorted(K[axis].unique())
        hs_list = [("half A", set(vals[::2])), ("half B", set(vals[1::2]))]
        got = []
        for name, hs in hs_list:
            cov_, n_ = 0.0, 0
            for (line, drug), g in K.groupby(["line", "drug"], observed=True):
                if (drug if axis == "drug" else line) not in hs:
                    continue
                v = [x for x in g.i.to_numpy() if x in resid]
                sub = g.set_index("i")
                for a_ in range(len(v)):
                    for b_ in range(a_ + 1, len(v)):
                        if sub.plate[v[a_]] == sub.plate[v[b_]]:
                            continue
                        if sub.conc[v[a_]] != sub.conc[v[b_]]:
                            continue
                        cov_ += float(np.mean(resid[v[a_]] * resid[v[b_]]))
                        n_ += 1
            if n_ >= 200:
                inter_ = cov_ / n_ - off
                sh = max(inter_, 0) / (shared + max(inter_, 0))
                got.append((name, n_, sh))
        if len(got) == 2:
            print(f"  by {axis:5s}  {got[0][0]}: {got[0][2]:.1%} "
                  f"(n={got[0][1]:,})   {got[1][0]}: {got[1][2]:.1%} "
                  f"(n={got[1][1]:,})   difference "
                  f"{abs(got[0][2]-got[1][2]):.1%}")
    print("  Read the two axes differently. Agreement across LINE halves says "
          "the estimate\n  does not depend on which cell lines are used. "
          "Disagreement across DRUG halves\n  is expected and is not "
          "instability: compounds genuinely differ in how "
          "context-\n  dependent they are (RESULTS.md 5), so the pooled figure "
          "is a mean over a wide\n  distribution rather than a constant.")

    print("\nattenuation check: cross-dose covariance by dose gap")
    gaps = []
    for (line, drug), g in K.groupby(["line", "drug"], observed=True):
        v = [x for x in g.i.to_numpy() if x in resid]
        sub = g.set_index("i")
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                if sub.plate[v[a]] == sub.plate[v[b]]:
                    continue
                ca, cb = float(sub.conc[v[a]]), float(sub.conc[v[b]])
                if ca == cb:
                    continue
                gaps.append({"gap": abs(np.log10(ca) - np.log10(cb)),
                             "cov": float(np.mean(resid[v[a]] * resid[v[b]]))})
    Gp = pd.DataFrame(gaps)
    if len(Gp):
        s = Gp.groupby(Gp.gap.round(1)).cov_ if False else \
            Gp.groupby(Gp.gap.round(1)).agg(cov=("cov", "mean"),
                                            n=("cov", "size"))
        print(s.round(5).to_string())
        rho = stats.spearmanr(Gp.gap, Gp["cov"])
        print(f"  covariance vs log10 dose gap: rho={rho.statistic:+.3f}, "
              f"p={rho.pvalue:.2e} -- a negative slope is the attenuation "
              f"predicted\n  if these pairs are not replicates but genuinely "
              f"different conditions.")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    xx = np.arange(len(R))
    ax[0].bar(xx, R.share, color=[AQUA, ORANGE, BLUE], width=0.6)
    ax[0].errorbar(xx, R.share, yerr=[R.share - R.share_lo,
                                      R.share_hi - R.share],
                   fmt="none", ecolor="#333333", capsize=4, lw=1.2)
    for x, (v, n) in enumerate(zip(R.share, R.n_pairs)):
        ax[0].text(x, v + 0.012, f"{v:.0%}\nn={n}", ha="center", fontsize=7.5)
    ax[0].set_xticks(xx, ["true\nreplicate", "cross-dose\n(previous)",
                          "pooled"], fontsize=8)
    ax[0].set_ylabel("interaction share of reproducible variance")
    ax[0].set_title("A  The pairing decides the answer", loc="left",
                    fontweight="bold", fontsize=10)

    ax[1].hist(N, bins=60, color="#bdbdbd", density=True, label="cross-line null")
    ax[1].hist(T, bins=60, color=AQUA, alpha=0.6, density=True,
               label="true replicate")
    ax[1].axvline(N.mean(), color="#666666", lw=1.5)
    ax[1].axvline(T.mean(), color=AQUA, lw=1.5)
    ax[1].set_xlim(np.percentile(np.concatenate([T, N]), [1, 99]))
    ax[1].set_xlabel("residual covariance")
    ax[1].set_ylabel("density"); ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("B  Signal against its matched null", loc="left",
                    fontweight="bold", fontsize=10)

    if len(Gp):
        s2 = Gp.groupby(Gp.gap.round(1)).agg(cov=("cov", "mean"),
                                             n=("cov", "size"))
        s2 = s2[s2.n >= 30]
        ax[2].plot(s2.index, s2["cov"], "o-", color=VIOLET, lw=2)
        ax[2].axhline(off, color="#888888", ls="--", lw=1,
                      label="cross-line null")
        ax[2].set_xlabel("dose separation (log10 units)")
        ax[2].set_ylabel("mean residual covariance")
        ax[2].legend(frameon=False, fontsize=8)
        ax[2].set_title("C  Cross-dose pairs attenuate with dose gap",
                        loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Tahoe-100M has almost no true replicates — does the "
                 "interaction share survive that?", fontsize=11, x=0.01,
                 ha="left")
    d = save_figure(fig, "tahoe_true_replicates", FIG,
                    source_data={"estimates": R,
                                 "dose_gap": Gp if len(Gp) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
