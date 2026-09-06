#!/usr/bin/env python3
"""Is a context-specific response a different response, or the same one harder?

RESULTS.md sec.40 ranked mechanism classes by how often their (line, drug) pairs
show reproducible interaction: protein synthesis inhibitors 100%, HDAC 81%, MEK
26%, DNA synthesis/repair 8%, against 5.2% overall. That ranking was flagged as
partly a magnitude artefact, because cross-replicate agreement rises with
signal-to-noise and the top classes are the most cytotoxic compounds in the
atlas.

Adjusting for magnitude statistically would leave the interesting question
unanswered, and there is a biological one underneath it. An interaction can take
two forms:

  POTENCY   every line runs the SAME programme and differs only in how far.
            Under a cytotoxic drug this is what differential killing looks like:
            lines die at different rates, so the shared stress-and-arrest
            response is engaged to different degrees. It is a real interaction
            statistically, and it is the general-sensitivity axis re-entering
            through the dose-response curve rather than anything drug-specific.

  REWIRING  lines run DIFFERENT programmes. This is what a pathway dependency
            predicts: a MEK inhibitor should collapse ERK output in a
            MAPK-driven line and do something else in a line that does not
            depend on it.

The two are distinguishable without any prior knowledge. Collect a drug's
interaction residuals across lines into a lines x genes matrix. If every line
responds along one direction with a different amplitude, that matrix is **rank
one**. Rewiring needs more directions. Rank is estimated from the CROSS-REPLICATE
covariance between plates 6 and 14, so noise -- which is full-rank and would
otherwise manufacture apparent rewiring -- contributes nothing.

Two things follow:

  1. a magnitude-matched enrichment test, comparing each class against control
     conditions matched on residual magnitude, so the sec.40 ranking can be
     restated with that confound removed;
  2. a potency/rewiring split per drug, and the biological check that decides
     the interpretation -- is a drug's leading direction the shared
     stress programme (sec.33), or something of its own?

Outputs: results/tables/potency_vs_rewiring.csv
         figure bundle results/figures/00_manuscript/potency_vs_rewiring/
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
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"
MIN_LINES = 12


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=400)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk,
                                                    responsive_genes)
    X, cond = load_pseudobulk(ROOT / "data/processed/pseudobulk_full")
    G, DELTA = build_deltas(X, cond, keep_plate=True, exclude_plates=())
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {len(resp)} genes", flush=True)

    # interaction residual, both main effects out, as sec.31/sec.39
    resid = {}
    for (dr, cc), g in K.groupby(["drug", "conc"], observed=True):
        ii = g.i.to_numpy(); ln = g.line.to_numpy()
        if len(np.unique(ln)) < 2:
            continue
        tot = D[ii].sum(0)
        csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
        ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
        for i, c in zip(ii, ln):
            n_out = len(ii) - ccnt[c]
            if n_out >= 1:
                resid[i] = D[i] - (tot - csum[c]) / n_out
    dv = K.drug.to_numpy()
    alpha = {}
    for (ln, pl), g in K.groupby(["line", "plate"], observed=True):
        ii = [i for i in g.i.to_numpy() if i in resid]
        by = {}
        for i in ii:
            by.setdefault(dv[i], []).append(resid[i])
        by = {d: np.mean(v, axis=0) for d, v in by.items()}
        if len(by) < 3:
            continue
        tot_a, n_a = np.sum(list(by.values()), axis=0), len(by)
        for i in ii:
            alpha[i] = (tot_a - by[dv[i]]) / (n_a - 1)
    for (dr, cc, pl), g in K.groupby(["drug", "conc", "plate"], observed=True):
        ii = [i for i in g.i.to_numpy() if i in alpha]
        if len(ii) < 2:
            continue
        m = np.mean([alpha[i] for i in ii], axis=0)
        for i in ii:
            resid[i] = resid[i] - (alpha[i] - m)

    # ---------- 1. magnitude-matched enrichment ----------
    M = pd.read_csv(TAB / "sparse_tahoe_pairs.csv")
    mag = {}
    for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"], observed=True):
        ii = [i for i in g.i.to_numpy() if i in resid]
        if ii:
            mag[(ln, dr, cc)] = float(np.mean(
                [np.linalg.norm(resid[i]) for i in ii]))
    M["mag"] = [mag.get((r.line, r.drug, r.conc), np.nan)
                for r in M.itertuples()]
    M = M[M.mag.notna()].copy()
    M["magbin"] = pd.qcut(M.mag, 10, labels=False, duplicates="drop")
    dm = pd.read_csv(TAB / "drug_metadata.csv")
    M["moa"] = M.drug.astype(str).map(dict(zip(dm.drug.astype(str),
                                               dm["moa-fine"].astype(str))))
    M["flag"] = M.q < 0.05
    print(f"\n1. MAGNITUDE-MATCHED ENRICHMENT ({len(M)} conditions, "
          f"{M.flag.sum()} flagged)", flush=True)
    print(f"   flagged rate by residual-magnitude decile: "
          f"{[round(v, 3) for v in M.groupby('magbin').flag.mean().tolist()]}")
    rng = np.random.default_rng(0)
    rows = []
    for moa, g in M.groupby("moa", observed=True):
        if moa in ("nan", "unclear", "") or len(g) < 25:
            continue
        obs = float(g.flag.mean())
        # controls drawn from the SAME magnitude deciles, in the same
        # proportions, but from other mechanism classes
        want = g.magbin.value_counts()
        null = []
        pool = M[M.moa != moa]
        for _ in range(args.n_perm):
            take = []
            for b, k in want.items():
                cand = pool.index[pool.magbin == b]
                if len(cand):
                    take.extend(rng.choice(cand, min(k, len(cand)),
                                           replace=len(cand) < k))
            if take:
                null.append(float(M.loc[take].flag.mean()))
        null = np.array(null)
        pv = float(((null >= obs).sum() + 1) / (len(null) + 1))
        rows.append({"moa": moa, "n": len(g), "rate": obs,
                     "matched_rate": float(null.mean()), "p": pv})
    E = pd.DataFrame(rows).sort_values("rate", ascending=False)
    E["q"] = np.minimum(E.p * len(E), 1.0)
    print("   class                          raw   magnitude-matched control   p")
    for r in E.head(9).itertuples():
        print(f"   {r.moa[:28]:28s} {r.rate:5.0%}   {r.matched_rate:5.0%}"
              f"   {r.p:.4f}{'  *' if r.q < 0.05 else ''}")
    print("   A class whose raw rate matches its control is enriched only "
          "because its\n   responses are large; one that exceeds it is "
          "specifically context-dependent.")

    # ---------- 2. potency vs rewiring, per drug ----------
    print(f"\n2. POTENCY vs REWIRING — rank of the between-line interaction",
          flush=True)
    rows = []
    for (dr, cc), g in K.groupby(["drug", "conc"], observed=True):
        per = {}
        for ln, gl in g.groupby("line", observed=True):
            ii = [i for i in gl.i.to_numpy() if i in resid]
            pl = [K.plate.iloc[i] for i in ii]
            if len(set(pl)) < 2:
                continue
            p0 = pl[0]
            a = [i for i, q in zip(ii, pl) if q == p0]
            b = [i for i, q in zip(ii, pl) if q != p0]
            per[ln] = (np.mean([resid[i] for i in a], axis=0),
                       np.mean([resid[i] for i in b], axis=0))
        if len(per) < MIN_LINES:
            continue
        lines = sorted(per)
        Ga = np.stack([per[l][0] for l in lines])
        Gb = np.stack([per[l][1] for l in lines])
        # cross-replicate covariance between LINES; noise is independent
        # between plates so it contributes nothing to this matrix
        C = 0.5 * (Ga @ Gb.T + Gb @ Ga.T)
        w = np.linalg.eigvalsh(C)[::-1]
        pos = w[w > 0]
        if not len(pos) or pos.sum() <= 0:
            continue
        rows.append({"drug": dr, "conc": cc, "n_lines": len(lines),
                     "reproducible": float(pos.sum()),
                     "frac_rank1": float(pos[0] / pos.sum()),
                     "n_dim": float((pos.sum() ** 2) / (pos ** 2).sum())})
    R = pd.DataFrame(rows)
    R["moa"] = R.drug.astype(str).map(dict(zip(dm.drug.astype(str),
                                               dm["moa-fine"].astype(str))))
    print(f"   {len(R)} (drug, dose) combinations with >= {MIN_LINES} "
          f"replicated lines")
    print(f"   median share of the interaction in ONE direction: "
          f"{R.frac_rank1.median():.1%}")
    print(f"   median effective number of directions: {R.n_dim.median():.2f}")
    byc = (R[R.moa.notna() & ~R.moa.isin(["nan", "unclear", ""])]
           .groupby("moa").agg(n=("frac_rank1", "size"),
                               rank1=("frac_rank1", "median"),
                               ndim=("n_dim", "median")).query("n>=4")
           .sort_values("rank1", ascending=False))
    # print the classes ONCE, ordered. A first version printed head(6) and
    # tail(6) as "most potency-like" and "most rewiring-like" without noticing
    # that only 6 classes clear the filter, so the two tables were identical.
    print(f"\n   all {len(byc)} mechanism classes with >=4 (drug, dose) "
          f"combinations, ordered\n   from most potency-like to most "
          f"rewiring-like:")
    print(byc.round(3).to_string())
    dropped = (R[R.moa.notna() & ~R.moa.isin(["nan", "unclear", ""])]
               .groupby("moa").size())
    print(f"   ({int((dropped < 4).sum())} further classes had <4 and are not "
          f"shown; protein synthesis\n   inhibitors are among them, so the "
          f"class with the highest raw flagged rate\n   cannot be placed on "
          f"this axis.)")
    R.to_csv(TAB / "potency_vs_rewiring.csv", index=False)
    E.to_csv(TAB / "magnitude_matched_enrichment.csv", index=False)

    print("\n   Interpretation. A drug whose between-line interaction is "
          "rank-one is not\n   acting differently in different lines -- it is "
          "acting the SAME way to a\n   different degree, which for a cytotoxic "
          "compound is differential killing.\n   That is the "
          "general-sensitivity axis returning through the dose-response\n   "
          "curve, and it is why a magnitude-matched control is the right test "
          "for those\n   classes.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
    top = E.head(9)
    yy = np.arange(len(top))[::-1]
    w_ = 0.38
    ax[0].barh(yy + w_ / 2, top.rate, w_, color=ORANGE, label="observed")
    ax[0].barh(yy - w_ / 2, top.matched_rate, w_, color=GREY,
               label="magnitude-matched control")
    ax[0].set_yticks(yy, [f"{m[:26]} ({n})" for m, n in zip(top.moa, top.n)],
                     fontsize=6.2)
    ax[0].set_xlabel("fraction of pairs flagged")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a  Enrichment with magnitude held fixed", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].hist(R.frac_rank1.dropna(), bins=30, color=VIOLET, alpha=0.85)
    ax[1].axvline(R.frac_rank1.median(), color=ORANGE, lw=2.2,
                  label=f"median {R.frac_rank1.median():.0%}")
    ax[1].set_xlabel("share of the interaction in a single direction")
    ax[1].set_ylabel("(drug, dose) combinations")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].text(0.03, 0.95, "1.0 = every line runs the same\nprogramme, "
               "differing only in degree", transform=ax[1].transAxes,
               va="top", fontsize=6.8, color="#444")
    ax[1].set_title("b  Potency or rewiring?", loc="left", fontweight="bold",
                    fontsize=9.5)

    if len(byc):
        sel = pd.concat([byc.head(5), byc.tail(5)])
        yy2 = np.arange(len(sel))[::-1]
        ax[2].barh(yy2, sel.rank1,
                   color=[ORANGE if v > byc.rank1.median() else AQUA
                          for v in sel.rank1], height=0.68)
        ax[2].axvline(float(byc.rank1.median()), ls="--", color="#555", lw=1.3)
        ax[2].set_yticks(yy2, [f"{m[:28]} ({int(n)})" for m, n in
                               zip(sel.index, sel.n)], fontsize=6.2)
        ax[2].set_xlabel("share in one direction")
    ax[2].set_title("c  By mechanism", loc="left", fontweight="bold",
                    fontsize=9.5)
    fig.suptitle("Separating a different response from the same response "
                 "harder", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "potency_vs_rewiring", FIG,
                    source_data={"per_drug": R, "enrichment": E,
                                 "by_moa": byc.reset_index()}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
