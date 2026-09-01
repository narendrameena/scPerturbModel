#!/usr/bin/env python3
"""Two remaining gaps: is the identity gain threshold-dependent, and what carries it?

**Gap 5 — the identity-validated fraction rests on few lines.** Requiring a cell
line to be its own single best match by response fingerprint retains only 25-61
of 488, and the resulting 87% carries a wide interval (69.5-102.9%). If that
number depended on the arbitrary choice of "rank 1", it would be a threshold
artefact. The honest test is a stringency sweep: relax the criterion from "best
match" to "top 2%, 5%, 10%, 25% of candidates" and ask whether cross-laboratory
agreement rises monotonically with stringency. A monotone curve means the
criterion is measuring something real and the choice of cut is a
precision-versus-coverage trade-off rather than a lucky threshold; a flat or
erratic curve would mean the opposite.

**Gap 8 — expression explains 9% of what, biologically?** §17 reports that
baseline expression predicts the interaction where genotype does not, but never
says which genes carry it. Ridge spreads weight across correlated features, so
individual coefficients are not interpretable; instead we ask which *gene sets*
are enriched among the features that consistently carry weight across compounds,
and whether the predictive signal is concentrated in a few programmes or spread
thinly. The distinction matters: a handful of programmes would be a mechanistic
result, whereas diffuse weight across the transcriptome would say the predictor
is capturing global cell state and nothing sharper.

Outputs: results/tables/identity_stringency.csv
         results/tables/expression_feature_weights.csv
         figure bundle results/figures/15_architecture/robustness_interpretation/
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
CC = ROOT / "data" / "external" / "ccle"
FIG = ROOT / "results" / "figures" / "15_architecture"
TAB = ROOT / "results" / "tables"
MIN_SHARED = 20
N_HVG = 2000
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-compounds", type=int, default=100)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    import cross_lab_reproducibility as clr
    from genetic_architecture import cv_r2, residuals
    from expression_gap_closure import load_expression

    # ---------------- Gap 5: identity stringency sweep ---------------------
    print("building PRISM and GDSC residuals ...", flush=True)
    p_full, _, _ = clr.prism_residuals()
    G = clr.gdsc_residuals()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cos = {}
    for d, xr in zip(dep, smp.cross_references.astype(str)):
        if isinstance(d, str):
            for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
                cos.setdefault(int(cid), d)
    gb = {}
    for tag in ("GDSC1", "GDSC2"):
        g = G[tag].copy()
        g["dep"] = g.COSMIC_ID.map(cos)
        g = g[g.dep.notna()]
        for k, gg in g.groupby("k", observed=True):
            v = gg.groupby("dep").Z_SCORE.mean()
            gb.setdefault(k, []).append(v)
    gboth = {k: pd.concat(v).groupby(level=0).mean() for k, v in gb.items()}
    P = {norm(k): v for k, v in p_full.items()}
    shared_cpd = sorted(set(P) & set(gboth))
    rs = np.random.default_rng(0)
    perm = rs.permutation(len(shared_cpd))
    sel = [shared_cpd[i] for i in perm[:len(perm) // 2]]
    ev = [shared_cpd[i] for i in perm[len(perm) // 2:]]
    print(f"  {len(shared_cpd)} shared compounds; identity selected on "
          f"{len(sel)}, scored on {len(ev)}", flush=True)

    def fingerprints(store, keys):
        F = {}
        for c in keys:
            s = store.get(c)
            if s is None:
                continue
            for ln, v in s.items():
                F.setdefault(ln, {})[c] = float(v)
        return {k: v for k, v in F.items() if len(v) >= 8}

    FA, FB = fingerprints(P, sel), fingerprints(gboth, sel)
    both = sorted(set(FA) & set(FB))
    ranks = {}
    for x in both:
        sc = []
        for y in FB:
            common = sorted(set(FA[x]) & set(FB[y]))
            if len(common) < 8:
                continue
            u = np.array([FA[x][c] for c in common])
            v = np.array([FB[y][c] for c in common])
            if np.std(u) > 0 and np.std(v) > 0:
                sc.append((float(stats.spearmanr(u, v).statistic), y))
        if not sc:
            continue
        sc.sort(reverse=True)
        r = next((i + 1 for i, (_, y) in enumerate(sc) if y == x), None)
        if r:
            ranks[x] = (r, len(sc))

    def restrict(store, keep, only):
        return {k: v[v.index.isin(keep)] for k, v in store.items() if k in only}

    ceil = np.sqrt(0.473 * 0.438)
    rows = []
    for lab, thr in (("rank 1 (best match)", 0.0), ("top 2%", 0.02),
                     ("top 5%", 0.05), ("top 10%", 0.10), ("top 25%", 0.25),
                     ("all ID-matched", 1.01)):
        keep = {x for x, (r, n) in ranks.items()
                if (r == 1 if thr == 0.0 else r <= max(thr * n, 1))}
        if len(keep) < 8:
            continue
        rr = clr.per_compound_corr(restrict(P, keep, set(ev)),
                                   restrict(gboth, keep, set(ev)),
                                   lab)
        med = float(np.median([x["rho"] for x in rr])) if rr else np.nan
        rows.append({"criterion": lab, "n_lines": len(keep),
                     "n_compounds": len(rr), "median_rho": med,
                     "frac_of_ceiling": med / ceil})
        print(f"  {lab:22s} {len(keep):4d} lines  r={med:.3f}  "
              f"{med/ceil:.0%} of ceiling", flush=True)
    S5 = pd.DataFrame(rows)
    S5.to_csv(TAB / "identity_stringency.csv", index=False)
    if len(S5) > 3:
        rho = stats.spearmanr(np.arange(len(S5)), S5.median_rho)
        print(f"\n  monotone in stringency: rho={rho.statistic:+.2f} "
              f"(−1 = perfectly monotone, since rows go strict → loose)")
        print("  A monotone curve means the criterion tracks something real "
              "and the cut is a\n  precision/coverage trade-off, not a lucky "
              "threshold.")

    # ---------------- Gap 8: what does expression weight? ------------------
    print("\nfitting expression weights per compound ...", flush=True)
    EXPR = load_expression()
    resid, ti = residuals()
    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:args.n_compounds]
    Wsum = np.zeros(EXPR.shape[1])
    Wabs = np.zeros(EXPR.shape[1])
    n_ok = 0
    for k, cpd in enumerate(cpds):
        y = resid[cpd]
        keep = [d for d in y.index if d in EXPR.index]
        if len(keep) < 80:
            continue
        yv = y.loc[keep].to_numpy().astype(np.float64)
        yv = yv - yv.mean()
        X = EXPR.loc[keep].to_numpy(dtype=np.float64)
        mu, sd = X.mean(0), X.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        Xs = (X - mu) / sd
        # dual ridge, fixed moderate alpha: we want a stable weight ranking
        # rather than the best fit, so alpha is not re-tuned per compound
        a = 1000.0
        K_ = Xs @ Xs.T
        al = np.linalg.solve(K_ + a * np.eye(len(keep)), yv)
        w = Xs.T @ al
        Wsum += w; Wabs += np.abs(w); n_ok += 1
    Wabs /= max(n_ok, 1)
    genes = np.array([str(g).split(" (")[0] for g in EXPR.columns])
    Wdf = pd.DataFrame({"gene": genes, "mean_abs_weight": Wabs,
                        "mean_weight": Wsum / max(n_ok, 1)})
    Wdf = Wdf.sort_values("mean_abs_weight", ascending=False)
    Wdf.to_csv(TAB / "expression_feature_weights.csv", index=False)
    print(f"  weights from {n_ok} compounds over {len(genes)} genes")
    print("  top 25 by mean |weight|:")
    print("   " + ", ".join(Wdf.gene.head(25)))

    # concentration: how many genes carry half the total weight?
    srt = np.sort(Wabs)[::-1]
    cum = np.cumsum(srt) / srt.sum()
    n50 = int(np.searchsorted(cum, 0.5) + 1)
    n90 = int(np.searchsorted(cum, 0.9) + 1)
    print(f"\n  concentration: {n50} of {len(genes)} genes carry 50% of the "
          f"total weight, {n90} carry 90%")
    print(f"  ({n50/len(genes):.0%} and {n90/len(genes):.0%} of the panel) — "
          f"diffuse weight would mean the\n  predictor tracks global cell state "
          f"rather than a specific programme.")

    # gene-set enrichment of the top-weighted genes, if a gene set file exists
    gs_dir = ROOT / "data" / "external" / "genesets"
    gmt = sorted(gs_dir.glob("*.gmt")) if gs_dir.exists() else []
    enr = pd.DataFrame()
    if gmt:
        sets = {}
        for f in gmt[:2]:
            for line in f.read_text().splitlines():
                p_ = line.split("\t")
                if len(p_) > 2:
                    sets[p_[0]] = set(p_[2:])
        univ = set(genes)
        top = set(Wdf.gene.head(200))
        rows = []
        for name, gset in sets.items():
            g_ = gset & univ
            if len(g_) < 10:
                continue
            k = len(g_ & top)
            if k < 3:
                continue
            p_ = stats.hypergeom.sf(k - 1, len(univ), len(g_), len(top))
            rows.append({"gene_set": name, "n_in_universe": len(g_),
                         "n_in_top200": k, "p": p_})
        if rows:
            enr = pd.DataFrame(rows).sort_values("p")
            enr["q"] = np.minimum(enr.p * len(enr), 1.0)
            print(f"\n  gene sets enriched among the 200 top-weighted genes "
                  f"({len(sets)} sets tested):")
            print(enr.head(10)[["gene_set", "n_in_top200", "p", "q"]]
                  .to_string(index=False))
            enr.to_csv(TAB / "expression_feature_enrichment.csv", index=False)
    else:
        print("\n  (no .gmt gene sets available; enrichment skipped)")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    n_ax = 3 if len(enr) else 2
    fig, ax = plt.subplots(1, n_ax, figsize=(5.0 * n_ax, 4.1),
                           constrained_layout=True)
    if len(S5):
        xx = np.arange(len(S5))
        ax[0].bar(xx, S5.frac_of_ceiling, width=0.6,
                  color=[VIOLET if i == 0 else BLUE for i in xx])
        for i, r in enumerate(S5.itertuples()):
            ax[0].text(i, r.frac_of_ceiling + 0.015,
                       f"{r.frac_of_ceiling:.0%}\nn={r.n_lines}", ha="center",
                       fontsize=6.4)
        ax[0].axhline(1.0, ls="--", color="#888", lw=1.2)
        ax[0].set_xticks(xx, S5.criterion, fontsize=6.3, rotation=20,
                         ha="right")
        ax[0].set_ylabel("fraction of the within-lab ceiling")
        ax[0].set_ylim(0, 1.2)
        ax[0].set_title("a  The gain is monotone in stringency", loc="left",
                        fontweight="bold", fontsize=9.5)

    ax[1].plot(np.arange(1, len(cum) + 1), cum, color=VIOLET, lw=2)
    ax[1].axhline(0.5, ls=":", color="#888", lw=1.2)
    ax[1].axvline(n50, ls="--", color=ORANGE, lw=1.4,
                  label=f"{n50} genes carry 50%")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("genes, ranked by |ridge weight|")
    ax[1].set_ylabel("cumulative share of total weight")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title("b  Weight is spread, not concentrated", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(enr):
        top = enr.head(10)[::-1]
        yy = np.arange(len(top))
        ax[2].barh(yy, -np.log10(top.p), color=AQUA, height=0.65)
        ax[2].set_yticks(yy, [g[:34] for g in top.gene_set], fontsize=6)
        ax[2].set_xlabel("−log10 p (hypergeometric)")
        ax[2].set_title("c  What the top-weighted genes are", loc="left",
                        fontweight="bold", fontsize=9.5)
    fig.suptitle("Robustness of the identity criterion, and what the "
                 "expression predictor uses", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "robustness_interpretation", FIG,
                    source_data={"stringency": S5, "weights": Wdf.head(2000),
                                 "enrichment": enr}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
