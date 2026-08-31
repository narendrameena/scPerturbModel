#!/usr/bin/env python3
"""Why do Tahoe-100M and LINCS disagree about which mechanisms are context-specific?

Phase II reproduced the mechanism *effect* but not the *ranking* (rho = -0.26
over six shared classes; RAF inhibitors most context-dependent in LINCS, among
the least in Tahoe). Three explanations were on the table. One is already
excluded: both assay at 24 h, so timing is not it. The other two are testable by
recomputing the Tahoe index under restrictions that make it resemble LINCS:

  measurement   LINCS measures 978 landmark genes directly; we use the
                transcriptome. 946 landmarks are present in Tahoe, so the index
                can be recomputed on exactly those genes.
  dose          LINCS is dosed mainly at 10 uM; Tahoe at 0.05/0.5/5 uM. The
                closest comparison is Tahoe's top dose alone.
  panel         only 3 Tahoe lines appear in LINCS Phase II, too few to match by
                restriction. Instead we ask whether the ranking is *stable* to
                panel choice at all, by recomputing it on random subsets of ~30
                lines. If it is not, panel differences alone can explain the
                disagreement without either dataset being wrong.

Agreement is measured per COMPOUND rather than per mechanism, since the panels
share far more compounds than mechanism classes, giving a much better-powered
comparison than the six-class rank correlation.

Outputs: results/tables/lincs_discrepancy_investigation.csv
         figure bundle results/figures/12_external/lincs_discrepancy/
"""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import build_deltas, load_pseudobulk
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
LIN = ROOT / "data" / "external" / "lincs"
FIG = ROOT / "results" / "figures" / "12_external"
TAB = ROOT / "results" / "tables"
MIN_PAIRS = 10
N_SUBSAMPLE = 6
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def cdi_table(D, K, gene_idx=None, lines=None, doses=None):
    """Per-compound context-dependence index under a restriction."""
    sub = K
    if lines is not None:
        sub = sub[sub.line.isin(lines)]
    if doses is not None:
        sub = sub[sub.conc.isin(doses)]
    if len(sub) == 0:
        return pd.DataFrame(columns=["drug", "cdi"])
    Dg = D if gene_idx is None else D[:, gene_idx]
    out = []
    for drug, gd in sub.groupby("drug", observed=True):
        if gd.line.nunique() < 5:
            continue
        cons, resid = 0.0, {}
        for _, gc in gd.groupby("conc", observed=True):
            ii = gc.i.to_numpy()
            mu = Dg[ii].mean(0)
            cons += float(np.mean(mu ** 2)) * len(ii)
            for i in ii:
                resid[i] = Dg[i] - mu
        cons /= len(gd)
        cov, n = 0.0, 0
        for _, gc in gd.groupby("line", observed=True):
            v = gc.i.to_numpy(); pl = gc.plate.to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if pl[a] != pl[b]:
                        cov += float(np.mean(resid[v[a]] * resid[v[b]])); n += 1
        if n < MIN_PAIRS:
            continue
        ctx = max(cov / n, 0.0)
        den = cons + ctx
        out.append({"drug": drug, "cdi": ctx / den if den > 0 else np.nan,
                    "effect": den, "n_lines": gd.line.nunique(), "n_pairs": n})
    return pd.DataFrame(out).dropna()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--lincs", default=str(TAB / "lincs_mechanism_ranking.csv"))
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    genes = pd.read_csv(Path(ROOT / args.pb_dir) / "genes.csv")
    gsym = genes.gene_symbol.astype(str).str.upper().to_numpy()
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {DELTA.shape[1]} genes", flush=True)

    gi = pd.read_csv(LIN / "GSE70138_Broad_LINCS_gene_info.txt.gz", sep="\t")
    lm = set(gi[gi.pr_is_lm == 1].pr_gene_symbol.astype(str).str.upper())
    lm_idx = np.where(np.isin(gsym, list(lm)))[0]
    # a size-matched random gene set, so any landmark effect is not just "fewer genes"
    rand_idx = np.sort(rng.choice(DELTA.shape[1], len(lm_idx), replace=False))
    top_dose = sorted(K.conc.unique())[-1]
    all_lines = sorted(K.line.unique())
    print(f"{len(lm_idx)} landmark genes in Tahoe; top dose {top_dose}", flush=True)

    lincs = pd.read_csv(args.lincs)
    lincs["key"] = lincs.pert_iname.map(norm)
    lincs = lincs[lincs.estimable].drop_duplicates("key")[["key", "index"]]

    variants = {
        "full (all genes, all doses, 47 lines)": dict(),
        "landmark genes only": dict(gene_idx=lm_idx),
        "random genes, size-matched": dict(gene_idx=rand_idx),
        f"top dose only ({top_dose} uM)": dict(doses=[top_dose]),
        "landmark genes + top dose": dict(gene_idx=lm_idx, doses=[top_dose]),
    }
    for s in range(N_SUBSAMPLE):
        variants[f"random 30-line panel #{s + 1}"] = dict(
            lines=list(rng.choice(all_lines, 30, replace=False)))

    recs, tabs = [], {}
    for name, kw in variants.items():
        t = cdi_table(DELTA, K, **kw)
        if not len(t):
            continue
        t["key"] = t.drug.map(norm)
        tabs[name] = t
        m = t.merge(lincs, on="key")
        rho = stats.spearmanr(m.cdi, m["index"]) if len(m) >= 8 else None
        recs.append({"variant": name, "n_compounds": len(t),
                     "median_cdi": float(t.cdi.median()),
                     "n_shared_with_lincs": len(m),
                     "rho_vs_lincs": float(rho.statistic) if rho else np.nan,
                     "p_vs_lincs": float(rho.pvalue) if rho else np.nan})
        print(f"  {name:40s} n={len(t):4d}  median={t.cdi.median():.3f}  "
              f"shared={len(m):3d}  rho_vs_LINCS="
              f"{rho.statistic:+.3f}" if rho else f"  {name}: too few shared",
              flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(TAB / "lincs_discrepancy_investigation.csv", index=False)

    # panel stability: agreement of Tahoe with ITSELF across random 30-line panels
    subs = [k for k in tabs if k.startswith("random 30-line")]
    stab = []
    for a in range(len(subs)):
        for b in range(a + 1, len(subs)):
            m = tabs[subs[a]].merge(tabs[subs[b]], on="key",
                                    suffixes=("_a", "_b"))
            if len(m) >= 20:
                stab.append(stats.spearmanr(m.cdi_a, m.cdi_b).statistic)
    base = tabs["full (all genes, all doses, 47 lines)"]
    print("\n--- interpretation ---")
    if stab:
        print(f"panel stability: Tahoe vs Tahoe across random 30-line panels, "
              f"median rho = {np.median(stab):+.3f} (n={len(stab)} pairs)")
    for name in ("landmark genes only", "random genes, size-matched"):
        if name in tabs:
            m = base.merge(tabs[name], on="key", suffixes=("_full", "_v"))
            print(f"{name:32s} vs full Tahoe: rho="
                  f"{stats.spearmanr(m.cdi_full, m.cdi_v).statistic:+.3f}")

    # ---------------- confounding effects ----------------
    # CDI is a ratio whose numerator and denominator both scale with response
    # magnitude, and the two panels differ in dose (10 uM vs 5 uM). Anything
    # that changes effect size can therefore move the index without any change
    # in context-dependence, so we check it explicitly.
    lp = pd.read_csv(args.lincs)
    lp["key"] = lp.pert_iname.map(norm)
    lp = lp[lp.estimable].drop_duplicates("key")
    m = base.merge(lp[["key", "index", "shared", "context_specific", "n_pairs"]],
                   on="key", suffixes=("_t", "_l"))
    m["effect_l"] = m.shared + m.context_specific
    print("\n--- confounding effects ---")
    conf = []
    for lab, a, b in (("Tahoe CDI vs its own effect size", "cdi", "effect"),
                      ("Tahoe CDI vs n lines", "cdi", "n_lines"),
                      ("Tahoe CDI vs n replicate pairs", "cdi", "n_pairs_t"),
                      ("LINCS index vs its own effect size", "index", "effect_l")):
        if a in m and b in m:
            r = stats.spearmanr(m[a], m[b])
            conf.append({"test": lab, "rho": r.statistic, "p": r.pvalue})
            print(f"  {lab:38s} rho={r.statistic:+.3f}  p={r.pvalue:.1e}")
    # is the disagreement explained by differing effect size?
    if len(m) >= 12:
        m["d_index"] = stats.zscore(m["index"]) - stats.zscore(m.cdi)
        m["d_effect"] = stats.zscore(np.log1p(m.effect_l)) - stats.zscore(np.log1p(m.effect))
        r = stats.spearmanr(m.d_index, m.d_effect)
        conf.append({"test": "disagreement vs effect-size difference",
                     "rho": r.statistic, "p": r.pvalue})
        print(f"  {'disagreement vs effect-size difference':38s} "
              f"rho={r.statistic:+.3f}  p={r.pvalue:.1e}  (n={len(m)})")
        # partial correlation of the two indices, controlling for effect sizes
        def resid(y, xs):
            A = np.column_stack([np.ones(len(y))] + [stats.zscore(x) for x in xs])
            return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
        ra = resid(stats.zscore(m.cdi.to_numpy()),
                   [np.log1p(m.effect.to_numpy()), np.log1p(m.effect_l.to_numpy())])
        rb = resid(stats.zscore(m["index"].to_numpy()),
                   [np.log1p(m.effect.to_numpy()), np.log1p(m.effect_l.to_numpy())])
        rp = stats.spearmanr(ra, rb)
        conf.append({"test": "Tahoe vs LINCS, effect size partialled out",
                     "rho": rp.statistic, "p": rp.pvalue})
        print(f"  {'Tahoe vs LINCS, effect partialled out':38s} "
              f"rho={rp.statistic:+.3f}  p={rp.pvalue:.1e}")
    pd.DataFrame(conf).to_csv(TAB / "lincs_discrepancy_confounds.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    r = res.dropna(subset=["rho_vs_lincs"]).sort_values("rho_vs_lincs")
    yy = np.arange(len(r))
    ax[0].barh(yy, r.rho_vs_lincs,
               color=[ORANGE if "random 30-line" in v else BLUE for v in r.variant],
               height=0.65)
    ax[0].set_yticks(yy, [f"{v[:38]} (n={n})" for v, n in
                          zip(r.variant, r.n_shared_with_lincs)], fontsize=7)
    ax[0].axvline(0, color="#888888", lw=0.9)
    ax[0].set_xlabel("per-compound agreement with LINCS (Spearman rho)")
    ax[0].set_title("A  Does any restriction reconcile Tahoe with LINCS?",
                    loc="left", fontweight="bold", fontsize=10)

    if stab:
        ax[1].hist(stab, bins=12, color=AQUA, alpha=0.85)
        ax[1].axvline(np.median(stab), color="#333333", ls="--", lw=1.2,
                      label=f"median {np.median(stab):+.2f}")
        rl = r.rho_vs_lincs.max()
        ax[1].axvline(rl, color=ORANGE, ls=":", lw=1.6,
                      label=f"best Tahoe-vs-LINCS {rl:+.2f}")
        ax[1].set_xlabel("Spearman rho, Tahoe vs Tahoe (different 30-line panels)")
        ax[1].set_ylabel("panel pairs")
        ax[1].set_title("B  How stable is the index to panel choice?",
                        loc="left", fontweight="bold", fontsize=10)
        ax[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Separating measurement, dose and panel as causes of the "
                 "Tahoe/LINCS disagreement", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "lincs_discrepancy", FIG,
                    source_data={"variants": res,
                                 "panel_stability": pd.DataFrame({"rho": stab})},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
