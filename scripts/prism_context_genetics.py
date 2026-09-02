#!/usr/bin/env python3
"""Is context-dependent drug response genotype-linked? A well-powered test in PRISM.

Every transcriptional attempt to link the context x drug interaction to genotype
failed, but all of them ran at 43-50 cell lines, where a line-level association
needs a large effect to clear multiple testing. Those negatives may simply be a
power ceiling (RESULTS.md §10), and we cannot tell which from Tahoe alone.

PRISM Repurposing changes the phenotype to buy the power: viability rather than
transcription, but ~740 cell lines instead of 47, ~1,500 compounds at 8 doses,
on replicate detection plates. The decomposition carries over unchanged -- the
response is a scalar per (line, compound, dose, replicate) instead of a gene
vector -- so we can ask, at fifteen times the context count:

  1. does the architecture hold for a viability phenotype?
  2. does mechanism structure context-dependence here too?
  3. **is the line x compound interaction associated with driver mutations?**
  4. **at how many contexts does that association become detectable?**

Question 4 is what makes 3 interpretable. Recovering BRAF -> vemurafenib at 740
lines is a positive control, not a discovery; the informative quantity is the
number of contexts at which it disappears, because that is the number Tahoe
would need. We subsample PRISM to 20, 47, 100, ... lines and measure it.

The interaction estimator is the project's standard one: covariance of per-line
residuals between INDEPENDENT REPLICATES (X1/X2/X3 detection plates) at matched
compound and dose, so uncorrelated noise cancels in expectation. Cross-dose
agreement is deliberately not used -- doses are not replicates, and the many
no-effect low doses would enter as pure noise.

Outputs: results/tables/prism_decomposition.csv
         results/tables/prism_mechanism_ranking.csv
         results/tables/prism_genotype_scan.csv
         results/tables/prism_genotype_power.csv
         figure bundle results/figures/13_prism/prism_context_genetics/
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

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_LINES = 100
MIN_MUT = 15
MIN_REP_PAIRS = 3
SUB_N = [20, 30, 47, 75, 100, 150, 250, 400, 600]
N_BOOT = 40
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.empty(n); prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def mwu_z(ranks, Mmat, n_mut, n_tot):
    """Vectorised Mann-Whitney z for many gene sets at once.

    ranks: (n_lines,) rank-transformed residuals. Mmat: (n_genes, n_lines)
    membership. Returns z per gene under the normal approximation, which is what
    scipy uses at these sample sizes anyway -- but as one matmul rather than
    n_genes separate calls, so the scan and its subsampling loop stay affordable.
    """
    n_wt = n_tot - n_mut
    Rsum = Mmat @ ranks
    mu = n_mut * (n_tot + 1) / 2.0
    sd = np.sqrt(np.maximum(n_mut * n_wt * (n_tot + 1) / 12.0, 1e-9))
    return (Rsum - mu) / sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-lines", type=int, default=MIN_LINES)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    # rows are pool_line (PR300_ACH-..., PR500_ACH-...); a line assayed in both
    # pools appears twice, and the two pools screen disjoint compound sets, so
    # averaging across pools (skipping NaN) recovers one row per cell line
    # PRISM row names are pool_line, but 8 are pool_line_FAILED_STR. The old
    # split("_")[-1] returned the literal "STR" for all eight, and
    # groupby.mean() then averaged eight different cell lines into one
    # fabricated line that carried data in 32,230 of 36,076 profiles and was
    # published in our own tables. Parse the ACH id, and DROP the STR failures
    # rather than average an unauthenticated culture into an authenticated one:
    # all eight lines also appear as clean rows, so nothing is lost. 770 rows
    # -> 737 distinct cell lines.
    keep = [not i.endswith("_FAILED_STR") for i in lfc.index]
    lfc = lfc[keep]
    lfc.index = [re.search(r"(ACH-\d+)", i).group(1) for i in lfc.index]        # -> DepMap IDs
    n_raw = len(lfc)
    lfc = lfc.groupby(level=0).mean()
    print(f"{n_raw} pool-rows -> {len(lfc)} unique cell lines", flush=True)

    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["rep"] = ti.detection_plate.astype(str).str.extract(r"_(X\d)")[0]
    ti = ti[ti.rep.notna()]
    ti["dose_s"] = ti.dose.round(4).astype(str)
    print(f"{len(lfc)} lines x {lfc.shape[1]} columns; {ti.name.nunique()} "
          f"compounds, replicates {sorted(ti.rep.unique())}", flush=True)

    L = lfc.to_numpy(dtype=np.float32)
    lines = np.array(lfc.index)
    colpos = {c: i for i, c in enumerate(lfc.columns)}

    # collapse to one column per (compound, dose, replicate)
    keys, mats = [], []
    for (cpd, dose, rep), g in ti.groupby(["name", "dose_s", "rep"],
                                          observed=True):
        idx = [colpos[c] for c in g.column_name]
        if len(idx) > 1:
            with np.errstate(invalid="ignore"):
                v = np.nanmean(L[:, idx], axis=1)
        else:
            v = L[:, idx[0]]
        keys.append((cpd, dose, rep)); mats.append(v)
    K = pd.DataFrame(keys, columns=["compound", "dose", "rep"])
    R = np.stack(mats, axis=1)                       # lines x (cpd,dose,rep)
    print(f"{R.shape[1]} (compound, dose, replicate) profiles", flush=True)

    # ---- replicate-validated decomposition ----
    recs, resid_store = [], {}
    for cpd, gc in K.groupby("compound", observed=True):
        add_num, add_den = 0.0, 0
        cov, npair = 0.0, 0
        res_by_dose = []
        for dose, gd in gc.groupby("dose", observed=True):
            cols = gd.index.to_numpy()
            if len(cols) < 2:
                continue
            sub = R[:, cols]                          # lines x replicates
            with np.errstate(invalid="ignore"):
                line_mean = np.nanmean(sub, axis=1)
            valid = np.isfinite(line_mean)
            n = int(valid.sum())
            if n < args.min_lines:
                continue
            tot = line_mean[valid].sum()
            loo = (tot - line_mean[valid]) / (n - 1)  # leave-one-line-out prior
            add_num += float(np.mean(loo ** 2)) * n; add_den += n
            res = np.full(sub.shape, np.nan, dtype=np.float32)
            res[valid] = sub[valid] - loo[:, None]
            for a in range(len(cols)):
                for b in range(a + 1, len(cols)):
                    m = np.isfinite(res[:, a]) & np.isfinite(res[:, b])
                    if m.sum() >= args.min_lines:
                        cov += float(np.mean(res[m, a] * res[m, b])); npair += 1
            with np.errstate(invalid="ignore"):
                res_by_dose.append(pd.Series(np.nanmean(res, axis=1),
                                             index=lines, name=dose))
        if npair < MIN_REP_PAIRS or add_den == 0:
            continue
        add = add_num / add_den
        inter = max(cov / npair, 0.0)
        den = add + inter
        resid_store[cpd] = pd.concat(res_by_dose, axis=1).mean(axis=1).dropna()
        recs.append({"compound": cpd, "n_lines": int(len(resid_store[cpd])),
                     "n_rep_pairs": npair, "additive": add, "interaction": inter,
                     "cdi": inter / den if den > 0 else np.nan})
    P = pd.DataFrame(recs).dropna(subset=["cdi"])
    ta, tii = P.additive.median(), P.interaction.median()
    print(f"\n{len(P)} compounds decomposed (>= {args.min_lines} lines, "
          f">= {MIN_REP_PAIRS} cross-replicate pairs)")
    print(f"  additive {ta:.4f} ({ta/(ta+tii):.0%} of reproducible), "
          f"interaction {tii:.4f} ({tii/(ta+tii):.0%})")
    print(f"  median CDI {P.cdi.median():.3f}  "
          f"[IQR {P.cdi.quantile(.25):.3f}-{P.cdi.quantile(.75):.3f}]")
    # 95% interval on the pooled split, resampling COMPOUNDS -- the shares were
    # previously reported as bare point estimates
    _b = np.array([[P.additive.sample(len(P), replace=True,
                                      random_state=s).median(),
                    P.interaction.sample(len(P), replace=True,
                                         random_state=s).median()]
                   for s in range(1000)])
    _sh = _b[:, 1] / (_b[:, 0] + _b[:, 1])
    print(f"  interaction share {tii/(ta+tii):.1%} "
          f"[{np.percentile(_sh, 2.5):.1%}-{np.percentile(_sh, 97.5):.1%}] "
          f"(bootstrap over {len(P)} compounds)")

    moa = ti.drop_duplicates("name").set_index("name").moa
    P["moa"] = P.compound.map(moa)
    P.to_csv(TAB / "prism_decomposition.csv", index=False)
    known = P[P.moa.notna() & (P.moa.astype(str) != "nan")]
    grp = (known.groupby("moa").cdi.agg(["median", "size"]).query("size >= 5")
           .sort_values("median", ascending=False))
    grp.to_csv(TAB / "prism_mechanism_ranking.csv")
    print(f"\nmechanism ranking ({len(known)} annotated, {len(grp)} classes n>=5):")
    print(grp.head(15).round(3).to_string())
    if len(grp) > 2:
        h, p = stats.kruskal(*[known[known.moa == m].cdi.to_numpy()
                               for m in grp.index])
        print(f"Kruskal-Wallis: H={h:.1f}, p={p:.2e}")

    # ---------------- genotype scan ----------------
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    dmg = mut[mut.depmap.notna() & mut.Variant_Classification.astype(str)
              .str.contains("Missense|Nonsense|Frame_Shift|Splice|In_Frame",
                            na=False)]
    gc_ = dmg.groupby("Gene_symbol").depmap.nunique()
    genes = sorted(gc_[gc_ >= 40].index)
    gsets = {g: np.array(sorted(set(d.depmap)))
             for g, d in dmg[dmg.Gene_symbol.isin(genes)].groupby("Gene_symbol")}
    print(f"\ngenotype scan: {dmg.depmap.nunique()} lines with calls, "
          f"{len(genes)} genes mutated in >=40 lines", flush=True)

    top_cpds = P.nlargest(300, "n_lines").compound.tolist()
    scan = []
    for cpd in top_cpds:
        r = resid_store[cpd]
        idx = r.index.to_numpy(); vals = r.to_numpy()
        M = np.stack([np.isin(idx, gsets[g]) for g in genes])
        nm = M.sum(1)
        keep = (nm >= MIN_MUT) & ((len(idx) - nm) >= MIN_MUT)
        if not keep.any():
            continue
        ranks = stats.rankdata(vals)
        z = mwu_z(ranks, M[keep].astype(np.float64), nm[keep], len(idx))
        pv = 2 * stats.norm.sf(np.abs(z))
        wk = np.where(keep)[0]
        med = np.array([np.median(vals[M[i]]) - np.median(vals[~M[i]])
                        for i in wk])
        scan.append(pd.DataFrame({"compound": cpd,
                                  "gene": [genes[i] for i in wk],
                                  "n_mut": nm[keep], "n_lines": len(idx),
                                  "effect": med, "z": z, "p": pv}))
    S = pd.concat(scan, ignore_index=True)
    S["q"] = bh(S.p.to_numpy())
    S = S.sort_values("p")
    # 1.5M rows: compress, and rasterise the matching scatter below -- a
    # vector point per test makes a 230 MB SVG that no viewer will open
    S.to_csv(TAB / "prism_genotype_scan.csv.gz", index=False,
             compression="gzip")
    hits = S[S.q < 0.05]
    print(f"{len(S)} (compound, gene) tests; {len(hits)} at FDR<0.05, "
          f"{(S.q<0.10).sum()} at 0.10; {hits.compound.nunique()} distinct "
          f"compounds, {hits.gene.nunique()} distinct genes")
    print("\ntop genotype x compound associations:")
    print(S.head(15)[["compound", "gene", "n_mut", "n_lines", "effect", "p", "q"]]
          .to_string(index=False))

    # ---------------- how many contexts does linkage need? ----------------
    ref = hits.drop_duplicates("compound").head(12)[["compound", "gene"]] \
              .values.tolist()
    # a scan at n contexts still pays a multiple-testing price; charge the same
    # per-test threshold the full scan used so the curve is comparable
    thr = 0.05 / len(S) * len(top_cpds)
    print(f"\npower curve on {len(ref)} confirmed associations, "
          f"{args.n_boot} subsamples/point, alpha={thr:.2e}", flush=True)
    pw = []
    for n_sub in SUB_N:
        det = []
        for cpd, gene in ref:
            r = resid_store[cpd]
            idx = r.index.to_numpy(); vals = r.to_numpy()
            memb = np.isin(idx, gsets[gene])
            if len(idx) < n_sub:
                continue
            k = 0
            for _ in range(args.n_boot):
                s = rng.choice(len(idx), n_sub, replace=False)
                m = memb[s]
                if m.sum() < 3 or (~m).sum() < 3:
                    continue
                _, p = stats.mannwhitneyu(vals[s][m], vals[s][~m],
                                          alternative="two-sided")
                k += p < thr
            det.append(k / args.n_boot)
        if det:
            pw.append({"n_contexts": n_sub, "detection_rate": float(np.mean(det)),
                       "n_assoc": len(det)})
    PW = pd.DataFrame(pw)
    PW.to_csv(TAB / "prism_genotype_power.csv", index=False)
    print(PW.round(3).to_string(index=False))

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    sh = [ta / (ta + tii), tii / (ta + tii)]
    ax[0].bar([0, 1], sh, width=0.6, color=[BLUE, ORANGE])
    for x, v in enumerate(sh):
        ax[0].text(x, v + 0.012, f"{v:.0%}", ha="center", fontsize=11)
    ax[0].set_xticks([0, 1], ["compound\n(additive)",
                              "line × compound\n(interaction)"])
    ax[0].set_ylabel("share of reproducible variance"); ax[0].set_ylim(0, 1.05)
    ax[0].set_title(f"A  PRISM viability, {P.n_lines.median():.0f} lines/compound",
                    loc="left", fontweight="bold", fontsize=10)

    ax[1].scatter(S.effect, -np.log10(S.p), s=5, alpha=0.25, color="#9e9e9e",
                  edgecolors="none", rasterized=True)
    ax[1].scatter(hits.effect, -np.log10(hits.p), s=18, color=ORANGE,
                  edgecolors="none", label=f"FDR<0.05 (n={len(hits)})")
    for r_ in hits.head(7).itertuples():
        ax[1].annotate(f"{r_.gene}·{r_.compound[:11]}",
                       (r_.effect, -np.log10(r_.p)), fontsize=6,
                       xytext=(4, 2), textcoords="offset points")
    ax[1].set_xlabel("median residual, mutant − wild-type")
    ax[1].set_ylabel("−log10 p")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title(f"B  Genotype × compound at {P.n_lines.median():.0f} lines",
                    loc="left", fontweight="bold", fontsize=10)

    if len(PW):
        ax[2].plot(PW.n_contexts, PW.detection_rate, "o-", color=VIOLET, lw=2,
                   ms=6)
        ax[2].axvline(47, color=ORANGE, ls="--", lw=1.4)
        ax[2].annotate("Tahoe-100M\n(47 lines)", (47, 0.85), fontsize=8,
                       color=ORANGE, ha="left", xytext=(6, 0),
                       textcoords="offset points")
        ax[2].axhline(0.8, color="#888888", ls=":", lw=1)
        ax[2].set_xscale("log")
        ax[2].set_xlabel("cell lines sampled")
        ax[2].set_ylabel("fraction of true associations recovered")
        ax[2].set_ylim(-0.02, 1.02)
        ax[2].set_title("C  How many contexts genotype linkage needs",
                        loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Is context-dependent drug response genotype-linked, and at "
                 "what context count is that visible?", fontsize=11, x=0.01,
                 ha="left")
    d = save_figure(fig, "prism_context_genetics", FIG,
                    source_data={"per_compound": P,
                                 "by_mechanism": grp.reset_index(),
                                 "genotype_scan": S[S.q < 0.25],
                                 "power_curve": PW},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
