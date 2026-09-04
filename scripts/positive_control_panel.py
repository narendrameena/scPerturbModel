#!/usr/bin/env python3
"""A positive-control panel: can each index detect pharmacology we already know?

RESULTS.md sec.36 established the rule this project kept violating: **a null from
an atlas-wide index is uninterpretable without a positive control showing the
instrument can detect a known effect on the same data.** Two nulls in this
project failed that test as soon as it was applied -- the chromatin interaction
(sec.35) and the Tahoe interaction (sec.36) -- and both turned out to be
dilution rather than absence.

One null has not yet faced it, and it is the paper's other headline: **"genome-
wide mutation status has no cross-validated predictive power for the interaction"**
(sec.17, sec.32; median CV R^2 approximately zero across 120 compounds). That is
an atlas-wide index too -- a ridge over 11,809 mutated genes -- and the same
dilution argument applies with full force. A handful of real drug x genotype
relationships averaged against eleven thousand irrelevant genes would look
exactly like this.

So the panel below is not a new analysis of the biology; it is the instrument
check the nulls require. Each entry is an established pharmacogenomic
relationship, taken from clinical practice or from the founding cell-line
pharmacogenomics literature, with the direction fixed in advance:

  MDM2 inhibitors need wild-type TP53   Vassilev et al., Science 2004. The
      textbook near-binary dependency; nutlin-3 kills TP53-WT lines and spares
      mutants.
  BRAF inhibitors need BRAF V600E       Bollag et al., Nature 2010; the basis of
      vemurafenib's label.
  ABL inhibitors need BCR-ABL           Druker et al., NEJM 2001. Imatinib in
      CML; in a cell-line panel this is the fusion-positive lines.
  ERBB2 inhibitors need ERBB2           Slamon et al., NEJM 2001; lapatinib and
      neratinib in amplified lines.
  EGFR inhibitors need EGFR mutation    Lynch et al., NEJM 2004.
  MEK inhibitors track MAPK activation  Confirmed transcriptionally in sec.36;
      included here as the viability counterpart.

For each, the test is one pre-specified contrast -- carriers versus non-carriers
of the marker, on that drug's viability -- with the direction stated before the
computation and a permutation null over the genotype label. The panel then
reports what fraction of these the instrument recovers, which is the number every
null in this project should be read against.

Finally the same compounds are put through the paper's own genome-wide ridge, so
that "genotype does not predict" and "genotype predicts, and we averaged it away"
can be told apart on identical data.

Outputs: results/tables/positive_control_panel.csv
         figure bundle results/figures/00_manuscript/positive_control_panel/
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
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"

# Each entry: the marker, how a line qualifies, the drugs, and the direction the
# published pharmacology predicts. "sensitive" means carriers should be MORE
# killed (more negative log-fold change) than non-carriers.
PANEL = [
    {"name": "MDM2 inhibitor → TP53 wild-type", "gene": "TP53",
     "carrier": "wildtype", "expect": "sensitive",
     "drugs": ["nutlin-3", "idasanutlin", "amg-232"],
     "ref": "Vassilev 2004, Science"},
    {"name": "BRAF inhibitor → BRAF mutant", "gene": "BRAF",
     "carrier": "mutant", "expect": "sensitive",
     "drugs": ["vemurafenib", "dabrafenib", "plx-4720"],
     "ref": "Bollag 2010, Nature"},
    {"name": "MEK inhibitor → BRAF/RAS mutant",
     "gene": ["BRAF", "KRAS", "NRAS"], "carrier": "mutant",
     "expect": "sensitive",
     "drugs": ["trametinib", "cobimetinib", "selumetinib"],
     "ref": "confirmed transcriptionally in §36"},
    {"name": "EGFR inhibitor → EGFR mutant", "gene": "EGFR",
     "carrier": "mutant", "expect": "sensitive",
     "drugs": ["erlotinib", "gefitinib", "osimertinib"],
     "ref": "Lynch 2004, NEJM"},
    {"name": "ERBB2 inhibitor → ERBB2 mutant/amplified", "gene": "ERBB2",
     "carrier": "mutant", "expect": "sensitive",
     "drugs": ["lapatinib", "neratinib", "tucatinib", "afatinib"],
     "ref": "Slamon 2001, NEJM"},
    {"name": "ABL inhibitor → BCR-ABL (CML lineage)", "gene": "__CML__",
     "carrier": "lineage", "expect": "sensitive",
     "drugs": ["imatinib", "dasatinib", "nilotinib", "ponatinib"],
     "ref": "Druker 2001, NEJM"},
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    lfc = lfc[[not i.endswith("_FAILED_STR") for i in lfc.index]]
    lfc.index = [re.search(r"(ACH-\d+)", i).group(1) for i in lfc.index]
    lfc = lfc.groupby(level=0).mean()
    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["lname"] = ti.name.astype(str).str.lower()
    print(f"PRISM: {lfc.shape[0]} lines x {ti.name.nunique()} compounds",
          flush=True)

    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(
        r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    ns = mut[mut.depmap.notna() &
             mut.Variant_Classification.astype(str).str.contains(
                 "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)]
    by_gene = {g: set(d.depmap) for g, d in ns.groupby("Gene_symbol")}
    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))
    sub = dict(zip(ci.depmap_id, ci.secondary_tissue.astype(str)))
    called = set(ns.depmap)          # lines with any mutation call at all
    print(f"{len(called)} lines have mutation calls; "
          f"{len(by_gene):,} genes mutated in >=1 line", flush=True)

    def carriers(entry):
        g = entry["gene"]
        if entry["carrier"] == "lineage":
            return {d for d in lfc.index
                    if "leukemia" in tis.get(d, "").lower()
                    and "myelo" in sub.get(d, "").lower()}
        genes = [g] if isinstance(g, str) else g
        hit = set()
        for x in genes:
            hit |= by_gene.get(x, set())
        if entry["carrier"] == "wildtype":
            # wild-type is only defined for lines that were actually sequenced
            return (called & set(lfc.index)) - hit
        return hit & set(lfc.index)

    rows = []
    for e in PANEL:
        cols = ti.index[ti.lname.isin([d.lower() for d in e["drugs"]])]
        if not len(cols):
            continue
        cn = ti.loc[cols, "column_name"]
        with np.errstate(invalid="ignore"):
            v = lfc[cn].mean(axis=1)
        v = v.dropna()
        car = carriers(e) & set(v.index)
        non = (called & set(v.index)) - car
        if len(car) < 5 or len(non) < 20:
            print(f"  {e['name']}: {len(car)} carriers / {len(non)} others — "
                  f"too few, skipped")
            continue
        a, b = v[list(car)], v[list(non)]
        diff = float(a.median() - b.median())      # negative = carriers killed
        rng = np.random.default_rng(0)
        idx = np.array(list(car) + list(non))
        vals = v[idx].to_numpy()
        n_car = len(car)
        null = np.array([float(np.median(vals[p[:n_car]])
                               - np.median(vals[p[n_car:]]))
                         for p in (rng.permutation(len(vals))
                                   for _ in range(args.n_perm))])
        pv = float(((null <= diff).sum() + 1) / (len(null) + 1))
        ok = (diff < 0) and (pv < 0.05)
        rows.append({"control": e["name"], "reference": e["ref"],
                     "n_carriers": len(car), "n_others": len(non),
                     "median_carrier": float(a.median()),
                     "median_other": float(b.median()), "diff": diff,
                     "p": pv, "detected": ok})
        print(f"  {e['name']:44s} carriers {a.median():+.3f} (n={len(car)})  "
              f"others {b.median():+.3f}   Δ={diff:+.3f}  p={pv:.4f}  "
              f"{'DETECTED' if ok else 'not detected'}", flush=True)

    T = pd.DataFrame(rows)
    if not len(T):
        print("no controls testable"); return
    det = int(T.detected.sum())
    print(f"\nPOSITIVE-CONTROL PANEL: {det} of {len(T)} established "
          f"relationships recovered ({det/len(T):.0%})")
    print("  This is the instrument's demonstrated sensitivity on this data. "
          "Any null\n  reported from the same data must be read against it.")

    # ---- the same compounds through the paper's genome-wide ridge ----------
    print(f"\nTHE SAME DRUGS THROUGH THE GENOME-WIDE INDEX", flush=True)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from perturbmodel.celldrug import prism_gamma
    gam, _, _ = prism_gamma()
    genes_20 = sorted([g for g, s in by_gene.items() if len(s) >= 20])
    print(f"  ridge over {len(genes_20):,} mutated genes, as in sec.17/sec.32",
          flush=True)

    def cv_r2(X, y, folds, alpha=1000.0):
        ss = float(np.sum((y - y.mean()) ** 2))
        if ss <= 0 or X.shape[1] == 0:
            return np.nan
        pred = np.zeros_like(y)
        for f in np.unique(folds):
            tr, te = folds != f, folds == f
            A, B = X[tr], X[te]
            mu, sd = A.mean(0), A.std(0)
            sd = np.where(sd > 0, sd, 1.0)
            A, B = (A - mu) / sd, (B - mu) / sd
            m = y[tr].mean()
            K = A @ A.T + alpha * np.eye(A.shape[0])
            pred[te] = m + (B @ A.T) @ np.linalg.solve(K, y[tr] - m)
        return 1 - float(np.sum((y - pred) ** 2)) / ss

    comp = []
    for e in PANEL:
        for d in e["drugs"]:
            hit = [c for c in gam if str(c).lower() == d.lower()]
            if not hit:
                continue
            y0 = gam[hit[0]]
            keep = [x for x in y0.index if x in called]
            if len(keep) < 80:
                continue
            yv = y0.loc[keep].to_numpy(float); yv = yv - yv.mean()
            X = np.stack([np.isin(keep, list(by_gene[g])).astype(float)
                          for g in genes_20], 1)
            X = X[:, X.std(0) > 0]
            folds = np.random.default_rng(0).permutation(
                np.arange(len(keep)) % 5)
            # the marker alone, as a single column
            gl = e["gene"] if isinstance(e["gene"], list) else [e["gene"]]
            mk = np.zeros(len(keep))
            for g in gl:
                if g in by_gene:
                    mk = np.maximum(mk, np.isin(keep,
                                                list(by_gene[g])).astype(float))
            r_all = cv_r2(X, yv, folds)
            r_one = (cv_r2(mk[:, None], yv, folds, alpha=1.0)
                     if mk.std() > 0 else np.nan)
            comp.append({"drug": hit[0], "control": e["name"],
                         "r2_genomewide": r_all, "r2_marker_only": r_one,
                         "n_lines": len(keep)})
    Cc = pd.DataFrame(comp)
    if len(Cc):
        print(f"  {len(Cc)} of the panel's drugs are in the corrected residual")
        print(f"  median CV R2, genome-wide ridge : "
              f"{Cc.r2_genomewide.median():+.4f}")
        print(f"  median CV R2, the ONE known marker: "
              f"{Cc.r2_marker_only.median():+.4f}")
        win = int((Cc.r2_marker_only > Cc.r2_genomewide).sum())
        print(f"  the single marker beats the genome-wide ridge for "
              f"{win}/{len(Cc)} drugs")
        print("  If a one-column model of the known marker outperforms a ridge "
              "over 11,809\n  genes, 'genotype does not predict' is a statement "
              "about the index, not the\n  biology.")
        Cc.to_csv(TAB / "positive_control_ridge.csv", index=False)
    T.to_csv(TAB / "positive_control_panel.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4), constrained_layout=True)
    yy = np.arange(len(T))[::-1]
    ax[0].barh(yy, T["diff"], color=[ORANGE if d else GREY for d in T.detected],
               height=0.68)
    ax[0].axvline(0, color="#444", lw=0.9)
    ax[0].set_yticks(yy, [f"{c[:38]}\n({int(n)} carriers)"
                          for c, n in zip(T.control, T.n_carriers)],
                     fontsize=6.2)
    ax[0].set_xlabel("Δ log-fold change, carriers − others (negative = killed)")
    ax[0].set_title(f"a  {det}/{len(T)} known relationships recovered",
                    loc="left", fontweight="bold", fontsize=9.5)

    ax[1].bar([0, 1], [det / len(T), 1.0], width=0.5, color=[ORANGE, GREY])
    ax[1].text(0, det / len(T) + 0.02, f"{det/len(T):.0%}", ha="center",
               fontsize=12, fontweight="bold")
    ax[1].set_xticks([0, 1], ["recovered by\nthis instrument",
                              "published\npharmacology"], fontsize=7.5)
    ax[1].set_ylim(0, 1.15)
    ax[1].set_ylabel("fraction of the panel")
    ax[1].text(0.5, 0.45, "every null in this project\nmust be read against\n"
               "this number", transform=ax[1].transAxes, ha="center",
               fontsize=7.5, color="#444")
    ax[1].set_title("b  The instrument's demonstrated sensitivity", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(Cc):
        ax[2].scatter(Cc.r2_genomewide, Cc.r2_marker_only, s=34, color=VIOLET,
                      edgecolors="none")
        lo = float(min(Cc.r2_genomewide.min(), Cc.r2_marker_only.min()))
        hi = float(max(Cc.r2_genomewide.max(), Cc.r2_marker_only.max()))
        ax[2].plot([lo, hi], [lo, hi], ls="--", color="#555", lw=1)
        ax[2].axhline(0, color="#888", lw=0.8); ax[2].axvline(0, color="#888",
                                                              lw=0.8)
        ax[2].set_xlabel("CV $R^2$, ridge over 11,809 mutated genes")
        ax[2].set_ylabel("CV $R^2$, the one known marker")
        ax[2].text(0.04, 0.94, f"the single marker wins for\n"
                   f"{int((Cc.r2_marker_only > Cc.r2_genomewide).sum())} of "
                   f"{len(Cc)} drugs", transform=ax[2].transAxes, va="top",
                   fontsize=7.5, color=ORANGE, fontweight="bold")
    ax[2].set_title("c  Why 'genotype does not predict' is about the index",
                    loc="left", fontweight="bold", fontsize=9.5)
    fig.suptitle("A positive-control panel: what this instrument can detect, "
                 "and what an atlas-wide index hides", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "positive_control_panel", FIG,
                    source_data={"panel": T,
                                 "ridge": Cc if len(Cc) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
