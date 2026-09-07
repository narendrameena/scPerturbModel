#!/usr/bin/env python3
"""The design calculation against every atlas in scPerturb, not five.

RESULTS.md sec.46 showed that a design calculation predicts, from three integers
per atlas, which of five published perturbation atlases could resolve a context x
perturbation interaction. Five is few, and they were the five this project
happened to analyse. A reviewer is right to call that anecdotal.

scPerturb (Peidli et al., *Nature Methods* 2024) harmonises the field's
perturbation datasets into one schema with common ``obs`` fields, which makes
every one of them a test case: read the metadata, count contexts, perturbations
and replicates, and ask what the calculation says each design can detect. No
expression values are needed -- the calculation uses only the design.

Two things come out of it:

  1. **A benchmark.** The distribution of minimum detectable interaction share
     across the field, and how many datasets can resolve an effect of the size
     this project measures (0.5-70% depending on atlas). If most cannot, that is
     a statement about the field's designs rather than about any one study.
  2. **The binding constraint, per dataset.** For each, whether contexts,
     perturbations or replicates is what limits it, and what would have to change.

The prediction that keeps this honest: the calculation is applied uniformly, with
the same shared noise constant used in sec.46 and no per-dataset tuning. Datasets
this project has already measured (Tahoe, LINCS, OP3, sci-Plex, Spear-ATAC) are
scored the same way as the rest and their known outcomes are checked, not fitted.

Outputs: results/tables/atlas_design_benchmark.csv
         figure bundle results/figures/00_manuscript/atlas_benchmark/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.atlas_meta import (CTX, MISSING, PERT, REP, describe, nlev,
                                     pick, read_obs)
from perturbmodel.design import (min_detectable_share, n_pairs,
                                 required_replicates)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
SP = ROOT / "data" / "external" / "scperturb_all"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"
SNR = 0.20                      # the shared constant of sec.46, unchanged

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.05,
                    help="interaction share the field would want to detect")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    files = sorted(SP.glob("*.h5ad"))
    print(f"{len(files)} scPerturb datasets on disk", flush=True)
    rows = []
    for i, f in enumerate(files):
        d = describe(f)
        rows.append(d)
        if "error" in d:
            print(f"  [{i+1}/{len(files)}] {f.stem[:44]:44s} skipped "
                  f"({d['error']})", flush=True)
        else:
            print(f"  [{i+1}/{len(files)}] {f.stem[:44]:44s} "
                  f"{d['n_contexts']:4d} ctx x {d['n_perturbations']:5d} pert "
                  f"x {d['n_replicates']} rep", flush=True)
    D = pd.DataFrame(rows)
    ok = D[D.get("error").isna()].copy() if "error" in D else D.copy()
    if not len(ok):
        print("nothing usable"); return

    ok["n_pairs"] = [n_pairs(r.n_contexts, r.n_perturbations, r.n_replicates)
                     for r in ok.itertuples()]
    ok["mds"] = [min_detectable_share(r.n_contexts, r.n_perturbations,
                                      r.n_replicates, SNR)
                 for r in ok.itertuples()]
    ok["resolves_target"] = ok.mds <= args.target
    ok["reps_needed"] = [required_replicates(args.target, r.n_contexts,
                                             r.n_perturbations, SNR)
                         for r in ok.itertuples()]
    ok = ok.sort_values("mds")
    ok.to_csv(TAB / "atlas_design_benchmark.csv", index=False)

    def binding(r):
        if r.n_contexts < 2:
            return "single context"
        if r.n_replicates < 2:
            return "no annotated replicates"
        if r.mds > args.target:
            return "underpowered"
        return "adequate"
    ok["binding"] = [binding(r) for r in ok.itertuples()]

    n_norep = int((ok.n_replicates < 2).sum())
    n_1ctx = int((ok.n_contexts < 2).sum())
    n_rc = int(ok.rep_col.notna().sum())
    print(f"\nBENCHMARK over {len(ok)} datasets")
    print(f"  {n_rc}/{len(ok)} carry any replicate-like annotation at all "
          f"(replicate/batch/\n     sample/donor); {n_norep} "
          f"({n_norep/len(ok):.0%}) have no CONDITION measured twice.")
    print(f"  {n_1ctx} ({n_1ctx/len(ok):.0%}) have a single context, so a "
          f"context x perturbation\n     interaction is undefined in them, not "
          f"merely underpowered.")
    est = ok[(ok.n_replicates >= 2) & (ok.n_contexts >= 2)]
    print(f"  {len(est)} ({len(est)/len(ok):.0%}) can support the estimator at "
          f"all.")
    if len(est):
        n_ok = int(est.resolves_target.sum())
        print(f"     of those, {n_ok} ({n_ok/len(est):.0%}) can detect a "
              f"{args.target:.0%} interaction; median smallest\n     "
              f"detectable share {est.mds.median():.3f}")
        print("\n  the datasets that CAN answer the question:")
        print(est.head(10)[["dataset", "n_contexts", "n_perturbations",
                            "n_replicates", "mds"]].round(4).to_string(
                                index=False))

    print("\n  binding constraint:")
    print(ok.binding.value_counts().to_string())

    print("\n  by perturbation type — a CRISPR screen in one line is "
          "single-context by\n  design, so the interesting row is the drug "
          "screens:")
    bt = pd.crosstab(ok.pert_type, ok.binding)
    print(bt.to_string())
    drug = ok[ok.pert_type.str.contains("drug|compound|chem", case=False,
                                        na=False)]
    if len(drug):
        n_multi = int((drug.n_contexts >= 2).sum())
        print(f"\n  {len(drug)} drug/compound datasets: {n_multi} use more than "
              f"one context,\n  {int(((drug.n_contexts >= 2) & (drug.n_replicates >= 2)).sum())}"
              f" of those also replicate a condition.")
    print(f"""
  CAVEAT, and it bounds the claim. This counts replicate structure RECOVERABLE
  FROM THE DEPOSITED METADATA, not what was run at the bench. An experiment may
  have had replicates that harmonisation dropped. That distinction does not
  weaken the practical point: an index of context-dependence computed by a
  reanalyst — or by the atlas's own authors from the released object — can only
  use annotation that is present, and where it is absent the cross-replicate
  covariance has nothing to average and the usual substitute is to split one
  well's cells, which sec.38 shows returns 0.500 from data containing no
  interaction whatever.

  (the five atlases of sec.46 are scored by this same rule and the same shared
  noise constant; their known outcomes are checked there, not refitted here)""")

    # atlases of sec.46, which are not in scPerturb (Tahoe is newer, LINCS is
    # bulk) and are shown as reference points
    REF = [("Tahoe-100M", 48, 95, 2), ("LINCS phase 1", 71, 831, 3),
           ("OP3", 6, 147, 3), ("Spear-ATAC", 3, 41, 5)]

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.5), constrained_layout=True)

    # (a) the design space, with the region where the estimator is defined
    a0 = ax[0]
    # the estimator needs BOTH >1 context and >1 replicate, so the region is a
    # corner and not a band
    a0.add_patch(plt.Rectangle((1.7, 1.7), 1e4, 1e4, color=AQUA, alpha=0.10,
                               zorder=0, lw=0))
    a0.axvline(1.7, color="#888", lw=1.0, ls=":", zorder=1)
    a0.axhline(1.7, color="#888", lw=1.0, ls=":", zorder=1)
    cols = {"CRISPR": BLUE, "drug": ORANGE, "cytokine": VIOLET}
    rng = np.random.default_rng(0)
    for t, g in ok.groupby("pert_type"):
        j = np.exp(rng.normal(0, 0.035, len(g)))
        a0.scatter(g.n_contexts * j, g.n_replicates * np.exp(
            rng.normal(0, 0.045, len(g))), s=34,
            color=cols.get(t, GREY), alpha=0.8, label=f"{t} ({len(g)})",
            edgecolor="white", lw=0.6, zorder=3)
    OFF = {"Tahoe-100M": (8, -11), "LINCS phase 1": (8, 4),
           "OP3": (8, 4), "Spear-ATAC": (8, 4)}
    for nm, c, p_, r in REF:
        a0.scatter([c], [r], marker="*", s=190, color="#111", zorder=4)
        a0.annotate(nm, (c, r), textcoords="offset points",
                    xytext=OFF.get(nm, (8, 4)), fontsize=7, color="#111",
                    zorder=5)
    a0.set_xscale("log"); a0.set_yscale("log")
    a0.set_xlim(0.72, 420); a0.set_ylim(0.72, 230)
    a0.set_xlabel("contexts (cell lines or cell types)")
    a0.set_ylabel("replicates per condition")
    a0.set_title("a  Where the field's designs sit", loc="left",
                 fontweight="bold", fontsize=9.5)
    a0.legend(frameon=False, fontsize=7.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=3, handletextpad=0.3,
              columnspacing=1.4)
    a0.text(0.97, 0.96, "estimator defined\nonly in this corner",
            transform=a0.transAxes, fontsize=7.5, color="#1b7a58",
            fontweight="bold", va="top", ha="right")

    # (b) the funnel
    a1 = ax[1]
    steps = [("all datasets", len(ok)),
             ("> 1 context", int((ok.n_contexts >= 2).sum())),
             ("+ a replicated\ncondition",
              int(((ok.n_contexts >= 2) & (ok.n_replicates >= 2)).sum())),
             (f"+ can detect\n{args.target:.0%}",
              int(((ok.n_contexts >= 2) & (ok.n_replicates >= 2)
                   & ok.resolves_target).sum()))]
    y = np.arange(len(steps))[::-1]
    a1.barh(y, [s[1] for s in steps], color=[BLUE, BLUE, ORANGE, AQUA],
            height=0.6)
    for yy, (lab, v) in zip(y, steps):
        a1.text(v + 0.6, yy, str(v), va="center", fontsize=9,
                fontweight="bold")
    a1.set_yticks(y, [s[0] for s in steps], fontsize=8)
    a1.set_xlabel("datasets")
    a1.set_xlim(0, len(ok) * 1.15)
    a1.set_title("b  How many can answer the question", loc="left",
                 fontweight="bold", fontsize=9.5)

    # (c) drug datasets only -- CRISPR screens are single-context by design
    a2 = ax[2]
    if len(drug):
        dd = drug.sort_values("n_contexts")
        yy = np.arange(len(dd))
        a2.barh(yy, dd.n_contexts, color=[AQUA if (c >= 2 and r >= 2) else ORANGE
                                          for c, r in zip(dd.n_contexts,
                                                          dd.n_replicates)],
                height=0.62)
        a2.set_yticks(yy, list(dd.dataset), fontsize=6)
        a2.axvline(2, color="#555", ls="--", lw=1.2)
        a2.set_xlabel("contexts")
        a2.set_xscale("log")
        a2.text(0.97, 0.06, "green = also has a\nreplicated condition",
                transform=a2.transAxes, ha="right", fontsize=7,
                color="#1b7a58", fontweight="bold")
    a2.set_title("c  Drug screens only", loc="left", fontweight="bold",
                 fontsize=9.5)

    fig.suptitle(f"The design calculation applied to all {len(ok)} scPerturb "
                 "datasets, uniformly and without tuning", fontsize=10.5,
                 x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "atlas_benchmark", FIG,
                    source_data={"designs": ok, "all": D}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
