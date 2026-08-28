#!/usr/bin/env python3
"""Which genes carry context-dependent drug response? A variance decomposition.

Our modelling established that line-specific drug response is a faint line x
drug INTERACTION. This asks the biological follow-up: the interaction is faint
*on average over genes* — is it concentrated in an identifiable set of genes
and pathways?

Per gene g, over all (line, drug, dose) conditions:
  delta            = observed log1p-CPM response vs plate-matched control
  prior            = additive drug x dose main effect (mean over other lines)
  residual  R      = delta - prior                      (line-specific part)
  var_total(g)     = Var[delta]                         total response variance
  var_resid(g)     = Var[R]                             line-specific + noise
  var_interact(g)  = Cov[R_dose_i, R_dose_j] within the same (line, drug)
                     -> the REPRODUCIBLE part of R, i.e. real interaction,
                        since independent noise does not covary across doses
  noise(g)         = var_resid - var_interact

Reported per gene: interaction_fraction = var_interact / var_total, plus the
lines/drugs driving it. Genes are then ranked and handed to pathway analysis.

Outputs: results/tables/response_variance_decomposition.csv
         figure bundle results/figures/06_diagnostics/response_variance/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                load_pseudobulk)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "06_diagnostics"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev47")
    ap.add_argument("--tag", default="dev47")
    ap.add_argument("--min-expr", type=float, default=0.05,
                    help="min mean log1p CPM in controls to keep a gene")
    args = ap.parse_args()

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    allm = np.ones(len(G), dtype=bool)
    PRIOR = additive_prior(G, DELTA, allm, loo=True)
    R = DELTA - PRIOR
    genes = pd.read_csv(Path(ROOT / args.pb_dir) / "genes.csv")
    print(f"{len(G)} conditions x {DELTA.shape[1]} genes", flush=True)

    # keep genes actually expressed in controls (avoid variance from zeros)
    ctrl_rows = cond[(cond.drug == "DMSO_TF")].row.to_numpy()
    ctrl_mean = np.log1p(X[ctrl_rows]).mean(0)
    keep = np.where(ctrl_mean > args.min_expr)[0]
    print(f"{len(keep)} genes pass the expression filter", flush=True)

    var_total = DELTA[:, keep].var(axis=0)
    var_resid = R[:, keep].var(axis=0)

    # reproducible (interaction) variance: covariance of residuals between
    # different doses of the SAME (line, drug); independent noise cancels.
    pair_i, pair_j = [], []
    idx = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                        "drug": G.drug})
    for _, grp in idx.groupby(["line", "drug"], observed=True):
        v = grp.i.to_numpy()
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                pair_i.append(v[a]); pair_j.append(v[b])
    pair_i, pair_j = np.array(pair_i), np.array(pair_j)
    print(f"{len(pair_i)} within-(line,drug) dose pairs", flush=True)

    A = R[pair_i][:, keep]
    B = R[pair_j][:, keep]
    var_interact = ((A - A.mean(0)) * (B - B.mean(0))).mean(0)

    out = pd.DataFrame({
        "gene_symbol": genes.gene_symbol.to_numpy()[keep],
        "ensembl_id": genes.ensembl_id.to_numpy()[keep],
        "ctrl_mean_expr": ctrl_mean[keep],
        "var_total": var_total,
        "var_residual": var_resid,
        "var_interaction": var_interact,
        "noise": var_resid - var_interact,
    })
    out["interaction_fraction"] = (out.var_interaction / out.var_total).clip(0, 1)
    out["additive_fraction"] = ((out.var_total - out.var_residual)
                                / out.var_total).clip(0, 1)
    out["reproducibility"] = (out.var_interaction / out.var_residual).clip(0, 1)
    out = out.sort_values("var_interaction", ascending=False)
    out.to_csv(TAB / f"response_variance_decomposition_{args.tag}.csv",
               index=False)

    med_add = out.additive_fraction.median()
    med_int = out.interaction_fraction.median()
    print(f"\nmedian variance shares: additive {med_add:.1%}, "
          f"interaction {med_int:.1%}, noise {1 - med_add - med_int:.1%}")
    print("\ntop 25 context-dependent genes (highest interaction variance):")
    print(out.head(25)[["gene_symbol", "var_total", "var_interaction",
                        "interaction_fraction"]].to_string(index=False))

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)

    shares = [out.additive_fraction.median(), out.interaction_fraction.median()]
    shares.append(max(0.0, 1 - sum(shares)))
    axes[0].bar(range(3), shares, width=0.6, color=[BLUE, ORANGE, "#bdbdbd"])
    for x, v in enumerate(shares):
        axes[0].text(x, v + 0.01, f"{v:.0%}", ha="center", fontsize=10)
    axes[0].set_xticks(range(3), ["drug x dose\n(additive)",
                                  "line x drug\n(interaction)", "noise"])
    axes[0].set_ylabel("median share of response variance")
    axes[0].set_title("A  Where response variance lives", loc="left",
                      fontweight="bold", fontsize=10)

    axes[1].scatter(out.var_total, out.var_interaction, s=4, alpha=0.25,
                    color=BLUE, edgecolors="none")
    top = out[~out.gene_symbol.str.startswith("ENSG")].head(10)
    axes[1].scatter(top.var_total, top.var_interaction, s=20, color=ORANGE,
                    edgecolors="none")
    # stagger labels alternately to stop the top genes overprinting
    for n, r in enumerate(top.itertuples()):
        dx, dy = (6, 4) if n % 2 == 0 else (-6, -9)
        axes[1].annotate(r.gene_symbol, (r.var_total, r.var_interaction),
                         fontsize=7, color="#333333",
                         ha="left" if n % 2 == 0 else "right",
                         xytext=(dx, dy), textcoords="offset points")
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].set_xlabel("total response variance")
    axes[1].set_ylabel("interaction variance")
    axes[1].set_title("B  Context-dependent genes", loc="left",
                      fontweight="bold", fontsize=10)

    axes[2].hist(out.interaction_fraction, bins=60, color=AQUA, alpha=0.75)
    axes[2].axvline(med_int, color="#333333", ls="--", lw=1.2,
                    label=f"median {med_int:.1%}")
    axes[2].set_xlabel("interaction fraction of response variance")
    axes[2].set_ylabel("genes")
    axes[2].set_title("C  A long tail of context-dependent genes", loc="left",
                      fontweight="bold", fontsize=10)
    axes[2].legend(frameon=False, fontsize=9)

    fig.suptitle("Decomposition of drug-response variance across "
                 f"{G.cell_line_id.nunique()} cell lines and "
                 f"{G.drug.nunique()} drugs", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"response_variance_{args.tag}", FIG,
                    source_data=out, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
