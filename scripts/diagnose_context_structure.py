#!/usr/bin/env python3
"""Why does every unseen-line approach fail? Test the structure of context.

Define the context residual of a condition as
    R(line, drug, dose) = delta(line, drug, dose) - additive_prior(drug, dose)
i.e. what a line does differently from the average line for that treatment.

If line-specific modulation were a LINE property, R would be correlated across
different drugs within the same line — measure drug A on a new line and you
learn something about its response to drug B. If instead it is a per-(line,
drug) INTERACTION, cross-drug correlation is ~0, and no line-level descriptor
(covariates, baseline expression, a probe-fitted embedding) can ever transfer.

Reported per dataset:
  within_line_across_drug  mean r between R of drug pairs within a line
  within_drug_across_line  mean r between R of line pairs for the same drug
                           (positive control: shared drug-specific structure)
  same_drug_across_dose    mean r between doses of the same (line, drug)
                           (positive control: R is reproducible signal, not noise)
  cross_line_same_drug_pairs, permuted null for each.

Outputs: results/tables/context_structure.csv + figure bundle
         results/figures/06_diagnostics/context_structure/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "06_diagnostics"
TAB = ROOT / "results" / "tables"
MAX_PAIRS = 4000
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#9e9e9e"


def corr_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise Pearson r between paired row sets."""
    A = A - A.mean(1, keepdims=True)
    B = B - B.mean(1, keepdims=True)
    num = (A * B).sum(1)
    den = np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)
    return num / np.where(den > 0, den, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev47")
    ap.add_argument("--tag", default="dev47")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    allm = np.ones(len(G), dtype=bool)
    resp = responsive_genes(DELTA, allm)
    PRIOR = additive_prior(G, DELTA, allm, loo=True)   # leave-one-out prior
    R = (DELTA - PRIOR)[:, resp]                       # context residual
    print(f"{len(G)} conditions, residual over {len(resp)} responsive genes")

    lines = G.cell_line_id.to_numpy()
    drugs = G.drug.to_numpy()
    concs = G.conc.to_numpy()

    def sample_pairs(match_cols, differ_cols, n=MAX_PAIRS):
        """Pairs of conditions equal on match_cols and different on differ_cols."""
        key_m = pd.Series(list(zip(*[{"line": lines, "drug": drugs,
                                      "conc": concs}[c] for c in match_cols])))
        key_d = pd.Series(list(zip(*[{"line": lines, "drug": drugs,
                                      "conc": concs}[c] for c in differ_cols])))
        idx_by_m = {k: g.to_numpy() for k, g in
                    pd.Series(np.arange(len(G))).groupby(key_m)}
        out = []
        for k, idx in idx_by_m.items():
            if len(idx) < 2:
                continue
            for _ in range(min(len(idx), 40)):
                i, j = rng.choice(idx, 2, replace=False)
                if key_d[i] != key_d[j]:
                    out.append((i, j))
        if len(out) > n:
            out = [out[t] for t in rng.choice(len(out), n, replace=False)]
        return np.array(out)

    tests = {
        "within_line_across_drug": (["line"], ["drug"]),
        "within_drug_across_line": (["drug"], ["line"]),
        "same_line_drug_across_dose": (["line", "drug"], ["conc"]),
    }
    recs = []
    for name, (m, d) in tests.items():
        pairs = sample_pairs(m, d)
        if len(pairs) == 0:
            continue
        r = corr_rows(R[pairs[:, 0]], R[pairs[:, 1]])
        # permuted null: same pair count, random condition pairs
        pi = rng.integers(0, len(G), len(pairs))
        pj = rng.integers(0, len(G), len(pairs))
        rn = corr_rows(R[pi], R[pj])
        recs += [{"test": name, "kind": "observed", "r": v} for v in r]
        recs += [{"test": name, "kind": "null", "r": v} for v in rn]
        print(f"{name:28s} observed median r = {np.nanmedian(r):+.3f} "
              f"(null {np.nanmedian(rn):+.3f}, n={len(r)})")

    res = pd.DataFrame(recs)
    res.to_csv(TAB / f"context_structure_{args.tag}.csv", index=False)

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(8.5, 4.4), constrained_layout=True)
    order = [t for t in tests if (res.test == t).any()]
    labels = {"within_line_across_drug": "same line,\ndifferent drugs",
              "within_drug_across_line": "same drug,\ndifferent lines",
              "same_line_drug_across_dose": "same line+drug,\ndifferent doses"}
    colors = {"within_line_across_drug": BLUE,
              "within_drug_across_line": ORANGE,
              "same_line_drug_across_dose": AQUA}
    pos = 0
    ticks, ticklabels = [], []
    for t in order:
        for kind, alpha, col in (("observed", 0.65, colors[t]),
                                 ("null", 0.35, GREY)):
            v = res[(res.test == t) & (res.kind == kind)].r.dropna()
            vp = ax.violinplot([v], positions=[pos], widths=0.75,
                               showmedians=True, showextrema=False)
            vp["bodies"][0].set_facecolor(col)
            vp["bodies"][0].set_alpha(alpha)
            vp["bodies"][0].set_edgecolor("none")
            vp["cmedians"].set_color("#333333")
            ax.text(pos, 1.02, f"{np.nanmedian(v):+.2f}", ha="center",
                    fontsize=8, color="#555555")
            pos += 1
        ticks.append(pos - 1.5)
        ticklabels.append(labels[t])
        pos += 0.6
    ax.axhline(0, color="#888888", lw=0.9, ls="--", zorder=0)
    ax.set_xticks(ticks, ticklabels)
    ax.set_ylabel("Pearson r between context residuals")
    ax.set_ylim(-1.05, 1.12)
    ax.set_title("Is line-specific drug response a line property or a "
                 "line x drug interaction?", loc="left", fontweight="bold",
                 fontsize=10)
    ax.text(0.99, 0.03, "solid = observed pairs, grey = permuted null",
            transform=ax.transAxes, ha="right", fontsize=8, color="#666666")
    d = save_figure(fig, "context_structure", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
