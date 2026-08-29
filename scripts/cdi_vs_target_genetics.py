#!/usr/bin/env python3
"""Does a drug's context-dependence follow the genetics of its TARGET?

Our genotype scans asked whether a cell line's mutations predict how it
responds — and failed, because with ~50 lines a driver x mechanism scan is
underpowered. This inverts the question to the well-powered direction: we have
367 drugs with context-dependence indices (CDI) and annotated targets, so we
can ask whether a drug is context-dependent BECAUSE its target varies across
the panel.

Three target-level predictors, each computed over the atlas cell lines:

  target_mut_freq   fraction of lines carrying a driver alteration in any
                    annotated target of the drug
  target_expr_var   variance across lines of the target's baseline (DMSO)
                    expression, averaged over the drug's targets
  target_expr_mean  mean baseline expression of the targets (a control: a drug
                    whose target is simply not expressed should be inert, not
                    context-dependent)

Hypothesis: a drug whose target is mutated in some lines and not others, or
whose target's abundance varies, should show a more line-specific response.
Nuclear receptors would be the counter-case — NR3C1 is rarely mutated, so a
high CDI for glucocorticoids would have to come from downstream chromatin
context rather than target genetics, which is itself informative.

Outputs: results/tables/cdi_vs_target_genetics.csv
         figure bundle results/figures/10_drugs/cdi_vs_target_genetics/
"""
import argparse
import ast
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import load_pseudobulk
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
FIG = ROOT / "results" / "figures" / "10_drugs"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def parse_targets(x):
    if isinstance(x, (list, np.ndarray)):
        return [str(t).strip().upper() for t in x if str(t).strip()]
    if not isinstance(x, str) or not x.strip() or x.strip().lower() == "none":
        return []
    try:
        v = ast.literal_eval(x)
        if isinstance(v, (list, tuple)):
            return [str(t).strip().upper() for t in v if str(t).strip()]
    except Exception:
        pass
    return [t.strip().upper() for t in x.replace(";", ",").split(",")
            if t.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()

    cdi = pd.read_csv(TAB / "drug_context_dependence.csv")
    dm = pd.read_parquet(META / "drug_metadata.parquet")
    cdi = cdi.merge(dm[["drug", "targets"]], on="drug", how="left",
                    suffixes=("", "_dm"))
    tcol = "targets_dm" if "targets_dm" in cdi.columns else "targets"
    cdi["target_list"] = cdi[tcol].map(parse_targets)

    # --- baseline expression + mutation status across the atlas lines ---
    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    genes = pd.read_csv(Path(ROOT / args.pb_dir) / "genes.csv")
    dmso = cond[(cond.drug == "DMSO_TF") & (cond.plate != "plate14")]
    base = np.stack([np.average(X[g.row.to_numpy()], axis=0,
                                weights=g.n_cells)
                     for _, g in dmso.groupby("cell_line_id", observed=True)])
    lines = sorted(dmso.cell_line_id.unique())
    cpm = np.log1p(base / (base.sum(1, keepdims=True) + 1e-9) * 1e6)
    sym_pos = {}
    for i, s in enumerate(genes.gene_symbol.astype(str).str.upper()):
        sym_pos.setdefault(s, i)
    expr_var = cpm.var(0)
    expr_mean = cpm.mean(0)

    cl = pd.read_parquet(META / "cell_line_metadata.parquet")
    cl = cl[cl.Cell_ID_Cellosaur.isin(lines)]
    mut_lines = cl.groupby("Driver_Gene_Symbol").Cell_ID_Cellosaur.nunique()
    n_lines = len(lines)
    print(f"{n_lines} lines, {len(cdi)} drugs, "
          f"{cdi.target_list.map(len).gt(0).sum()} with annotated targets")

    rows = []
    for r in cdi.itertuples():
        ts = [t for t in r.target_list if t in sym_pos]
        if not ts:
            continue
        rows.append({
            "drug": r.drug, "cdi": r.cdi, "moa": getattr(r, "_5", None),
            "n_targets": len(ts),
            "target_mut_freq": float(np.mean(
                [mut_lines.get(t, 0) / n_lines for t in ts])),
            "target_expr_var": float(np.mean([expr_var[sym_pos[t]] for t in ts])),
            "target_expr_mean": float(np.mean([expr_mean[sym_pos[t]] for t in ts])),
        })
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "cdi_vs_target_genetics.csv", index=False)
    print(f"{len(df)} drugs with mappable targets\n")

    for col, label in (("target_mut_freq", "target mutation frequency"),
                       ("target_expr_var", "target expression variance"),
                       ("target_expr_mean", "target expression level")):
        rho = stats.spearmanr(df[col], df.cdi, nan_policy="omit")
        print(f"CDI vs {label:32s} rho={rho.statistic:+.3f}  p={rho.pvalue:.2e}")

    # nuclear receptors as the informative counter-case
    nr = df[df.drug.str.contains(
        "sone|solone|nide|Clobetasol|Triamcinolone|Betamethasone|tretinoin|"
        "Retino", case=False, na=False)]
    if len(nr):
        print(f"\nnuclear-receptor-like drugs (n={len(nr)}): "
              f"CDI {nr.cdi.median():.3f}, target mutation freq "
              f"{nr.target_mut_freq.median():.3f} "
              f"(all drugs: {df.cdi.median():.3f}, "
              f"{df.target_mut_freq.median():.3f})")

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    for ax, col, lab, c in ((axes[0], "target_mut_freq",
                             "fraction of lines with a driver mutation\nin the drug's target", BLUE),
                            (axes[1], "target_expr_var",
                             "variance of target expression across lines", AQUA),
                            (axes[2], "target_expr_mean",
                             "mean target expression (control)", ORANGE)):
        ax.scatter(df[col], df.cdi, s=12, alpha=0.5, color=c, edgecolors="none")
        rho = stats.spearmanr(df[col], df.cdi, nan_policy="omit")
        ax.set_xlabel(lab)
        ax.set_ylabel("context-dependence index")
        ax.set_title(f"rho={rho.statistic:+.2f}, p={rho.pvalue:.1e}",
                     loc="left", fontweight="bold", fontsize=9.5)
    fig.suptitle("Does a drug's context-dependence follow its target's genetics?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cdi_vs_target_genetics", FIG, source_data=df,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
