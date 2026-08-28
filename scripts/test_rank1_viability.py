#!/usr/bin/env python3
"""Is the context residual a rank-1 viability effect, or genuinely idiosyncratic?

McFarland et al. (MIX-Seq, Nat Commun 2020) implicitly model the line-specific
part of a drug response as rank-1:

    residual(line, drug) ~ s(line, drug) * v          (v shared across all
                                                       lines and drugs)

with s the pair's sensitivity and v a common viability / stress / cell-cycle
direction. If true, that MECHANISTICALLY EXPLAINS our interaction term and
predicts its failure to transfer (sensitivity is a pair property). If false,
the interaction is high-dimensional and idiosyncratic — a stronger claim.

Noise is ~90% of the residual, so a naive SVD would mostly decompose noise.
Instead every component is validated by SPLIT-HALF over plates: for each
(line, drug) measured on >=2 plates we form two independent estimates A and B
from disjoint plate sets, fit components on A, and ask whether the per-pair
loadings reproduce on B. Only reproducible components are interpretable.

Reported:
  - reproducibility r of each component (A-loadings vs B-loadings)
  - share of REPRODUCIBLE residual variance carried by component 1
  - what component 1 is: top genes, and its correlation with an independent
    viability readout (cell-cycle arrest measured from the single-cell data)

Outputs: results/tables/rank1_viability_<tag>.csv
         figure bundle results/figures/06_diagnostics/rank1_viability_<tag>/
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
N_COMP = 10
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev47")
    ap.add_argument("--tag", default="dev47")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    allm = np.ones(len(G), dtype=bool)
    resp = responsive_genes(DELTA, allm)
    PRIOR = additive_prior(G, DELTA, allm, loo=True)
    R = (DELTA - PRIOR)[:, resp]
    genes = pd.read_csv(Path(ROOT / args.pb_dir) / "genes.csv")
    sym = genes.gene_symbol.to_numpy()[resp]
    print(f"{len(G)} plate-level conditions, {len(resp)} responsive genes")

    # ---------------- split-half over plates within (line, drug) ----------
    A_rows, B_rows, keys = [], [], []
    idx = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                        "drug": G.drug, "plate": G.plate, "conc": G.conc})
    for (ln, dr), grp in idx.groupby(["line", "drug"], observed=True):
        pls = grp.plate.unique()
        if len(pls) < 2:
            continue
        rng.shuffle(pls)
        half = max(1, len(pls) // 2)
        a = grp[grp.plate.isin(pls[:half])].i.to_numpy()
        b = grp[grp.plate.isin(pls[half:])].i.to_numpy()
        A_rows.append(R[a].mean(0))
        B_rows.append(R[b].mean(0))
        keys.append((ln, dr))
    A = np.stack(A_rows)
    B = np.stack(B_rows)
    keys = pd.DataFrame(keys, columns=["cell_line_id", "drug"])
    print(f"{len(A)} (line, drug) pairs with >=2 plates for split-half")

    Ac = A - A.mean(0, keepdims=True)
    Bc = B - B.mean(0, keepdims=True)

    # components fitted on half A only
    U, S, Vt = np.linalg.svd(Ac, full_matrices=False)
    comps = Vt[:N_COMP]                     # gene-space directions
    loadA = Ac @ comps.T
    loadB = Bc @ comps.T

    recs = []
    for k in range(N_COMP):
        r = np.corrcoef(loadA[:, k], loadB[:, k])[0, 1]
        var_share = float(S[k] ** 2 / (S ** 2).sum())
        # reproducible variance carried: covariance of the two halves on this axis
        cov_ab = float(np.mean((loadA[:, k] - loadA[:, k].mean())
                               * (loadB[:, k] - loadB[:, k].mean())))
        recs.append({"component": k + 1, "var_share_halfA": var_share,
                     "split_half_r": float(r), "reproducible_cov": cov_ab})
    comp_df = pd.DataFrame(recs)
    tot_repro = comp_df.reproducible_cov.clip(lower=0).sum()
    comp_df["share_of_reproducible"] = (comp_df.reproducible_cov.clip(lower=0)
                                        / max(tot_repro, 1e-12))
    print("\ncomponent reproducibility (split-half over plates):")
    print(comp_df.round(4).to_string(index=False))

    # total reproducible variance in the residual, for context
    total_cov = float(np.mean(np.sum(Ac * Bc, axis=1)))
    print(f"\ncomponent 1 carries {comp_df.share_of_reproducible[0]:.1%} of the "
          f"reproducible residual; top-3 {comp_df.share_of_reproducible[:3].sum():.1%}")

    # ---------------- what is component 1? ----------------
    v1 = comps[0]
    order = np.argsort(v1)
    top_pos = [sym[i] for i in order[::-1][:20]]
    top_neg = [sym[i] for i in order[:20]]
    print(f"\ncomponent 1 (+): {', '.join(top_pos[:12])}")
    print(f"component 1 (−): {', '.join(top_neg[:12])}")

    # independent viability readout: cell-cycle arrest from the single cells
    viab = None
    cc_path = ROOT / "results/tables/cell_cycle_log2or.csv"
    if cc_path.exists():
        cc = pd.read_csv(cc_path)
        cpc = pd.read_csv(ROOT / "results/tables/cells_per_condition.csv")
        id_of = dict(cpc.drop_duplicates("cell_name")[["cell_name", "cell_line"]]
                     .itertuples(index=False))
        cc["cell_line_id"] = cc.cell_name.map(id_of)
        g2m = (cc[cc.phase == "G2M"].groupby(["cell_line_id", "drug"],
                                             observed=True).log2_or.median())
        g1 = (cc[cc.phase == "G1"].groupby(["cell_line_id", "drug"],
                                           observed=True).log2_or.median())
        key_idx = pd.MultiIndex.from_frame(keys)
        viab = pd.DataFrame({
            "g2m": g2m.reindex(key_idx).to_numpy(),
            "g1": g1.reindex(key_idx).to_numpy(),
            "load1": loadA[:, 0], "load2": loadA[:, 1]})
        ok = viab.dropna()
        if len(ok) > 30:
            r_g2m = np.corrcoef(ok.load1, ok.g2m)[0, 1]
            r_g1 = np.corrcoef(ok.load1, ok.g1)[0, 1]
            print(f"\ncomponent-1 loading vs cell-cycle arrest (n={len(ok)}): "
                  f"G2M r={r_g2m:+.3f}, G1 r={r_g1:+.3f}")
            comp_df.attrs["r_g2m"] = r_g2m

    out = comp_df.copy()
    out.to_csv(TAB / f"rank1_viability_{args.tag}.csv", index=False)
    pd.DataFrame({"component1_gene": sym, "loading": v1}).to_csv(
        TAB / f"rank1_component1_genes_{args.tag}.csv", index=False)

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)

    axes[0].bar(comp_df.component, comp_df.share_of_reproducible, width=0.65,
                color=[ORANGE if c == 1 else BLUE for c in comp_df.component])
    axes[0].set_xlabel("component")
    axes[0].set_ylabel("share of reproducible residual variance")
    axes[0].set_title("A  Is the residual rank-1?", loc="left",
                      fontweight="bold", fontsize=10)
    axes[0].text(0.97, 0.93, f"comp 1: {comp_df.share_of_reproducible[0]:.0%}",
                 transform=axes[0].transAxes, ha="right", fontsize=9,
                 color="#444444")

    axes[1].bar(comp_df.component, comp_df.split_half_r, width=0.65, color=AQUA)
    axes[1].axhline(0, color="#888888", lw=0.9)
    axes[1].set_xlabel("component")
    axes[1].set_ylabel("split-half reproducibility (r)")
    axes[1].set_title("B  Which components are real?", loc="left",
                      fontweight="bold", fontsize=10)

    if viab is not None and len(viab.dropna()) > 30:
        ok = viab.dropna()
        axes[2].scatter(ok.load1, ok.g2m, s=9, alpha=0.4, color=VIOLET,
                        edgecolors="none")
        axes[2].set_xlabel("component-1 loading (per line × drug)")
        axes[2].set_ylabel("G2M arrest, log2 odds ratio")
        rr = np.corrcoef(ok.load1, ok.g2m)[0, 1]
        axes[2].set_title(f"C  Component 1 vs viability (r={rr:+.2f})",
                          loc="left", fontweight="bold", fontsize=10)
    else:
        axes[2].axis("off")

    fig.suptitle("Testing the rank-1 viability model of the context residual",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"rank1_viability_{args.tag}", FIG,
                    source_data={"components": comp_df,
                                 "component1_genes": pd.DataFrame(
                                     {"gene": sym, "loading": v1})},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
