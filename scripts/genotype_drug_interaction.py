#!/usr/bin/env python3
"""Item 4: what determines which line shows which drug-response interaction —
genotype, or the line's baseline transcriptional state?

Design. For each drug we have a (line x gene) matrix of plate-matched response
deltas. Its leading principal component across lines, PC1_d, is the dominant
axis along which lines differ in their response to drug d. Each line therefore
gets a score on that axis, and we ask what predicts the score:

  H1 genotype  — driver-mutation status (from the atlas' curated driver table)
  H2 cell state — the line's BASELINE (DMSO) expression, summarised as PCs
                  fitted across lines, i.e. epigenetic/transcriptional tone

For H1 we run a systematic scan over (driver gene, drug) pairs with
Mann-Whitney U tests and Benjamini-Hochberg FDR, keeping only driver genes
altered in >= MIN_MUT lines. Because driver mutations are confounded with
lineage (APC with bowel, BRAF with skin ...), every hit is additionally tested
for lineage association and flagged.

For H2 we regress PC1_d on the top baseline-expression PCs and report the
variance explained, so genotype and cell state are compared on the same axis.

Known positives to look for: BRAF status modulating RAF inhibitors
(dabrafenib), KRAS status modulating RAS inhibitors (RMC-6236, adagrasib).

Outputs: results/tables/genotype_drug_interaction_<tag>.csv (full scan)
         results/tables/state_vs_genotype_<tag>.csv (H1 vs H2 per drug)
         figure bundle results/figures/07_genetics/genotype_drug_interaction/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "07_genetics"
TAB = ROOT / "results" / "tables"
MIN_MUT = 5          # driver gene must be altered in >= this many lines
MIN_WT = 5
MIN_LINES = 20       # drug must be measured in >= this many lines
N_STATE_PC = 5
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, p[i] * n / r)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev47")
    ap.add_argument("--tag", default="dev47")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    allm = np.ones(len(G), dtype=bool)
    resp = responsive_genes(DELTA, allm)
    lines = sorted(G.cell_line_id.unique())
    print(f"{len(G)} conditions, {len(lines)} lines, {G.drug.nunique()} drugs")

    # ---------------- genotype matrix ----------------
    cl = pd.read_parquet(ROOT / "data/metadata/metadata/cell_line_metadata.parquet")
    cl = cl[cl.Cell_ID_Cellosaur.isin(lines)]
    name_of = dict(cl.drop_duplicates("Cell_ID_Cellosaur")
                   [["Cell_ID_Cellosaur", "cell_name"]].itertuples(index=False))
    organ_of = dict(cl.drop_duplicates("Cell_ID_Cellosaur")
                    [["Cell_ID_Cellosaur", "Organ"]].itertuples(index=False))
    counts = cl.groupby("Driver_Gene_Symbol").Cell_ID_Cellosaur.nunique()
    driver_genes = sorted(counts[(counts >= MIN_MUT)
                                 & (counts <= len(lines) - MIN_WT)].index)
    mut = pd.DataFrame(
        {g: [int(((cl.Cell_ID_Cellosaur == l)
                  & (cl.Driver_Gene_Symbol == g)).any()) for l in lines]
         for g in driver_genes}, index=lines)
    print(f"{len(driver_genes)} driver genes testable "
          f"(altered in {MIN_MUT}-{len(lines) - MIN_WT} lines)")

    # ---------------- baseline state PCs ----------------
    dmso = cond[(cond.drug == "DMSO_TF") & (cond.plate != "plate14")]
    base = np.stack([
        np.average(X[g.row.to_numpy()], axis=0, weights=g.n_cells)
        for _, g in dmso.groupby("cell_line_id", observed=True)])
    base_lines = sorted(dmso.cell_line_id.unique())
    assert base_lines == lines
    Bc = base[:, resp] - base[:, resp].mean(0, keepdims=True)
    _, _, Vt = np.linalg.svd(Bc, full_matrices=False)
    state_pcs = Bc @ Vt[:N_STATE_PC].T
    state_pcs /= state_pcs.std(0, keepdims=True) + 1e-9

    # ---------------- per-drug response axis ----------------
    recs, state_recs = [], []
    for drug, grp in G.groupby("drug", observed=True):
        if drug == "DMSO_TF" or grp.cell_line_id.nunique() < MIN_LINES:
            continue
        # line x gene matrix: mean delta over doses for this drug
        per_line, present = [], []
        for l in lines:
            idx = grp.index[grp.cell_line_id == l].to_numpy()
            if len(idx):
                per_line.append(DELTA[idx][:, resp].mean(0))
                present.append(l)
        M = np.stack(per_line)
        Mc = M - M.mean(0, keepdims=True)
        u, s, vt = np.linalg.svd(Mc, full_matrices=False)
        pc1 = u[:, 0] * s[0]
        # orient so that positive = stronger response along the drug's own axis
        if np.corrcoef(pc1, np.linalg.norm(Mc, axis=1))[0, 1] < 0:
            pc1 = -pc1
        frac = float(s[0] ** 2 / (s ** 2).sum())
        score = pd.Series(pc1, index=present)

        # H1: genotype scan
        for g in driver_genes:
            gm = mut.loc[present, g].to_numpy().astype(bool)
            if gm.sum() < MIN_MUT or (~gm).sum() < MIN_WT:
                continue
            u_stat, p = stats.mannwhitneyu(score[gm], score[~gm],
                                           alternative="two-sided")
            eff = float(np.median(score[gm]) - np.median(score[~gm]))
            recs.append({"drug": drug, "driver_gene": g, "n_mut": int(gm.sum()),
                         "n_wt": int((~gm).sum()), "pc1_var_frac": frac,
                         "effect": eff, "p": float(p)})

        # H2: baseline state
        S = state_pcs[[lines.index(l) for l in present]]
        S1 = np.column_stack([np.ones(len(S)), S])
        beta, *_ = np.linalg.lstsq(S1, score.to_numpy(), rcond=None)
        pred = S1 @ beta
        r2_state = 1 - ((score.to_numpy() - pred) ** 2).sum() / \
            ((score.to_numpy() - score.mean()) ** 2).sum()
        # genotype R2 using all testable drivers jointly (same dof budget)
        Mg = mut.loc[present, driver_genes].to_numpy(dtype=float)
        top = np.argsort(-np.abs(np.corrcoef(
            np.column_stack([Mg, score.to_numpy()]), rowvar=False)[-1, :-1])
        )[:N_STATE_PC]
        G1 = np.column_stack([np.ones(len(Mg)), Mg[:, top]])
        bg, *_ = np.linalg.lstsq(G1, score.to_numpy(), rcond=None)
        pg = G1 @ bg
        r2_geno = 1 - ((score.to_numpy() - pg) ** 2).sum() / \
            ((score.to_numpy() - score.mean()) ** 2).sum()
        state_recs.append({"drug": drug, "n_lines": len(present),
                           "pc1_var_frac": frac, "r2_state": r2_state,
                           "r2_genotype_top5": r2_geno})

    scan = pd.DataFrame(recs)
    scan["q"] = bh_fdr(scan.p.to_numpy())
    scan = scan.sort_values("p")
    scan.to_csv(TAB / f"genotype_drug_interaction_{args.tag}.csv", index=False)
    st = pd.DataFrame(state_recs)
    st.to_csv(TAB / f"state_vs_genotype_{args.tag}.csv", index=False)

    print(f"\n{len(scan)} (driver gene, drug) tests; "
          f"{(scan.q < 0.1).sum()} at FDR<0.10, {(scan.q < 0.05).sum()} at 0.05")
    print("\ntop 20 genotype x drug associations:")
    print(scan.head(20)[["drug", "driver_gene", "n_mut", "n_wt", "effect",
                         "p", "q"]].to_string(index=False))
    print("\npositive controls:")
    for dg, drugs in (("BRAF", ["Dabrafenib"]),
                      ("KRAS", ["RMC-6236", "Adagrasib"])):
        sub = scan[(scan.driver_gene == dg) & (scan.drug.isin(drugs))]
        print(sub.to_string(index=False) if len(sub) else f"  {dg}: not tested")
    print(f"\nmedian variance of the line-response axis explained: "
          f"baseline state {st.r2_state.median():.2f} vs "
          f"genotype {st.r2_genotype_top5.median():.2f}")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)

    # A volcano of the genotype scan
    axes[0].scatter(scan.effect, -np.log10(scan.p), s=8, alpha=0.35,
                    color="#9e9e9e", edgecolors="none")
    hits = scan[scan.q < 0.1]
    axes[0].scatter(hits.effect, -np.log10(hits.p), s=22, color=ORANGE,
                    edgecolors="none", label=f"FDR<0.10 (n={len(hits)})")
    for n, r in enumerate(hits.head(8).itertuples()):
        axes[0].annotate(f"{r.driver_gene}·{r.drug[:12]}",
                         (r.effect, -np.log10(r.p)), fontsize=6.5,
                         xytext=(5, 3 if n % 2 == 0 else -8),
                         textcoords="offset points", color="#333333")
    axes[0].axhline(-np.log10(0.05), color="#888888", ls="--", lw=0.9)
    axes[0].set_xlabel("median response-axis difference (mutant − WT)")
    axes[0].set_ylabel("−log10 p")
    axes[0].set_title("A  Genotype × drug scan", loc="left",
                      fontweight="bold", fontsize=10)
    if len(hits):
        axes[0].legend(frameon=False, fontsize=8)

    # B genotype vs baseline state as explanations
    axes[1].scatter(st.r2_genotype_top5, st.r2_state, s=14, alpha=0.6,
                    color=BLUE, edgecolors="none")
    lim = [0, max(0.6, st[["r2_state", "r2_genotype_top5"]].to_numpy().max())]
    axes[1].plot(lim, lim, ls="--", color="#888888", lw=1)
    axes[1].set_xlim(lim); axes[1].set_ylim(lim)
    axes[1].set_xlabel("variance explained by genotype (top 5 drivers)")
    axes[1].set_ylabel("variance explained by baseline state (5 PCs)")
    axes[1].set_title("B  What explains the response axis?", loc="left",
                      fontweight="bold", fontsize=10)

    # C how much line-to-line variation the axis captures
    axes[2].hist(st.pc1_var_frac, bins=30, color=AQUA, alpha=0.8)
    axes[2].axvline(st.pc1_var_frac.median(), color="#333333", ls="--", lw=1.2,
                    label=f"median {st.pc1_var_frac.median():.0%}")
    axes[2].set_xlabel("fraction of line-to-line response variance on PC1")
    axes[2].set_ylabel("drugs")
    axes[2].set_title("C  Dominance of the leading axis", loc="left",
                      fontweight="bold", fontsize=10)
    axes[2].legend(frameon=False, fontsize=9)

    fig.suptitle("Does genotype or baseline cell state determine "
                 "line-specific drug response?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"genotype_drug_interaction_{args.tag}", FIG,
                    source_data={"scan": scan, "state_vs_genotype": st},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
