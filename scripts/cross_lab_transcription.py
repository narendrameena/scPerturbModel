#!/usr/bin/env python3
"""Does the cross-laboratory result hold for transcription, or only viability?

The viability arm (cross_lab_reproducibility.py) found that ~56% of the
reproducible line-specific drug response survives the move from PRISM (Broad,
pooled-barcode viability) to GDSC (Sanger, fitted IC50). That comparison changes
laboratory AND assay at once, and it measures a scalar phenotype. If the number
is a property of laboratories rather than of the killing assay, it should recur
for a transcriptional readout.

Two comparisons, matched in construction:

  within-lab   LINCS phase 1 (GSE92742) vs phase 2 (GSE70138) — both Broad, both
               L1000, different experiments years apart. 1,034 shared compounds,
               18 shared cell lines. This is the transcriptional ceiling.
  cross-lab    Tahoe-100M (Vevo/Arc, single-cell) vs LINCS (Broad, bulk L1000).
               ~133 shared compounds.

**The metric differs from the viability arm for a reason.** There the residual
was a scalar and could be correlated ACROSS cell lines, needing >=20 shared lines
per compound. Here only a handful of lines are shared, so correlating across
lines would be hopeless — but each (line, compound) carries a whole residual
GENE VECTOR, so the profile itself can be correlated over hundreds of genes. That
asks the same question ("does this line's compound-specific deviation reproduce")
with the power placed on the gene axis instead of the line axis.

Residuals are the project's standard construction: the line's profile minus the
leave-one-line-out mean profile for that compound, so the compound's shared
effect is removed and only the context-specific part is compared.

Outputs: results/tables/cross_lab_transcription.csv
         figure bundle results/figures/16_crosslab/cross_lab_transcription/
"""
import argparse
import gzip
import re
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
LIN = ROOT / "data" / "external" / "lincs"
FIG = ROOT / "results" / "figures" / "16_crosslab"
TAB = ROOT / "results" / "tables"
MIN_GENES = 100
MIN_LINES_PER_CPD = 3
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load_lincs(gctx, inst_info, gene_info, tag):
    """Per (cell line, compound) residual profile on the landmark genes."""
    gi = pd.read_csv(gene_info, sep="\t")
    lm = gi[gi.pr_is_lm == 1] if "pr_is_lm" in gi.columns else gi
    landmark = lm.pr_gene_id.astype(str).tolist()
    sym = dict(zip(gi.pr_gene_id.astype(str), gi.pr_gene_symbol.astype(str)))
    inst = pd.read_csv(inst_info, sep="\t", low_memory=False)
    cp = inst[inst.pert_type == "trt_cp"].copy()
    cp["k"] = cp.pert_iname.map(norm)
    keep = cp.groupby("k").cell_id.nunique().loc[lambda s: s >= 2].index
    cp = cp[cp.k.isin(keep)]
    print(f"  {tag}: {len(cp)} instances, {cp.k.nunique()} compounds, "
          f"{cp.cell_id.nunique()} lines", flush=True)
    from cmapPy.pandasGEXpress.parse_gctx import parse
    g = parse(str(gctx), rid=landmark, cid=cp.inst_id.tolist())
    M = g.data_df.T
    meta = cp.set_index("inst_id").loc[M.index]
    genes = np.array([sym.get(c, c) for c in M.columns])
    X = M.to_numpy(dtype=np.float32)
    out = {}
    for k, gg in meta.groupby("k", observed=True):
        pos = {c: i for i, c in enumerate(M.index)}
        by_line = {}
        for cid, gl in gg.groupby("cell_id", observed=True):
            ii = [pos[x] for x in gl.index]
            by_line[cid] = X[ii].mean(0)
        if len(by_line) < 2:
            continue
        names = list(by_line)
        A = np.stack([by_line[c] for c in names])
        tot = A.sum(0)
        for j, c in enumerate(names):
            out[(c, k)] = A[j] - (tot - A[j]) / (len(names) - 1)
    return out, genes


def load_tahoe(pb_dir):
    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk)
    X, cond = load_pseudobulk(ROOT / pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    genes = np.array(cond["genes"]) if isinstance(cond, dict) and "genes" in cond \
        else None
    cl = pd.read_parquet(ROOT / "data/metadata/metadata/cell_line_metadata.parquet")
    cvcl2name = dict(zip(cl.Cell_ID_Cellosaur.astype(str),
                         cl.cell_name.astype(str)))
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug})
    K["cname"] = K.line.astype(str).map(cvcl2name).fillna(K.line.astype(str))
    K["k"] = K.drug.map(norm)
    out = {}
    for k, gg in K.groupby("k", observed=True):
        by_line = {}
        for cn, gl in gg.groupby("cname", observed=True):
            by_line[norm(cn)] = DELTA[gl.i.to_numpy()].mean(0)
        if len(by_line) < 2:
            continue
        names = list(by_line)
        A = np.stack([by_line[c] for c in names])
        tot = A.sum(0)
        for j, c in enumerate(names):
            out[(c, k)] = A[j] - (tot - A[j]) / (len(names) - 1)
    return out, G


def compare(a, ga, b, gb, label):
    """Correlate residual profiles on the shared genes, per (line, compound)."""
    ia = {g: i for i, g in enumerate(ga)}
    shared = [g for g in gb if g in ia]
    if len(shared) < MIN_GENES:
        print(f"  {label}: only {len(shared)} shared genes — skipped")
        return []
    ai = np.array([ia[g] for g in shared])
    ib = {g: i for i, g in enumerate(gb)}
    bi = np.array([ib[g] for g in shared])
    keys = set(a) & set(b)
    rows = []
    for (line, cpd) in keys:
        x = a[(line, cpd)][ai]
        y = b[(line, cpd)][bi]
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        r = stats.pearsonr(x, y)
        rows.append({"comparison": label, "line": line, "compound": cpd,
                     "n_genes": len(shared), "r": float(r.statistic),
                     "p": float(r.pvalue)})
    print(f"  {label}: {len(rows)} (line, compound) pairs over "
          f"{len(shared)} shared genes", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    for gz, gc in ((LIN / "level4.gctx.gz", LIN / "level4.gctx"),
                   (LIN / "p1_level4.gctx.gz", LIN / "p1_level4.gctx")):
        if not gc.exists() and gz.exists():
            print(f"decompressing {gz.name} ...", flush=True)
            with gzip.open(gz, "rb") as f, open(gc, "wb") as o:
                shutil.copyfileobj(f, o, length=1 << 24)

    print("loading LINCS phase 1 ...", flush=True)
    p1, g1 = load_lincs(LIN / "p1_level4.gctx", LIN / "p1_inst_info.txt.gz",
                        LIN / "p1_gene_info.txt.gz", "phase1")
    print("loading LINCS phase 2 ...", flush=True)
    p2, g2 = load_lincs(LIN / "level4.gctx", LIN / "inst_info.txt.gz",
                        LIN / "GSE70138_Broad_LINCS_gene_info.txt.gz", "phase2")

    rows = compare(p1, g1, p2, g2, "LINCS p1 vs p2 (within-lab, Broad)")

    print("loading Tahoe ...", flush=True)
    try:
        th, G = load_tahoe(args.pb_dir)
        gt = np.array(G.columns) if hasattr(G, "columns") else None
        import anndata  # gene names come from the pseudobulk var index
        gt = np.array(pd.read_parquet(
            ROOT / "data/metadata/metadata/gene_metadata.parquet")
            .gene_symbol.astype(str)) if gt is None else gt
        rows += compare(th, gt, p1, g1, "Tahoe vs LINCS p1 (CROSS-LAB)")
        rows += compare(th, gt, p2, g2, "Tahoe vs LINCS p2 (CROSS-LAB)")
    except Exception as e:
        print(f"  Tahoe arm unavailable: {e}")

    R = pd.DataFrame(rows)
    if not len(R):
        print("no comparisons produced"); return
    R.to_csv(TAB / "cross_lab_transcription.csv", index=False)
    S = R.groupby("comparison").r.agg(["median", "mean", "size",
                                       lambda s: float((s > 0).mean())])
    S.columns = ["median_r", "mean_r", "n_pairs", "frac_positive"]
    print("\n=== residual-profile correlation per (line, compound) ===")
    print(S.round(3).to_string())

    within = R[R.comparison.str.contains("within-lab")].r.median()
    cross = R[R.comparison.str.contains("CROSS-LAB")].r.median()
    if np.isfinite(within) and np.isfinite(cross) and within > 0:
        print(f"\nwithin-lab (Broad, p1 vs p2): {within:.3f}")
        print(f"cross-lab (Tahoe vs LINCS):   {cross:.3f}")
        print(f"REPRODUCIBLE FRACTION (transcription) = {cross/within:.1%}")
        print("  Compare the viability arm's 56%. Agreement between the two "
              "readouts would\n  make this a property of laboratories rather "
              "than of a particular assay.")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    order = list(S.index)
    cols = [AQUA if "within" in o else ORANGE for o in order]
    data = [R[R.comparison == o].r.dropna() for o in order]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax[0].set_xticks(range(1, len(order) + 1),
                     [o.split("(")[0].strip()[:18] for o in order],
                     fontsize=7, rotation=15, ha="right")
    ax[0].axhline(0, color="#444444", lw=0.9)
    ax[0].set_ylabel("residual-profile r")
    for i, d in enumerate(data):
        if len(d):
            ax[0].text(i + 1, np.median(d) + 0.01, f"{np.median(d):.2f}",
                       ha="center", fontsize=9, fontweight="bold")
    ax[0].set_title("A  Transcriptional residual reproducibility", loc="left",
                    fontweight="bold", fontsize=10)
    for o, c in zip(order, cols):
        d = R[R.comparison == o].r.dropna()
        if len(d):
            ax[1].hist(d, bins=40, histtype="step", lw=2, color=c, density=True,
                       label=o.split("(")[0].strip()[:22])
    ax[1].axvline(0, color="#444444", lw=0.9)
    ax[1].set_xlabel("residual-profile r"); ax[1].set_ylabel("density")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("B  Distributions", loc="left", fontweight="bold",
                    fontsize=10)
    fig.suptitle("Does the line-specific TRANSCRIPTIONAL response reproduce "
                 "across laboratories?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cross_lab_transcription", FIG,
                    source_data={"per_pair": R, "summary": S.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
