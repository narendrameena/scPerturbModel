#!/usr/bin/env python3
"""Cell-type audit of the Tahoe-100M cell lines against Tabula Sapiens.

For every (tissue, cell_type) with >=MIN_CELLS cells in the independent
reference, build a mean log1p-CPM signature (raw counts from .raw when
present). Each Tahoe line's DMSO-control pseudobulk (dev47) is then ranked
against all signatures by Spearman correlation over the genes most variable
ACROSS signatures. Audit question: does each line's transcriptome match cell
types of its annotated organ of origin?

Outputs: results/tables/cell_type_audit.csv (+ signature matches),
         figure bundle results/figures/cell_type_audit/
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TS = ROOT / "data" / "external" / "tabula_sapiens"
PB = ROOT / "data" / "processed" / "pseudobulk_dev47"
FIG = ROOT / "results" / "figures" / "05_audit"
TAB = ROOT / "results" / "tables"
MIN_CELLS = 50
MAX_CELLS = 400
N_GENES = 3000
RNG = np.random.default_rng(0)


def tissue_signatures(path: Path):
    """Return (sig_df [types x genes], ensembl index) for one tissue h5ad."""
    import anndata as ad
    A = ad.read_h5ad(path, backed="r")
    raw = A.raw if A.raw is not None else A
    genes = pd.Index(raw.var_names.str.replace(r"\.\d+$", "", regex=True))
    out, names, ncells = [], [], []
    for ct, grp in A.obs.groupby("cell_type", observed=True):
        idx = np.where(A.obs.index.isin(grp.index))[0]
        if len(idx) < MIN_CELLS:
            continue
        if len(idx) > MAX_CELLS:
            idx = np.sort(RNG.choice(idx, MAX_CELLS, replace=False))
        X = raw.X[idx]
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        lib = X.sum(1, keepdims=True)
        sig = np.log1p(X / (lib + 1e-8) * 1e4).mean(0)
        out.append(sig)
        names.append(str(ct))
        ncells.append(len(idx))
    A.file.close()
    return pd.DataFrame(np.stack(out), index=names, columns=genes), ncells


def main():
    manifest = json.load(open(TS / "manifest.json"))
    tissue_organs = {t: m["organs"] for t, m in manifest.items()}

    # ---------------- reference signatures ----------------
    sigs, meta = [], []
    common = None
    per_tissue = {}
    for t in sorted(manifest):
        df, ncells = tissue_signatures(TS / manifest[t]["file"])
        per_tissue[t] = df
        common = df.columns if common is None else common.intersection(df.columns)
        print(f"{t}: {len(df)} cell types", flush=True)
    for t, df in per_tissue.items():
        sigs.append(df[common].to_numpy())
        meta += [(t, ct) for ct in df.index]
    S = np.concatenate(sigs)                       # (n_sigs, n_common_genes)
    meta = pd.DataFrame(meta, columns=["tissue", "cell_type"])
    print(f"{len(S)} signatures over {len(common)} common genes")

    # ---------------- Tahoe DMSO pseudobulk per line ----------------
    counts = np.load(PB / "pseudobulk_counts.npz")["counts"]
    cond = pd.read_csv(PB / "conditions.csv")
    gm = pd.read_parquet(ROOT / "data/metadata/metadata/gene_metadata.parquet"
                         ).sort_values("token_id").reset_index(drop=True)
    cpc = pd.read_csv(ROOT / "results/tables/cells_per_condition.csv")
    name_of = dict(cpc.drop_duplicates("cell_line")[["cell_line", "cell_name"]]
                   .itertuples(index=False))
    cl = pd.read_parquet(ROOT / "data/metadata/metadata/cell_line_metadata.parquet")
    organ_of = dict(cl.drop_duplicates("Cell_ID_Cellosaur")
                    [["Cell_ID_Cellosaur", "Organ"]].itertuples(index=False))

    dmso = cond[cond.drug == "DMSO_TF"]
    lines = sorted(dmso.cell_line_id.unique())
    profiles = {}
    for ln in lines:
        rows = dmso[dmso.cell_line_id == ln]
        v = np.average(counts[rows.index], axis=0, weights=rows.n_cells)
        profiles[ln] = np.log1p(v / (v.sum() + 1e-9) * 1e4)

    # map Tahoe genes -> common reference genes
    tah_pos = pd.Series(np.arange(len(gm)), index=gm.ensembl_id)
    shared = [g for g in common if g in tah_pos.index]
    ref_col = {g: i for i, g in enumerate(common)}
    S_shared = S[:, [ref_col[g] for g in shared]]
    tah_idx = tah_pos[shared].to_numpy()
    print(f"{len(shared)} genes shared with Tahoe")

    top = np.argsort(S_shared.var(axis=0))[-N_GENES:]
    S_top = S_shared[:, top]
    S_rank = np.apply_along_axis(rankdata, 1, S_top)
    S_rank -= S_rank.mean(1, keepdims=True)
    S_norm = S_rank / np.linalg.norm(S_rank, axis=1, keepdims=True)

    recs, heat = [], np.zeros((len(lines), len(per_tissue)))
    tissues = sorted(per_tissue)
    for li, ln in enumerate(lines):
        p = profiles[ln][tah_idx][top]
        pr = rankdata(p)
        pr -= pr.mean()
        pr /= np.linalg.norm(pr)
        rho = S_norm @ pr                      # Spearman via ranked cosine
        order = np.argsort(rho)[::-1]
        organs = tissue_organs
        expected = [t for t in tissues if organ_of.get(ln) in organs[t]]
        for k, t in enumerate(tissues):
            mask = (meta.tissue == t).to_numpy()
            heat[li, k] = rho[mask].max() if mask.any() else np.nan
        top3 = [(meta.tissue[i], meta.cell_type[i], float(rho[i]))
                for i in order[:3]]
        recs.append({
            "cell_line_id": ln, "cell_name": name_of.get(ln, ln),
            "organ": organ_of.get(ln, "?"),
            "has_reference_tissue": bool(expected),
            "top1": f"{top3[0][0]}|{top3[0][1]}", "rho1": round(top3[0][2], 3),
            "top2": f"{top3[1][0]}|{top3[1][1]}", "rho2": round(top3[1][2], 3),
            "top3": f"{top3[2][0]}|{top3[2][1]}", "rho3": round(top3[2][2], 3),
            "match_top1": top3[0][0] in expected if expected else None,
            "match_top3": any(t in expected for t, _, _ in top3)
            if expected else None,
        })
    res = pd.DataFrame(recs)

    # column-normalized view: z-score each tissue column across lines, so a
    # generically-correlated tissue (proliferative epithelium) cannot dominate;
    # asks "is THIS line more similar to tissue T than other lines are?"
    Hdf = pd.DataFrame(heat, index=[r["cell_name"] for r in recs],
                       columns=tissues)
    Z = (Hdf - Hdf.mean(0)) / Hdf.std(0)
    zm1, zm3 = [], []
    for i, r in enumerate(recs):
        expected = [t for t in tissues
                    if organ_of.get(r["cell_line_id"]) in tissue_organs[t]]
        if not expected:
            zm1.append(None); zm3.append(None)
            continue
        order = Z.iloc[i].sort_values(ascending=False).index.tolist()
        zm1.append(order[0] in expected)
        zm3.append(any(t in expected for t in order[:3]))
    res["zmatch_top1"], res["zmatch_top3"] = zm1, zm3
    res.to_csv(TAB / "cell_type_audit.csv", index=False)

    with_ref = res[res.has_reference_tissue]
    print(f"\norgan concordance (n={len(with_ref)} auditable lines):")
    print(f"  raw          top1 {with_ref.match_top1.mean():.1%}, "
          f"top3 {with_ref.match_top3.mean():.1%}")
    print(f"  column-znorm top1 {with_ref.zmatch_top1.mean():.1%}, "
          f"top3 {with_ref.zmatch_top3.mean():.1%}")
    print("\ndiscordant after z-normalization (top3 miss):")
    print(with_ref[~with_ref.zmatch_top3.astype(bool)]
          [["cell_name", "organ", "top1", "rho1"]].to_string(index=False))

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 8, "axes.grid": False,
                         "figure.facecolor": "white"})
    order_lines = np.argsort([res.organ[i] for i in range(len(lines))])
    fig, axes = plt.subplots(1, 3, figsize=(19, 9), constrained_layout=True,
                             gridspec_kw={"width_ratios": [3, 3, 1]})

    def organ_boxes(ax):
        for row, i in enumerate(order_lines):
            for k, t in enumerate(tissues):
                if organ_of.get(lines[i]) in tissue_organs[t]:
                    ax.add_patch(plt.Rectangle((k - .5, row - .5), 1, 1,
                                               fill=False, edgecolor="#eb6834",
                                               lw=1.2))
        ax.set_xticks(range(len(tissues)), tissues, rotation=45, ha="right")

    im = axes[0].imshow(heat[order_lines], aspect="auto", cmap="Blues")
    organ_boxes(axes[0])
    axes[0].set_yticks(range(len(lines)),
                       [f"{res.cell_name[i]} ({res.organ[i]})"
                        for i in order_lines], fontsize=6.5)
    fig.colorbar(im, ax=axes[0], label="max Spearman rho vs tissue cell types",
                 shrink=0.5)
    axes[0].set_title("A  Raw similarity (orange = annotated organ)",
                      loc="left", fontweight="bold", fontsize=10)

    zmax = float(np.nanmax(np.abs(Z.to_numpy())))
    imz = axes[1].imshow(Z.to_numpy()[order_lines], aspect="auto",
                         cmap="RdBu_r", vmin=-zmax, vmax=zmax)
    organ_boxes(axes[1])
    axes[1].set_yticks([])
    fig.colorbar(imz, ax=axes[1],
                 label="z-score within tissue column (across lines)",
                 shrink=0.5)
    axes[1].set_title("B  Column-normalized (removes generic-epithelium pull)",
                      loc="left", fontweight="bold", fontsize=10)

    rates = [with_ref.match_top1.mean(), with_ref.match_top3.mean(),
             with_ref.zmatch_top1.mean(), with_ref.zmatch_top3.mean()]
    cols = ["#9ec4ea", "#9ec4ea", "#2a78d6", "#2a78d6"]
    axes[2].bar(range(4), rates, width=0.6, color=cols)
    for x, v in enumerate(rates):
        axes[2].text(x, v + 0.02, f"{v:.0%}", ha="center", fontsize=9)
    axes[2].set_xticks(range(4), ["raw\ntop-1", "raw\ntop-3",
                                  "znorm\ntop-1", "znorm\ntop-3"], fontsize=8)
    axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel("organ concordance")
    axes[2].set_title("C  Concordance", loc="left", fontweight="bold",
                      fontsize=10)
    fig.suptitle("Cell-type audit: Tahoe lines vs Tabula Sapiens "
                 f"({len(S)} signatures, {len(with_ref)} auditable lines)",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cell_type_audit", FIG,
                    source_data={"audit": res,
                                 "heatmap": Hdf.iloc[order_lines].reset_index(
                                     names="cell_name"),
                                 "heatmap_znorm": Z.iloc[order_lines]
                                 .reset_index(names="cell_name")},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
