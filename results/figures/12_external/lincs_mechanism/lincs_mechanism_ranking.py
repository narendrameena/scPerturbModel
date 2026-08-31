#!/usr/bin/env python3
"""Independent verification of the mechanism ranking in LINCS L1000.

Our finding that a drug's context-dependence is set by its mechanism — with
nuclear-receptor agonists least transferable — rests on a single atlas
(Tahoe-100M). LINCS L1000 is an independent test bed: a different platform
(direct 978-gene measurement), a different lab, bulk rather than single-cell,
30 cell lines, 1,796 compounds, and — critically — 99% of (compound, cell) pairs
measured on two or more plates, so the replicate-based interaction estimate is
well supported, unlike OP3 where it was not.

Two design choices matter:

* **Landmark genes only.** The 12,328-gene matrices are ~11,350 genes imputed as
  linear functions of the 978 measured landmarks; including them would
  manufacture correlation structure. We use the 978 directly measured genes.
* **Level 4, not Level 5.** Level 5 aggregates replicates away — only 5% of its
  conditions retain more than one signature. Level 4 keeps one profile per well,
  so plates provide the independent replicate axis.

Mechanism labels are taken from the Tahoe drug metadata by compound-name match,
so both datasets are annotated identically and the rankings are comparable.

Outputs: results/tables/lincs_mechanism_ranking.csv
         results/tables/lincs_vs_tahoe_cdi.csv
         figure bundle results/figures/12_external/lincs_mechanism/
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

from perturbmodel.decompose import decompose
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
LIN = ROOT / "data" / "external" / "lincs"
FIG = ROOT / "results" / "figures" / "12_external"
TAB = ROOT / "results" / "tables"
MIN_CELL_LINES = 10
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gctx", default=str(LIN / "level4.gctx"))
    ap.add_argument("--min-lines", type=int, default=MIN_CELL_LINES)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    gz = LIN / "level4.gctx.gz"
    gctx = Path(args.gctx)
    if not gctx.exists() and gz.exists():
        print("decompressing Level 4 ...", flush=True)
        with gzip.open(gz, "rb") as f, open(gctx, "wb") as o:
            shutil.copyfileobj(f, o, length=1 << 24)

    gi = pd.read_csv(LIN / "GSE70138_Broad_LINCS_gene_info.txt.gz", sep="\t")
    lm_col = "pr_is_lm" if "pr_is_lm" in gi.columns else None
    landmark = (gi[gi[lm_col] == 1] if lm_col else gi).pr_gene_id.astype(str).tolist()
    sym = dict(zip(gi.pr_gene_id.astype(str), gi.pr_gene_symbol.astype(str)))
    print(f"{len(landmark)} landmark genes", flush=True)

    inst = pd.read_csv(LIN / "inst_info.txt.gz", sep="\t", low_memory=False)
    inst["plate"] = inst.det_plate.astype(str).str.rsplit("_", n=1).str[0]
    cp = inst[inst.pert_type == "trt_cp"]
    keep_cpd = (cp.groupby("pert_id").cell_id.nunique()
                .loc[lambda s: s >= args.min_lines].index)
    sel = inst[((inst.pert_type == "trt_cp") & inst.pert_id.isin(keep_cpd))
               | (inst.pert_type == "ctl_vehicle")].copy()
    print(f"{len(sel)} instances: {sel.pert_id.nunique()-1} compounds in "
          f">={args.min_lines} lines, {sel.cell_id.nunique()} cell lines",
          flush=True)

    from cmapPy.pandasGEXpress.parse_gctx import parse
    g = parse(str(gctx), rid=landmark, cid=sel.inst_id.tolist())
    M = g.data_df.T                                   # instances x genes
    print(f"loaded {M.shape[0]} x {M.shape[1]}", flush=True)

    import anndata as ad
    meta = sel.set_index("inst_id").loc[M.index]
    meta["pert"] = np.where(meta.pert_type == "ctl_vehicle", "DMSO",
                            meta.pert_id.astype(str))
    meta["dose"] = meta.pert_dose.round(2).astype(str)
    obs = meta[["cell_id", "pert", "dose", "plate", "pert_iname"]].copy()
    obs.index = obs.index.astype(str)
    A = ad.AnnData(X=M.to_numpy(dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=[sym.get(c, c) for c in M.columns]))

    res = decompose(A, context="cell_id", perturbation="pert", control="DMSO",
                    replicate="plate", dose="dose", n_genes=len(landmark))
    print("\n" + res.report(), flush=True)

    per = res.per_perturbation.copy()
    name = dict(zip(meta.pert_id.astype(str), meta.pert_iname.astype(str)))
    per["pert_iname"] = per.perturbation.map(name)
    per["key"] = per.pert_iname.map(norm)
    # Broad Repurposing Hub gives MOA for ~6,800 compounds; far better coverage
    # for LINCS than name-matching to the Tahoe panel (23 of 273 matched)
    hub_f = LIN / "repurposing_drugs.txt"
    if hub_f.exists():
        hub = pd.read_csv(hub_f, sep="\t", comment="!", low_memory=False)
        hub["key"] = hub.pert_iname.map(norm)
        hub = hub.drop_duplicates("key")[["key", "moa"]].rename(
            columns={"moa": "moa-fine"})
        per = per.merge(hub, on="key", how="left")
    else:
        dm = pd.read_parquet(ROOT / "data/metadata/metadata/drug_metadata.parquet")
        dm["key"] = dm.drug.map(norm)
        per = per.merge(dm[["key", "moa-fine"]], on="key", how="left")
    per.to_csv(TAB / "lincs_mechanism_ranking.csv", index=False)

    ok = per[per.estimable]
    print(f"\n{len(ok)}/{len(per)} compounds estimable "
          f"(>=25 cross-plate pairs); median index {ok['index'].median():.3f}")
    known = ok[ok["moa-fine"].notna() & (ok["moa-fine"] != "unclear")]
    grp = (known.groupby("moa-fine")["index"].agg(["median", "size"])
           .query("size >= 3").sort_values("median", ascending=False))
    print(f"\nmechanism ranking in LINCS ({len(known)} annotated compounds):")
    print(grp.round(3).to_string())
    if len(grp) > 2:
        h, p = stats.kruskal(*[known[known["moa-fine"] == m]["index"].to_numpy()
                               for m in grp.index])
        print(f"Kruskal-Wallis across {len(grp)} mechanisms: H={h:.1f}, p={p:.2e}")

    # rank agreement with Tahoe
    # the two panels use different MOA vocabularies (Tahoe title-case, Hub
    # lower-case, e.g. "Glucocorticoid receptor agonist" vs "glucocorticoid
    # receptor agonist"), so normalise before joining mechanism names
    tah = pd.read_csv(TAB / "drug_context_dependence.csv")
    tah = tah[tah["moa-fine"].notna() & (tah["moa-fine"] != "unclear")].copy()
    tah["m"] = tah["moa-fine"].str.lower().str.strip()
    tg = tah.groupby("m").cdi.median()
    gl = grp.copy(); gl.index = gl.index.str.lower().str.strip()
    both = pd.DataFrame({"lincs": gl["median"].groupby(level=0).median(),
                         "tahoe": tg}).dropna()
    if len(both) >= 5:
        rho = stats.spearmanr(both.lincs, both.tahoe)
        print(f"\nmechanism-ranking agreement LINCS vs Tahoe: "
              f"Spearman rho={rho.statistic:+.3f} (p={rho.pvalue:.3f}, "
              f"n={len(both)} mechanisms)")
        both.to_csv(TAB / "lincs_vs_tahoe_cdi.csv")
        print(both.round(3).sort_values("lincs", ascending=False).to_string())

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    n = 2 + (len(both) >= 5)
    fig, ax = plt.subplots(1, n, figsize=(5.2 * n, 4.4),
                           constrained_layout=True, squeeze=False)
    ax = ax[0]
    top = grp.head(14)
    yy = np.arange(len(top))[::-1]
    ax[0].barh(yy, top["median"], color=BLUE, height=0.65)
    ax[0].set_yticks(yy, [f"{m[:32]} (n={int(k)})"
                          for m, k in zip(top.index, top["size"])], fontsize=7.5)
    ax[0].axvline(ok["index"].median(), color="#888888", ls="--", lw=0.9)
    ax[0].set_xlabel("context-dependence index (median)")
    ax[0].set_title("A  Mechanism ranking in LINCS L1000", loc="left",
                    fontweight="bold", fontsize=10)

    ax[1].hist(ok["index"], bins=40, color=AQUA, alpha=0.85)
    ax[1].set_xlabel("context-dependence index")
    ax[1].set_ylabel("compounds")
    ax[1].set_title(f"B  {len(ok)} estimable compounds", loc="left",
                    fontweight="bold", fontsize=10)

    if len(both) >= 5:
        ax[2].scatter(both.tahoe, both.lincs, s=34, color=ORANGE,
                      edgecolors="none")
        for m, r in both.iterrows():
            ax[2].annotate(m[:16], (r.tahoe, r.lincs), fontsize=6,
                           xytext=(3, 2), textcoords="offset points")
        rho = stats.spearmanr(both.lincs, both.tahoe)
        ax[2].set_xlabel("mechanism median CDI, Tahoe-100M")
        ax[2].set_ylabel("mechanism median index, LINCS")
        ax[2].set_title(f"C  Cross-platform agreement "
                        f"(rho={rho.statistic:+.2f})", loc="left",
                        fontweight="bold", fontsize=10)

    fig.suptitle("Does the mechanism ranking replicate in an independent "
                 "platform?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "lincs_mechanism", FIG,
                    source_data={"per_compound": per,
                                 "by_mechanism": grp.reset_index(),
                                 "cross_dataset": both.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
