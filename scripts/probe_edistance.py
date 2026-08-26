#!/usr/bin/env python3
"""Phase 1 sanity probe: E-distance dose trend on the dev subset (paper Fig 4B-C).

Per cell line: sample <=CAP cells per (drug, conc, plate), featurize as log1p CPM
over the line's top-2000 HVGs, PCA to 10-d, then E-distance of each treated
population vs the same-plate DMSO_TF population. Median across plates then lines.

Expected: E-distance increases with dose; proteasome/HDAC/PI3K-AKT inhibitors and
homoharringtonine/dinaciclib among the strongest.

Outputs: results/tables/edistance_dev.csv + figure bundle results/figures/edistance_dose_moa/
Pass --replot to regenerate the figure from the saved table (skips the long sampling).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
import yaml

from perturbmodel.evaluation import e_distance
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "interim" / "dev_subset"
META = ROOT / "data" / "metadata" / "metadata"
FIG = ROOT / "results" / "figures" / "01_paper_replication"
TAB = ROOT / "results" / "tables"
CAP = 300          # cells per (line,drug,conc,plate) group
MIN_CELLS = 50
N_HVG = 2000
N_PC = 10
BLUE = "#2a78d6"


def compute_table() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    lines = [c["cellosaurus"] for c in
             yaml.safe_load(open(SRC / "selection.yaml"))["cell_lines"]]
    gm = pd.read_parquet(META / "gene_metadata.parquet").sort_values("token_id")
    lut = np.full(int(gm.token_id.max()) + 1, -1, dtype=np.int64)
    lut[gm.token_id.to_numpy()] = np.arange(len(gm))
    n_genes = len(gm)

    # ---------- pass 1: reservoir-sample cells per group ----------
    store: dict[tuple, list] = {}   # key -> list of (genes_idx, counts)
    seen: dict[tuple, int] = {}
    shards = sorted((SRC / "data").glob("train-*.parquet"))
    for si, shard in enumerate(shards, 1):
        df = pq.read_table(shard).to_pandas()
        for r in df.itertuples():
            g, e = np.asarray(r.genes), np.asarray(r.expressions, dtype=np.float32)
            if e[0] < 0:
                g, e = g[1:], e[1:]
            key = (r.cell_line_id, r.drug, r.conc, r.plate)
            seen[key] = seen.get(key, 0) + 1
            bucket = store.setdefault(key, [])
            if len(bucket) < CAP:
                bucket.append((lut[g], e))
            else:  # reservoir replacement keeps an unbiased sample
                j = rng.integers(0, seen[key])
                if j < CAP:
                    bucket[j] = (lut[g], e)
        if si % 200 == 0 or si == len(shards):
            print(f"[sample {si}/{len(shards)}] {len(store)} groups", flush=True)

    # ---------- pass 2: per line -> HVG/PCA -> E-distance ----------
    rows = []
    for line in lines:
        keys = [k for k in store if k[0] == line and seen[k] >= MIN_CELLS]
        if not keys:
            continue
        cells, owners = [], []
        for k in keys:
            for gi, e in store[k]:
                cells.append((gi, e)); owners.append(k)
        M = np.zeros((len(cells), n_genes), dtype=np.float32)
        for i, (gi, e) in enumerate(cells):
            ok = gi >= 0
            M[i, gi[ok]] = e[ok]
        M = np.log1p(M / (M.sum(axis=1, keepdims=True) + 1e-9) * 1e4)
        hvg = np.argsort(M.var(axis=0))[-N_HVG:]
        Mh = M[:, hvg]
        Mh -= Mh.mean(axis=0, keepdims=True)
        U, S, _ = np.linalg.svd(Mh, full_matrices=False)
        Z = torch.from_numpy((U[:, :N_PC] * S[:N_PC]).astype(np.float32))
        owners = np.array(["\x1f".join(map(str, k)) for k in owners])
        by_group = {k: Z[owners == "\x1f".join(map(str, k))] for k in keys}
        for k in keys:
            cl, drug, conc, plate = k
            if drug == "DMSO_TF":
                continue
            ctrl = by_group.get((cl, "DMSO_TF", 0.0, plate))
            if ctrl is None or len(ctrl) < MIN_CELLS:
                continue
            rows.append((cl, drug, float(conc), plate, seen[k],
                         e_distance(by_group[k], ctrl)))
        print(f"{line}: {sum(1 for r in rows if r[0] == line)} conditions",
              flush=True)
        del M, Mh

    ed = pd.DataFrame(rows, columns=["cell_line_id", "drug", "conc", "plate",
                                     "n_cells", "e_distance"])
    dr = pd.read_parquet(META / "drug_metadata.parquet")[["drug", "moa-fine"]]
    ed = ed.merge(dr, on="drug", how="left")
    ed.to_csv(TAB / "edistance_dev.csv", index=False)
    print(f"table -> {TAB / 'edistance_dev.csv'} ({len(ed)} rows)")
    return ed


def make_figure(ed: pd.DataFrame):
    ed = ed.copy()
    ed["conc"] = ed.conc.round(3)   # float32 storage artifacts -> clean doses
    agg = (ed.groupby(["drug", "moa-fine", "conc", "cell_line_id"], observed=True)
           .e_distance.median().groupby(level=[0, 1, 2]).median()
           .rename("med_e").reset_index())

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)

    doses = sorted(agg.conc.unique())
    data = [np.log1p(agg[agg.conc == d].med_e.clip(lower=0)) for d in doses]
    vp = axes[0].violinplot(data, positions=range(len(doses)), widths=0.7,
                            showmedians=True, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(BLUE); b.set_alpha(0.55); b.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[0].set_xticks(range(len(doses)), [f"{d:g} uM" for d in doses])
    axes[0].set_ylabel("log1p E-distance (median over lines)")
    axes[0].set_title("A  E-distance rises with dose", loc="left",
                      fontweight="bold")

    top = (agg[agg.conc == max(doses)].groupby("moa-fine", observed=True)
           .med_e.median().sort_values(ascending=False).head(12))
    axes[1].barh(range(len(top))[::-1], np.log1p(top.values), height=0.62,
                 color=BLUE)
    axes[1].set_yticks(range(len(top))[::-1], top.index, fontsize=8)
    axes[1].set_xlabel("log1p median E-distance at top dose")
    axes[1].set_title("B  Strongest MOAs at 5 uM", loc="left", fontweight="bold")

    fig.suptitle("E-distance probes, dev subset (cf. paper Fig 4B-C)",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "edistance_dose_moa", FIG,
                    source_data={"edistance_per_condition": ed,
                                 "edistance_aggregated": agg},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--replot", action="store_true",
                    help="regenerate figure from results/tables/edistance_dev.csv")
    args = ap.parse_args()
    table = (pd.read_csv(TAB / "edistance_dev.csv") if args.replot
             else compute_table())
    make_figure(table)
