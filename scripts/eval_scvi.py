#!/usr/bin/env python3
"""Phase 2: evaluate the scVI latent-shift baseline on held-out conditions.

Method (scGen-style latent arithmetic, all within the 4k-HVG space):
  1. Per TRAIN (line, drug, conc, plate): latent shift dz = mean z(treated)
     - mean z(plate-matched DMSO of that line); average across plates, then
     across train lines -> dz(drug, conc).
  2. For a TEST condition (L, D, C): take L's DMSO cells (plate-matched),
     predict z' = z + dz(D, C), decode both z and z' through the scVI decoder;
     predicted delta = log1p(1e4*decoded(z')) - log1p(1e4*decoded(z)), averaged
     over control cells.
  3. Compare to the true pseudobulk delta (log1p CPM-1e4 on the same genes),
     with the SAME metrics and split as eval_baselines: r on the 2,000 most
     delta-variable genes (train-derived) and on the condition's top-100 DE
     genes. The additive-shift baseline is recomputed in the 4k space so the
     comparison is apples-to-apples.

Requires: dev_subset_hvg.h5ad, scvi_latent_dev.npz, results/checkpoints/scvi_dev.
Outputs: results/tables/scvi_eval.csv + figure bundle results/figures/scvi_vs_additive/
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from perturbmodel.evaluation import held_out_condition_triples
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures" / "02_baselines"
TAB = ROOT / "results" / "tables"
MIN_CELLS = 50
N_RESP = 2000          # "most responsive" genes for r_hvg
MAX_CTRL_CELLS = 500   # control cells decoded per test condition
BLUE, ORANGE = "#2a78d6", "#eb6834"


def pearson(a, b):
    a, b = a - a.mean(), b - b.mean()
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / den) if den > 0 else np.nan


def main():
    import anndata as ad
    import scvi

    adata = ad.read_h5ad(PROC / "dev_subset_hvg.h5ad")
    Z = np.load(PROC / "scvi_latent_dev.npz")["latent"]
    assert len(Z) == adata.shape[0], "latent/adata row mismatch"
    obs = adata.obs.reset_index(drop=True)
    obs["conc"] = obs.conc.astype(float)
    X = adata.X.tocsr()

    model = scvi.model.SCVI.load(str(ROOT / "results" / "checkpoints" / "scvi_dev"),
                                 adata=adata)
    model.module.eval()

    # ---------------- condition bookkeeping ----------------
    obs["cond4"] = obs.groupby(["cell_line_id", "drug", "conc", "plate"],
                               observed=True).ngroup()
    c4 = (obs.groupby("cond4", observed=True)
          .agg(cell_line_id=("cell_line_id", "first"), drug=("drug", "first"),
               conc=("conc", "first"), plate=("plate", "first"),
               n=("cond4", "size")))
    c4 = c4[c4.plate != "plate14"]
    rows_of = {g: np.where(obs.cond4.to_numpy() == g)[0] for g in c4.index}

    ctrl4 = c4[c4.drug == "DMSO_TF"]
    ctrl_of = {(r.cell_line_id, r.plate): g for g, r in ctrl4.iterrows()}

    # true pseudobulk log1p CPM profiles per cond4 (4k space)
    def pb_profile(rows):
        sub = X[rows]
        cpm = sub.sum(axis=0).A1
        return np.log1p(cpm / (cpm.sum() + 1e-9) * 1e4)

    # ---------------- true deltas per (L,D,C), plate-matched ----------------
    trip_rows = {}
    for g, r in c4[(c4.drug != "DMSO_TF") & (c4.n >= MIN_CELLS)].iterrows():
        cg = ctrl_of.get((r.cell_line_id, r.plate))
        if cg is None or c4.loc[cg].n < MIN_CELLS:
            continue
        key = (r.cell_line_id, r.drug, r.conc)
        trip_rows.setdefault(key, []).append((g, cg, r.n))

    triples = pd.DataFrame([k for k in trip_rows],
                           columns=["cell_line_id", "drug", "conc"])
    true_delta = {}
    for key, entries in trip_rows.items():
        num, den = 0.0, 0.0
        for g, cg, n in entries:
            num = num + n * (pb_profile(rows_of[g]) - pb_profile(rows_of[cg]))
            den += n
        true_delta[key] = num / den

    # ---------------- split (same module/seed as eval_baselines) ----------------
    test_triples = held_out_condition_triples(triples, seed=0, test_frac=0.2)
    train_keys = [k for k in trip_rows if k not in test_triples]
    test_keys = [k for k in trip_rows if k in test_triples]
    print(f"{len(train_keys)} train / {len(test_keys)} test conditions")

    D_train = np.stack([true_delta[k] for k in train_keys])
    resp = np.sort(np.argsort(D_train.var(axis=0))[-N_RESP:])

    # ---------------- additive baseline in 4k space ----------------
    add_shift = {}
    for k in train_keys:
        add_shift.setdefault(k[1:], []).append(true_delta[k])
    add_shift = {p: np.mean(v, axis=0) for p, v in add_shift.items()}

    # ---------------- scVI latent shifts from TRAIN ----------------
    dz_pool = {}
    for key in train_keys:
        num, den = 0.0, 0.0
        for g, cg, n in trip_rows[key]:
            num = num + n * (Z[rows_of[g]].mean(0) - Z[rows_of[cg]].mean(0))
            den += n
        dz_pool.setdefault(key[1:], []).append(num / den)
    dz = {p: np.mean(v, axis=0) for p, v in dz_pool.items()}

    # ---------------- decode helper ----------------
    @torch.no_grad()
    def decode(z_np):
        z = torch.as_tensor(z_np, dtype=torch.float32)
        lib = torch.full((len(z), 1), float(np.log(1e4)))
        out = model.module.generative(z=z, library=lib,
                                      batch_index=torch.zeros(len(z), 1,
                                                              dtype=torch.long))
        px = out["px"]
        scale = px.scale if hasattr(px, "scale") else px["scale"]
        return scale.numpy()  # per-cell normalized expression, sums to 1

    rng = np.random.default_rng(0)
    recs = []
    for key in test_keys:
        L, D, C = key
        if (D, C) not in dz:
            continue
        # control cells of this line (any train plate), capped
        ctrl_rows = np.concatenate(
            [rows_of[g] for (lc, pl), g in ctrl_of.items() if lc == L])
        if len(ctrl_rows) > MAX_CTRL_CELLS:
            ctrl_rows = rng.choice(ctrl_rows, MAX_CTRL_CELLS, replace=False)
        z0 = Z[ctrl_rows]
        base = decode(z0).mean(0)
        pert = decode(z0 + dz[(D, C)]).mean(0)
        pred_scvi = np.log1p(1e4 * pert) - np.log1p(1e4 * base)

        t = true_delta[key]
        de100 = np.argsort(np.abs(t))[-100:]
        for name, pred in (("scvi_latent_shift", pred_scvi),
                           ("additive_shift", add_shift[(D, C)])):
            recs.append({
                "cell_line_id": L, "drug": D, "conc": C, "model": name,
                "r_hvg": pearson(pred[resp], t[resp]),
                "r_de100": pearson(pred[de100], t[de100]),
                "rmse_hvg": float(np.sqrt(np.mean((pred[resp] - t[resp]) ** 2))),
            })
    res = pd.DataFrame(recs)
    res.to_csv(TAB / "scvi_eval.csv", index=False)
    print(res.groupby("model")[["r_hvg", "r_de100", "rmse_hvg"]]
          .agg(["mean", "median"]).round(3).to_string())

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for ax, metric, title in ((axes[0], "r_hvg", "A  r on 2,000 most responsive genes"),
                              (axes[1], "r_de100", "B  r on condition's top-100 DE genes")):
        data = [res[res.model == m][metric].dropna()
                for m in ("additive_shift", "scvi_latent_shift")]
        vp = ax.violinplot(data, positions=[0, 1], widths=0.7, showmedians=True,
                           showextrema=False)
        for body, c in zip(vp["bodies"], (ORANGE, BLUE)):
            body.set_facecolor(c); body.set_alpha(0.55); body.set_edgecolor("none")
        vp["cmedians"].set_color("#333333")
        ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
        ax.set_xticks([0, 1], ["additive shift", "scVI latent shift"])
        ax.set_ylabel("Pearson r (predicted vs true delta)")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("scVI latent-shift vs additive baseline, held-out conditions",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "scvi_vs_additive", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
