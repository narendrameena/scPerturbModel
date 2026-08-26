#!/usr/bin/env python3
"""Phase 3 model 2 evaluation: does the conditional VAE generate the right
cell populations for held-out conditions?

Per test (line, drug, conc) triple (N_CELLS real/generated cells each):
  population level (10-d PCA per triple, fit on real cells, log1p CPM-1e4):
    e_noise     E(real_treated half1, half2)        noise floor
    e_gen_real  E(generated_treated, real_treated)  model quality (lower better)
    e_ctrl_real E(real_control, real_treated)       'no-change' reference
  pseudobulk level (4k HVG space):
    pred_delta = pb(gen_treated) - pb(gen_control)  vs
    true_delta = pb(real_treated) - pb(real_control)
    -> r_de100 / r_hvg (responsive genes = model-1 resp set restricted to HVGs)

Outputs: results/tables/cvae_eval.csv + figure bundle results/figures/cvae_eval/
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from perturbmodel.evaluation import e_distance, held_out_condition_triples
from perturbmodel.evaluation.delta_eval import pearson
from perturbmodel.models.cond_vae import CondVAE
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures" / "03_models"
TAB = ROOT / "results" / "tables"
N_CELLS = 300
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RNG = np.random.default_rng(0)


def lognorm(dense: np.ndarray) -> np.ndarray:
    lib = dense.sum(1, keepdims=True)
    return np.log1p(dense / (lib + 1e-8) * 1e4)


def pca_embed(reference: np.ndarray, others: list[np.ndarray], k: int = 10):
    mu = reference.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(reference - mu, full_matrices=False)
    P = Vt[:k].T
    return [(reference - mu) @ P] + [(o - mu) @ P for o in others]


def main():
    import anndata as ad

    ck = torch.load(ROOT / "results/checkpoints/cond_vae/model.pt",
                    weights_only=False)
    model = CondVAE(n_genes=ck["n_genes"], n_lines=len(ck["lines"]),
                    fp_dim=ck["fp_dim"], n_latent=ck["n_latent"])
    model.load_state_dict(ck["state_dict"])
    model.eval()
    lines, mu_d, sd_d = ck["lines"], ck["dose_mu"], ck["dose_sd"]

    ecfp = np.load(PROC / "drug_ecfp.npz", allow_pickle=True)
    fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))

    adata = ad.read_h5ad(PROC / "dev_subset_hvg.h5ad", backed="r")
    obs = adata.obs.reset_index(drop=True)
    obs["conc"] = obs.conc.astype(float)
    keep = ~obs.is_replicate_plate

    triples = obs.loc[keep, ["cell_line_id", "drug", "conc"]].drop_duplicates()
    test_triples = held_out_condition_triples(triples, seed=0, test_frac=0.2)
    print(f"{len(test_triples)} test triples")

    # rows needed: test-triple cells + per-line controls
    sub = obs[keep.to_numpy()]
    grp_rows = {k: g.index.to_numpy() for k, g in
                sub.groupby(["cell_line_id", "drug", "conc"], observed=True)}
    want_rows, ctrl_rows_line = [], {}
    for t in test_triples:
        rows = grp_rows.get(t, np.array([], dtype=int))
        if len(rows) > N_CELLS * 2:
            rows = RNG.choice(rows, N_CELLS * 2, replace=False)
        want_rows.append((t, np.sort(rows)))
    for ln in lines:
        rows = np.where(keep.to_numpy() & (obs.cell_line_id == ln).to_numpy()
                        & (obs.drug == "DMSO_TF").to_numpy())[0]
        ctrl_rows_line[ln] = np.sort(RNG.choice(rows, min(len(rows), 1000),
                                                replace=False))
    all_rows = np.unique(np.concatenate(
        [r for _, r in want_rows] + list(ctrl_rows_line.values())))
    Xall = adata[all_rows].to_memory().X.toarray().astype(np.float32)
    pos = {r: i for i, r in enumerate(all_rows)}
    print(f"loaded {len(all_rows):,} cells")

    # responsive genes: model-1 resp (full-gene positions) restricted to HVGs
    m1 = torch.load(ROOT / "results/checkpoints/phase3_delta/full.pt",
                    weights_only=False)
    gm = pd.read_parquet(
        ROOT / "data/metadata/metadata/gene_metadata.parquet"
    ).sort_values("token_id").reset_index(drop=True)
    hvg = pd.read_csv(PROC / "hvg_genes.csv")
    tok_to_pos = pd.Series(np.arange(len(gm)), index=gm.token_id)
    hvg_full_pos = tok_to_pos[hvg.token_id].to_numpy()
    resp4k = np.where(np.isin(hvg_full_pos, m1["resp"]))[0]
    print(f"{len(resp4k)} responsive genes within the 4k HVG space")

    recs = []
    for (ln, drug, conc), rows in want_rows:
        if drug not in fp_of or len(rows) < 100:
            continue
        real = Xall[[pos[r] for r in rows]]
        ctrl = Xall[[pos[r] for r in ctrl_rows_line[ln]]]
        n = min(N_CELLS, len(real) // 2, len(ctrl))
        RNG.shuffle(real)
        half1, half2 = real[:n], real[len(real) - n:]
        ctrl_s = ctrl[RNG.choice(len(ctrl), n, replace=False)]

        li = torch.tensor(lines.index(ln))
        fp = torch.as_tensor(fp_of[drug])
        dose_t = torch.tensor([(np.log10(conc) - mu_d) / sd_d, 0.0],
                              dtype=torch.float32)
        dose_c = torch.tensor([0.0, 1.0], dtype=torch.float32)
        libs = torch.as_tensor(ctrl_s.sum(1))
        with torch.no_grad():
            gen_t = model.generate(n, li, fp, dose_t, libs).numpy()
            gen_c = model.generate(n, li, torch.zeros_like(fp), dose_c,
                                   libs).numpy()

        ln_real1, ln_real2, ln_ctrl, ln_gen = (lognorm(half1), lognorm(half2),
                                               lognorm(ctrl_s), lognorm(gen_t))
        ref = np.concatenate([ln_real1, ln_real2, ln_ctrl])
        emb = pca_embed(ref, [ln_gen])
        r1 = emb[0][:n]
        r2 = emb[0][n:2 * n]
        rc = emb[0][2 * n:]
        rg = emb[1]

        true_delta = lognorm(real).mean(0) - lognorm(ctrl_s).mean(0)
        pred_delta = ln_gen.mean(0) - lognorm(gen_c).mean(0)
        de100 = np.argsort(np.abs(true_delta))[-100:]
        recs.append({
            "cell_line_id": ln, "drug": drug, "conc": conc, "n": n,
            "e_noise": e_distance(torch.as_tensor(r1), torch.as_tensor(r2)),
            "e_gen_real": e_distance(torch.as_tensor(rg), torch.as_tensor(
                np.concatenate([r1, r2]))),
            "e_ctrl_real": e_distance(torch.as_tensor(rc), torch.as_tensor(
                np.concatenate([r1, r2]))),
            "r_de100": pearson(pred_delta[de100], true_delta[de100]),
            "r_hvg": pearson(pred_delta[resp4k], true_delta[resp4k]),
        })
    res = pd.DataFrame(recs)
    res.to_csv(TAB / "cvae_eval.csv", index=False)
    med = res[["e_noise", "e_gen_real", "e_ctrl_real", "r_de100", "r_hvg"]].median()
    closer = (res.e_gen_real < res.e_ctrl_real).mean()
    print(med.round(3).to_string())
    print(f"generated population closer to real than control is: "
          f"{closer:.1%} of {len(res)} conditions")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    cols = [("e_noise", "real vs real\n(noise floor)", AQUA),
            ("e_gen_real", "generated\nvs real", BLUE),
            ("e_ctrl_real", "control vs real\n(no-change)", ORANGE)]
    data = [np.log1p(res[c].clip(lower=0)) for c, _, _ in cols]
    vp = axes[0].violinplot(data, positions=range(3), widths=0.7,
                            showmedians=True, showextrema=False)
    for body, (_, _, col) in zip(vp["bodies"], cols):
        body.set_facecolor(col); body.set_alpha(0.6); body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[0].set_xticks(range(3), [l for _, l, _ in cols])
    axes[0].set_ylabel("log1p E-distance")
    axes[0].set_title("A  Population match (held-out conditions)",
                      loc="left", fontweight="bold", fontsize=10)

    vp = axes[1].violinplot([res.r_hvg.dropna(), res.r_de100.dropna()],
                            positions=[0, 1], widths=0.7, showmedians=True,
                            showextrema=False)
    for body in vp["bodies"]:
        body.set_facecolor(BLUE); body.set_alpha(0.6); body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[1].axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
    axes[1].set_xticks([0, 1], ["r (responsive genes)", "r (top-100 DE)"])
    axes[1].set_ylabel("Pearson r (predicted vs true delta)")
    axes[1].set_title("B  Pseudobulk delta accuracy", loc="left",
                      fontweight="bold", fontsize=10)
    fig.suptitle("Phase 3 model 2: conditional VAE generation quality",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cvae_eval", FIG, source_data=res, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
