#!/usr/bin/env python3
"""Phase 3 model 1: train + evaluate the context-residual delta model.

pred(L,D,C) = additive_prior(D,C) + residual(line_emb, ECFP, dose), residual
zero-initialized and trained only on the 2,000 most responsive genes (all other
genes pass the prior through unchanged, so the model can never be worse than
additive outside the responsive set).

Anti-leakage: TRAIN rows use leave-one-line-out priors (their own delta never
appears in their input); TEST rows use the standard all-train-lines prior,
identical to the additive baseline.

Variants trained: 'full' and 'no_line' (context ablation).
Outputs: results/tables/phase3_delta_eval.csv,
         results/figures/phase3_delta_eval/ (bundle),
         results/checkpoints/phase3_delta/{full,no_line}.pt
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from perturbmodel.evaluation import held_out_condition_triples
from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                delta_metrics, load_pseudobulk,
                                                responsive_genes)
from perturbmodel.models.context_residual import ContextResidualDelta
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
CKPT = ROOT / "results" / "checkpoints" / "phase3_delta"
COLORS = {"additive": "#eb6834", "no_line": "#eda100", "full": "#2a78d6"}
SEED = 0

torch.manual_seed(SEED)
np.random.seed(SEED)

# ---------------- data ----------------
X, cond = load_pseudobulk(ROOT / "data" / "processed" / "pseudobulk_dev")
G, DELTA = build_deltas(X, cond)
test_triples = held_out_condition_triples(G, seed=0, test_frac=0.2)
is_test = np.array([(l, d, c) in test_triples
                    for l, d, c in zip(G.cell_line_id, G.drug, G.conc)])
train_mask = ~is_test
resp = responsive_genes(DELTA, train_mask)
PRIOR = additive_prior(G, DELTA, train_mask, loo=True)  # LOO on train rows only
print(f"{train_mask.sum()} train / {is_test.sum()} test triples; "
      f"{len(resp)} responsive genes")

lines = sorted(G.cell_line_id.unique())
line_idx = np.array([lines.index(l) for l in G.cell_line_id])
ecfp = np.load(ROOT / "data" / "processed" / "drug_ecfp.npz", allow_pickle=True)
fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))
missing = sorted(set(G.drug) - set(fp_of))
assert not missing, f"drugs without fingerprints: {missing}"
FP = np.stack([fp_of[d] for d in G.drug]).astype(np.float32)
logdose = np.log10(G.conc.to_numpy(dtype=np.float64))
mu, sd = logdose[train_mask].mean(), logdose[train_mask].std()
LOGDOSE = ((logdose - mu) / sd).astype(np.float32)[:, None]

# fit/val split of TRAIN triples for early stopping
rng = np.random.default_rng(1)
val_of_train = rng.random(len(G)) < 0.1
fit_mask = train_mask & ~val_of_train
val_mask = train_mask & val_of_train

T = lambda a: torch.as_tensor(a)
Y = T(DELTA[:, resp])           # residual target space
P = T(PRIOR[:, resp])
li, fp_t, ld = T(line_idx), T(FP), T(LOGDOSE)


def train_variant(use_line: bool, tag: str):
    model = ContextResidualDelta(n_lines=len(lines), n_genes=len(resp),
                                 use_line=use_line)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    fit_idx = np.where(fit_mask)[0]
    best_val, best_state, patience = np.inf, None, 0
    for epoch in range(3000):
        model.train()
        perm = np.random.permutation(fit_idx)
        for s in range(0, len(perm), 128):
            b = perm[s:s + 128]
            opt.zero_grad()
            pred = model(li[b], fp_t[b], ld[b], P[b])
            loss = torch.mean((pred - Y[b]) ** 2)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vb = np.where(val_mask)[0]
            vloss = float(torch.mean(
                (model(li[vb], fp_t[vb], ld[vb], P[vb]) - Y[vb]) ** 2))
        if vloss < best_val - 1e-6:
            best_val, patience = vloss, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            patience += 1
            if patience >= 200:
                break
    model.load_state_dict(best_state)
    CKPT.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": best_state, "resp": resp, "lines": lines,
                "dose_mu": mu, "dose_sd": sd, "use_line": use_line},
               CKPT / f"{tag}.pt")
    print(f"[{tag}] stopped epoch {epoch}, val MSE {best_val:.5f}")
    return model


results = {"additive": None}
results["no_line"] = train_variant(False, "no_line")
results["full"] = train_variant(True, "full")

# ---------------- evaluate on test triples ----------------
recs = []
test_idx = np.where(is_test)[0]
with torch.no_grad():
    preds_resp = {"additive": P.numpy()}
    for tag in ("no_line", "full"):
        m = results[tag]
        m.eval()
        preds_resp[tag] = m(li, fp_t, ld, P).numpy()
for i in test_idx:
    true = DELTA[i]
    for tag, pr in preds_resp.items():
        pred = PRIOR[i].copy()
        pred[resp] = pr[i]
        recs.append({"cell_line_id": G.cell_line_id[i], "drug": G.drug[i],
                     "conc": G.conc[i], "model": tag,
                     **delta_metrics(pred, true, resp)})
res = pd.DataFrame(recs)
res.to_csv(TAB / "phase3_delta_eval.csv", index=False)
order = ["additive", "no_line", "full"]
print(res.groupby("model")[["r_hvg", "r_de100", "rmse_hvg"]]
      .agg(["mean", "median"]).round(3).reindex(order).to_string())

# paired improvement over additive per condition
piv = res.pivot_table(index=["cell_line_id", "drug", "conc"], columns="model",
                      values="r_de100", observed=True)
print(f"\nfull beats additive on {int((piv['full'] > piv['additive']).sum())}"
      f"/{len(piv)} test conditions "
      f"(median paired gain {float((piv['full'] - piv['additive']).median()):+.4f})")

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.facecolor": "white"})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
labels = {"additive": "additive\nshift", "no_line": "residual\n(no context)",
          "full": "residual\n(+line context)"}
for ax, metric, title in ((axes[0], "r_hvg", "A  r on 2,000 most responsive genes"),
                          (axes[1], "r_de100", "B  r on condition's top-100 DE genes")):
    data = [res[res.model == m][metric].dropna() for m in order]
    vp = ax.violinplot(data, positions=range(len(order)), widths=0.7,
                       showmedians=True, showextrema=False)
    for body, m in zip(vp["bodies"], order):
        body.set_facecolor(COLORS[m]); body.set_alpha(0.6)
        body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    ax.set_xticks(range(len(order)), [labels[m] for m in order])
    ax.set_ylabel("Pearson r (predicted vs true delta)")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
fig.suptitle("Phase 3: context-residual model vs additive baseline "
             "(held-out conditions)", fontsize=11, x=0.01, ha="left")
d = save_figure(fig, "phase3_delta_eval", FIG, source_data=res, script=__file__)
print(f"figure bundle -> {d}")
