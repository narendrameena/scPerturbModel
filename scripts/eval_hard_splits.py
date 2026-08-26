#!/usr/bin/env python3
"""Phase 3: hard generalization splits for the residual delta models.

MODE line — leave-one-LINE-out (8 folds): can a line's covariates (driver
  mutations + organ) supply context for a NEVER-SEEN line? Compare
  additive-shift (mean of the other 7 lines) vs CovariateResidualDelta.

MODE drug — drug-grouped 5-fold CV: all conditions of held-out DRUGS are
  unseen, so the additive prior is zero and prediction must come from chemistry
  (ECFP) + context. ContextResidualDelta is trained with prior-dropout (30% of
  train rows see a zero prior) so it learns to produce full deltas without a
  prior. Baseline: nearest training drug by ECFP cosine, same dose.

Outputs: results/tables/hard_splits_eval.csv + bundle results/figures/hard_splits_eval/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                delta_metrics, load_pseudobulk,
                                                responsive_genes)
from perturbmodel.models.context_residual import (ContextResidualDelta,
                                                  CovariateResidualDelta)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
COLORS = {"additive": "#eb6834", "covariate": "#2a78d6",
          "nn_drug": "#eda100", "ecfp_model": "#2a78d6"}
SEED = 0
T = lambda a: torch.as_tensor(a)


def fit(model, ctx, fp, ld, P, Y, fit_idx, val_idx, prior_dropout=0.0,
        max_epochs=3000, patience=200):
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    rng = np.random.default_rng(SEED)
    best, best_state, wait = np.inf, None, 0
    for _ in range(max_epochs):
        model.train()
        perm = rng.permutation(fit_idx)
        for s in range(0, len(perm), 128):
            b = perm[s:s + 128]
            Pb = P[b]
            if prior_dropout > 0:
                mask = T((rng.random(len(b)) > prior_dropout)
                         .astype(np.float32))[:, None]
                Pb = Pb * mask
            opt.zero_grad()
            loss = torch.mean((model(ctx[b], fp[b], ld[b], Pb) - Y[b]) ** 2)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = float(torch.mean(
                (model(ctx[val_idx], fp[val_idx], ld[val_idx], P[val_idx])
                 - Y[val_idx]) ** 2))
        if vl < best - 1e-6:
            best, wait = vl, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            wait += 1
            if wait >= patience:
                break
    model.load_state_dict(best_state)
    return model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["line", "drug", "both"], default="both")
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev")
    ap.add_argument("--tag", default="", help="suffix for outputs, e.g. dev47")
    args = ap.parse_args()
    suf = f"_{args.tag}" if args.tag else ""
    torch.manual_seed(SEED)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    lf = np.load(ROOT / "data/processed/line_features.npz", allow_pickle=True)
    line_feat_of = dict(zip(lf["lines"].tolist(), lf["feat"]))
    ecfp = np.load(ROOT / "data/processed/drug_ecfp.npz", allow_pickle=True)
    fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))

    LF = T(np.stack([line_feat_of[l] for l in G.cell_line_id]).astype(np.float32))
    FP_np = np.stack([fp_of[d] for d in G.drug]).astype(np.float32)
    FP = T(FP_np)
    logdose = np.log10(G.conc.to_numpy(dtype=np.float64))
    lines = sorted(G.cell_line_id.unique())
    line_idx = T(np.array([lines.index(l) for l in G.cell_line_id]))
    rng = np.random.default_rng(1)
    recs = []

    def eval_rows(test_idx, preds, resp, mode, fold):
        for i in test_idx:
            for name, P_full in preds.items():
                recs.append({"mode": mode, "fold": fold, "model": name,
                             "cell_line_id": G.cell_line_id[i],
                             "drug": G.drug[i], "conc": G.conc[i],
                             **delta_metrics(P_full[i], DELTA[i], resp)})

    # ---------------- MODE line ----------------
    if args.mode in ("line", "both"):
        for ln in lines:
            test_mask = (G.cell_line_id == ln).to_numpy()
            train_mask = ~test_mask
            PRIOR = additive_prior(G, DELTA, train_mask, loo=True)
            resp = responsive_genes(DELTA, train_mask)
            mu, sd = logdose[train_mask].mean(), logdose[train_mask].std()
            LD = T(((logdose - mu) / sd).astype(np.float32)[:, None])
            P_r, Y_r = T(PRIOR[:, resp]), T(DELTA[:, resp])
            val = rng.random(len(G)) < 0.1
            model = fit(CovariateResidualDelta(LF.shape[1], len(resp)),
                        LF, FP, LD, P_r, Y_r,
                        np.where(train_mask & ~val)[0],
                        np.where(train_mask & val)[0])
            with torch.no_grad():
                pr = model(LF, FP, LD, P_r).numpy()
            pred_cov = PRIOR.copy()
            pred_cov[:, resp] = pr
            eval_rows(np.where(test_mask)[0],
                      {"additive": PRIOR, "covariate": pred_cov},
                      resp, "line", ln)
            print(f"[line fold {ln}] done", flush=True)

    # ---------------- MODE drug ----------------
    if args.mode in ("drug", "both"):
        drugs = sorted(G.drug.unique())
        rng2 = np.random.default_rng(2)
        order = rng2.permutation(len(drugs))
        folds = [sorted(np.array(drugs)[order[k::5]]) for k in range(5)]
        fpn = FP_np / (np.linalg.norm(FP_np, axis=1, keepdims=True) + 1e-9)
        for k, held in enumerate(folds):
            test_mask = G.drug.isin(held).to_numpy()
            train_mask = ~test_mask
            PRIOR = additive_prior(G, DELTA, train_mask, loo=True)  # 0 for test
            resp = responsive_genes(DELTA, train_mask)
            mu, sd = logdose[train_mask].mean(), logdose[train_mask].std()
            LD = T(((logdose - mu) / sd).astype(np.float32)[:, None])
            P_r, Y_r = T(PRIOR[:, resp]), T(DELTA[:, resp])
            val = rng.random(len(G)) < 0.1
            model = fit(ContextResidualDelta(len(lines), len(resp)),
                        line_idx, FP, LD, P_r, Y_r,
                        np.where(train_mask & ~val)[0],
                        np.where(train_mask & val)[0], prior_dropout=0.3)
            with torch.no_grad():
                pr = model(line_idx, FP, LD, P_r).numpy()
            pred_m = PRIOR.copy()
            pred_m[:, resp] = pr

            # nearest-train-drug baseline (ECFP cosine, same conc)
            tr_drugs = sorted(set(G.drug[train_mask]))
            tr_fp = np.stack([fpn[G.drug.tolist().index(d)] for d in tr_drugs])
            shift = {}
            for i in np.where(train_mask)[0]:
                shift.setdefault((G.drug[i], G.conc[i]), []).append(DELTA[i])
            shift = {kk: np.mean(v, axis=0) for kk, v in shift.items()}
            pred_nn = np.zeros_like(DELTA)
            for i in np.where(test_mask)[0]:
                sims = tr_fp @ fpn[i]
                for j in np.argsort(sims)[::-1]:
                    kk = (tr_drugs[j], G.conc[i])
                    if kk in shift:
                        pred_nn[i] = shift[kk]
                        break
            eval_rows(np.where(test_mask)[0],
                      {"nn_drug": pred_nn, "ecfp_model": pred_m},
                      resp, "drug", f"fold{k}")
            print(f"[drug fold {k}] {len(held)} drugs held out", flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(TAB / f"hard_splits_eval{suf}.csv", index=False)
    for mode in res["mode"].unique():
        print(f"\n== {mode} ==")
        print(res[res["mode"] == mode].groupby("model")
              [["r_hvg", "r_de100", "rmse_hvg"]]
              .agg(["mean", "median"]).round(3).to_string())

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    panels = [("line", ["additive", "covariate"],
               "A  Held-out cell line (8-fold)",
               ["additive shift\n(7 other lines)", "+ mutation/organ\ncontext"]),
              ("drug", ["nn_drug", "ecfp_model"],
               "B  Held-out drug (5-fold)",
               ["nearest drug\nby ECFP", "ECFP model\n(prior-dropout)"])]
    for ax, (mode, models, title, labels) in zip(axes, panels):
        sub = res[res["mode"] == mode]
        if not len(sub):
            continue
        data = [sub[sub.model == m].r_de100.dropna() for m in models]
        vp = ax.violinplot(data, positions=range(len(models)), widths=0.7,
                           showmedians=True, showextrema=False)
        for body, m in zip(vp["bodies"], models):
            body.set_facecolor(COLORS[m]); body.set_alpha(0.6)
            body.set_edgecolor("none")
        vp["cmedians"].set_color("#333333")
        ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
        ax.set_xticks(range(len(models)), labels)
        ax.set_ylabel("Pearson r, top-100 DE genes")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Phase 3: hard generalization splits (dev subset)",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"hard_splits_eval{suf}", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
