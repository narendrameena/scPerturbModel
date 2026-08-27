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
FIG = ROOT / "results" / "figures" / "04_generalization"
TAB = ROOT / "results" / "tables"
COLORS = {"additive": "#eb6834", "covariate": "#2a78d6",
          "expr_context": "#1baf7a",
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
    ap.add_argument("--context", choices=["covariate", "expression"],
                    default="covariate",
                    help="line-mode context: static covariates, or the line's "
                         "own DMSO control transcriptome (per-fold PCA, fit on "
                         "training lines only; defined for unseen lines)")
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--replot", action="store_true",
                    help="regenerate figure from the saved eval CSV")
    args = ap.parse_args()
    suf = f"_{args.tag}" if args.tag else ""
    torch.manual_seed(SEED)

    if args.replot:
        res = pd.read_csv(TAB / f"hard_splits_eval{suf}.csv")
        make_figure(res, suf)
        return

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
        # control (DMSO) log1p-CPM profile per line, plate-weighted
        dmso = cond[(cond.drug == "DMSO_TF") & (cond.plate != "plate14")]
        ctrl_prof = np.stack([
            np.average(X[grp.row.to_numpy()], axis=0, weights=grp.n_cells)
            for _, grp in dmso.groupby("cell_line_id", observed=True)])
        ctrl_lines = sorted(dmso.cell_line_id.unique())
        assert ctrl_lines == lines, "control profiles must cover all lines"
        row_line = np.array([lines.index(l) for l in G.cell_line_id])
        model_name = "covariate" if args.context == "covariate" else "expr_context"

        for ln in lines:
            test_mask = (G.cell_line_id == ln).to_numpy()
            train_mask = ~test_mask
            PRIOR = additive_prior(G, DELTA, train_mask, loo=True)
            resp = responsive_genes(DELTA, train_mask)
            mu, sd = logdose[train_mask].mean(), logdose[train_mask].std()
            LD = T(((logdose - mu) / sd).astype(np.float32)[:, None])
            P_r, Y_r = T(PRIOR[:, resp]), T(DELTA[:, resp])
            val = rng.random(len(G)) < 0.1

            if args.context == "expression":
                # per-fold PCA of control profiles, TRAIN lines only
                tr_l = [i for i, l in enumerate(lines) if l != ln]
                gvar = ctrl_prof[tr_l].var(axis=0)
                gsel = np.argsort(gvar)[-5000:]
                M = ctrl_prof[:, gsel]
                mu_c = M[tr_l].mean(axis=0, keepdims=True)
                _, _, Vt = np.linalg.svd(M[tr_l] - mu_c, full_matrices=False)
                k = min(args.n_pcs, Vt.shape[0])
                F = (M - mu_c) @ Vt[:k].T
                F = F / (F[tr_l].std(axis=0, keepdims=True) + 1e-8)
                # held-out line can project far outside the training range;
                # unclipped, the ReLU trunk extrapolates it into huge residuals
                F = np.clip(F, -4.0, 4.0)
                CTX = T(F[row_line].astype(np.float32))
            else:
                CTX = LF

            model = fit(CovariateResidualDelta(CTX.shape[1], len(resp)),
                        CTX, FP, LD, P_r, Y_r,
                        np.where(train_mask & ~val)[0],
                        np.where(train_mask & val)[0])
            with torch.no_grad():
                pr = model(CTX, FP, LD, P_r).numpy()
            pred_ctx = PRIOR.copy()
            pred_ctx[:, resp] = pr
            eval_rows(np.where(test_mask)[0],
                      {"additive": PRIOR, model_name: pred_ctx},
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

    make_figure(res, suf)


def make_figure(res: pd.DataFrame, suf: str):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    nfl = res.loc[res["mode"] == "line", "fold"].nunique()
    nfd = res.loc[res["mode"] == "drug", "fold"].nunique()
    lbl = {"additive": f"additive shift\n({max(nfl - 1, 0)} other lines)",
           "covariate": "+ mutation/organ\ncontext",
           "expr_context": "+ control-expression\ncontext (PCA)",
           "nn_drug": "nearest drug\nby ECFP",
           "ecfp_model": "ECFP model\n(prior-dropout)"}
    line_models = [m for m in ("additive", "covariate", "expr_context")
                   if ((res["mode"] == "line") & (res.model == m)).any()]
    drug_models = [m for m in ("nn_drug", "ecfp_model")
                   if ((res["mode"] == "drug") & (res.model == m)).any()]
    panels = [p for p in
              [("line", line_models, f"A  Held-out cell line ({nfl}-fold)",
                [lbl[m] for m in line_models]),
               ("drug", drug_models, f"B  Held-out drug ({nfd}-fold)",
                [lbl[m] for m in drug_models])] if p[1]]
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 4.2),
                             constrained_layout=True, squeeze=False)
    for ax, (mode, models, title, labels) in zip(axes[0], panels):
        sub = res[res["mode"] == mode]
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
    fig.suptitle("Phase 3: hard generalization splits"
                 f"{' (dev47: 47 lines)' if suf else ' (dev subset)'}",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"hard_splits_eval{suf}", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
