#!/usr/bin/env python3
"""Phase 3: few-shot transfer to an unseen cell line.

Motivation: our unseen-line experiments produced a triple negative — learned
embeddings are undefined for a new line, and neither static covariates
(mutations/organ) nor the line's baseline transcriptome beat the additive
baseline. The realistic deployment middle ground: measure k probe drugs on the
new line, infer its embedding from those, predict everything else.

Protocol per held-out line L (leave-one-line-out):
  1. Train the context-residual model on the other 46 lines (their embeddings).
  2. Reserve a FIXED evaluation drug set for L (never used for probing), so the
     evaluation set is identical across all k — the k-curve is apples-to-apples.
  3. For each k: draw k probe drugs from the remaining pool, take L's observed
     deltas for those (all doses) as the "measurement", and fit ONLY L's
     embedding vector by gradient descent (all other weights frozen), starting
     from the mean of the training-line embeddings with an L2 pull toward it.
  4. Predict L's evaluation conditions; compare to the additive baseline.

k=0 means the mean training embedding (no measurement) — the honest zero-shot
model. Probe strategies: 'random' (k drugs at random) and 'high_effect' (drugs
with the largest mean |delta| across training lines — does panel design matter?).

Few-shot hyperparameters are fixed a priori (no validation data exists in the
few-shot setting by definition): 300 Adam steps, lr 1e-2, L2 0.01 toward mean.

Outputs: results/tables/few_shot_eval<suf>.csv
         figure bundle results/figures/04_generalization/few_shot_eval<suf>/
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
from perturbmodel.models.context_residual import ContextResidualDelta
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "04_generalization"
TAB = ROOT / "results" / "tables"
SEED = 0
K_VALUES = [0, 1, 2, 5, 10, 20]
N_EVAL_DRUGS = 20          # fixed evaluation panel per (fold, seed)
FEWSHOT_STEPS = 300
FEWSHOT_LR = 1e-2
FEWSHOT_L2 = 0.01
COLORS = {"random": "#2a78d6", "high_effect": "#1baf7a", "additive": "#eb6834"}
T = lambda a: torch.as_tensor(a)


def fit_base(model, li, fp, ld, P, Y, fit_idx, val_idx,
             max_epochs=3000, patience=200):
    """Train the full model on training-line rows (same recipe as model 1)."""
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    rng = np.random.default_rng(SEED)
    best, best_state, wait = np.inf, None, 0
    for _ in range(max_epochs):
        model.train()
        perm = rng.permutation(fit_idx)
        for s in range(0, len(perm), 128):
            b = perm[s:s + 128]
            opt.zero_grad()
            torch.mean((model(li[b], fp[b], ld[b], P[b]) - Y[b]) ** 2).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = float(torch.mean((model(li[val_idx], fp[val_idx], ld[val_idx],
                                         P[val_idx]) - Y[val_idx]) ** 2))
        if vl < best - 1e-6:
            best, wait, best_state = vl, 0, {k: v.clone()
                                             for k, v in model.state_dict().items()}
        else:
            wait += 1
            if wait >= patience:
                break
    model.load_state_dict(best_state)
    return model


def finetune_with_probes(model, held_idx, mean_emb, probe_idx, train_idx,
                         li, fp, ld, P, Y, epochs=60, lr=3e-4, probe_frac=0.15):
    """Adaptation variant: fine-tune the WHOLE model on training rows plus the
    new line's probe rows, so the residual head can co-adapt to the new
    embedding instead of treating it as an unseen point. Probe rows are
    oversampled to `probe_frac` of each epoch or they are drowned out."""
    with torch.no_grad():
        model.line_emb.weight[held_idx] = mean_emb
    if len(probe_idx) == 0:
        return 0.0
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    rng = np.random.default_rng(SEED)
    reps = max(1, int(probe_frac * len(train_idx) / max(len(probe_idx), 1)))
    pool = np.concatenate([train_idx, np.tile(probe_idx, reps)])
    model.train()
    for _ in range(epochs):
        perm = rng.permutation(pool)
        for s in range(0, len(perm), 128):
            b = perm[s:s + 128]
            opt.zero_grad()
            torch.mean((model(li[b], fp[b], ld[b], P[b]) - Y[b]) ** 2).backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        return float(torch.norm(model.line_emb.weight[held_idx] - mean_emb))


def fit_embedding(model, held_idx, mean_emb, probe_idx, li, fp, ld, P, Y):
    """Fit ONLY the held-out line's embedding on the probe conditions.

    Returns the L2 displacement of the fitted embedding from the mean —
    a diagnostic that separates 'the fit did nothing' from 'the fit moved
    the embedding but predictions barely changed' (i.e. no headroom)."""
    with torch.no_grad():
        model.line_emb.weight[held_idx] = mean_emb
    if len(probe_idx) == 0:
        return 0.0
    for p in model.parameters():
        p.requires_grad_(False)
    model.line_emb.weight.requires_grad_(True)
    opt = torch.optim.Adam([model.line_emb.weight], lr=FEWSHOT_LR)
    model.eval()                      # keep dropout off: k is tiny
    for _ in range(FEWSHOT_STEPS):
        opt.zero_grad()
        pred = model(li[probe_idx], fp[probe_idx], ld[probe_idx], P[probe_idx])
        loss = torch.mean((pred - Y[probe_idx]) ** 2) + FEWSHOT_L2 * torch.sum(
            (model.line_emb.weight[held_idx] - mean_emb) ** 2)
        loss.backward()
        # gradient reaches only the held-out row; zero the rest defensively
        g = model.line_emb.weight.grad
        mask = torch.zeros_like(g)
        mask[held_idx] = 1.0
        g.mul_(mask)
        opt.step()
    for p in model.parameters():
        p.requires_grad_(True)
    with torch.no_grad():
        return float(torch.norm(model.line_emb.weight[held_idx] - mean_emb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_dev47")
    ap.add_argument("--tag", default="dev47")
    ap.add_argument("--folds", type=int, default=0,
                    help="pilot: only the first N lines (0 = all)")
    ap.add_argument("--seeds", type=int, default=3,
                    help="probe draws per fold")
    ap.add_argument("--adapt", choices=["embedding", "finetune"],
                    default="embedding",
                    help="'embedding': fit only the new line's embedding "
                         "(frozen head). 'finetune': retrain the whole model "
                         "with the probe rows included (head co-adapts).")
    ap.add_argument("--k-values", type=int, nargs="*", default=None,
                    help="override the k sweep (e.g. --k-values 1 5 20)")
    ap.add_argument("--strategies", nargs="*", default=["random", "high_effect"])
    ap.add_argument("--replot", action="store_true")
    args = ap.parse_args()
    k_values = args.k_values if args.k_values is not None else K_VALUES
    suf = f"_{args.tag}" if args.tag else ""
    torch.manual_seed(SEED)

    if args.replot:
        make_figure(pd.read_csv(TAB / f"few_shot_eval{suf}.csv"), suf)
        return

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    ecfp = np.load(ROOT / "data/processed/drug_ecfp.npz", allow_pickle=True)
    fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))
    FP = T(np.stack([fp_of[d] for d in G.drug]).astype(np.float32))
    logdose = np.log10(G.conc.to_numpy(dtype=np.float64))
    lines = sorted(G.cell_line_id.unique())
    line_idx_np = np.array([lines.index(l) for l in G.cell_line_id])
    li = T(line_idx_np)
    fold_lines = lines[: args.folds] if args.folds else lines
    print(f"{len(G)} conditions, {len(lines)} lines, "
          f"{len(fold_lines)} folds", flush=True)

    recs = []
    for fi, ln in enumerate(fold_lines):
        held_idx = lines.index(ln)
        test_mask = (G.cell_line_id == ln).to_numpy()
        train_mask = ~test_mask
        PRIOR = additive_prior(G, DELTA, train_mask, loo=True)
        resp = responsive_genes(DELTA, train_mask)
        mu, sd = logdose[train_mask].mean(), logdose[train_mask].std()
        ld = T(((logdose - mu) / sd).astype(np.float32)[:, None])
        P, Y = T(PRIOR[:, resp]), T(DELTA[:, resp])

        rng_v = np.random.default_rng(SEED + fi)
        val = rng_v.random(len(G)) < 0.1
        model = fit_base(ContextResidualDelta(len(lines), len(resp)),
                         li, FP, ld, P, Y,
                         np.where(train_mask & ~val)[0],
                         np.where(train_mask & val)[0])
        base_state = {k: v.clone() for k, v in model.state_dict().items()}
        with torch.no_grad():
            others = [i for i in range(len(lines)) if i != held_idx]
            mean_emb = model.line_emb.weight[others].mean(0).clone()

        # drug-level effect sizes on TRAIN lines (for the high_effect strategy)
        eff = (pd.DataFrame({"drug": G.drug[train_mask],
                             "m": np.abs(DELTA[train_mask][:, resp]).mean(1)})
               .groupby("drug").m.mean().sort_values(ascending=False))
        test_drugs = sorted(set(G.drug[test_mask]))

        for seed in range(args.seeds):
            rng = np.random.default_rng(1000 * seed + fi)
            perm = rng.permutation(test_drugs)
            eval_drugs = set(perm[:N_EVAL_DRUGS])
            pool = [d for d in perm[N_EVAL_DRUGS:]]
            eval_idx = np.where(test_mask & G.drug.isin(eval_drugs).to_numpy())[0]
            if len(eval_idx) == 0:
                continue
            # additive baseline on the same evaluation conditions
            for i in eval_idx:
                recs.append({"line": ln, "seed": seed, "k": None,
                             "strategy": "additive", "drug": G.drug[i],
                             "conc": G.conc[i],
                             **delta_metrics(PRIOR[i], DELTA[i], resp)})

            pool_by_effect = [d for d in eff.index if d in pool]
            # 'oracle' = fit the embedding on the ENTIRE probe pool: the
            # architecture's ceiling for this line, so a flat k-curve can be
            # attributed to missing headroom rather than to few-shot data.
            plan = [(s, k) for s in args.strategies for k in k_values]
            plan.append(("oracle", -1))
            for strategy, k in plan:
                if strategy == "random":
                    probes = list(rng.permutation(pool)[:k])
                elif strategy == "high_effect":
                    if k == 0:
                        continue          # identical to random k=0
                    probes = pool_by_effect[:k]
                else:
                    probes = list(pool)
                if k > 0 and len(probes) < k:
                    continue
                probe_idx = np.where(
                    test_mask & G.drug.isin(probes).to_numpy())[0] \
                    if probes else np.array([], dtype=int)
                model.load_state_dict(base_state)
                if args.adapt == "embedding":
                    shift = fit_embedding(model, held_idx, mean_emb, probe_idx,
                                          li, FP, ld, P, Y)
                else:
                    shift = finetune_with_probes(
                        model, held_idx, mean_emb, probe_idx,
                        np.where(train_mask & ~val)[0], li, FP, ld, P, Y)
                with torch.no_grad():
                    pr = model(li[eval_idx], FP[eval_idx], ld[eval_idx],
                               P[eval_idx]).numpy()
                for j, i in enumerate(eval_idx):
                    pred = PRIOR[i].copy()
                    pred[resp] = pr[j]
                    recs.append({"line": ln, "seed": seed, "k": k,
                                 "strategy": strategy, "n_probe_drugs":
                                 len(probes), "emb_shift": shift,
                                 "drug": G.drug[i], "conc": G.conc[i],
                                 **delta_metrics(pred, DELTA[i], resp)})
        print(f"[fold {fi + 1}/{len(fold_lines)} {ln}] done", flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(TAB / f"few_shot_eval{suf}.csv", index=False)
    summ = (res.groupby(["strategy", "k"], dropna=False)
            .agg(r_mean=("r_de100", "mean"), r_median=("r_de100", "median"),
                 emb_shift=("emb_shift", "median"), n=("r_de100", "size"))
            .round(4))
    print(summ.to_string())
    base = res[res.strategy == "additive"].r_de100.median()
    orac = res[res.strategy == "oracle"].r_de100.median()
    print(f"\nadapt={args.adapt}; headroom (oracle - additive): {orac - base:+.4f}")
    for k in k_values:
        s = res[(res.strategy == "random") & (res.k == k)].r_de100.median()
        frac = (s - base) / (orac - base) if orac > base else float("nan")
        print(f"  k={k:>2}: {s - base:+.4f} ({frac:.0%} of headroom)")
    make_figure(res, suf)


def make_figure(res: pd.DataFrame, suf: str):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    add = res[res.strategy == "additive"]
    base_med = add.r_de100.median()
    orac = res[res.strategy == "oracle"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    # A: k-curve. Band = IQR of PER-FOLD medians (uncertainty in the estimate),
    # not spread across conditions, which is far larger and would hide the
    # effect entirely.
    for strat in ("random", "high_effect"):
        sub = res[res.strategy == strat]
        if not len(sub):
            continue
        ks = sorted(sub.k.dropna().unique())
        fold_med = sub.groupby(["k", "line"]).r_de100.median()
        med = [fold_med[k].median() for k in ks]
        lo = [np.percentile(fold_med[k], 25) for k in ks]
        hi = [np.percentile(fold_med[k], 75) for k in ks]
        axes[0].plot(ks, med, "-o", color=COLORS[strat], lw=2, ms=6,
                     label=f"{strat.replace('_', ' ')} probes")
        axes[0].fill_between(ks, lo, hi, color=COLORS[strat], alpha=0.15,
                             lw=0)
    axes[0].axhline(base_med, color=COLORS["additive"], ls="--", lw=1.5,
                    label=f"additive baseline ({base_med:.3f})")
    if len(orac):
        axes[0].axhline(orac.r_de100.median(), color="#4a3aa7", ls=":", lw=1.5,
                        label=f"oracle: all probes ({orac.r_de100.median():.3f})")
    axes[0].set_xlabel("k probe drugs measured on the unseen line")
    axes[0].set_ylabel("Pearson r, top-100 DE genes")
    axes[0].set_title("A  Few-shot transfer curve "
                      "(band = IQR of per-fold medians)", loc="left",
                      fontweight="bold", fontsize=9.5)
    axes[0].legend(frameon=False, fontsize=8.5)

    # B: paired per-condition gain over additive, random strategy
    key = ["line", "seed", "drug", "conc"]
    a = add.set_index(key).r_de100
    ks = sorted(res[res.strategy == "random"].k.dropna().unique())
    data, labels = [], []
    for k in ks:
        s = res[(res.strategy == "random") & (res.k == k)].set_index(key).r_de100
        d = (s - a).dropna()
        data.append(d.to_numpy())
        labels.append(f"k={int(k)}")
    vp = axes[1].violinplot(data, positions=range(len(data)), widths=0.7,
                            showmedians=True, showextrema=False)
    for body in vp["bodies"]:
        body.set_facecolor(COLORS["random"]); body.set_alpha(0.55)
        body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[1].axhline(0, color="#888888", lw=0.9, ls="--", zorder=0)
    axes[1].set_xticks(range(len(labels)), labels)
    axes[1].set_ylabel("paired gain in r vs additive")
    axes[1].set_title("B  Per-condition gain (random probes)", loc="left",
                      fontweight="bold", fontsize=10)
    fig.suptitle("Few-shot transfer to an unseen cell line "
                 "(leave-one-line-out, fixed evaluation panel)",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"few_shot_eval{suf}", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
