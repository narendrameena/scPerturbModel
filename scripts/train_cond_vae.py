#!/usr/bin/env python3
"""Phase 3 model 2: train the conditional VAE on single cells (CPU/SLURM).

Data: dev_subset_hvg.h5ad, split=='train' cells, plate14 excluded, capped at
CAP cells per (line, drug, conc, plate) group so DMSO controls stay represented.
Controls train with fp=0, is_control=1. KL warmup over the first 10 epochs;
early stopping on val loss (split=='val' cells, same capping).

Output: results/checkpoints/cond_vae/model.pt (+ training curve CSV)
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from perturbmodel.models.cond_vae import CondVAE

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
CKPT = ROOT / "results" / "checkpoints" / "cond_vae"
CAP = 400
SEED = 0


def select_rows(obs: pd.DataFrame, split: str, cap: int,
                rng: np.random.Generator) -> np.ndarray:
    m = (obs.split == split) & (~obs.is_replicate_plate)
    idx = np.where(m)[0]
    sub = obs.iloc[idx]
    picked = []
    for _, grp in sub.groupby(["cell_line_id", "drug", "conc", "plate"],
                              observed=True):
        g = grp.index.to_numpy()
        picked.append(g if len(g) <= cap else rng.choice(g, cap, replace=False))
    return np.sort(np.concatenate(picked))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--n-latent", type=int, default=16)
    ap.add_argument("--kl-warmup", type=int, default=10)
    ap.add_argument("--patience", type=int, default=5)
    args = ap.parse_args()

    torch.manual_seed(SEED)
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK",
                                             os.cpu_count() or 8)))
    rng = np.random.default_rng(SEED)

    adata = ad.read_h5ad(PROC / "dev_subset_hvg.h5ad", backed="r")
    obs = adata.obs.reset_index(drop=True)
    obs["conc"] = obs.conc.astype(float)

    tr_rows = select_rows(obs, "train", CAP, rng)
    va_rows = select_rows(obs, "val", 100, rng)
    print(f"train cells {len(tr_rows):,}  val cells {len(va_rows):,}", flush=True)

    X_tr = adata[tr_rows].to_memory().X.tocsr()
    X_va = adata[va_rows].to_memory().X.tocsr()
    obs_tr, obs_va = obs.iloc[tr_rows], obs.iloc[va_rows]

    lines = sorted(obs.cell_line_id.unique())
    ecfp = np.load(PROC / "drug_ecfp.npz", allow_pickle=True)
    fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))
    fp_dim = len(next(iter(fp_of.values())))

    treated = obs_tr.conc > 0
    logdose = np.log10(obs_tr.conc[treated])
    mu, sd = float(logdose.mean()), float(logdose.std())

    def features(o: pd.DataFrame):
        li = torch.as_tensor([lines.index(l) for l in o.cell_line_id])
        fp = torch.as_tensor(np.stack([
            fp_of[d] if d in fp_of else np.zeros(fp_dim, np.float32)
            for d in o.drug]))
        isc = (o.conc.to_numpy() == 0)
        ld = np.where(isc, 0.0, (np.log10(np.where(isc, 1.0, o.conc)) - mu) / sd)
        dose = torch.as_tensor(
            np.stack([ld, isc.astype(float)], axis=1).astype(np.float32))
        return li, fp, dose

    li_tr, fp_tr, dose_tr = features(obs_tr)
    li_va, fp_va, dose_va = features(obs_va)

    model = CondVAE(n_genes=adata.shape[1], n_lines=len(lines),
                    fp_dim=fp_dim, n_latent=args.n_latent)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

    def batch_tensors(X, rows_local):
        dense = torch.as_tensor(X[rows_local].toarray())
        lib = dense.sum(1, keepdim=True)
        logn = torch.log1p(dense / (lib + 1e-8) * 1e4)
        return logn, dense, lib

    def val_loss():
        model.eval()
        tot, n = 0.0, 0
        with torch.no_grad():
            for s in range(0, X_va.shape[0], 1024):
                sl = np.arange(s, min(s + 1024, X_va.shape[0]))
                logn, dense, lib = batch_tensors(X_va, sl)
                nll, kl = model(logn, dense, lib, li_va[sl], fp_va[sl],
                                dose_va[sl])
                tot += float(nll + kl) * len(sl)
                n += len(sl)
        return tot / n

    hist, best, patience, best_state = [], np.inf, 0, None
    n_tr = X_tr.shape[0]
    for epoch in range(args.epochs):
        model.train()
        beta = min(1.0, (epoch + 1) / args.kl_warmup)
        perm = rng.permutation(n_tr)
        run_nll = run_kl = 0.0
        for s in range(0, n_tr, args.batch_size):
            b = perm[s:s + args.batch_size]
            logn, dense, lib = batch_tensors(X_tr, b)
            opt.zero_grad()
            nll, kl = model(logn, dense, lib, li_tr[b], fp_tr[b], dose_tr[b])
            (nll + beta * kl).backward()
            opt.step()
            run_nll += float(nll) * len(b)
            run_kl += float(kl) * len(b)
        vl = val_loss()
        hist.append({"epoch": epoch, "beta": beta, "train_nll": run_nll / n_tr,
                     "train_kl": run_kl / n_tr, "val_loss": vl})
        print(f"epoch {epoch:3d} beta {beta:.2f} "
              f"nll {run_nll/n_tr:9.2f} kl {run_kl/n_tr:7.2f} val {vl:9.2f}",
              flush=True)
        if beta >= 1.0:                      # only judge after warmup
            if vl < best - 1e-3:
                best, patience = vl, 0
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
            else:
                patience += 1
                if patience >= args.patience:
                    print("early stop"); break

    CKPT.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": best_state or model.state_dict(),
                "lines": lines, "dose_mu": mu, "dose_sd": sd,
                "n_genes": adata.shape[1], "fp_dim": fp_dim,
                "n_latent": args.n_latent},
               CKPT / "model.pt")
    pd.DataFrame(hist).to_csv(CKPT / "training_curve.csv", index=False)
    print(f"model -> {CKPT/'model.pt'}")


if __name__ == "__main__":
    main()
