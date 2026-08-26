#!/usr/bin/env python3
"""Phase 2: train the scVI baseline on the dev subset (CPU cluster).

Trains on split=='train' cells only (held-out conditions and plate14 never seen),
then writes the 10-d latent for ALL cells — downstream scripts use it for
E-distance evaluation and the latent-shift baseline.

Requires: data/processed/dev_subset_hvg.h5ad (scripts/make_training_matrix.py).

Outputs:
  results/checkpoints/scvi_dev/         scvi-tools model
  data/processed/scvi_latent_dev.npz    latent (all cells, training order of h5ad)

Usage: python scripts/train_scvi_dev.py [--max-cells 1500000] [--epochs 10]
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
H5AD = ROOT / "data" / "processed" / "dev_subset_hvg.h5ad"
CKPT = ROOT / "results" / "checkpoints" / "scvi_dev"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-cells", type=int, default=1_500_000,
                    help="subsample of train cells (CPU budget)")
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--n-latent", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=4096)
    ap.add_argument("--litlogger", action="store_true",
                    help="log metrics to Lightning.ai via litlogger "
                         "(needs a Lightning account / API key)")
    args = ap.parse_args()

    import scvi  # local import: heavy

    torch.set_num_threads(min(32, os.cpu_count() or 8))
    adata = ad.read_h5ad(H5AD)
    print(f"loaded {adata.shape[0]:,} cells x {adata.shape[1]} genes")

    train_mask = (adata.obs.split == "train").to_numpy()
    idx = np.where(train_mask)[0]
    if len(idx) > args.max_cells:
        idx = np.random.default_rng(0).choice(idx, args.max_cells, replace=False)
    sub = adata[np.sort(idx)].copy()
    print(f"training on {sub.shape[0]:,} cells")

    trainer_kwargs = {}
    if args.litlogger:
        try:
            from litlogger import LightningLogger
            trainer_kwargs["logger"] = LightningLogger(name="scvi_dev")
            print("litlogger enabled -> Lightning.ai experiments")
        except Exception as exc:  # no account/API key or offline
            print(f"litlogger unavailable ({exc}); continuing without it")

    scvi.model.SCVI.setup_anndata(sub)
    model = scvi.model.SCVI(sub, n_latent=args.n_latent, gene_likelihood="nb")
    model.train(max_epochs=args.epochs, batch_size=args.batch_size,
                early_stopping=True, enable_progress_bar=False,
                **trainer_kwargs)
    CKPT.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(CKPT), overwrite=True)
    print(f"model -> {CKPT}")

    # latent for ALL cells (batched through the trained model)
    full = ad.read_h5ad(H5AD)
    scvi.model.SCVI.setup_anndata(full)
    z = model.get_latent_representation(full, batch_size=8192)
    np.savez_compressed(ROOT / "data" / "processed" / "scvi_latent_dev.npz",
                        latent=z.astype(np.float32))
    print(f"latent {z.shape} -> data/processed/scvi_latent_dev.npz")


if __name__ == "__main__":
    main()
