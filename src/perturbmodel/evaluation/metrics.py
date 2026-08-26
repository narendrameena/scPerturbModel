"""Evaluation metrics for perturbation prediction.

E-distance follows Peidli et al. 2024 (scPerturb), as used in the Tahoe-100M
paper: computed between perturbed and control cell populations in a latent
space (the paper used a 10-d scVI latent).
"""
from __future__ import annotations

import torch


def _mean_pairwise_dist(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    return torch.cdist(x, y, p=2).mean()


def e_distance(x: torch.Tensor, y: torch.Tensor, unbiased: bool = True) -> float:
    """Energy distance between two cell populations x (n, d) and y (m, d).

    E = 2 * E||x - y|| - E||x - x'|| - E||y - y'||
    With `unbiased=True` the within-population terms exclude self-distances.
    """
    x = torch.as_tensor(x, dtype=torch.float32)
    y = torch.as_tensor(y, dtype=torch.float32)
    delta = 2.0 * _mean_pairwise_dist(x, y)
    for pop in (x, y):
        n = pop.shape[0]
        d = torch.cdist(pop, pop, p=2)
        if unbiased and n > 1:
            within = d.sum() / (n * (n - 1))  # exclude zero diagonal
        else:
            within = d.mean()
        delta = delta - within
    return float(delta)
