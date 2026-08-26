"""Common interface every perturbation model implements, so training and
evaluation code stay model-agnostic (baselines and deep models alike)."""
from __future__ import annotations

from abc import abstractmethod

import torch
from torch import nn


class PerturbationModel(nn.Module):
    """Predict the perturbed transcriptomic state of a cell.

    Contract: given a control/basal cell state and a perturbation description
    (drug embedding + dose + cellular context), produce parameters of the
    predicted expression distribution (e.g. negative binomial mean/dispersion),
    not just a point estimate — heterogeneity is the point of this dataset.
    """

    @abstractmethod
    def forward(
        self,
        counts: torch.Tensor,        # (B, n_genes) control-state counts
        drug: torch.Tensor,          # (B, d_drug) drug embedding (onehot/ECFP/...)
        dose: torch.Tensor,          # (B, 1) log-dose
        context: torch.Tensor,       # (B, d_ctx) cell line / mutation context
    ) -> dict[str, torch.Tensor]:
        """Return {'mean': ..., 'dispersion': ..., ...} for the perturbed state."""
        raise NotImplementedError

    @torch.no_grad()
    def predict_pseudobulk_delta(self, *args, **kwargs) -> torch.Tensor:
        """Convenience: expected perturbed-minus-control mean expression;
        default derives it from forward(); override for baselines."""
        raise NotImplementedError
