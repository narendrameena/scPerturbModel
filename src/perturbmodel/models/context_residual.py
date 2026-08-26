"""Phase 3 model 1: context-residual delta prediction.

pred_delta(L, D, C) = additive_prior(D, C) + residual(line_emb, ECFP(D), dose)

The residual head is zero-initialized, so the model starts exactly at the
additive-shift baseline; training can only move away from it where the data
supports context-dependence. Drug identity enters purely through chemistry
(ECFP), so unseen drugs are representable.
"""
from __future__ import annotations

import torch
from torch import nn


class ContextResidualDelta(nn.Module):
    def __init__(self, n_lines: int, n_genes: int, fp_dim: int = 1024,
                 d_line: int = 32, d_drug: int = 128, d_dose: int = 16,
                 hidden: int = 512, dropout: float = 0.1,
                 use_line: bool = True):
        super().__init__()
        self.use_line = use_line
        self.line_emb = nn.Embedding(n_lines, d_line)
        self.drug_net = nn.Sequential(nn.Linear(fp_dim, d_drug), nn.ReLU())
        self.dose_net = nn.Sequential(nn.Linear(1, d_dose), nn.ReLU())
        self.trunk = nn.Sequential(
            nn.Linear(d_line + d_drug + d_dose, hidden), nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.ReLU(), nn.Dropout(dropout),
        )
        self.out = nn.Linear(hidden, n_genes)
        nn.init.zeros_(self.out.weight)
        nn.init.zeros_(self.out.bias)

    def forward(self, line_idx: torch.Tensor, fp: torch.Tensor,
                logdose: torch.Tensor, prior: torch.Tensor) -> torch.Tensor:
        le = self.line_emb(line_idx)
        if not self.use_line:            # ablation: no cellular context
            le = torch.zeros_like(le)
        h = torch.cat([le, self.drug_net(fp), self.dose_net(logdose)], dim=-1)
        return prior + self.out(self.trunk(h))
