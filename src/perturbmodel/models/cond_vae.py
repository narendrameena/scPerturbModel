"""Phase 3 model 2: conditional VAE over single-cell HVG counts.

Direct conditioning (the Phase 2 lesson): the condition vector
c = [line_emb, drug_net(ECFP), dose_net(std_logdose, is_control)] enters BOTH
the encoder and the decoder — no post-hoc latent arithmetic. Likelihood is
negative binomial with per-gene dispersion (scVI-style), so generated cells are
counts, comparable to real cells after identical normalization.

Generation for a condition: z ~ N(0, I), decode(z, c) -> px_scale (sums to 1
over genes); counts ~ NB(mean = px_scale * library, theta) with library sizes
resampled from real cells of the same line.
"""
from __future__ import annotations

import torch
from torch import nn
from torch.distributions import NegativeBinomial


class CondVAE(nn.Module):
    def __init__(self, n_genes: int, n_lines: int, fp_dim: int = 1024,
                 d_line: int = 32, d_drug: int = 128, d_dose: int = 16,
                 n_latent: int = 16, hidden: int = 512, dropout: float = 0.1):
        super().__init__()
        self.line_emb = nn.Embedding(n_lines, d_line)
        self.drug_net = nn.Sequential(nn.Linear(fp_dim, d_drug), nn.ReLU())
        self.dose_net = nn.Sequential(nn.Linear(2, d_dose), nn.ReLU())
        d_cond = d_line + d_drug + d_dose

        self.encoder = nn.Sequential(
            nn.Linear(n_genes + d_cond, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.ReLU(),
        )
        self.z_mu = nn.Linear(hidden, n_latent)
        self.z_logvar = nn.Linear(hidden, n_latent)

        self.decoder = nn.Sequential(
            nn.Linear(n_latent + d_cond, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.ReLU(),
        )
        self.px_logits = nn.Linear(hidden, n_genes)   # softmax -> px_scale
        self.log_theta = nn.Parameter(torch.zeros(n_genes))
        self.n_latent = n_latent

    def condition(self, line_idx, fp, dose_feat):
        return torch.cat([self.line_emb(line_idx), self.drug_net(fp),
                          self.dose_net(dose_feat)], dim=-1)

    def forward(self, x_lognorm, counts, library, line_idx, fp, dose_feat):
        """Return (nb_nll, kl) per batch (means)."""
        c = self.condition(line_idx, fp, dose_feat)
        h = self.encoder(torch.cat([x_lognorm, c], dim=-1))
        mu, logvar = self.z_mu(h), self.z_logvar(h).clamp(-8, 8)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        px_scale = torch.softmax(self.px_logits(
            self.decoder(torch.cat([z, c], dim=-1))), dim=-1)
        px_mu = px_scale * library
        theta = torch.exp(self.log_theta).clamp(1e-3, 1e4)
        nb = NegativeBinomial(total_count=theta,
                              logits=(px_mu + 1e-8).log() - theta.log())
        nll = -nb.log_prob(counts).sum(-1).mean()
        kl = 0.5 * (mu.pow(2) + logvar.exp() - 1 - logvar).sum(-1).mean()
        return nll, kl

    @torch.no_grad()
    def generate(self, n: int, line_idx, fp, dose_feat, library,
                 sample_counts: bool = True):
        """Generate n cells for ONE condition. line_idx: scalar tensor;
        fp: (fp_dim,); dose_feat: (2,); library: (n,) sampled real libraries."""
        c = self.condition(line_idx.expand(n), fp.expand(n, -1),
                           dose_feat.expand(n, -1))
        z = torch.randn(n, self.n_latent)
        px_scale = torch.softmax(self.px_logits(
            self.decoder(torch.cat([z, c], dim=-1))), dim=-1)
        px_mu = px_scale * library[:, None]
        if not sample_counts:
            return px_mu
        theta = torch.exp(self.log_theta).clamp(1e-3, 1e4)
        nb = NegativeBinomial(total_count=theta,
                              logits=(px_mu + 1e-8).log() - theta.log())
        return nb.sample()
