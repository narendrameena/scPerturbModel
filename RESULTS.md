# Results

Modelling drug-perturbation response in **Tahoe-100M** (Zhang et al. 2025,
[doi:10.1101/2025.02.20.639398](https://doi.org/10.1101/2025.02.20.639398)).
All numbers below are reproducible from this repository; every figure is a
bundle (PNG+SVG+PDF + source CSVs + generating script) under
`results/figures/`, and every heavy job ran through SLURM (`jobs/`).

Two working subsets are used throughout:

| Subset | Lines | Drugs | Cells | Purpose |
|---|---|---|---|---|
| `dev8` | 8 (coverage-picked, mutation-stratified) | 50 | 5.88 M | fast iteration, largest context signal |
| `dev47` | all 47 analysed lines | 50 | 20.47 M | scale, cross-line generalization |

Metrics are computed on plate-matched pseudobulk deltas (log1p CPM):
`r_hvg` = Pearson r over the 2,000 most response-variable genes (selected on
training data only), `r_de100` = Pearson r over each condition's own top-100
differentially expressed genes.

---

## 1. The atlas reproduces, and so do its published analyses
`results/figures/01_paper_replication/`

- **Cell-cycle pharmacology** (paper Fig 6C–E) recovered across all 47 lines
  from ~100 M cells: palbociclib → G1 arrest (log2 OR +0.54), dinaciclib →
  S/G2M (+1.76/+1.19), panobinostat and belinostat → G2M (+0.77/+0.62),
  tucidinostat negative (−0.09), carbamazepine inert (+0.12), and paclitaxel
  minimal (−0.01) — including every reported exception.
- **E-distance** rises with dose (median log1p 0.49 → 0.51 → 0.75) and ranks
  mechanisms as published: protein-synthesis inhibitors ≫ microtubule >
  proteasome > PI3K/AKT > HDAC.
- **Replicate plates** (6 vs 14) agree at r = 0.989 for matched conditions vs
  0.929 unmatched — slightly cleaner than the published 0.975/0.915.

## 2. The additive baseline is very strong, and latent arithmetic loses to it
`results/figures/02_baselines/`

Held-out (line, drug, dose) conditions, dev8:

| Predictor | r_hvg | r_de100 | RMSE |
|---|---|---|---|
| no change | — | — | 0.318 |
| random train delta (same line, wrong drug) | 0.159 | 0.343 | 0.458 |
| **additive shift** (same drug+dose, other lines) | **0.570** | **0.777** | **0.268** |

An scVI latent-shift model (scGen-style arithmetic in a 10-d latent space,
trained on 1.5 M cells) reaches r_de100 0.623 against 0.901 for the additive
baseline in the same gene space — the bottleneck discards drug-specific
structure. This independently reproduces the field's central critique
([Nature Methods 2025](https://www.nature.com/articles/s41592-025-02772-6)) on
chemical rather than genetic perturbations.

## 3. A context-residual model beats the baseline, and the gain is real
`results/figures/03_models/`

`prediction = additive_prior(drug, dose) + residual(line embedding, ECFP
fingerprint, log dose)`, with the residual zero-initialised (so the model
*starts* at the baseline) and leave-one-line-out priors on training rows.

| Split | additive r_hvg | model r_hvg | conditions improved | median paired gain (r_de100) |
|---|---|---|---|---|
| dev8, held-out condition | 0.570 | 0.729 | 234/238 (98%) | +0.037 |
| dev47, held-out condition | 0.498 | 0.544 | 1288/1464 (88%) | +0.006 |
| **dev47, held-out (line,drug) pair** | 0.498 | **0.564** | 1369/1495 (92%) | **+0.010** |

The no-context ablation lands *exactly* on the additive baseline in every run,
so the gain is attributable to context by construction. Crucially the gain
**survives the strict pair split**, which removes all other doses of a
line–drug pair from training — it is not dose leakage.

A conditional VAE over single cells (direct drug/dose/line conditioning, NB
likelihood) generates populations closer to the real perturbed cells than the
unperturbed control on 69.6% of held-out conditions (E-distance 1.24 vs 1.66,
noise floor ≈ 0), with delta fidelity r_de100 0.701.

## 4. Unseen cell lines: a triple negative, and its explanation
`results/figures/04_generalization/`, `results/figures/06_diagnostics/`

Nothing we tried lets the model predict a **never-seen** line's response
modulation:

| Context source for an unseen line | r_de100 | vs additive |
|---|---|---|
| driver mutations + organ (dev8, 8-fold) | 0.764 | 0.768 |
| driver mutations + organ (dev47, 47-fold) | 0.692 | 0.692 |
| the line's own baseline (DMSO) transcriptome, PCA (dev47) | 0.692 | 0.692 |
| probe-fitted embedding, frozen head (dev8, 8-fold) | 0.765 | 0.767 |
| probe-fitted embedding, frozen head (dev47, 47-fold × 3 seeds, 8,420 conditions) | 0.692 | 0.691 |

The frozen-head result holds at both scales and even with an **oracle** that
fits the embedding on all ~30 available probe drugs (dev8 −0.0008, dev47
+0.0009 — i.e. nothing), while the embedding demonstrably moves. The failure is
therefore not a shortage of probe data.

The structure diagnostic explains all four at once. Writing a condition's
**context residual** as `delta − additive_prior`, and comparing conditions
measured on **different plates** (see the batch caveat below), its correlation
is:

- **+0.067** across *doses* of the same (line, drug) — the residual is real,
  reproducible signal;
- **+0.013** across *drugs* within the same line — almost no line-level transfer;
- **−0.009** across *lines* for the same drug;
- ≈ 0.00 for permuted nulls in all three cases.

> **Batch caveat, and why it matters.** Restricting to cross-plate pairs is
> essential. The same comparisons *within* a plate give +0.48 and +0.40 — about
> sevenfold higher — because residual plate structure survives plate-matched
> DMSO normalisation. Any analysis of this atlas that pools within- and
> cross-plate comparisons will substantially overstate reproducibility, and an
> earlier version of this section did exactly that.

So line-specific drug response is a **line × drug interaction, not a line
property**. There is no line-level quantity for a covariate, a baseline
profile, or a probe-fitted embedding to recover — which is exactly why all
three fail, and why the measured gains are small wherever the effect is real.

**Variance decomposition** (per gene, cross-plate dose pairs only, 47 lines ×
50 drugs): additive drug×dose effect **6.9%** of response variance, line×drug
interaction **2.7%**, noise 90.4%. Of the *reproducible* signal, roughly
**72% is the drug's average effect and 28% is context interaction**. The genes
carrying the interaction are biologically coherent rather than random —
acute-phase (SAA1, SAA2, SAA2-SAA4), imprinted/oncofetal (H19), secreted
protease inhibitors (SERPINB2, SERPINB3) and MAPK feedback (SPRY4) — none of
them canonical drug-target genes.

## 5. What does work for a new line: measure ~20 compounds and fine-tune
`results/figures/04_generalization/few_shot_eval_dev8ft/`

Adapting to a new line succeeds only if the residual head is allowed to
co-adapt with the new embedding. Gain in r_de100 over the additive baseline,
leave-one-line-out with a fixed evaluation panel reserved before probing:

| probe drugs measured on the new line | frozen head (dev47) | fine-tuned (dev8) | **fine-tuned (dev47)** |
|---|---|---|---|
| 1 | +0.0002 | +0.0022 | +0.0001 |
| 5 | +0.0006 | +0.0086 | **+0.0199** |
| 20 | +0.0009 | +0.0102 | **+0.0343** |
| all ~30 (oracle) | +0.0009 | — | +0.0350 |

At full scale the effect is large and monotone: **20 probe compounds recover
98% of the oracle ceiling** (r_de100 0.691 → 0.725), and even 5 compounds
recover 57%. The identical probe data yields *nothing* when the residual head
is frozen (n = 8,420 conditions), so the barrier is architectural, not
informational: the learned context space behaves like a lookup over training
lines rather than a space one can interpolate into. Practically, a new cell
model needs neither its genotype nor its baseline profile — it needs a handful
of measured compounds and a fine-tuning pass.

**Unseen drugs**, by contrast, generalize through chemistry: an ECFP-conditioned
model trained with prior-dropout reaches r_de100 0.388 vs 0.267 for a
nearest-drug-by-fingerprint baseline (dev47, 5-fold drug-grouped CV).

## 6. Independent audit: the atlas's cell-line identities hold up
`results/figures/05_audit/`

47 lines scored against **184 (tissue, cell-type) signatures** built from 12
Tabula Sapiens tissues over 60,602 shared genes. Organ concordance is 17%
top-1 and 43% top-3 after per-tissue normalisation, with clear identity blocks:
all nine colon lines strongest at large intestine, plus stomach, liver, breast,
bladder and uterine matches, and a neuroendocrine positive control (SHP-77, a
small-cell lung cancer line, matching a neuronal signature). Lines that miss
collapse onto a *generic proliferative-epithelium* signature rather than a
confident wrong tissue — culture dedifferentiation, not mislabelling. This
also motivates using learned line embeddings over tissue-of-origin labels.

---

## Reproducing

```bash
python scripts/download_tahoe100m.py --mode full     # 429 GB from HuggingFace
python scripts/select_dev_subset.py && sbatch jobs/dev47_extract.sbatch
sbatch jobs/dev47_pseudobulk.sbatch
sbatch jobs/dev47_delta.sbatch                       # model 1
sbatch jobs/phase3_cvae_train.sbatch                 # model 2
sbatch jobs/dev47_hardsplits.sbatch                  # generalization
sbatch jobs/few_shot_finetune_dev8.sbatch            # few-shot transfer
sbatch jobs/cell_type_audit.sbatch                   # identity audit
```

See `docs/related_work_perturbation_models.md` for how these results sit
against published methods, and `docs/data_notes.md` for verified dataset facts.
