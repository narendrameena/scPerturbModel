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
| probe-fitted embedding, frozen head — *including an oracle using all ~30 probe drugs* | 0.764 | 0.765 |

The structure diagnostic explains all four at once. Writing a condition's
**context residual** as `delta − additive_prior`, its correlation is:

- **+0.13** across *doses* of the same (line, drug) — the residual is real,
  reproducible signal;
- **+0.04** across *drugs* within the same line — almost no line-level transfer;
- **−0.01** across *lines* for the same drug;
- ≈ 0.00 for permuted nulls in all three cases.

So line-specific drug response is a **line × drug interaction, not a line
property**. There is no line-level quantity for a covariate, a baseline
profile, or a probe-fitted embedding to recover — which is exactly why all
three fail, and why the measured gains are small wherever the effect is real.

## 5. What does work for a new line: measure ~20 compounds and fine-tune
`results/figures/04_generalization/few_shot_eval_dev8ft/`

Adapting to a new line succeeds only if the residual head is allowed to
co-adapt with the new embedding (dev8, 8-fold, fixed evaluation panel, paired
per condition against the additive baseline):

| probe drugs measured on the new line | frozen head | fine-tuned |
|---|---|---|
| 1 | 0.000 | +0.0022 |
| 5 | 0.000 | +0.0086 |
| 20 | 0.000 | **+0.0102** |
| all ~30 (oracle) | −0.0008 | — |

≈20 probe compounds recover approximately the entire available line-level gain
(+0.0101 measured for *seen* lines under the strict pair split). The learned
context space therefore behaves like a lookup over training lines rather than
a space one can interpolate into — a concrete constraint for virtual-cell
models, and a concrete protocol for screening a new cell model.

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
