# Results

Modelling drug-perturbation response in **Tahoe-100M** (Zhang et al. 2025,
[doi:10.1101/2025.02.20.639398](https://doi.org/10.1101/2025.02.20.639398)).
All numbers below are reproducible from this repository; every figure is a
bundle (PNG+SVG+PDF + source CSVs + generating script) under
`results/figures/`, and every heavy job ran through SLURM (`jobs/`).

Three working subsets are used throughout:

| Subset | Lines | Drugs | Cells | Purpose |
|---|---|---|---|---|
| `dev8` | 8 (coverage-picked, mutation-stratified) | 50 | 5.88 M | fast iteration, largest context signal |
| `dev47` | all 47 analysed lines | 50 | 20.47 M | scale, cross-line generalization |
| `full` | all 50 lines | **379** | 95.6 M | every drug; mechanism-level analyses |

Metrics are computed on plate-matched pseudobulk deltas (log1p CPM):
`r_hvg` = Pearson r over the 2,000 most response-variable genes (selected on
training data only), `r_de100` = Pearson r over each condition's own top-100
differentially expressed genes.

---

## Headline findings

**Methodological.**
1. **A mean baseline is only as good as the number of contexts it averages.**
   Going from 5 to 45 cell lines improves the additive baseline by **+10.6%**
   (r 0.649 → 0.718) with everything else held fixed — the same order as the
   improvement recent models report over baselines estimated from six lines.
   Benchmarks that subsample contexts for tractability will overstate headroom
   (§2).
2. **Plate structure survives plate-matched normalisation.** Within-plate
   comparisons in this atlas agree ~7× better than cross-plate ones
   (r 0.45 vs 0.06); pooling them overstates reproducibility (§4).

**Architecture of drug response.**
3. Of the *reproducible* transcriptional response, **~59% is the drug's average
   effect and ~41% is line×drug interaction** (full atlas, 379 drugs) (§4).
4. That interaction is a **line×drug interaction, not a line property**: it
   reproduces across doses (r = 0.062) but barely transfers across drugs
   (0.019) or lines (−0.002). Hence no line-level descriptor can predict it (§4).
5. It is **not a viability artefact** — the rank-1 model is falsified; at least
   six components reproduce on held-out plates, and the leading one is
   uncorrelated with cell-cycle arrest (§4).
6. At the cell level it is a **uniform population shift, not a pre-existing
   subpopulation** (20.5M cells; median shape excess 0.91× the noise floor) (§4).

**Predictability.**
7. **Context-dependence is set by the drug's mechanism** (Kruskal–Wallis
   p = 2×10⁻³): nuclear-receptor agonists are the most context-specific
   (glucocorticoid 0.363), metabolic and RAF inhibitors the most conserved
   (0.087–0.119). Chemistry alone predicts it only weakly (r = 0.19) (§5).
8. For a **new cell line**, genotype, organ and baseline expression all fail —
   but measuring ~20 **arbitrary** compounds and fine-tuning recovers 98% of the
   achievable gain. Probe-panel design does not help (§4, §6).

**Clinical.**
9. In 10,921 TCGA tumours the programs track driver genotype (2,128 associations
   at FDR<0.05) but their **survival associations do not survive** adjustment for
   stage, grade and infiltration (§7).

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

**The baseline's strength depends on how many contexts it averages — which may
explain why the field disagrees about it.** The additive prior is a mean, so it
should sharpen as more cell lines contribute. Holding the evaluation conditions
fixed (leave-one-line-out on the full atlas) and varying only the number of
other lines averaged:

| lines averaged | additive r (top-100 DE) |
|---|---|
| 5 | 0.649 |
| 10 | 0.690 |
| 20 | 0.709 |
| 45 | **0.718** |

That is a **+10.6%** improvement from 5 to 45 lines, attributable to context
count alone. For comparison, MAP ([Nat Mach Intell 2026](https://doi.org/10.1038/s42256-026-01286-w))
reports +12.3% for its model over the best baseline while restricting
Tahoe-100M to **six** cell lines. The gap a sparse benchmark opens up in the
baseline is therefore of the same order as the gap such benchmarks attribute to
the model. This does not show that any particular model's gains are artefactual
— a well-estimated baseline may still be beaten — but it does mean **a mean
baseline must be estimated from all available contexts before headroom is
claimed**, and that studies using few contexts will systematically overstate it.

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
**72% is the drug's average effect and 28% is context interaction**.

Repeating this on the **full atlas** — all 379 drugs, 65,918 conditions
aggregated from every one of the 95.6 M cells — replicates the structure and
shifts the balance toward context: additive **4.8%**, interaction **3.3%**,
noise 91.9%, i.e. **59% / 41%** of the reproducible signal. The three residual
correlations reproduce almost exactly (dose +0.062, drug +0.019, line −0.002),
as does the plate confound (same-plate +0.447 vs cross-plate +0.062). More
mechanistically diverse compounds yield proportionally more context-specific
response, which section 5 explains.

The genes carrying the interaction are biologically coherent rather than random —
acute-phase (SAA1, SAA2, SAA2-SAA4), imprinted/oncofetal (H19), secreted
protease inhibitors (SERPINB2, SERPINB3) and MAPK feedback (SPRY4) — none of
them canonical drug-target genes.

**The interaction is not a viability artefact.** The obvious deflationary
explanation is the rank-1 model implicit in MIX-Seq (McFarland et al. 2020):
residual ≈ (that pair's sensitivity) × (one shared death/stress direction),
which would explain both the interaction and its failure to transfer. Testing
it with split-half validation over plates (2,247 line×drug pairs measured on
≥2 plates, components fitted on one half and scored on the other):

- component 1 carries only **34%** of reproducible residual variance (top 3:
  49%) — not rank-1;
- **at least six components reproduce** on held-out plates (split-half
  r = 0.33–0.48), while components 3, 4 and 8 do not (r < 0), confirming the
  procedure discriminates signal from noise;
- component 1 is **uncorrelated with cell-cycle arrest** measured independently
  from the single-cell data (G2M r = −0.05, G1 r = +0.05).

Component 1 instead loads on stress and stemness/adhesion genes — HMOX1
(oxidative stress), HSPA1B (proteotoxic stress), LGR5 (Wnt/stem), COL17A1,
SEMA7A, PLA2R1. So context-specific drug response is a **multi-dimensional set
of reproducible programs**, not a single "some lines die more" axis.

## 5. What *is* predictable: context-dependence is set by the drug's mechanism
`results/figures/10_drugs/`

Sections 3–4 ask what makes a *cell line* respond differently. The
complementary question has a positive answer. For each of 367 drugs we split
its reproducible effect into a **conserved** part (its dose-matched mean
response across lines) and a **context** part (line-specific variance,
estimated from cross-plate, cross-dose residual covariance), and define

    CDI = context / (conserved + context)

so 0 means the drug does the same thing in every line and 1 means its
reproducible effect is entirely line-specific. Median CDI is 0.160, and
**mechanism explains a significant share of it** (Kruskal–Wallis H = 46.9,
p = 2.3 × 10⁻³ across 24 classes):

| most context-dependent | CDI | most conserved | CDI |
|---|---|---|---|
| Glucocorticoid receptor agonist | 0.363 | Glucose transporter inhibitor | 0.087 |
| Proteasome inhibitor | 0.334 | Other MAPK inhibitor | 0.109 |
| MEK inhibitor | 0.324 | RAF inhibitor | 0.119 |
| Retinoic receptor agonist | 0.232 | Other TK inhibitor | 0.122 |

The leading class is mechanistically expected, which is what makes it both a
positive control and a result: glucocorticoid receptor agonists act through a
nuclear receptor whose transcriptional output is set by each cell's own
enhancer landscape, so their effects *should* be maximally context-specific.
Four of the twelve most context-dependent individual compounds are
corticosteroids (budesonide 0.428, dexamethasone 0.376, betamethasone 0.363,
triamcinolone 0.342), and retinoic receptor agonists — also nuclear receptors —
rank fourth by class. Clobetasol propionate ranks fifth overall (0.359) while
the atlas's GPT-derived annotation calls it "unclear", so the analysis
recovers a corticosteroid the metadata mislabels.

The **MEK versus RAF contrast** is the sharpest pharmacological observation
here: MEK inhibition is among the most context-dependent mechanisms (0.324)
while RAF inhibition is among the most conserved (0.119), one step apart in the
same pathway.

Two constraints on the claim. Chemical structure alone predicts CDI only weakly
(cross-validated r = 0.19, ridge on ECFP4), so context-dependence is a property
of the target and mechanism rather than of the molecule — it cannot be read off
a compound before choosing its target. And CDI correlates modestly with effect
size (Spearman ρ = +0.38), so potent drugs are somewhat more context-dependent,
though that does not account for the mechanism ranking.

**Practical consequence.** How many cell models a compound must be screened in
is predictable from its mechanism: nuclear-receptor, proteasome and MEK
pharmacology needs broad panels, whereas RAF, metabolic and most TK inhibitors
transfer.

**It is not a target-genetics phenomenon.** For the 254 drugs with mappable
annotated targets, CDI is *uncorrelated* with how often the target is mutated
across the panel (Spearman ρ = +0.10, p = 0.10) and with target expression
level (ρ = +0.08, p = 0.21); it is weakly *negatively* correlated with target
expression variance (ρ = −0.21, p = 8×10⁻⁴). The clearest case is the leading
class itself: nuclear-receptor drugs have a median CDI of 0.342 versus 0.154
overall, yet their targets are mutated in **zero** atlas lines.

Taken with the failure of the genotype scans (§4, 0/825 mechanism × driver
tests surviving FDR) and the identity of the response programs themselves
(§4: epithelial-vs-neuronal, EMT, chromatin/senescence), this points to a
single conclusion: **context-dependence in drug response is set by the cell's
regulatory and chromatin state, not by mutations in the drug's target or by
driver genotype.** That is why nuclear receptors — whose output is read out
through each cell's enhancer landscape — sit at the top of the ranking, and why
no mutation-level descriptor predicts the interaction.

## 6. What does work for a new line: measure ~20 compounds and fine-tune
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

**And the panel needs no design.** Comparing five ways of choosing those
compounds over 47 folds and 2,809 held-out conditions (median paired gain over
additive, win rate): at k = 5, chemically diverse +0.018 (76%) and **random
+0.017 (76%)** lead, while MOA-diverse +0.010 (64%), most-discriminative +0.006
(60%) and strongest-effect +0.005 (57%) trail. By **k = 20 every strategy
converges** on the oracle (0.721 vs 0.7217). Selecting the most potent or the
most line-discriminating compounds is therefore *worse* than sampling at
random, presumably because both pick atypical compounds while random sampling
matches the distribution the model is scored on. The protocol is simply: measure
~20 arbitrary compounds and fine-tune.

**Unseen drugs**, by contrast, generalize through chemistry: an ECFP-conditioned
model trained with prior-dropout reaches r_de100 0.388 vs 0.267 for a
nearest-drug-by-fingerprint baseline (dev47, 5-fold drug-grouped CV).

## 7. Patient anchoring (TCGA): genotype yes, survival no

The six components were scored across **10,921 TCGA tumours** (31–33 cancer
types), using each component's top ±100 loading genes and z-scoring **within**
cancer type so lineage cannot manufacture the signal.

- **Coherence — no.** Component genes are no more co-expressed in tumours than
  size-matched random gene sets (z −0.83 to +0.63). These are *response*
  programs, apparent when cells are perturbed, not standing co-expression
  modules. A real constraint on how far they generalise.
- **Genotype — yes.** 2,128 of 106,896 (driver × component) tests reach
  FDR < 0.05, led by KRAS×C7, BRAF×C7 (both MAPK, same component), IDH1×C10,
  TP53×C1, CTNNB1×C1. Where 49 cell lines showed nothing, ~8,700 patients
  resolve the association — though effect sizes are small (0.05–0.11 z).
- **Survival — does not survive adjustment.** Cox models fitted within each
  cancer type and combined across them give C6 z = +3.63 (p = 3e-4) and
  C7 z = −2.19 unadjusted; after adjusting for age, sex, stage, grade and
  stromal/immune infiltration, **no component remains significant** once the
  six tests are accounted for (largest: C10 z = +2.32, p = 0.02, which does not
  survive correction across components). The apparent unadjusted associations
  are largely explained by disease burden and sample composition.

> **Correction.** An earlier median-split log-rank version of this analysis
> reported C1 as predicting worse and C7 better survival. Its sign convention
> was computed from median survival among events only via `(x or 0)`, which
> returns NaN rather than 0 because NaN is truthy, so the signs were unreliable.
> The Cox analysis supersedes it and that claim is withdrawn.

Caveat on purity: the TCGA ABSOLUTE table is not publicly downloadable from the
endpoints available here, so purity is proxied by hallmark stromal (EMT) and
immune infiltration scores. This captures the dominant non-tumour axis but is
weaker than ABSOLUTE.

## 8. Independent audit: the atlas's cell-line identities hold up
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
