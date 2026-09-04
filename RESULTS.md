# Results

> ## ⚠ STATUS (2026-09-03): every interaction number recomputed under a three-way decomposition
>
> A fourth defect, larger than the ~60 the three audits found, was located on
> 2026-09-03: **every residual in this project removed the compound main effect
> but not each cell line's general sensitivity** — its response to *every*
> compound, set by growth rate, seeding density and drug metabolism. That term
> reproduces across disjoint compound halves at **r = 0.989** (737 lines) and is
> shared between replicate detection plates exactly as a genuine interaction is,
> so it entered the pair covariance and was counted as context-dependence.
>
> Correct apportionment (PRISM, 737 lines × 1,324 compounds): **drug effect
> 94.1%, cell property 2.0%, cell–drug relation 3.9%.** Prior interaction
> numbers were roughly **1.5× too large**. Ben-David et al. (Nature 2018) name
> the same axis: their most resistant MCF7 strains are resistant *in general*,
> through downregulated drug-metabolism pathways.
>
> All estimators now take both main effects out of fold — β from other lines,
> α from other compounds and **within the same plate** — via
> `perturbmodel.celldrug`. See §27. Any section not yet re-run is stale.
>
> The three earlier audits' repairs still stand: the laboratory-versus-assay
> apportionment inverts on a matched compound set (§18), the copy-number claim
> is withdrawn (§17), and two dose-mechanism tests sit at their permutation
> null (§22).

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
3. Of the *reproducible* transcriptional response, the context×drug interaction
   is **11.5% [10.5–12.5%] in Tahoe** measured on its designed replicate plate,
   and **57% in LINCS phase 1** (6.1M same-dose replicate pairs; different
   platform and normalisation, so magnitudes are not directly comparable).
   Tahoe's replicate plate 14 is a deliberate copy of plate 6 that training
   pipelines discard; without it the only available pairing is cross-dose, which
   inflates the estimate to 20.7% — **nearly double** (§4, §14).
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
   but measuring ~20 **arbitrary** compounds and fine-tuning recovers 94% of the
   achievable gain. Probe-panel design does not help (§4, §6). **That gain is
   not context-specific pharmacology:** adding the line's mean residual over the
   same probes to the additive prior — one arithmetic step, no model — recovers
   145% of the headroom and beats the fine-tuned model by 0.023 *r* [0.018,
   0.029] (§32).

**Clinical.**
9. In 10,921 TCGA tumours the programs track driver genotype (2,128 associations
   at FDR<0.05) but their **survival associations do not survive** adjustment for
   stage, grade and infiltration (§7).

**The binding constraint — now measured, not inferred.**
10. Tahoe-100M has 100 M cells but only **~47 usable contexts**. Every
    per-condition analysis here (thousands of conditions) is well powered and
    replicates; every line-level predictor we tried failed — plausibly one
    power ceiling hit four times rather than four separate biological
    negatives (§10).
11. **PRISM proves it.** At 738 cell lines the genotype x compound scan recovers
    the clinical biomarker set de novo (TP53/MDM2i, BRAF/vemurafenib,
    PIK3CA/alpelisib, KRAS/MEKi; 80 hits at FDR<0.05). Subsampling those same
    associations to **47 lines recovers 4% of them**; 80% power needs **~400
    contexts**, 95% needs ~600. The Tahoe genotype negatives are a power
    ceiling, not a biological absence (§11).

**Scope limits discovered by replication.**
12. **The mechanism ranking does not reproduce outside Tahoe.** It is solid
    *within* Tahoe (Kruskal–Wallis p = 6.0×10⁻⁴ across 24 mechanisms) but gives
    rho = +0.09 (n.s.) against LINCS phase 1 and −0.09 against PRISM viability
    over 67 classes — a well-powered null. An earlier +0.56 was an artefact of a
    biased estimator and is withdrawn (§12, §16).
13. **Context-dependence is dose-dependent.** It rises ~4x with dose and peaks
    just below lethality, then collapses when every line is dying (PRISM,
    rho = +0.33, p = 4x10⁻⁹²; peak 0.384 at ~2.5 uM vs 0.218 at 10 uM).
    The rising limb replicates in LINCS transcription (rho = +0.20, p = 1x10⁻⁵).
    A CDI is only comparable at matched dose (§13).

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

> **Variance shares here are upper bounds — see §14.** Tahoe replicates only 4.7%
> of its (line, drug, dose) combinations, so the pairing behind these shares is
> largely cross-dose rather than true replicate.
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

> **Scope narrowed by §16.** This ranking is significant *within* Tahoe
> (Kruskal–Wallis p = 6.0×10⁻⁴) but does not reproduce against LINCS phase 1 or
> PRISM. Read it as a within-atlas result.
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

**The principle replicates in LINCS L1000; the specific ranking does not.**
Applying the same analysis to LINCS Phase II (30 cell lines, 273 compounds in
≥10 lines, 894 plates, 1.45 million cross-replicate pairs, 978 directly measured
genes) reproduces the *existence* of the effect — mechanism significantly
structures context-dependence there too (Kruskal–Wallis H = 48.2,
P = 1.6 × 10⁻³ across 24 classes, 216 annotated compounds). But the **order does
not transfer**. Across the six mechanism classes present in both panels the rank
correlation is **ρ = −0.26 (P = 0.62)**, and the largest single disagreement is
stark: RAF inhibitors are the *most* context-dependent class in LINCS (0.588)
and among the *least* in Tahoe-100M (0.119).

Three differences could produce this and we cannot presently separate them: the
cell panels (LINCS' compounds in ≥10 lines are dominated by a small core panel
of oncology models, Tahoe's by 47 diverse lines), the measurement (978 landmark
genes versus the transcriptome), and dose and timing. Note also that LINCS'
well-replicated compounds are almost entirely kinase inhibitors, so the
nuclear-receptor result that leads our Tahoe ranking is **not testable** there.

We therefore report the mechanism *effect* as replicated and the mechanism
*ranking* as Tahoe-specific pending a panel that shares more classes. The
ranking should not be used prospectively across platforms on this evidence.

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

## 9. The decomposition replicates in two independent atlases
`results/figures/12_external/`

Every result above comes from Tahoe-100M, so the decomposition was re-run
unchanged on two external atlases, mapping their columns onto the same three
roles (context / perturbation / independent replicate). The replicate axis is
what makes each estimate trustworthy: reproducible interaction is the covariance
of residuals between *independent replicates* of the same context × perturbation.

| | **Tahoe-100M** | **OP3** | **sciPlex3** |
|---|---|---|---|
| contexts | 47 cancer lines | **6 immune cell types** | 3 cancer lines |
| perturbations | 379 | 147 | 189 |
| replicate axis | plates | **donors** | experimental replicates |
| platform / lab | Parse GigaLab, Tahoe | 10x, NeurIPS 2023 | sci-Plex, Trapnell lab |
| additive share | 59% | 67% | 70% |
| **interaction share** | **41%** | **33%** | **30%** |
| residual across *replicates* | +0.062 | +0.135 | +0.048 |
| residual across *perturbations* | +0.019 | +0.007 | +0.026 |
| residual across *contexts* | −0.002 | −0.100 | −0.040 |

The central structural claim holds in all three: the residual is **reproducible
across independent replicates** of the same context × perturbation, yet
**barely transfers across perturbations** within a context. Context-dependence
is a context × perturbation interaction rather than a context property — now
shown in primary human immune cells as well as cancer cell lines, across three
platforms, three labs, and three different kinds of replicate (plate, donor,
experimental repeat).

**The per-drug ranking does not replicate — for a quantifiable reason.** We also
recomputed the per-compound context-dependence index (§5) in OP3. It fails, and
the failure is one of power rather than of biology: a compound in OP3 has a
median of **18** cross-replicate pairs from which to estimate its own
interaction (6 cell types × 3 donors), against **142** in Tahoe-100M. The
consequence is stark — only **23%** of OP3 compounds yield a positive
interaction estimate at all, versus **96%** in Tahoe; the rest come out negative
and clamp to zero, which is what an unbiased but noisy covariance estimator does
when the true value is small. Only 15 compounds are shared by name between the
two panels, too few for a cross-system correlation (Spearman ρ = +0.06,
p = 0.84). Of the three OP3 compounds with enough signal to estimate, two agree
closely with their Tahoe values (belinostat 0.270 vs 0.215; LY2090314 0.144 vs
0.160) and one does not (dabrafenib 0.276 vs 0.119) — an anecdote, not evidence.

So the **aggregate** decomposition replicates because it pools 10,317 pairs,
while the **per-compound** ranking needs per-compound replication that OP3 does
not have. This is the same lesson as §10 applied to a different axis: pooled
estimates are robust, stratified ones are limited by the replication depth
available within each stratum.

Two caveats worth stating. The across-context row is partly forced negative by
the leave-one-out prior (expected ≈ −1/(k−1), i.e. −0.20 at k = 6 and −0.022 at
k = 47), so the informative contrast is replicates versus perturbations, not
that row. And the additive share rises as contexts fall (59% → 67% → 70% at 47,
6 and 3 contexts): with few contexts the prior is noisier and its mean square is
inflated, the same effect quantified directly in §2.

## 10. The binding constraint is context count, not cell count

Four independent attempts to predict a line's context-deviance from line-level
features have failed:

| line-level predictor | result |
|---|---|
| driver mutations × mechanism (825 tests, full atlas) | 0 at FDR<0.10 |
| baseline (DMSO) transcriptome, 47-fold LOO | 0.692 vs additive 0.692 |
| drug-target mutation frequency and abundance (264 drugs) | ρ = +0.10, p = 0.10 |
| DNA methylation level and heterogeneity, 20 features, 43 lines | 0 at FDR<0.10 |

The DNA-methylation panel is worth spelling out because it was the strongest
available test of the regulatory hypothesis: CCLE bisulfite profiles of 928
lines (45 of our 50), including epigenetic *heterogeneity* metrics — proportion
of discordant reads, methylation entropy, haplotype load, local pairwise
discordance — not merely methylation level. Best feature: promoter PDR,
ρ = −0.25, p = 0.10; all twenty trending weakly negative; unchanged after
adjusting for sequencing depth, which is itself uncorrelated with the index
(ρ = −0.23, p = 0.14). ATAC-seq would be the definitive test but ENCODE covers
only 3 of our lines and DepMap is not scriptable.

**These four negatives should probably be read as one.** With 43–50 lines, a
line-level correlation needs |ρ| > 0.30 to clear nominal significance, so any
effect of realistic size is invisible. The atlas contains 100 million cells but
only about **47 usable contexts**, and for line-level questions that is what
determines power — it is effectively an *n = 47 study* however many cells it
holds.

This explains the pattern of everything above. Analyses resting on thousands of
*conditions* — the variance decomposition (§4), the drug-mechanism ranking
(§5, 367 drugs), the few-shot protocol (§6), the baseline-scaling result (§2) —
are well powered and replicate across subsets. Analyses resting on ~47 *lines*
all fail. That is a statement about experimental design, not about biology: we
cannot conclude that genotype or chromatin state is irrelevant to
context-dependence, only that it is not detectable at this context count.

The implication for the field is direct. Effort is currently going into cell
count and model capacity, but for the context-transfer problem that these
atlases are explicitly built to solve, **the binding constraint is the number
of distinct cellular contexts profiled**. A follow-up atlas with 500 lines and
one-tenth the cells per condition would answer the questions this one cannot.

**§11 tests that recommendation rather than leaving it as a hunch.** Repeating
the genotype scan in PRISM at 738 cell lines recovers the effects that fail
here, and subsampling it back down shows they become invisible at 47 — putting
the required context count at ~400 for 80% power.

## 11. PRISM settles the power question: genotype linkage is real, and needs ~400 lines

§10 argued that four failed line-level predictors were probably one power
ceiling hit four times, but that argument could not be closed from Tahoe alone —
"undetectable at n=47" and "absent" look identical at n=47. Settling it requires
the same test at a context count where a true effect *must* show.

**PRISM Repurposing** (Corsello et al. 2020) supplies it: pooled viability for
**738 cell lines x ~1,500 compounds x 8 doses** on replicate detection plates,
with DepMap mutation calls for 1,257 of those lines. The phenotype is a scalar
rather than a transcriptome, but the decomposition is unchanged — the
leave-one-line-out shared response, and an interaction estimated as the
covariance of per-line residuals between **independent replicate plates** at
matched compound and dose.

**The architecture holds.** 1,444 compounds decomposed over >=100 lines:
**78% of reproducible variance is the compound's shared effect, 22% is
line x compound interaction** (median CDI 0.199, IQR 0.115-0.287) — the same
shape as Tahoe's 72/28 for transcription, in a different lab, assay and
phenotype.

**Genotype linkage appears, and it is textbook.** Scanning 300 compounds against
5,048 genes mutated in >=40 lines (1,504,298 tests) returns **80 associations at
FDR<0.05** across 24 compounds — and the top of the list is the clinical
biomarker set, recovered de novo:

| compound | gene | n mutant / n lines | effect | q |
|---|---|---|---|---|
| idasanutlin (MDM2i) | TP53 | 497 / 734 | +0.75 (mutant resistant) | 6e-36 |
| CGM097 (MDM2i) | TP53 | 353 / 479 | +0.59 (mutant resistant) | 4e-24 |
| alpelisib (PI3Kai) | PIK3CA | 102 / 732 | −0.20 (mutant sensitive) | 9e-09 |
| vemurafenib (RAFi) | BRAF | 104 / 714 | −0.15 (mutant sensitive) | 4e-08 |
| dabrafenib (RAFi) | BRAF | 104 / 734 | −0.22 (mutant sensitive) | 4e-07 |
| ipatasertib (AKTi) | PTEN | 66 / 479 | −0.26 (mutant sensitive) | 5e-06 |
| copanlisib (PI3Ki) | PIK3CA | 102 / 711 | −0.37 (mutant sensitive) | 2e-05 |
| pimasertib (MEKi) | KRAS | 135 / 714 | −0.31 (mutant sensitive) | 2e-04 |

Every sign is the pharmacologically correct one: MDM2 inhibitors need wild-type
p53, and BRAF/PIK3CA/PTEN/KRAS lesions sensitise to the matched pathway
inhibitor. These are not discoveries — they are a **positive control that the
estimator and the genotype join are sound**.

**The informative quantity is where they disappear.** Subsampling PRISM to fewer
cell lines and re-testing the 12 confirmed associations at the same per-test
threshold gives the recovery curve:

| cell lines | 20 | 30 | **47** | 75 | 100 | 150 | 250 | 400 | 600 |
|---|---|---|---|---|---|---|---|---|---|
| associations recovered | 0% | 1% | **4%** | 9% | 15% | 24% | 45% | **72%** | **96%** |

**At Tahoe's 47 contexts, effects this large are recovered 4% of the time.**
That is a direct, quantitative vindication of §10: the genotype negatives in
this atlas are a power ceiling, not a biological absence. It also converts §10's
closing recommendation from a hunch into a number — **~400 contexts for 80%
power, ~600 for 95%** — for effects at the top of the real effect-size
distribution. Weaker, more typical effects need more.

The caveat is honest and worth stating: PRISM's phenotype is viability, and a
compound whose *killing* is genotype-linked need not have genotype-linked
*transcription*. The power curve transfers between readouts only to the extent
that effect sizes do. What it establishes firmly is that n=47 cannot see effects
this size in any readout.

Figure bundle: `results/figures/13_prism/prism_context_genetics/`.
Tables: `prism_decomposition.csv`, `prism_genotype_scan.csv`,
`prism_genotype_power.csv`.

## 12. Which drugs transfer depends on what you measure

> **Superseded in part by §16 and §18.** The cross-dataset comparisons here were
> computed before the estimator correction (§14) and before cross-laboratory
> reproducibility was measured (§18). The readout-decoupling result stands; the
> transcription-to-transcription agreement does not.

The mechanism ranking (§5) was tested against LINCS phase 2 and looked shaky,
but that comparison shared only **six** mechanism classes with Tahoe — too few
to tell "different biology" from "too few points". Two better arms now exist.

**LINCS phase 1 (GSE92742) is consistent with it, but does not confirm it.** Phase 1 is an order of magnitude
larger than phase 2: 373,198 instances, **2,834 compounds, 70 lines, 1,960
replicate plates**, all estimable, 6.1M cross-replicate pairs. The architecture
replicates (47% shared / 53% interaction on 978 landmark genes), and the
mechanism ranking trends with Tahoe over ten shared classes:
**rho = +0.56, p = 0.093 — a positive trend that does not reach significance**,
and ten classes is still thin. MEK inhibitors (0.464) and glucocorticoid receptor
agonists (0.452, n=31) sit at the context-specific end in both.

**PRISM does not.** Against the same ranking, viability gives
**rho = −0.18 to −0.23 (n.s.)** versus Tahoe, and — with 67 shared classes, so
this is *not* an underpowered comparison — **rho = −0.09 (n.s.)** versus LINCS
phase 1. A consensus transcriptional rank built from Tahoe and LINCS correlates
with viability at rho = −0.19 (p = 0.13, n = 69).

| comparison | readouts | shared classes | rho |
|---|---|---|---|
| Tahoe vs LINCS-1 | transcription vs transcription | 10 | **+0.56** (p=0.09) |
| Tahoe vs PRISM | transcription vs viability | 11 | −0.18 (n.s.) |
| LINCS-1 vs PRISM | transcription vs viability | 67 | −0.09 (n.s.) |

The reading is that **mechanism sets transcriptional context-dependence, and
that ordering does not carry to viability**. Transcriptional rewiring and
differential killing are decoupled: a drug whose expression response is highly
line-specific is not thereby a drug whose killing is line-specific. This scopes
§5 — the claim is about transcriptional response, and should be stated that way
— and it is a substantive result in its own right, since the two readouts are
routinely used interchangeably as "drug response".

One confound is controlled: CDI is a ratio, so an inert compound has a
near-zero numerator *and* denominator and lands at 0 by default rather than
because it transfers. LINCS phase 1 is full of these (35 classes — antihistamines,
adrenergics — sit at exactly 0.000). Restricting to compounds with a detectable
shared effect changes none of the three correlations above.

Figure bundle: `results/figures/13_prism/three_platform_synthesis/`.

## 13. Context-dependence is dose-dependent, and peaks just below lethality

Every result so far pooled doses, treating transferability as one number per
compound. It is not — it is a curve, and the curve has a shape.

**PRISM (viability, 1,443 compounds, 8 doses, ~740 lines).** Within compound,
CDI rises with dose: median Spearman **rho = +0.33, 74% of compounds positive,
p = 4x10⁻⁹²**. The medians by dose position:

| dose position (median uM) | 1 (0.0006) | 2 (0.002) | 3 (0.010) | 4 (0.039) | 5 (0.16) | 6 (0.63) | 7 (2.5) | 8 (10) |
|---|---|---|---|---|---|---|---|---|
| CDI | 0.098 | 0.106 | 0.145 | 0.168 | 0.163 | 0.289 | **0.384** | 0.218 |
| mean effect (log2) | +0.31 | +0.30 | +0.26 | +0.23 | +0.28 | +0.13 | −0.25 | **−1.67** |

Context-dependence climbs roughly fourfold to a **peak at ~2.5 uM, then
collapses at the top dose** (0.391 vs 0.218 paired within compound,
p = 1x10⁻⁵⁹, n = 1,398). The collapse is a ceiling effect, and the mean-effect
row shows it directly: at 10 uM the average line has lost 1.67 log2 of viability
— everything is dying, so there is little left to differ about and the shared
component dominates.

**LINCS phase 1 (transcription, 298 compounds, 0.08–10 uM) shows the same rising
limb**: CDI **rho = +0.20, 60% positive, p = 1x10⁻⁵**, rising 0.455 → 0.505
across dose quartiles. Decomposed, the interaction grows faster than the shared
response (0.213 → 0.370, +74%, versus 0.269 → 0.366, +36%) — which is *why* the
ratio rises. LINCS shows no collapse, consistent with the collapse being a
consequence of mass killing rather than of dose as such.

**Tahoe cannot resolve this** and should not be quoted on it — a point since
strengthened: with the replicate plate restored (§14) the within-drug alignment
trend loses significance entirely (p = 0.11, was p = 9×10⁻⁹ when computed on the
plate-14-excluded build), confirming that the apparent Tahoe dose signal was an
artefact of the pairing rather than a weak real effect. Its three
concentrations (0.05 / 0.5 / 5 uM) produce almost identical shared-response
magnitude (0.047 / 0.048 / 0.054), so there is no gradient to read; splitting by
concentration also leaves too few cross-plate replicate pairs per (drug, dose)
to estimate the interaction at all.

Two consequences. **Scientifically**, the interaction is largest where lines
differ most in *where they sit on the dose-response curve* and vanishes once
dose saturates that difference — which is a mechanistic statement about what
the interaction is. **Practically**, a context-dependence index is only
comparable at matched dose; CDI values quoted across studies at different doses
are not the same quantity, and probe-panel and benchmark designs (§6) should
specify dose near the peak, not at the maximum.

Figure bundles: `results/figures/13_prism/dose_vs_context/`,
`results/figures/13_prism/lincs_dose_vs_context/`.

---

## 14. CORRECTED: Tahoe *can* measure the interaction — on the replicate plate everyone discards

> **This section replaces an earlier version that claimed the opposite.** That
> version reported 4.7% of conditions replicated and an undetectable interaction
> (0.0%, p = 0.97). Both figures were artefacts of our own pipeline, which
> inherited a `plate14` exclusion. The corrected analysis is below; the error and
> its cause are kept because the cause is the finding.

**Plate 14 is a deliberate biological replicate of plate 6.** The Tahoe authors
state it explicitly — it was included "to highlight the reproducibility of the
Mosaic platform" — and it carries 6.2M cells across 50 lines and 95 drugs, with
**100% of its 4,746 (line, drug, dose) triples also present on plate 6**. They
then excluded it from *training* and reserved it for validation, which is correct
for fitting a model and exactly wrong for measuring replicate agreement. Our
`build_deltas` inherited that exclusion as a default, and every downstream
estimate silently lost the atlas's only source of same-dose replicates.

**Replication counted from the released metadata**, with no processing of ours
(`scripts/verify_tahoe_replicates.py`, 100.6M cell records):

| filter | (line, drug, dose) triples | replicated | % |
|---|---|---|---|
| all cells | 56,879 | 7,714 | **13.56%** |
| atlas's own `pass_filter == full` | 56,877 | 7,691 | **13.52%** |
| ≥10 cells per plate | 56,395 | 7,366 | **13.06%** |
| *excluding plate 14* | 56,877 | 3,045 | *5.35%* |

The three inclusive filters agree, so 13.5% is a property of the experimental
design rather than of any threshold. Dropping plate 14 more than halves it.

**Measured on the designed replicates**, against a matched cross-context null:

| pairing | pairs | interaction share | *P* vs null |
|---|---|---|---|
| **true replicate** (same line, drug **and dose**) | 11,492 | **11.5% [10.5–12.5%]** | 1.7×10⁻⁷ |
| cross-dose (same line and drug, different dose) | 67,744 | 20.7% [20.5–20.9%] | 1.1×10⁻³² |
| pooled | 79,236 | 19.5% [19.3–19.8%] | 7.0×10⁻²⁸ |

So the interaction in Tahoe is real, well determined, and **11.5%** — while the
cross-dose pairing, which is what remains once plate 14 is dropped, inflates it
to **20.7%, nearly twice the true value**. Cross-dose covariance also decays with
dose separation (ρ = −0.020, p = 1.5×10⁻⁷), the attenuation expected if those
pairs are different conditions rather than repeats.

**What this means, and it is more useful than the claim it replaces.** The atlas
was designed with the replication needed to measure context-dependence. The
difficulty is that the replicate plate is routinely discarded — the atlas authors
withhold it from training, and any pipeline that copies that convention, as ours
did, is left with a cross-dose pairing that silently doubles the estimate. The
recommendation is therefore not "build a better atlas" but **"analyse the one you
have with the replicate plate in, and never treat doses as replicates"**.

Where genuine replicates are abundant the same estimator gives larger values —
LINCS phase 1, with 6.1M same-dose cross-plate pairs, gives 43% shared / 57%
interaction — but that is a different platform with different normalisation
(plate-wise z-scores rather than log-CPM deltas) and the two shares are not
directly comparable in magnitude.

Figure bundle: `results/figures/04_decomposition/tahoe_true_replicates/`.
Tables: `tahoe_replicate_structure.csv`, `tahoe_true_replicates.csv`.

## 15. No novel allele survives independent validation in GDSC

§11 recovers the known biomarkers; the question worth asking is what is left once
they are removed. Screening 111,589 allele × compound tests and setting aside
hits matching known pharmacogenomics (BRAF→RAF/MEK/EGFR, TP53→MDM2,
PIK3CA/PTEN→PI3K/AKT, KRAS→MEK and 25 further gene→mechanism rules) left 45
unknown hits, of which 4 were also allele-specific and passed lineage and MSI
controls. They were then tested in **GDSC** — a different laboratory, a different
assay (fitted IC50 rather than pooled barcode viability), 972 of 978 cell lines
mapped through Cellosaurus.

**The validation pipeline is demonstrably able to detect real effects.** Known
biomarkers were run first as positive controls and **11 of 11 replicated with the
correct sign**:

| control | drugs | result |
|---|---|---|
| BRAF p.V600E | dabrafenib, trametinib, selumetinib, PLX-4720, SB590885 | all 5, Δz −1.47 to −2.57, p 3×10⁻³¹ to 1×10⁻¹⁹ |
| PIK3CA p.E545K | alpelisib, pictilisib, taselisib | all 3, p 6×10⁻⁶ to 2×10⁻⁴ |
| PIK3CA p.H1047R | alpelisib, pictilisib, taselisib | all 3, p 1×10⁻⁶ to 4×10⁻² |

**None of the candidates survives.**

| candidate | outcome in GDSC |
|---|---|
| MUC16 p.T12415A × lapatinib | **does not replicate** (Δz −0.155, p = 0.15, n = 31) |
| HMCN1 p.K2374fs × EGFR class | 6 of 16 EGFR inhibitors raw, **0 of 16 after lineage adjustment** |
| MUC16 p.T12415A × EGFR class | 6 raw, 4 after lineage adjustment (chance ≈ 0.4) — but see below |
| ZMYM5 p.K463fs × mozavaptan | untestable; no vasopressin-receptor drug in GDSC |

**The one partial survivor is a germline polymorphism, not a somatic variant.**
CCLE calls variants without matched normals, so common SNPs pass into the
somatic table. Three signatures identify MUC16 p.T12415A as one:

1. All 40 carriers share a **single genomic position** (g.chr19:9021080T>C).
2. Carriers show **no excess mutational burden** (median 352 vs 358 atlas-wide),
   so it is not selected and not a burden proxy.
3. It heads a list of recurrent MUC16 variants whose other members are
   **synonymous** — p.T13039T, p.V13913V, p.R13041R. Synonymous changes recurring
   across unrelated lines are germline by definition, and no somatic hotspot
   keeps that company.

The same test flags HMCN1 p.K2374fs and ZMYM5 p.K463fs as germline too. An
association between a germline SNP and drug response is population structure —
ancestry correlates with cell-line provenance, which the coarse TCGA lineage
label only partly absorbs — not a somatic biomarker.

**Verdict: no credible novel allele-level discovery.** The validation is
informative rather than merely null, because the same pipeline recovers every
known biomarker at overwhelming significance in the same data.

The allele-resolution result that *is* solid is methodological: BRAF V600E is
significant for 11 compounds whose gene-level BRAF test fails, in both directions
(sensitising to RAF/MEK inhibitors, resistance-conferring to EGFR inhibitors),
because only 50 of 104 BRAF-mutant lines carry V600E and the rest dilute it.
That is a real argument that gene-level indicators lose associations, but it is a
statement about method that re-derives known biology. **The Nature Genetics route
through allele-level discovery is closed on this data.**

Two lessons are recorded rather than buried. The earlier split-sample validation
(94–99%) was **circular** — candidates were selected using every line, so
re-splitting those lines measures effect size, not generalisation. And the first
class-level test matched drug targets by substring, where `"her"` matches
`"other"`, inflating the EGFR class to 170 drugs including 5-fluorouracil;
word-boundary matching on `PUTATIVE_TARGET` gives the correct 16.

## 16. Mechanism ranking: significant within Tahoe, not reproducible across datasets

With every dataset on the same matched-null estimator, the mechanism claim of §5
narrows sharply.

| comparison | shared classes | rho |
|---|---|---|
| within Tahoe (Kruskal–Wallis across mechanisms) | 24 | H = 50.2, **p = 6.0×10⁻⁴** |
| Tahoe vs LINCS phase 1 | 10 | **+0.09** (p = 0.80) |
| Tahoe vs PRISM viability | 11 | −0.18 (n.s.) |
| LINCS phase 1 vs PRISM viability | 67 | −0.09 (n.s.) |

Mechanism structures context-dependence *within* Tahoe, robustly. It does not
transfer to another transcriptional platform, and it does not transfer to
viability. The earlier ρ = +0.56 reported against LINCS was produced by the
biased per-perturbation estimator and disappears once both sides use the same
corrected one — it should not be cited.

## 17. Genotype has no cross-validated predictive power; lineage has a little

§11 showed that at 738 cell lines the genotype scan recovers the clinical
biomarker set, and reported that a single allele explains a median 5.5% of a
compound's interaction. **That 5.5% figure was in-sample and is withdrawn** — it
was a point-biserial R² for an allele chosen on the same data, so selection
inflated it. Cross-validated, it is approximately zero.

The proper question is how much of the line × compound interaction *any* genotype
block predicts out of sample. Testing four blocks per compound with 5-fold
cross-validated ridge (alpha chosen inside each training fold), across 150
compounds at ~740 lines:

| predictor block | median CV R² | compounds with R² > 0 | p |
|---|---|---|---|
| **lineage** (primary tissue) | **+0.0145** | **72.7%** | 4.9×10⁻¹⁵ |
| mutational burden | −0.0023 | 32.0% | 1.4×10⁻² |
| nonsynonymous variants | −0.0058 | 33.3% | 7.3×10⁻⁷ |
| synonymous variants | −0.0064 | 22.0% | 6.9×10⁻¹³ |
| everything together | −0.0029 | 42.7% | 1.9×10⁻² |

A negative cross-validated R² means the block predicts *worse than the mean*.
**Genome-wide mutation status carries no generalisable information about which
compounds a cell line responds to unusually.** Lineage does, but explains only
about 1.5% of the interaction.

**The synonymous control isolates the mechanistic part.** A synonymous variant
cannot change a protein, so it cannot cause a response difference — but it
carries the same ancestry, lineage and germline-contamination structure as a
nonsynonymous variant in the same gene. Restricting both blocks to the 3,435
genes carrying each class in ≥20 lines makes them the same size and the same
genes, so only protein-changing capacity differs. The excess is:

- nonsynonymous − synonymous, marginal: **+0.0025** (57.3% positive, p = 0.028)
- the same, each added over lineage: **+0.0036** (58.0% positive, p = 0.025)

So a genuine mechanistic component exists and is statistically detectable, but it
is **~0.3% of the interaction variance**. Everything else that looked genetic is
population structure and lineage.

**One caveat runs in our favour.** Ben-David et al. (2018, Nature) show germline
variant calls are far more reproducible across laboratories than somatic ones
(allelic-fraction r = 0.95 vs 0.86; a median 19% of non-silent mutations appear
in only one of CCLE/GDSC). Synonymous variants are largely germline, so the
synonymous block carries *less* measurement error than the nonsynonymous somatic
block. The comparison is therefore conservative — the true mechanistic excess
could be somewhat larger than 0.3%, though not by enough to change the reading.

### Expression, not genotype, is the predictor — §17 as first written was incomplete

The blocks above omitted the one the field considers strongest. Adding CCLE
baseline expression (top 2,000 variable genes) and baseline protein (RPPA, 214
antibodies) under the identical cross-validated ridge, over 120 compounds:

| predictor block | median CV *R*² | compounds positive | p |
|---|---|---|---|
| **baseline expression** | **+0.0927** | **92.5%** | 2.7×10⁻²⁰ |
| **baseline protein (RPPA)** | **+0.0731** | 88.3% | 1.3×10⁻¹⁹ |
| lineage | +0.0201 | 74.2% | 8.2×10⁻¹⁵ |
| copy number | +0.0051 | 57.5% | 6.9×10⁻³ |
| nonsynonymous variants | +0.0002 | 50.8% | 0.74 |
| expression + lineage | +0.0929 | 92.5% | 2.8×10⁻²⁰ |
| copy number + expression | +0.0768 | 87.5% | 1.2×10⁻¹⁷ |

**Copy number confirms Schlüter & Schönhuth (2025) in relative terms and refutes
them in absolute ones.** On the 120 compounds where all blocks are measurable on
identical lines, copy number does beat mutations (+0.0064, p = 5.7×10⁻⁴) — their
claim holds — but it trails expression by a wide margin (−0.0880, p = 3.3×10⁻²¹).
Copy number also *adds nothing over expression*: the joint block (+0.0768) scores
**below** expression alone (+0.0927), the 2,000 extra largely uninformative
predictors diluting the fit even under cross-validation. The ordering is
expression > protein > lineage > copy number > mutations.

Baseline expression predicts the interaction **4.6× better than lineage**
(difference +0.068, p = 4.7×10⁻²⁰), protein nearly as well, and **lineage adds
nothing once expression is included** — it was a coarse proxy for expression
state all along. Genotype remains at zero.

The corrected statement is therefore not "nothing predicts the interaction" but
**"molecular state predicts it and genotype does not"**. The effect is modest in
absolute terms (≈9% of variance) but consistent, appearing in 92.5% of compounds,
The comparison with Schlüter & Schönhuth (2025) is set out above.

**Context.** In UK Biobank, genome-wide common variation explains 9–17% of drug
response (Sadowski et al. 2024, *Cell Genomics*: statin–LDL 9%, statin–A1c 10%,
statin–glucose 11%, metformin–BMI 17%). Our cell-line system yields far less.
Two explanations are live and we cannot separate them here: cell lines genuinely
have less genetic determination of response than patients, or the cross-dataset
genotype–phenotype mismatch Ben-David quantifies destroys it before we measure.

Figure bundle: `results/figures/15_architecture/genetic_architecture/`.

## 18. Cross-laboratory reproducibility, and why cell-line identity is the limiting step

Ben-David et al. (2018, *Nature*) showed that "the same" cell line differs between
laboratories — a median 19% of non-silent mutations appear in only one of
CCLE/GDSC, and 48 of 55 compounds active against one MCF7 strain were completely
inactive against another. They demonstrated it in a single cell line, at one
dose, with viability. Scaled up, it sets a hard ceiling on every context-transfer
model trained on one atlas: **of the line × compound interaction we attribute to
cell-line identity, how much reproduces elsewhere?**

The design separates laboratory divergence from assay noise, because a bare
cross-lab correlation conflates them:

| comparison | what varies | compounds | median Spearman r | % of ceiling |
|---|---|---|---|---|
| PRISM replicate plates | nothing | 1,435 | 0.473 | — |
| GDSC1 vs GDSC2 | experiment, assay version (both Sanger) | 123 | 0.438 | — |
| **PRISM vs GDSC** | **laboratory + assay** | 187 | **0.255** | **56%** |
| PRISM vs GDSC, **identity-validated** | as above, lines checked from data | 68 | **0.397** | **87%** |

**The ceiling is itself low**: even repeating the measurement within one
laboratory, the line's compound-specific deviation agrees at only r ≈ 0.45. Half
the "irreproducibility" usually blamed on other labs is present before you leave
the building.

### Identity, not protocol, is what breaks

Every cross-atlas comparison begins by asserting that a line in one atlas is the
same as a line in another, from an identifier that is then never checked. We
checked it: each line gets a **response fingerprint** — its residual across the
compounds both atlases share — and must be its own best match among all
candidates.

It usually is not. Of 488 testable COSMIC→DepMap-matched lines, **only 5–12% are
their own best match**, with the identifier-matched partner ranking a median 82nd
of 971 candidates. That rank is far better than chance (~486), so identity does
carry real information — it is simply not unique.

**Restricting to lines whose identity the data confirms raises cross-lab
agreement from 56% to 87% of the within-lab ceiling.** This is *not* circular:
identity is validated on one random half of the shared compounds and the
agreement is measured on the disjoint other half (all-lines control on the same
held-out compounds: 59%).

> **Interpretation corrected (2026-08-31).** We first read reciprocal-best-hit
> failures as evidence of cell-line divergence between laboratories, following
> Ben-David. **CCLE baseline expression does not support that reading.** For the
> 423 lines whose identifier match was outranked, the outranking line is a
> near-random expression neighbour — median expression rank **306 of 798**
> against a chance expectation of 399 — and only **7%** share the primary tissue.
> Genuine divergence would put the better match among the line's close relatives;
> a near-random partner instead indicates the best-hit is largely **fingerprint
> noise**. The identifier match remains far above tissue-matched random (0.235 vs
> 0.048), so identity carries real signal — it is simply not resolved sharply
> enough for "is it its own best match" to diagnose divergence.
>
> Unchanged: the held-out gain from 59% to 87% is real and predictive. Changed:
> its **mechanism**. Validation selects lines whose cross-laboratory signal is
> strong and self-consistent — a sound practical filter — but we cannot claim it
> specifically detects divergent cultures. The 5–12% reciprocal-best-hit rate is
> a property of the metric's resolution, not an estimate of how many cell lines
> are wrong. This section's title should be read as "cell-line identity is where
> the recoverable loss sits", not as a demonstration that cultures have diverged.

### What each matching strategy actually buys

Running the same comparison under a ladder of matching rules (`matching_strategy.py`):

| strategy | median fingerprint r | n |
|---|---|---|
| random pairing | −0.002 | 9,753 |
| same tissue, random line | 0.048 | 9,740 |
| **identifier (standard practice)** | **0.235** | 489 |
| identifier + fingerprint-validated | 0.403 | 61 |
| best available hit (upper bound) | 0.419 | 489 |

The last two rungs are selected on the outcome and are shown to bound the room
available, not as independent estimates — the honest, non-circular version of the
same gain is the held-out-compound result above (59% → 87%).

Identifier matching beats same-tissue pairing decisively (p = 1.5×10⁻¹²⁴), so
cell-line identity is carrying genuine information well beyond lineage. But the
best available partner reaches 0.419 — **nearly twice what the identifier finds**
— so a large amount of cross-atlas similarity is being left on the table by
trusting the label.

### Transfer scales steeply with signal strength

ρ = +0.54, p = 1.8×10⁻²⁹, n = 364:

| interaction strength quartile | cross-lab r | fraction of ceiling |
|---|---|---|
| Q1 (weakest) | 0.165 | 36% |
| Q2 | 0.230 | 51% |
| Q3 | 0.378 | 83% |
| Q4 (strongest) | 0.455 | **100%** |

The strongest quartile transfers at the within-lab ceiling — perfectly. Weak
compounds do not transfer at all, and modelling them across atlases is chasing
noise.

### It is a property of laboratories, not of the killing assay

The same analysis on transcription (`cross_lab_transcription.py`), with LINCS
phase 1 vs phase 2 as the within-lab control:

| comparison | median residual-profile r | pairs |
|---|---|---|
| LINCS p1 vs p2 (within-lab, Broad) | 0.061 | 5,803 |
| Tahoe vs LINCS (cross-lab) | 0.032 | 489 |
| **reproducible fraction** | **46%** |

The transcription figure was reported as 52% until the sciPlex3 arm was added
(§23), whose rows were then pooled into the headline by accident along with the
identity-validated rows the viability arm excludes. Rebuilt the same way as the
viability number — non-validated rows, Tahoe versus LINCS only — it is **46%**,
against viability's 56%. The audit in §25 now covers this quantity so the same
slip cannot recur silently. |

46% for transcription against 56% for viability, from different institutions and
wholly different assays. The number is a property of laboratories rather than of
a particular readout. The identity check behaves the same way: within the Broad
(LINCS p1 vs p2) **16 of 16** name-matched lines are their own best match; across
labs (Tahoe vs LINCS) only **3 of 6**.

A caveat specific to this arm: the transcriptional within-lab ceiling is itself
only r ≈ 0.06, so the line-specific transcriptional residual is barely
reproducible even within one laboratory, and the 46% is a ratio of two small
numbers.

### Consequences

This reframes two of our own results. The cross-dataset non-replication of the
mechanism ranking (§12, §16) needs no biological explanation — a ~56%
reproducible fraction concentrated in strong-interaction compounds suffices. And
the genotype negatives (§10–11, §17) gain a second cause alongside power: we join
CCLE genotype to PRISM and GDSC response, and for 88% of lines the two atlases
are not measuring behaviourally identical cultures.

**Practical recommendation.** Any cross-atlas analysis should validate cell-line
identity from the data rather than trusting the identifier, and should weight or
filter compounds by the strength of their line-specific component. Doing both
recovers agreement from 56% to 87% of what the assay can achieve with itself.

### Separating laboratory from assay

PRISM vs GDSC changes institution *and* assay at once, so on its own it cannot
apportion the loss. GDSC1 and GDSC2 are separate screening versions run at the
same institution over different concentration ranges (GDSC1 0.0004–16 µM, GDSC2
0.00001–8 µM) with different assay chemistry, so that comparison changes the
**assay with the laboratory held fixed**. The ladder therefore apportions it:

| step | r | cost |
|---|---|---|
| same lab, same assay, repeat plates | 0.473 | — |
| same lab, **different assay version** | 0.438 | −0.035 (**16%** of the drop) |
| **different laboratory** and assay | 0.255 | −0.183 (**84%** of the drop) |

**Changing the assay within one laboratory costs 16% of the total loss; changing
laboratory costs the remaining 84%.** The limitation is now partial rather than
total: this is two points on a ladder, not a factorial design, and a cross-lab
comparison with the assay held exactly fixed would settle it. CTRP would have
provided one, but the NCI CTD2 data portal has been retired (every historic URL
now redirects to `studycatalog.cancer.gov`), so it could not be obtained.

### Remaining limitations

*Identity verification.* Baseline expression is now available (legacy CCLE
distribution) and has been used — the interpretation note above is what it
changed. It remains a **single snapshot from one institution**, so it is an
independent molecular reference but cannot compare two laboratories' cultures
against each other; a per-laboratory molecular profile would be needed and does
not exist for these atlases. Fingerprinting also cannot separate culture divergence from
misidentification, and genuinely related lines are similar, so the 5–12%
reciprocal-best-hit rate is an **upper bound on the rate of true identity
problems** — the median rank of 82 of 971 shows identity is informative but not
unique.

*Transcriptional arm.* Its within-lab ceiling is only r ≈ 0.06, so the 46%
figure is a ratio of two small numbers. It agrees with the viability arm, which
is the substantive point, but should not be quoted precisely.

*Selection.* The identity-validated rungs in the matching ladder are selected on
the outcome and bound the available room rather than estimating it; the
non-circular value is the held-out-compound result (59% → 87%).

Figure bundles: `results/figures/16_crosslab/cross_lab_reproducibility/`,
`.../matching_strategy/`, `.../cross_lab_transcription/`.

## 19. The estimator recovers a known interaction share; the shortcuts do not

Every number in this document rests on one estimator, and until now it was
defended by argument plus the observation that the alternatives disagree with it.
Disagreement is not evidence of correctness. `scripts/simulate_estimator.py`
generates data whose interaction share is fixed by construction —

    d[c,p,q,r,w] = shared[p,q] + interaction[c,p]·a(q) + batch[r] + noise

over contexts *c*, perturbations *p*, doses *q*, replicate plates *r* and wells
*w* — and asks each estimator to return it.

**Recovery across the range** (mean of 3 seeds; 30 contexts, 20 perturbations,
3 doses, 2 plates × 2 wells, batch σ = 0.8, dose persistence 0.7):

| true share | 0.0 | 0.1 | 0.2 | 0.3 | 0.5 | mean bias |
|---|---|---|---|---|---|---|
| **recommended** | **0.011** | 0.086 | 0.161 | 0.235 | 0.381 | **−0.045** |
| pooled batch | **0.149** | 0.213 | 0.277 | 0.341 | 0.466 | +0.069 |
| in-sample prior | 0.001 | 0.070 | 0.142 | 0.214 | 0.358 | −0.063 |
| cross-dose as replicate | 0.005 | 0.043 | 0.086 | 0.131 | 0.233 | −0.120 |
| residual variance | **0.498** | 0.535 | 0.572 | 0.608 | 0.681 | +0.359 |

The decisive column is the first. **On data containing no interaction at all,
residual variance reports 50% and pooled batch reports 15%** — both are
manufacturing the quantity they claim to measure. Only the recommended estimator
returns approximately zero.

**Sensitivity to noise** (true share fixed at 0.2), which no estimator of a
*reproducible* component should show:

| noise σ | 0.5 | 1.0 | 2.0 | 4.0 |
|---|---|---|---|---|
| residual variance | 0.384 | 0.587 | 0.819 | **0.939** |
| recommended | 0.187 | 0.193 | 0.214 | 0.280 |

Residual variance tracks noise almost perfectly — at σ = 4 it reports 94%
interaction from data that is 20%. The recommended estimator is flat to within a
few points up to σ = 2 and inflates modestly beyond it.

**Two honest caveats, because the validation produces them.**

*The recommended estimator is conservative.* Its recovery is highly linear
(R² = 0.999) but the slope is **0.74**, so it underestimates by roughly a
quarter in this regime. Reported shares should be read as **lower bounds**,
which is the safe direction but should not be glossed.

*It is stable near Tahoe's context count and drifts outside it.* At a true share
of 0.2 the estimate is 0.241 at 10 contexts, 0.204 at 20, 0.192 at 47 and 0.182
at 100 — accurate in the 20–50 range that these atlases occupy, mildly optimistic
below it and mildly conservative above.

Figure bundle: `results/figures/00_manuscript/estimator_simulation/`.

## 20. Confidence intervals on the headline quantities

The pooled shares and the cross-laboratory quantities were reported as point
estimates. Bootstrapping them — over cross-replicate pairs for the decomposition
(reservoir-sampled, 200,000 pairs), over compounds for everything else — gives:

| quantity | estimate | 95% interval |
|---|---|---|
| LINCS phase 1 interaction share | 57% | 56.9–57.4% |
| LINCS phase 2 interaction share | 47% | 46.4–47.1% |
| PRISM interaction share | 22.1% | 20.2–23.5% |
| Tahoe, true replicates | 11.5% | 10.5–12.5% |
| Tahoe, cross-dose pairing | 20.7% | 20.5–20.9% |
| cross-lab reproducible fraction, all lines | 56.1% | 48.2–60.4% |
| cross-lab reproducible fraction, identity-validated | 87.3% | **69.5–102.9%** |
| **share of the loss due to laboratory** | 84% | **72–106%** |

The decomposition shares are tightly determined; their uncertainty is systematic,
not statistical. **The two derived cross-laboratory quantities are not**, and
both intervals cross or approach 100%, which changes how they should be stated:

- After identity validation, cross-laboratory agreement is **statistically
  indistinguishable from the within-laboratory ceiling** (87.3%, CI 69.5–102.9%).
  That is a cleaner claim than "87%".
- The laboratory accounts for 84% of the loss with a CI of 72–106%, so **the
  assay contribution is not distinguishable from zero**. The direction of the
  claim strengthens — protocol is not the problem — while the precise 84/16 split
  should not be quoted as if resolved.

A sampling bug is recorded here because it produced a plausible wrong answer: the
first version kept the *first* 200,000 cross-replicate pairs for the bootstrap,
but pairs arrive in group order, so the prefix covered a handful of contexts
rather than a sample of them, and the resulting interval (55.1–55.8%) did not
even contain the point estimate it was meant to bracket. Reservoir sampling fixed
it.

## 21. Against a standard mixed model, and why a tool is still warranted

The obvious objection to `pertdecomp` is that variance partitioning is solved:
fit `y ~ (1|perturbation) + (1|context) + (1|context:perturbation)` and read the
components off, as variancePartition, lme4 or GxEMM would. Fitting exactly that
with `statsmodels` on the simulated data of §19, where the truth is known:

| regime | true 0.1 | true 0.3 | mean bias, ours | mean bias, mixed model |
|---|---|---|---|---|
| replicated | 0.150 / 0.157 | 0.338 / 0.298 | +0.044 | +0.028 |
| batch-confounded | 0.106 / 0.154 | 0.261 / 0.357 | **−0.016** | +0.055 |
| unreplicated | **refuses** / 0.161 | **refuses** / 0.182 | — | −0.028 |

*(cells show ours / mixed model)*

Three things follow, and the first is not in our favour.

**Where both are identifiable they agree**, to 0.051 on average. The mixed model
is a perfectly valid alternative with replicates, and our estimator is not a
reinvention of it — it is the same quantity in closed form. We say so.

**Without replicates the mixed model returns a number anyway, and it is the same
number regardless of the truth**: 0.161 when the true share is 0.1, and 0.182
when it is 0.3. Context×perturbation and residual enter the likelihood
identically when each condition is measured once, so the split is set by the
optimiser rather than the data — and nothing in the output says so. This is the
practical case for a tool that **refuses**, and it is precisely the regime a
Tahoe-like atlas falls into once its replicate plate is dropped (§14).

**With batch structure ours is the less biased of the two** (−0.016 versus
+0.055), because same-batch pairs are excluded by construction rather than
requiring the user to specify a plate term.

Cost is a secondary but real consideration: 0.1 s for the closed form across all
genes against 0.18 s *per gene* for the mixed model, about 6 minutes for a
2,000-gene panel and proportionally worse for the 978-gene LINCS panels run
across four datasets.

Figure bundle: `results/figures/00_manuscript/variance_components/`.

## 22. Testing the dose explanation — supported away from the ceiling, not at it

§13 explained the dose curve by asserting that the interaction is largest where
lines differ most in *where they sit on the dose–response curve*. That was a
just-so story until tested. It predicts that the interaction should track the
between-line spread of response, dose by dose, within each compound.

| test | result |
|---|---|
| interaction vs between-line spread, within compound | median ρ = **+0.405**, 79% positive, p = 4×10⁻¹³⁵ (n = 1,442) |
| restricted to **non-saturating** doses (positions 1–6) | median ρ = **+0.486**, 82% positive, p = 6×10⁻¹⁴⁰ |
| dose position of peak interaction vs peak spread | ρ = +0.211, p = 6×10⁻¹⁶; identical for only **14%**, within one position for 45% |
| interaction peaks below the dose of maximal killing | **81%** of compounds |

**The explanation is supported, but not in the simple form it was stated.** The
interaction tracks the between-line spread strongly, and *more* strongly once the
saturating doses are removed (ρ = +0.49 versus +0.41) — which rules out the
competing account that the whole pattern is a ceiling artefact of the viability
assay, since a ceiling artefact cannot operate where there is no ceiling.

But the two peaks do **not** coincide: the spread peaks at the top dose (8/8)
while the interaction peaks at 6/8. So at saturating doses the lines still differ
in measured viability — the spread is real — yet that difference **stops
reproducing between replicate plates**. The refinement is therefore that
context-dependence tracks *reproducible* between-line spread, and at the lethal
dose the spread survives while its reproducibility does not.

Figure bundle: `results/figures/13_prism/dose_mechanism/`.

## 23. Remaining gaps closed, and one that could not be

### The identity criterion is not a threshold artefact (gap 5)

The 87% figure retains only 25 of 488 lines, so it could have been an artefact of
the arbitrary "must be its own single best match" cut. Relaxing that cut, with
identity still selected on one half of the compounds and scored on the disjoint
other half:

| criterion | lines retained | cross-lab r | % of ceiling |
|---|---|---|---|
| rank 1 (best match) | 25 | 0.397 | **87%** |
| top 2% | 139 | 0.337 | 74% |
| top 5% | 198 | 0.334 | 73% |
| top 10% | 262 | 0.312 | 69% |
| top 25% | 367 | 0.299 | 66% |
| all identifier-matched | 488 | 0.270 | 59% |

The relationship is **perfectly monotone** (ρ = −1.00 against stringency rank).
The criterion is therefore measuring something real and the cut is a
precision-versus-coverage trade-off, not a lucky threshold. It also yields a more
practical operating point than the one we first reported: **"top 5%" keeps 198
lines — eight times as many — and still recovers 73% of the ceiling.**

### What baseline expression actually uses (gap 8)

§17 shows expression predicts the interaction where genotype does not, but not
what carries the signal. Fitting ridge weights per compound and pooling across
100 compounds:

- **The weight is diffuse.** 751 of 2,000 genes carry half the total weight and
  1,680 carry 90% — 38% and 84% of the panel.
- **No gene set survives correction.** Across 5,457 GO, Reactome and Hallmark
  sets, the strongest enrichment among the 200 top-weighted genes is
  p = 0.009, q = 1.0. Ion transport and redox terms appear at the top but none
  is significant.

So the predictor is **tracking global transcriptional state, not a specific
programme**. That is a negative result about mechanism and a positive one about
interpretation: it explains why no compact line-level descriptor has ever worked
(§4, §10), since the predictive information is spread across the transcriptome
rather than concentrated in features one could name or measure cheaply.

### A third laboratory: attempted and not achieved (gap 6)

All cross-laboratory results rest on two institutions. Two routes to a third
were tried.

**NCI-60** could not be obtained: CellMiner's archive endpoints return empty or
403 responses without a browser session, and DTP's bulk-data page no longer
exposes the GI50 files. CTRP is likewise gone, its NCI portal retired.

**sciPlex3** (Srivatsan et al., Trapnell laboratory) is a genuine third
institution and third technology, and shares A549, K562 and MCF7 with LINCS. It
returns median r = −0.055 (p1) and −0.063 (p2) — but this is **not evidence about
laboratories** and should not be read as such. sciPlex3 has only **three cell
lines**, so its leave-one-context-out shared response is the mean of *two* lines,
while LINCS's is the mean of sixty-nine. The two residuals are not the same
quantity, and §19 independently shows the estimator drifting outside the 20–50
context range. The arm is ill-posed rather than informative.

**Gap 6 therefore remains open.** A third laboratory arm needs a dataset with
both an independent institution and enough contexts, and no public dataset we
could reach satisfies both. This is the single most valuable addition new data
could make, and it is why the laboratory-versus-assay apportionment is reported
with an interval that reaches 100% (§20) rather than as a resolved split.

Figure bundle: `results/figures/15_architecture/robustness_interpretation/`.

## 24. Benchmark context count moves apparent model gains by ~20 points

§2 showed the additive baseline improves as it averages more contexts. Extending
that sweep to the context counts published benchmarks actually use makes it
directly consequential for how model gains are read.

| contexts averaged | baseline r_de100 | handicap vs a 45-context baseline |
|---|---|---|
| 2 | 0.550 | **21.1%** |
| 3 | 0.592 | 12.4% |
| 4 | 0.606 | 9.8% |
| 5 | 0.618 | 7.7% |
| 6 | 0.625 | 6.4% |
| 10 | 0.643 | 3.6% |
| 18 | 0.654 | **1.7%** |
| 45 | 0.666 | 0% |

A model that is **exactly as good in every dataset** would appear ~21% stronger
against a baseline built from 2 contexts than against one built from 18, purely
because the baseline is worse. Nothing about the model changes.

**This bears directly on State** (Adduri et al., *Cell* 2026), the current state
of the art. Its zero-shot arm holds out one context at a time across five query
datasets with 3, 3, 4, 6 and 18 contexts, so its perturbation-mean baseline is
estimated from 2, 2, 3, 5 and 17 contexts respectively. The paper reports "more
than 19% improvement" on the 18-context dataset and "several-fold improvements"
on the smallest ones — the same ordering the table above produces from baseline
quality alone.

**This is a partial confound and is stated as such.** The handicap differential
between a 2-context and an 18-context baseline is about 19 points, which is real
and runs in the reported direction, but it does not account for several-fold
gaps. State is doing something the baseline is not. What follows is narrower and
still useful: **a portion of the spread in reported improvements across
benchmarks reflects how many contexts each benchmark's baseline was estimated
from, not how much better the model is there**, and gains are only comparable
across benchmarks at matched context counts.

The concrete recommendation is cheap: **benchmarks should report the number of
contexts behind their mean baseline**, alongside the gain. That single number
moves apparent performance by up to a fifth.

See `docs/paper_state_2026.md` for the full reading, including where State's
setting corroborates our few-shot result (§6) rather than competing with it.

## 25. Two audits: are the numbers current, and did we fall for the same trap twice?

### Every headline number regenerates from current code

After this many corrections — the estimator, the matched null, the
leave-one-context-out prior, the plate-14 exclusion, the withdrawn mechanism
claim, the reservoir-sampled bootstrap — the live risk is prose quoting a number
computed under a superseded code path. Prose does not recompute itself.
`scripts/audit_numbers.py` recomputes 28 headline quantities from the tables the
analyses wrote and checks each appears in RESULTS.md or MANUSCRIPT.md.

**28 of 28 now verify.** The first run flagged one: the baseline at 18 contexts
was written as 0.655 where the table gives 0.654 — a rounding slip, not a stale
result, and now corrected. The audit is cheap and should be run before any
submission and after any estimator change.

### The replicate-discarding trap: checked in LINCS and PRISM too

Having been caught once by Tahoe's plate 14 (§14), the obvious question is
whether the other datasets carry an analogous convention that silently removes
replicate structure.

**LINCS: the trap exists and we had already avoided it.** Level 5 signatures
aggregate replicates away — only ~5% of Level 5 conditions retain more than one
signature — which is why this project uses Level 4, where each well is a separate
profile and plates provide the replicate axis. Neither phase's `inst_info`
carries a QC flag that would drop replicates; the gold/quality flags live in
`sig_info`, which applies to Level 5 and is not used here.

**PRISM: no trap, and an unexpected finding.** PRISM ships
`passed_str_profiling` for every cell line, and 248 of 1,120 (22%) fail it.
Short-tandem-repeat profiling is the field's standard molecular test of cell-line
identity, and we had stated in §23 that no orthogonal identity check was
available across these datasets. That was wrong — but the check turns out to be
unusable for a different reason: **all 738 lines with response data pass STR
profiling.** PRISM excludes failures from the released screen entirely, so there
is no contrast to test.

That absence is itself informative, and it is not circular. Every cell line in
the cross-laboratory comparison has been **STR-authenticated by the Broad**, so
outright misidentification is ruled out as the explanation for the
reciprocal-best-hit failures of §18. Combined with the expression evidence in
§23 — the outranking line is a near-random expression neighbour — the remaining
explanations are metric noise and genuine culture drift *within* correctly
labelled lines, not mislabelling. This independently supports the correction
already made rather than reopening it.

Tables: `number_audit.csv`, `str_identity_check.csv`.

## 26. Audit repairs: what changed and what it cost

Three adversarial audits (statistics, data handling, claims-versus-evidence) were
commissioned against this work and returned roughly sixty defects. The repairs
are recorded here because several change published numbers, and two were bugs in
code this project releases.

### The estimator had two biases with one root

Both members of a replicate pair subtracted the **same** leave-one-context-out
prior, so that prior's estimation noise entered their covariance as signal; and
the denominator `mean(P²)` estimated `var(shared)` **plus** noise/n_out **plus**
batch/n_rep rather than `var(shared)` alone. Consequences, measured on simulated
data with a known share:

- the recovery slope was **0.741**, and the attenuation was not a fixed
  "conservative" factor but ranged from **0.47 to 1.01** with nuisance variance,
  so every published share was scaled by an unknown dataset-specific amount;
- at a true share of **zero** the estimator returned **0.036**, not zero — the
  claim in §19 that only this estimator returns ~0 on null data was false, and it
  looked true only at the single noise level and context count tested.

Both are fixed by a **split prior**: two disjoint leave-one-context-out priors per
condition, pairing residual A with residual B, and estimating `var(shared)` as
⟨P_A, P_B⟩. Recovery slope **0.911**, and **exactly 0.000** at a true share of
zero. The library was also pairing across doses — the error this paper is built
on criticising — and now groups by dose.

The simulation itself had to be fixed first: it applied a plate offset with no
control subtraction, contributing −(b₁−b₂)²/4 to every cross-plate pair, which
has no counterpart in plate-matched data. My first attempt at the split prior
looked *worse* because of it.

Two biases remain, stated rather than hidden: the estimate still drifts upward at
extreme noise, and **overestimates at low context counts** (0.294 against a true
0.200 at n=10), which makes the OP3 (6 contexts) and sciPlex3 (3) arms
unreliable.

### `decompose()` crashed on every real dataset

`find_replicate_batches` was passed the caller's dataframe while hard-coding the
internal column names, so the module's own documented example raised `KeyError`.
It was introduced by the commit that *added* the feature together with two tests
that could not catch it, because every test in the suite names its columns
`ctx`/`pert`. A regression test now renames them.

### A cell line that did not exist

Eight PRISM rows are `pool_line_FAILED_STR`. `split("_")[-1]` returned the
literal string `"STR"` for all eight, and `groupby.mean()` averaged eight
different cell lines into one — a chimera carrying data in 32,230 of 36,076
profiles that was **published in our own table**
(`cross_lab_identity_viability.csv`). Fixed in eight scripts by parsing the ACH
identifier and dropping STR failures rather than averaging an unauthenticated
culture into an authenticated one; all eight lines also appear as clean rows, so
nothing is lost. **The count is 737 lines, not 738.**

### Three claims that do not survive proper accounting

| claim | why it fails | corrected |
|---|---|---|
| laboratory accounts for 84% of the cross-lab loss | the three rungs were medians over 1,435 / 123 / 187 **different** compounds, while agreement depends strongly on the compound | on the 75 compounds all three share: **assay 68%, laboratory 32%** |
| identity validation recovers 87% of the ceiling | numerator on 25 selected lines, denominator on all 488 — a mismatch that returns 182% in a simulation where every line is identical | **71%** against a matched ceiling; a reliability-selected control gives 59%, so the effect is real but smaller |
| copy number beats mutations (p = 5.7×10⁻⁴) | 150 compounds treated as independent; their residuals correlate at mean r = +0.23, effective n ≈ 11 | cluster bootstrap CI **[−0.002, +0.015]** — withdrawn |

### Two claims that survive, one in weaker form

The **dose-mechanism Test 1** correlated a quantity with itself: `spread` was the
SD of the replicate mean, so `spread² = interaction + σ²/k` algebraically and the
two agreed at r = +0.98 by construction. Recomputed with spread measured on a
**disjoint replicate**, the relationship holds at ρ = **+0.317** (74% positive,
p = 4×10⁻⁹³) against the circular +0.405. But its companion tests do not: against
a label-permutation null preserving both marginals, "peaks identical" is **at
chance** (14.3% vs 15.2%) and "peak below maximal killing" clears its null by
**1.8 points** (80.4% vs 78.6%). Neither should be cited as support.

The **identity effect** survives its control. Selecting the same number of lines
by within-PRISM replicate agreement alone — using no cross-laboratory information
— gives 59.2%, indistinguishable from all lines (59.3%), while identity selection
gives 70.9%. So the criterion is not merely picking well-measured lines.

### A false claim withdrawn

§25 asserted that "all 738 lines with response data pass STR profiling". The join
behind it was degenerate: every one of the 248 STR-failing rows has
`depmap_id = NaN`, so the lookup returned 747 True and 1 False and the test could
not fail. Parsed correctly from `row_name`, **10 lines carry FAILED_STR
records**. Of the nine that appear in the identity table, none is its own best
fingerprint match — too few to test formally, but pointing the opposite way to
the claim it replaced.

---

## 27. The cell property that was being counted as a cell–drug relation

Every residual in this project was built by subtracting the leave-one-line-out
compound mean. That removes the drug's average effect and nothing else, so what
remained was **two** things added together: the line's *general sensitivity* —
how it responds to every compound, set by growth rate, seeding density and drug
metabolism — and the interaction specific to the pairing. Only the second is a
cell–drug relation.

**The general term is large and almost perfectly reproducible.** A line's mean
residual computed on one half of the compounds agrees with the value computed on
a disjoint half at **ρ = +0.984** (permutation null 0.000 [−0.070, +0.078],
p = 0.005), and **570 of 736 lines (77.4%)** carry a reproducible line-level
shift at FDR < 0.05. It is not noise, and because it is identical across
replicate detection plates it survives every noise-cancelling device in the
estimator and lands directly in the pair covariance.

Ben-David et al. (Nature 2018) name the same axis independently: their most
resistant MCF7 strains are resistant *in general*, through downregulated
drug-metabolism pathways, not through anything drug-specific.

### The three-way decomposition

    y(line, drug, dose) = β(drug, dose) + α(line) + γ(line, drug) + noise

with both main effects estimated out of fold — β from **other lines**, α from
**other compounds** — so neither can absorb the term it is meant to leave alone.
α is additionally estimated **within a replicate plate**: a log-fold-change is
taken against controls on its own plate, and an α pooled across plates carries a
share of each plate's control noise and subtracts it from the wrong side of a
cross-plate pair. In simulation that error removed var(control)/2 from the
covariance and drove the estimate to exactly **0.000**.

| component | share | variance |
|---|---:|---:|
| drug effect (identical in every line) | 94.1% | 1.6338 |
| **cell property** (general sensitivity) | **2.0%** | 0.0343 |
| **cell–drug relation** (interaction) | **3.9%** | 0.0683 |

PRISM, 737 lines × 1,324 compounds. Each component is measured as a covariance
between two independent estimates, so noise contributes zero.

**Consequence: prior interaction numbers were ~1.5× too large** (0.1025 counted
where 0.0683 belongs).

### What this does to the two central claims

**"A drug×cell relation rather than a cell property" — overstated, rewritten.**
A cell property plainly exists. The relation is **2.0× larger** than it, which
supports *"the larger of the two"* and not *"rather than."* Asking the
line-consistency question of the corrected residual returns ρ = −0.080, which
only confirms the subtraction ran; the claim is settled by the two variances
above, not by that test.

**"Transcriptional state rather than genotype" — survives, and strengthens.**
Re-run on the corrected residual, with a per-compound permutation null (cell-line
labels shuffled within compound, ridge kernel factorised once per fold and reused
across permutations) and FDR across 120 compounds:

| predictor | compounds beating their own null | median CV R² | median null |
|---|---:|---:|---:|
| baseline expression | **119 / 120 (99.2%)** | +0.0778 | −0.1279 |
| lineage | 104 / 120 (86.7%) | +0.0257 | −0.0098 |
| mutations | 22 / 120 (18.3%) | −0.0413 | −0.0754 |

97 compounds are expression-significant and mutation-non-significant; **zero** go
the other way (binomial p = 1.26 × 10⁻²⁹). Compound-cluster bootstrap over 150
compounds: expression **+0.0814 [+0.0716, +0.0925]**, lineage +0.0175
[+0.0135, +0.0205], mutations +0.0018 [−0.0009, +0.0041] — the mutation interval
includes zero.

Critically, expression predicts the **relation** (+0.0814) *better* than it
predicts the **property** (+0.0679, same lines, same folds), so the result is not
general sensitivity in disguise. Before the correction the same test gave
107/120; the correction raised it to 119/120.

### Scripts, guards and scope

`src/perturbmodel/celldrug.py` holds the decomposition and a pair query
(`--cell MCF7 --drug dabrafenib`), and `prism_gamma()` is the single cached
source of corrected residuals for every downstream script.
`tests/test_celldrug.py` (9 tests) plants a known general-sensitivity term and a
per-plate control offset and asserts neither is reported as interaction — the
per-plate case is the one that silently returned zero while looking correct in
review. The estimator simulation now sweeps a planted general response: at 0,
0.5 and 1.0 SD the corrected estimator returns **0.211, 0.212, 0.211** against a
truth of 0.20, where the uncorrected one returns 0.170, 0.347, 0.595.

Sixteen scripts and three library modules were changed. Two reporting defects
surfaced during the rerun and were fixed with it: the cross-laboratory
"reproducible fraction, identity-validated" was dividing by an all-lines ceiling
and read **110%** — impossible for a fraction of a ceiling — and is 98.0%
[47.0%, 116.2%] against the matched ceiling; and both the context-count and
general-sensitivity sweeps in the simulator were labelling the default estimator
"recommended" without passing `split_prior=True`, so neither had been testing the
recommended estimator at all.

---

## 28. Does the same omission inflate other groups' published numbers? No — a negative result

§27 showed our interaction shares were inflated ~1.5x by a cell property left
inside the residual. That is a defect in our code. It becomes a *finding* only if
the same term sits inside statistics other groups have published. It does not,
and this section records the attempt and its failure, because the failure is
informative and because the opposite claim would have been easy to make.

Three published context-specificity statistics were recomputed on the publishing
group's own data, under a strict rule: **reproduce the published number first, by
the paper's own definition.** A correction applied to a statistic we cannot
reproduce measures our reimplementation, not their result.

| study | statistic | published | our replication | corrected | moved |
|---|---|---:|---:|---:|---:|
| Srivatsan 2020 (sci-Plex 3) | % of responsive genes that are cell-type-dependent | 48% | **74.4%** | 73.9% | 1.01x |
| Subramanian 2017 (CMap/LINCS) | % of compounds with a panel-conserved signature | 26% | **5.0%** | 4.3% | 1.15x, wrong direction |
| Ben-David 2018 (rule applied to PRISM) | % of active compounds inactive in ≥1 context | 87% | 94.7% | 90.4% | 1.05x |

**Neither published statistic reproduced** — 74.4% against 48%, and 5.0% against
26%. Our operationalisations are not theirs: sci-Plex used a Monocle
likelihood-ratio DE test on cells, and Subramanian used Level 5 moderated
z-scores and a transcriptional activity score, not Level 4 signature
correlation. So those two arms are uninformative about those papers, and **no
inflation claim is made for either**. The Ben-David row is their counting rule
applied to PRISM lines, not their 27 MCF7 strains, and is labelled an analogue.

**Even within our own implementations the correction barely moves anything** —
1.01x to 1.15x, against 1.5x for our own interaction share. So the effect is not
merely unproven for other groups; it is absent from these statistics as we
computed them.

### Why, and the diagnostic that also failed

The natural explanation is that these pipelines already remove the context main
effect — LINCS Level 4 z-scores within plate, sci-Plex takes deltas against
vehicle in the same line. A diagnostic was added to test that: the context main
effect as a share of the residual each pipeline constructs. It **does not predict
the effect** — across the three it is anti-ordered with it (sci-Plex has the
largest share, 8.3%, and the smallest movement; LINCS the smallest, 1.4%, and the
largest).

The reason it fails is the real lesson. That diagnostic measures the context term
against *total* residual variance, which is dominated by measurement noise: in
PRISM it returns **2.8%**, while the replicate-covariance decomposition of §27
puts the same term at **~33% of the reproducible residual** (0.0343 of
0.0343 + 0.0683). Our estimator compares components that noise cannot inflate,
because each is a covariance between two independent estimates. A statistic
computed on a raw residual is mostly measuring noise, so removing a term worth a
few percent of it changes little.

**Prediction, since tested — see §29.** The inflation should appear in
*replicate-validated* context-specificity statistics and not in raw-residual
ones. TRADE's deposit was obtained and the prediction **holds**: on the same
data, the omission moves the noise-corrected split by 11.7 points and the raw one
by 3.4. §29 also shows why the diagnostic in this section fails — it was measured
without centring the correction across contexts, which is what makes it a context
property rather than a shared programme.

**Consequence for the manuscript.** The defect is reported as ours, quantified,
with the tool that prevents it and the simulation that shows what it costs. The
claim that the field's published numbers are inflated is **not** made.

---

## 29. The prediction tested on another group's data: it holds

§28 could not show the inflation in other groups' published numbers and proposed
one reason: those statistics are computed on raw, noise-dominated residuals,
whereas ours compares components noise cannot inflate. That was a prediction with
no data behind it. Nadig, Replogle, Pogson et al. (*Nature Genetics* 57:1228,
2025) supply the missing case — their TRADE deposit ships **lfcSE beside every
log2 fold change** for four cell lines given the same 2,052 essential-gene
perturbations (K562, RPE1, HepG2, Jurkat; 6,642 shared genes). The standard
errors make the noise-corrected statistic constructible, and four cell lines make
the consistent component estimable as a between-line covariance, which noise
cannot inflate because two lines are measured independently.

**The correction must be a contrast, and that is easy to get wrong.** A first
version subtracted each line's uncentred mean response across the other
perturbations. In this dataset every perturbation is an essential-gene knockdown,
which drives a common growth-arrest programme, so those four vectors correlate
**with each other at r = +0.51** and with the cross-line shared response at
**r = +0.80** — one shared programme, not four line properties. Subtracting it
removed signal the statistic is meant to keep, and the split moved the wrong way.
Centred across lines, as `perturbmodel.celldrug` does, the same vectors correlate
at **r = −0.53** between lines and −0.19 with the shared response, and behave as
line properties should.

| statistic | line effect left in | line effect removed | shift |
|---|---:|---:|---:|
| raw (noise included) | 61.5% cell-type-dependent | 58.1% | **+3.4 points** |
| **noise-corrected** (covariance + published SEs) | 66.3% | **54.6%** | **+11.7 points** |

The noise-corrected split moves **3.4× more**, and by 11.7 points — a 1.21×
inflation of the cell-type-dependent share. Sampling noise is 38% of the raw
per-line variance here, which is why the raw statistic is comparatively inert.

**What this does and does not license.** Our raw replication gives 39% consistent
against TRADE's published 56%, so it is **not reproduced** and nothing here
restates or corrects their number. What it establishes is that the effect is not
peculiar to us or to viability: on another group's data, in a different
perturbation modality (CRISPR rather than chemical), across four cell lines, a
context-specificity statistic of the same form is inflated by the same omission
whenever it is noise-corrected. That is a statement about the method.

---

## 30. Does drug exposure remodel the expression programme, or move cells along it?

"The drug changes gene expression" is not a finding — at 5 µM almost anything
does. The question with content is whether exposure moves cells along directions
they already had, or into directions no untreated cell occupies.

The baseline programme space is the span of principal components of **untreated**
Tahoe pseudobulk across cell lines. Each response is split into the part inside
that span and the part outside. Every term is an inner product between deltas
measured on **two independent plates** (plate 14 duplicates plate 6), because
noise is close to isotropic and lands outside any low-dimensional span — a purely
noisy response would otherwise read as pure remodelling. Measured from single
observations the shift share is 8.5%; measured across plates it is 10.1%.

**The reference that makes the number readable** is a held-out untreated cell
line: a signal that is a shift by construction, since it *is* line-to-line
variation. Over 5,345 replicated conditions, with 36 lines fitting the span and
12 held out:

| baseline directions used | captures of a drug response | captures of a held-out untreated line |
|---:|---:|---:|
| 1 | 1.7% | 9.7% |
| 5 | 6.2% | 22.1% |
| 10 | 9.0% | 25.0% |
| 20 | 12.2% | 29.6% |
| 35 | 15.7% | 32.8% |

Drug responses are captured **two to five times less** than a held-out untreated
line at every k. **Drug exposure does not move cells along the axes that separate
untreated lines — it adds a new one.**

Three controls decide what kind of remodelling this is:

- **Not death.** The outside-span share is 89%, 86% and 90% at 0.05, 0.5 and
  5 µM — flat, and already at its full size at the lowest dose, where cells are
  not dying. Cytotoxic collapse would appear only at the top dose.
- **Cell identity is preserved.** Between-line spread, treated against the same
  lines untreated, has ratio 1.002 / 1.000 / 1.007 across the three doses (all
  *P* > 0.5). Lines stay exactly as far apart as they were; the drug does not
  overwrite what makes them different.
- **Coordination is untouched.** The gene–gene correlation structure across lines
  is essentially unchanged by treatment (ρ = +0.998 over 400 genes; mean absolute
  change 0.0127 against a label-permutation null of 0.150 [0.128, 0.173]) — the
  observed change is far *smaller* than any random split of the same samples.

**So: yes, but additively.** Drug exposure writes a new, drug-specific axis that
baseline variation does not contain, while leaving cell identity and the
co-expression architecture intact. It is not reprogramming and not convergence
onto a common stress state — it is a new programme laid on top of an unchanged
one.

This is also the geometric reason the decomposition in §27 works so cleanly. The
drug axis is close to orthogonal to the axes on which cells differ, which is why
a drug main effect can account for 94.1% of the variance while the cell property
and the cell–drug relation stay small and separable.

---

## 31. Gap sweep: what the corrected estimator does to the rest of the paper

§27 corrected the estimator; §28–30 tested how far the correction generalises. A
systematic sweep of every claim against the current tables found four analyses
that had never been re-run under it, and three claims that were wrong
independently of it.

### Analyses that had not been re-run

| analysis | before | after correction |
|---|---:|---:|
| **Tahoe, same-dose replicates** | 11.5% [10.5–12.5%], *P* = 2×10⁻⁷ | **0.5% [0.0–1.5%], *P* = 0.10** |
| Tahoe, cross-dose pairing | 20.7% | 9.2% |
| OP3 | — | 33.1% |
| sci-Plex 3 | — | 30.2% |

**The Tahoe result is the consequential one.** `tahoe_true_replicates.py` carried
its own inline estimator and never imported the corrected library, so the paper's
headline empirical claim was still computed from a residual containing each
line's general response. Corrected, Tahoe shows **no detectable context ×
compound interaction at matched dose**. With 11,492 replicate pairs the interval
bounds it below 1.5%, so this is an informative null rather than an absence of
power. The cross-dose pairing does not "nearly double" the estimate as previously
written — it **manufactures eighteen-fold**, 0.5% → 9.2%.

The largest atlas built for this measurement cannot resolve, at 95.6 million
cells, the quantity it is quoted for. OP3 and sci-Plex 3 — smaller, but with more
replicates per condition — do resolve it, at 33.1% and 30.2%.

An internal check confirms the correction did what it should: in both external
datasets, residual agreement *across perturbations within a context* now sits at
+0.007 (OP3) and +0.026 (sci-Plex 3), against strongly positive values before.
That is the quantity a leftover context main effect inflates, and it is now at
its null.

### Claims wrong independently of the estimator

**The readout-decoupling claim was drastically underpowered.** The manuscript
reported "ρ = −0.09 against PRISM viability over **67 shared classes** — a
well-powered null". Only **11–12** mechanism classes clear the minimum class size
on both platforms; the true values are ρ = −0.18 (PRISM, n = 12) and +0.19
(LINCS phase 1, n = 10). At that n even ρ = 0.6 would not reach significance. The
class count was wrong six-fold and the power claim was unsupportable. Restated as
an absence of evidence.

**The laboratory-versus-assay apportionment was still inverted in three places** —
the Discussion, the Limitations and a Results paragraph — despite having been
corrected in the status banner. It is assay 81% / laboratory 19%, with the
laboratory CI [−6%, 38%] including zero.

**The withdrawn STR claim was still asserted in the Limitations**, stating that
all 738 screened lines are authenticated. Ten carry `FAILED_STR` and the count is
737.

### Resolved since

The few-shot re-run with the honest floor completed; see §32. The floor **beats**
the fine-tuned model at k = 5 and k = 20 and beats the oracle at k = 20, so the
few-shot gain is the line's general response rather than a learned interaction.
The LINCS decompositions also re-ran: phase 1 74% → 70%, phase 2 63% → 61%, the
same direction as PRISM and Tahoe.

---

## 32. The few-shot gain is the line's general response, not a learned interaction

§27 added an honest floor to the few-shot evaluation: the additive prior plus the
held-out line's **mean residual over the same *k* probe compounds** the model
sees. It is one arithmetic step with no fitting, and the probe compounds are
disjoint from the evaluation compounds, so nothing leaks. It changes what the
few-shot result means.

Metric *r*<sub>de100</sub>, 47 held-out lines, 5,615 evaluation conditions:

| | k = 1 | k = 5 | k = 20 |
|---|---:|---:|---:|
| fine-tuned model | 0.6650 | 0.6869 | 0.7093 |
| **additive + line's general response** | 0.6388 | **0.7104** | **0.7323** |
| model − floor | +0.026 | **−0.024** | **−0.023** |
| 95% CI (bootstrap over held-out lines) | [+0.017, +0.037] | [−0.030, −0.018] | [−0.029, −0.018] |
| model share of headroom | −3% | 45% | 94% |
| **floor share of headroom** | −61% | **97%** | **145%** |

Additive baseline 0.6665; oracle (all probes, fitted) 0.7120, so the headroom is
+0.0454.

**The floor beats the fine-tuned model at k = 5 and k = 20** (*P* = 6×10⁻²⁰⁵ and
1×10⁻²⁴⁵, paired over identical evaluation rows; the bootstrap resamples whole
held-out lines, the unit of generalisation). At k = 20 it also beats the
**oracle** — a model fitted on every available probe — recovering 145% of the
nominal headroom.

At k = 1 the ordering reverses: one compound estimates a line's general response
too noisily to help, and the floor falls below the additive baseline. That is the
control that matters, because it shows the floor is not winning by construction —
it wins only once it has enough probes to estimate the scalar it represents.

**What this changes.** The few-shot result stands as an empirical recipe: measure
about 20 arbitrary compounds in a new line and prediction improves substantially.
What does not stand is the interpretation. The improvement is almost entirely the
model discovering **how strongly the new line responds to everything** — a single
scalar shift — and not which drugs it responds to specifically. A mean over probe
wells captures that scalar better than fine-tuning does.

This is the same conclusion §31 reaches from measurement rather than modelling:
at matched dose Tahoe shows no detectable context × compound interaction
(0.5%, CI [0.0–1.5%]). If there is no interaction to learn, a model given probe
data can only learn the context main effect — and the floor learns it more
efficiently. Three independent lines of evidence now agree: the decomposition,
the replicate analysis, and the few-shot benchmark.

**Recommendation.** Any benchmark for context generalisation should report this
floor. Without it, an architecture is credited for recovering a quantity that a
mean over probe wells recovers better, and the credit grows with *k* precisely
because the scalar is estimated better — which looks exactly like a model
learning more from more data.

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
sbatch jobs/lincs_phase1.sbatch                      # LINCS phase 1 (GSE92742)
sbatch jobs/prism.sbatch                             # PRISM decomposition + genotype
python scripts/prism_vs_tahoe.py                     # cross-platform ranking
python scripts/three_platform_synthesis.py           # Tahoe / LINCS / PRISM
sbatch jobs/dose_ctx.sbatch                          # dose x context (PRISM, Tahoe)
sbatch jobs/lincs_dose.sbatch                        # dose x context (LINCS)
python scripts/validate_alleles_gdsc.py            # GDSC validation of alleles
```

See `docs/related_work_perturbation_models.md` for how these results sit
against published methods, and `docs/data_notes.md` for verified dataset facts.
