# Transcriptional drug response transfers between cellular contexts to a degree set by drug mechanism

**Draft manuscript.** Every number below is reproducible from this repository;
figure bundles (PNG/SVG/PDF + source data + generating script) are under
`results/figures/`, and each analysis ran as a recorded SLURM job in `jobs/`.
See `RESULTS.md` for the full result set and `docs/methodology_positioning.md`
for an audit of these methods against the published literature.

---

## Abstract

Giga-scale single-cell perturbation atlases are being built to predict how any
cell responds to any drug, yet how much of a drug's effect is shared between
cellular contexts — and what the unshared part consists of — has not been
measured. Decomposing 95.6 million single-cell transcriptomes spanning 47 cancer
cell lines and 379 compounds, we find that **59% of the reproducible
transcriptional response to a drug is its context-averaged effect and 41% is a
cell-line × drug interaction**. That interaction is real but is *not a property
of the cell line*: it reproduces across independent replicates of the same
line–drug pair (r = 0.062) yet barely transfers to other drugs in the same line
(r = 0.019). The same architecture holds in primary human immune cells and in a
second cancer-line atlas, across three platforms and three kinds of replicate.
It is not a viability artefact — at least six transcriptional programmes carry
it, and the leading one is uncorrelated with cell-cycle arrest — and at
single-cell resolution it is a uniform population shift rather than selection of
a pre-existing subpopulation. How much of a compound's effect is
context-specific is predictable from **mechanism** (P = 2×10⁻³ across 24
classes): nuclear-receptor agonists are the least transferable (index 0.363) and
metabolic and RAF inhibitors the most conserved (0.087–0.119), while chemical
structure alone predicts poorly (r = 0.19) and target genotype not at all. For a
cell line never seen during training, no descriptor we tested — driver
mutations, tissue of origin, baseline transcriptome, DNA methylation — recovers
its interaction; measuring roughly 20 *arbitrary* compounds and fine-tuning
recovers 98% of the achievable gain, and panel design does not help. Finally,
the apparent strength of the mean baseline that such models must beat depends on
benchmark design: averaging 45 contexts rather than 5 improves it by 10.6%,
comparable to gains recently attributed to models evaluated on six contexts.
These results quantify what is transferable in drug response, what is not, and
what it costs to obtain.

---

## Introduction

Single-cell perturbation atlases now measure the transcriptional consequences of
thousands of chemical interventions across dozens of cellular contexts, and are
explicitly motivated by the goal of predictive "virtual cell" models. The
central difficulty such models face is *context transfer*: predicting a
compound's effect in a cell state that was not profiled with it.

Progress is assessed almost entirely by prediction accuracy, and the picture is
contested. Several benchmarks report that deep models do not outperform simple
mean baselines, while individual method papers report substantial gains. What
has not been established is the quantity that determines what any model can
achieve: how much of a drug's transcriptional effect is actually shared between
contexts, how much is context-specific, and what the context-specific part is
made of. The additive main-effects model appears throughout this literature as a
*baseline* and nowhere as an *object of study*; its residual is not analysed.

We treat that residual as the measurement. Using Tahoe-100M — 95.6 million cells
across 47 analysable cancer cell lines and 379 compounds at three doses — we
decompose the response into a context-averaged component and a context ×
compound interaction, establish that the interaction is reproducible signal
rather than noise using the atlas's own replicate structure, characterise what
it consists of, test every cell-line descriptor we could obtain against it, and
measure what it costs to acquire for a new context.

---

## Results

### The additive component dominates, and the remainder is an interaction

For every (line, drug, dose) condition we computed a plate-matched response
delta and an additive prior — the mean delta for that drug and dose over *other*
lines. Per gene, the additive component accounts for 4.8% of response variance
and the line × drug interaction 3.3%, the remainder being noise (91.9%); of the
**reproducible** signal, **59% is additive and 41% interaction** (Fig. 1a).

That the interaction is reproducible rather than residual noise is established
by the atlas's replicate structure. Comparing context residuals between
independent measurements, correlation is **+0.062** between doses of the same
line–drug pair, **+0.019** between different drugs in the same line, and
**−0.002** between lines for the same drug, against ≈0 permuted nulls (Fig. 1b).
The residual therefore contains real, repeatable signal that is specific to the
*pair* and does not generalise to the line.

**A batch caveat that changes the numbers.** All of the above is restricted to
comparisons made on different plates. The identical comparisons *within* a plate
give +0.45 to +0.48 — roughly sevenfold higher — because residual plate
structure survives plate-matched DMSO normalisation. Analyses of this atlas that
pool within- and cross-plate comparisons will substantially overstate
reproducibility.

### The architecture replicates across atlases, platforms and replicate designs

We re-ran the decomposition unchanged on two external atlases by mapping their
columns onto the same three roles (context, perturbation, independent
replicate):

| | Tahoe-100M | OP3 | sciPlex3 |
|---|---|---|---|
| contexts | 47 cancer lines | 6 immune cell types | 3 cancer lines |
| replicate axis | plates | donors | experimental repeats |
| interaction share | 41% | 33% | 30% |
| residual across replicates | +0.062 | +0.135 | +0.048 |
| residual across perturbations | +0.019 | +0.007 | +0.026 |

The defining contrast — reproducible across replicates, non-transferable across
perturbations — holds in primary human immune cells as well as cancer lines
(Fig. 1c). The additive share rises as contexts fall (59% → 67% → 70% at 47, 6
and 3 contexts), consistent with a prior estimated from fewer contexts being
noisier; we quantify that effect directly below.

### The interaction is multi-dimensional and is not a viability artefact

The obvious deflationary explanation is that the residual is rank-1 — each
pair's sensitivity times one shared cytotoxic direction — which would explain
both its existence and its failure to transfer. Fitting components on one half
of each pair's plates and scoring them on the other half, the leading component
carries only **34%** of reproducible residual variance, at least **six**
components reproduce (split-half r = 0.33–0.48) while others return negative
values, and the leading component is **uncorrelated with cell-cycle arrest**
measured independently from the same cells (G2M r = −0.05) (Fig. 2a–c). By
pathway enrichment the reproducible components correspond to epithelial-versus-
neuronal identity, ERBB4/PI3K signalling, immune-versus-ECM, EMT, ribosome
biogenesis, and chromatin/senescence programmes (Fig. 2d).

At single-cell resolution, projecting 20.5 million cells onto these components
shows the typical context-specific response is a **uniform shift of the whole
population**: the median Kolmogorov–Smirnov distance between treated cells and
shifted controls is 0.91× the within-control noise floor, with only 11% of
comparisons exceeding twice it. The strongest decile does change shape (2.13×)
but with variance *decreasing* (0.87), i.e. convergence rather than a
subpopulation separating (Fig. 2e). Context-dependence at 24 h is therefore not
selection of a pre-existing responsive subpopulation.

### Context-dependence is set by drug mechanism

For each of 367 compounds we computed a context-dependence index (CDI) —
reproducible interaction as a fraction of total reproducible effect. Mechanism
explains a significant share of it (Kruskal–Wallis H = 46.9, P = 2.3×10⁻³ across
24 classes; Fig. 3a):

| most context-dependent | CDI | most conserved | CDI |
|---|---|---|---|
| Glucocorticoid receptor agonist | 0.363 | Glucose transporter inhibitor | 0.087 |
| Proteasome inhibitor | 0.334 | Other MAPK inhibitor | 0.109 |
| MEK inhibitor | 0.324 | RAF inhibitor | 0.119 |
| Retinoic receptor agonist | 0.232 | Other TK inhibitor | 0.122 |

The leading class is mechanistically expected: nuclear-receptor agonists signal
through receptors whose transcriptional output is set by each cell's own
enhancer landscape, so their effects should be maximally context-specific. Four
of the twelve most context-dependent individual compounds are corticosteroids,
and retinoic-receptor agonists rank fourth by class. MEK inhibition is among the
most context-dependent mechanisms while RAF inhibition is among the most
conserved, despite acting one step apart in the same pathway.

Chemical structure alone predicts CDI only weakly (cross-validated r = 0.19 on
ECFP4 fingerprints; Fig. 3b), so this is a property of the target and mechanism
rather than of the molecule.

### It is a regulatory, not a target-genetic, phenomenon

Three independent tests fail to link the interaction to genotype. A scan of 25
mechanism classes against 33 driver genes yields no association surviving
correction (0/825 at FDR < 0.10) once each line's overall responsiveness is
removed. Across 254 compounds with annotated targets, CDI is uncorrelated with
how often the target is mutated in the panel (ρ = +0.10, P = 0.10) and with
target expression level (ρ = +0.08). The clearest case is the leading class
itself: nuclear-receptor compounds have a median CDI of 0.342 against 0.154
overall, yet their targets are mutated in **zero** atlas lines.

Together with the identity of the response programmes — epithelial and
mesenchymal identity, chromatin and senescence — this indicates that
context-dependence is set by the cell's regulatory state rather than by
mutations in the drug's target.

### No cell-line descriptor predicts the interaction for an unseen line

For a line held out entirely, we tested every descriptor available: driver
mutations plus tissue of origin, the line's own baseline transcriptome
(per-fold PCA), and DNA methylation level and heterogeneity from CCLE bisulfite
profiles covering 45 of our lines. None improves on the additive baseline
(0.692 versus 0.692 for the first two; 0/20 methylation features at FDR < 0.10)
(Fig. 4a). This follows from the structure above: if the interaction is a
property of the *pair*, no line-level quantity exists for a descriptor to
recover.

### Adaptation is possible, but requires the model to change, not just its input

Fitting only a held-out line's embedding with the rest of the model frozen
yields **nothing** at any number of probe compounds — including an oracle using
every available probe drug (+0.0009 over 8,420 conditions) — even though the
embedding demonstrably moves. Retraining the model with the same probe
measurements included instead recovers the gain: +0.020 at five compounds and
**+0.034 at twenty**, or 98% of the achievable ceiling (Fig. 4b). The barrier is
therefore architectural, not informational: the learned context space cannot be
entered post hoc.

Panel design does not matter. Comparing five selection strategies over 47 folds,
choosing the most potent or the most line-discriminating compounds performs
*worse* than random at five compounds, and by twenty compounds every strategy
converges on the oracle (Fig. 4c). Practically, a new cell model needs neither
its genotype nor its baseline profile: roughly twenty arbitrary compounds and a
fine-tuning pass.

### The apparent headroom for models depends on benchmark design

The mean baseline that perturbation models are measured against is itself an
average, so its quality depends on how many contexts contribute. Holding
evaluation conditions fixed and varying only that number, the additive baseline
improves from r = 0.649 (5 lines) to 0.718 (45), **+10.6%** (Fig. 5). For
comparison, a recent model reports +12.3% over the best baseline while
restricting this atlas to six cell lines. We do not conclude that any specific
model's gains are artefactual — a well-estimated baseline may still be beaten —
but a mean baseline must be estimated from all available contexts before
headroom is claimed, and benchmarks that subsample contexts for tractability
will systematically overstate it.

### Context count, not cell count, is the binding constraint

Four independent line-level predictors failed above. At 43–50 lines a
line-level correlation requires |ρ| > 0.30 to reach nominal significance, so
effects of realistic size are invisible; these four negatives are plausibly one
power ceiling encountered four times. The atlas contains 95.6 million cells but
only ~47 usable contexts, and for line-level questions the latter determines
power. The same accounting explains a failed replication: the per-compound
mechanism ranking is not estimable in OP3, where a compound has a median of 18
cross-replicate pairs against 142 here, and only 23% of compounds yield a
positive interaction estimate against 96%.

---

## Discussion

We set out to measure rather than predict, and the measurements bound what
prediction can achieve. Most of a drug's reproducible transcriptional effect —
about 59% — is shared across cellular contexts, which is why simple additive
baselines are hard to beat and why reports of their weakness are partly an
artefact of evaluating on few contexts. The remaining 41% is real and
reproducible, but it is an interaction between a specific drug and a specific
cell rather than a property of either, which is a sufficient explanation for the
repeated failure of cell-line descriptors in this field: there is no line-level
quantity to predict.

Two results are directly actionable. First, how broadly a compound must be
screened is predictable from its mechanism before any experiment: nuclear
receptor, proteasome and MEK pharmacology require broad panels, whereas RAF,
metabolic and most tyrosine-kinase inhibitors transfer. Second, a new cell model
can be brought into a trained model with roughly twenty arbitrary compounds and
a fine-tuning pass, and effort spent designing that panel is wasted.

Our results also identify what these atlases should optimise. Effort is
currently directed at cell count and model capacity, but for context transfer —
the problem the atlases exist to solve — the binding constraint is the number of
distinct contexts profiled, and for per-compound statements the replication
depth available within each. An atlas with 500 lines and one-tenth the cells per
condition would answer questions this one cannot.

**Limitations.** These are computational analyses of existing data with no new
experiments and no prospective validation. Every result derives from 24 h
transcriptional responses; nothing here speaks to survival, persistence or
resistance over longer timescales, where pre-existing subpopulations are known
to matter. The mechanism ranking rests on a single atlas and on mechanism
annotations that are partly machine-generated and demonstrably imperfect — we
independently identified two mislabelled compounds. The cell-line descriptor
tests are underpowered as described. Finally, chromatin accessibility, the most
direct test of the regulatory hypothesis, is unavailable for these lines; DNA
methylation is an imperfect proxy.

---

## Methods (summary)

Full implementations are in `scripts/`, with the SLURM job for each analysis in
`jobs/` and per-figure source data in each figure bundle.

**Data.** Tahoe-100M (95.6M cells passing full filters; 50 lines, 379 compounds,
three doses, 14 plates) aggregated to 65,918 (line, drug, dose, plate)
pseudobulk profiles. External replication used OP3 (NeurIPS 2023; 6 immune cell
types, 147 compounds, 3 donors) and sciPlex3 (3 lines, 189 compounds, 2
replicates, 24 h subset).

**Decomposition.** Responses are plate-matched deltas in log1p CPM against DMSO
controls of the same line and plate. The additive prior for a condition is the
mean delta of the same drug and dose over other lines (leave-one-out).
Reproducible interaction is the covariance of residuals between *independent
replicates* — different plates, donors or experimental repeats — of the same
context × perturbation, so independent noise and shared batch structure cancel.
All reported reproducibility statistics exclude same-plate comparisons.

**Components.** Fitted by SVD on one half of each pair's plates and scored on
the other; only components whose per-pair loadings reproduce are interpreted.
Enrichment uses hypergeometric tests against the responsive-gene universe with
Benjamini–Hochberg correction.

**Models.** `prediction = additive_prior + residual(line embedding, ECFP4, log
dose)`, residual zero-initialised so the model begins at the additive baseline;
leave-one-line-out priors on training rows prevent target leakage. Splits are
defined on (line, drug, dose) triples and, for the strict variant, on whole
(line, drug) pairs.

**Statistics.** Nulls are matched to the confound in each case: expression- and
size-matched gene sets for target-abundance tests, permuted pairs for residual
correlations, within-cancer-type z-scoring before pooling TCGA, and sequencing
depth as an explicit covariate in the methylation analysis.

---

## Figures

1. **The architecture of drug response.** (a) Variance shares. (b) The three
   residual correlations with permuted nulls, and the within- versus cross-plate
   contrast. (c) Replication in OP3 and sciPlex3.
2. **What the interaction consists of.** (a) Component variance shares.
   (b) Split-half reproducibility. (c) Leading component versus cell-cycle
   arrest. (d) Pathway enrichment of reproducible components. (e) Uniform shift
   versus subpopulation at single-cell resolution.
3. **Mechanism sets context-dependence.** (a) CDI by mechanism class.
   (b) Predictability from chemical structure. (c) Independence from target
   genetics.
4. **Predicting an unseen context.** (a) Failure of every line descriptor.
   (b) Frozen versus fine-tuned adaptation across probe counts. (c) Probe-panel
   strategies.
5. **Benchmark design determines apparent headroom.** Additive baseline quality
   versus the number of contexts averaged.
