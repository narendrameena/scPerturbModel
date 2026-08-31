# Laboratory, not assay, limits the reproducibility of context-dependent drug response

**Draft manuscript.** Every number below is reproducible from this repository;
figure bundles (PNG/SVG/PDF + source data + generating script) are under
`results/figures/`, and each analysis ran as a recorded SLURM job in `jobs/`.
See `RESULTS.md` for the full result set, `docs/methodology_rationale.md` for why
each method was chosen over its alternatives, and `docs/pertdecomp.md` for the
tool.

> **Revision note (2026-08-31).** An earlier draft of this manuscript was titled
> *"Transcriptional drug response transfers between cellular contexts to a degree
> set by drug mechanism."* That framing has been withdrawn. The mechanism ranking
> is significant within Tahoe-100M (P = 8.7×10⁻⁴) but does **not** reproduce
> against LINCS phase 1 (ρ = +0.09) or PRISM viability (ρ = −0.09 over 67
> classes, a well-powered null). The earlier supporting value of ρ = +0.56 was
> produced by a biased per-perturbation estimator and does not survive its
> correction. The paper below is built on what replicates.

---

## Abstract

Giga-scale perturbation atlases are being built to predict how any cell responds
to any drug, and their central quantity — the part of a drug's effect that is
specific to a cellular context rather than shared across contexts — is widely
quoted but rarely measured properly. We show that measuring it requires
independent replicates at matched dose, that the flagship single-cell atlas
largely lacks them, and that the dominant obstacle to transferring the quantity
between studies is the laboratory rather than the assay. In LINCS phase 1, where
6.1 million genuine same-dose replicate pairs are available, **43% of the
reproducible transcriptional response is the compound's shared effect and 57% is
a context × compound interaction**. In Tahoe-100M only **4.7% of (line, compound,
dose) combinations are measured on more than one plate**; its apparent replicates
are different *doses*, and on true replicates the interaction is not detectable
(*P* = 0.97), so the atlas cannot resolve the quantity it is most used to report.
Estimating the interaction is subject to four failure modes that are invisible in
the output — using residual variance rather than replicate covariance, pooling
same-batch comparisons, taking an in-sample rather than leave-one-context-out
shared response, and treating doses as replicates — each of which we quantify and
guard against in a released tool. Applying the corrected estimator across four
atlases, the interaction is dose-dependent, rising roughly fourfold to a peak
just below lethality and collapsing where every line is dying (*P* = 4×10⁻⁹²),
and transcriptional and viability context-dependence are **decoupled**
(ρ = −0.09 across 67 mechanism classes). Genome-wide mutation status has **no
cross-validated predictive power** for the interaction (median CV *R*² = +0.0002,
50.8% of compounds positive), whereas **baseline expression predicts it in 92.5%
of compounds** (+0.093) and baseline protein nearly as well (+0.073); lineage
(+0.020) adds nothing once expression is included. Finally, of the line-specific response that
an assay reproduces with itself, only **56% survives transfer to another
laboratory** — but **87%** survives once cell-line identity is verified from the
data rather than assumed from an identifier, on held-out compounds. Only 5–12% of
identifier-matched lines are their own best match across atlases. Apportioning
the loss, changing assay within one laboratory costs 16% of it and changing
laboratory the remaining 84%. Context count and replication, not cell count, are
the binding constraints on what these atlases can answer.

---

## Introduction

Single-cell and pooled perturbation atlases now measure the consequences of
thousands of chemical interventions across dozens to hundreds of cellular
contexts, motivated by predictive "virtual cell" models. The difficulty such
models face is *context transfer*: predicting a compound's effect in a cell state
it was not profiled with. Progress is judged almost entirely by prediction
accuracy, and the picture is contested — benchmarks report that deep models do
not beat simple mean baselines, while method papers report substantial gains.

Underneath that dispute sits a quantity nobody has measured carefully: how much
of a drug's effect is actually shared between contexts, how much is
context-specific, and how reliably either can be estimated. The additive
main-effects model appears throughout this literature as a *baseline* and nowhere
as an *object of study*; its residual is not analysed, and the conditions under
which that residual can be estimated at all are not stated.

We treat the residual as the measurement, and we treat its measurability as the
first question rather than an assumption. Three findings follow, in order of
consequence: the estimator is fragile in ways that change conclusions; the
flagship atlas lacks the replicate structure the estimator requires; and the
quantity, where it can be estimated, transfers between laboratories far less well
than it repeats within one, and the loss is attributable to the laboratory rather
than the assay. Ben-David and colleagues showed cultures of one cell line diverge
between institutions; whether that divergence is what we are measuring is a
question our data can pose but, as we show, not settle.

---

## Results

### Estimating a context × compound interaction is fragile in four specific ways

The interaction is naturally estimated as the covariance of response residuals
between independent replicates of the same context × compound × dose. Four
apparently reasonable shortcuts each corrupt it, and none is visible in the
output. On our data:

| shortcut | why it fails | measured cost |
|---|---|---|
| residual **variance** instead of replicate covariance | noise does not cancel | 82% "interaction" against a true value near zero |
| pooling **same-batch** comparisons | batch state is shared signal | within-plate pairs agree ~7× better than cross-plate; 41% vs 28% |
| **in-sample** shared response | residuals sum to zero, forcing E[r_a·r_b] = −σ²/n | per-drug covariance driven negative for 21 of 24 drugs; clamping then yields exact zeros |
| treating **doses as replicates** | doses are different conditions | true replicates covary −0.0031 vs +0.0062 for cross-dose pairs (*P* = 3×10⁻¹²⁰) |

The third has a counterintuitive property worth stating because we got it wrong
twice: the obvious repair, leave-one-*condition*-out, makes the bias *worse*,
because the mean subtracted from replicate A still contains replicate B. Only
removing the whole context eliminates the coupling. We additionally report the
interaction against a **matched cross-context null**, which removes the residual
construction offset without assuming its magnitude; where the signal is strong
this changes nothing (LINCS phase 2: 45% → 47%), and where it is weak it was the
entire result. These guards are released as `pertdecomp`.

### The flagship single-cell atlas cannot measure the quantity it is used to report

Tahoe-100M comprises 95.6 million cells across 47 analysable lines and 379
compounds. Only **2,549 of its 53,881 (line, compound, dose) combinations — 4.7%
— appear on more than one plate**. The 96.4% of (line, compound) pairs that span
plates do so because the atlas places different *doses* on different plates.

Splitting the two pairings against a matched null:

| pairing | pairs | interaction share | *P* vs null |
|---|---|---|---|
| true replicate (same line, compound **and dose**) | 6,482 | **0.0%** | 0.97 |
| cross-dose (the pairing this implies) | 58,630 | 21.6% | 7×10⁻²⁵ |

Restricted to identical (line, compound) sets so only same-dose versus
different-dose differs, true replicates covary *less* than different doses do
(*P* = 3×10⁻¹²⁰) — which no reproducible interaction can produce. The
true-replicate arm covers only 25 of 379 compounds, and weak-effect ones, so it is
underpowered rather than decisive. **Neither number is the answer: the atlas
cannot resolve it**, and the value is bracketed between 0% and 21.6%.

Where genuine replicates exist the quantity is large and well determined. LINCS
phase 1 supplies 6.1 million same-dose cross-plate pairs across 2,834 compounds
and 70 lines, giving **43% shared / 57% interaction**.

### The interaction is dose-dependent and peaks just below lethality

Pooling doses, as all prior treatments do, hides a curve. Within compound, the
context-dependence index rises with dose (PRISM, 1,443 compounds, median
Spearman ρ = +0.33, 74% positive, *P* = 4×10⁻⁹²), climbing roughly fourfold to a
peak near 2.5 µM and then **collapsing at the top dose**, where the average line
has lost 1.67 log₂ of viability (0.391 vs 0.218 paired within compound,
*P* = 1×10⁻⁵⁹). The collapse is a ceiling effect: when every line is dying there
is little left to differ about. The rising limb replicates in LINCS transcription
(ρ = +0.20, *P* = 1×10⁻⁵), where the interaction grows faster than the shared
response (+74% vs +36%), which is why the ratio rises. A context-dependence index
is therefore only comparable at matched dose.

### Transcriptional and viability context-dependence are decoupled

With every dataset on one corrected estimator, mechanism structures
context-dependence *within* Tahoe (Kruskal–Wallis *P* = 8.7×10⁻⁴ across 24
classes) but the ordering does not transfer: ρ = +0.09 (n.s.) against LINCS phase
1 and **ρ = −0.09 against PRISM viability over 67 shared classes — a well-powered
null**. A drug whose transcriptional response is highly line-specific is not
thereby a drug whose killing is line-specific. Since the two readouts are used
interchangeably as "drug response", this is a substantive caution.

### Molecular state predicts the interaction; genotype does not

At 738 PRISM cell lines the genotype scan recovers the clinical biomarker set de
novo (TP53 with MDM2 inhibitors, BRAF with vemurafenib and dabrafenib, PIK3CA
with alpelisib, KRAS with a MEK inhibitor; 80 associations at FDR < 0.05 over 1.5
million tests), confirming the estimator and the genotype join. But out of
sample, across 150 compounds with 5-fold cross-validated ridge:

| predictor block | median CV *R*² | compounds positive |
|---|---|---|
| **baseline expression** (2,000 genes) | **+0.0927** | **92.5%** |
| **baseline protein** (RPPA, 214 antibodies) | **+0.0731** | 88.3% |
| lineage | +0.0201 | 74.2% |
| copy number (2,000 genes) | +0.0051 | 57.5% |
| nonsynonymous variants | +0.0002 | 50.8% |
| mutational burden | −0.0023 | 32.0% |
| synonymous variants | −0.0064 | 22.0% |

**Genome-wide mutation status carries no generalisable information; baseline
molecular state does.** Expression beats lineage by 4.6× (p = 4.7×10⁻²⁰) and
lineage adds nothing on top of expression, so lineage was acting as a coarse
proxy for expression state. Copy number beats mutations (+0.0064, p = 5.7×10⁻⁴),
as Schlüter & Schönhuth report, but trails expression by 0.088 (p = 3.3×10⁻²¹)
and adds nothing to it — the joint block scores below expression alone. The
absolute effect remains modest (≈9% of variance) but is highly consistent across
compounds. Using synonymous variants as a
control — silent changes cannot alter a protein but carry identical ancestry,
lineage and germline-contamination structure, and we match both blocks to the
same 3,435 genes — isolates a genuinely mechanistic excess of **+0.0036**
(*P* = 0.025), about 0.3% of the interaction variance.

Subsampling the confirmed PRISM associations quantifies why line-level tests fail
in smaller atlases: recovery is 4% at 47 cell lines, 45% at 250, 72% at 400 and
96% at 600. The genotype negatives reported in 47-line atlases are a power
ceiling, not a biological absence.

Screening 111,589 allele × compound tests for associations that are *not* known
pharmacogenomics and are invisible to a gene-level indicator yields four
candidates; none survives validation in GDSC, and all three testable ones are
identified as **germline** by a single shared genomic position, no excess
mutational burden, and neighbouring recurrent variants that are synonymous. The
same pipeline recovers 11 of 11 known biomarkers in GDSC (BRAF V600E ×
dabrafenib, *P* = 3×10⁻³¹), so the negative is informative.

### Cross-laboratory transfer is limited by the laboratory, not the assay

Against within-laboratory ceilings of r = 0.473 (PRISM replicate plates) and
0.438 (GDSC1 vs GDSC2), the line-specific response transfers between institutions
at r = 0.255 — **56% of what the assay achieves with itself**. The ceiling is
itself low: even repeating a measurement in one laboratory, the line's
compound-specific deviation agrees at only r ≈ 0.45.

Cross-atlas comparisons assert that a line in one atlas is the same as a line in
another, from an identifier that is never checked. Verifying it — each line must
be its own best match among all candidates by **response fingerprint**, its
residual across the compounds both atlases share — shows the assertion usually
fails: of 488 COSMIC→DepMap-matched lines, **only 5–12% are their own best
match**, ranking a median 82nd of 971 candidates. That is far better than chance
(~486), so identity carries real information; it is simply not unique.

**Verifying identity raises transfer from 56% to 87% of the ceiling.** This is
not circular: identity is validated on one random half of the shared compounds
and agreement measured on the disjoint other half (all-lines control on the same
held-out compounds: 59%). We initially attributed this to detecting divergent
cultures; baseline expression does not support that (below), so the gain is
reported as a practical filter whose mechanism is unresolved.

A matching ladder shows what each rule buys: random pairing r = −0.002,
same-tissue random pairing 0.048, identifier 0.235, best available partner 0.419.
Identifier matching beats same-tissue pairing decisively (*P* = 1.5×10⁻¹²⁴), so
identity carries information well beyond lineage — but the best available partner
reaches nearly twice what the identifier finds.

Apportioning the loss across the ladder — repeat plates (0.473) → different assay
version within one laboratory (0.438) → different laboratory and assay (0.255) —
**changing assay costs 16% of the total drop and changing laboratory the
remaining 84%.**

Transfer scales steeply with signal strength (ρ = +0.54, *P* = 1.8×10⁻²⁹): the
strongest quartile of compounds transfers *at* the within-laboratory ceiling,
the weakest at 36%.

The result is a property of laboratories rather than of a killing assay. Repeating
the analysis on transcription gives **52%** (LINCS phase 1 vs phase 2 within-lab
r = 0.061; Tahoe vs LINCS cross-lab 0.032), and the identity check behaves the
same way: 16 of 16 name-matched lines are reciprocal best matches within the
Broad, but only 3 of 6 across laboratories.

---

## Discussion

Three conclusions follow, in increasing order of consequence for how these
atlases are built and used.

**The estimator must be stated, because it decides the answer.** The same data
yield 82%, 41%, 28%, 21.6% or 0% interaction depending on choices that are rarely
reported: variance or covariance, pooled or cross-batch, in-sample or
leave-one-context-out prior, doses or true replicates as the replicate axis. Any
figure quoted for "the context-specific fraction" is uninterpretable without
them.

**Replication, not scale, is what a perturbation atlas most lacks.** Tahoe-100M
holds 100 million cells and cannot measure a context × compound interaction,
because 95.3% of its conditions are unreplicated at matched dose. LINCS phase 1,
far smaller in cells, answers the question cleanly with 6.1 million replicate
pairs. Cells were spent where replicates and contexts were needed. Together with
the power curve — 4% recovery at 47 contexts, 72% at 400 — this gives a concrete
design target: an atlas of roughly 400 lines, two plates per condition, and
one-tenth the cells per condition would answer what the present one cannot.

**Cross-atlas integration should verify identity, not assume it.** The
laboratory step costs 84% of the loss and the assay step 16%, so the problem is
not protocol; and much of the loss is recoverable: checking that a line is its own best
match by response fingerprint lifts transfer from 56% to 87% of the achievable
ceiling. Combined with weighting compounds by the strength of their line-specific
component, this is a cheap and immediate improvement to any analysis that pools
atlases. It also dissolves one puzzle of our own: the mechanism ranking's failure to
replicate across datasets needs no biological explanation, since a ~56%
reproducible fraction concentrated in strong-interaction compounds suffices. What
the loss consists of — divergent cultures, unmeasured protocol variables, or
simply a noisy line-specific signal — we can bound but not resolve.

### Limitations

*Design.* These are computational analyses of existing data, with no new
experiments and no prospective validation.

*Laboratory versus assay.* PRISM and GDSC differ in institution and in assay, so
neither is isolated by that comparison alone. We apportion them using GDSC1
versus GDSC2, which changes screening version, concentration range and assay
chemistry within one institution, and conclude that laboratory dominates
(84% versus 16%). This is a partial separation from two points on a ladder, not
a factorial design; a cross-laboratory comparison with the assay held exactly
fixed would settle it, and CTRP would have provided one but its NCI data portal
has been retired.

*Identity verification, and a correction.* We first read reciprocal-best-hit
failure as evidence of culture divergence. Baseline expression, now obtained from
the legacy CCLE distribution, does not support that: for the 423 lines whose
identifier match was outranked, the outranking line sits at a median expression
rank of **306 of 798** (chance 399) and shares the primary tissue only **7%** of
the time. Genuine divergence would place the better match among close relatives;
a near-random partner instead indicates the best-hit is largely fingerprint
noise. The held-out gain from 59% to 87% stands, but its mechanism does not — it
selects lines whose cross-laboratory signal is strong and self-consistent, which
is a useful filter, not a demonstrated detector of divergent cultures. The 5–12%
figure characterises the metric's resolution, not the rate of identity problems.
CCLE expression is also a single snapshot from one institution, so it cannot
compare two laboratories' cultures against each other.

*Transcriptional arm.* Its within-laboratory ceiling is only r ≈ 0.06, so the 52%
reproducible fraction is a ratio of two small numbers and is correspondingly
uncertain. It agrees with the viability arm, which is the substantive point, but
should not be quoted to two significant figures.

*Interaction magnitude in Tahoe.* We report a bracket (0–21.6%) rather than a
point estimate, and this propagates: any Tahoe-derived share in this work is an
upper bound obtained from a pairing the atlas's design forced on us.

*Withdrawn claims.* Three earlier results in this project did not survive: the
mechanism ranking outside Tahoe; a median 5.5% of interaction variance explained
by a single allele, which was in-sample and is ~0 cross-validated; and the
survival associations of §7, which did not survive adjustment for stage, grade
and purity. They are retained in `RESULTS.md` with the reasons.

*Scope.* All transcriptional results derive from short-timescale responses;
nothing here speaks to persistence or acquired resistance, where pre-existing
subpopulations are known to matter. Mechanism annotations are partly
machine-generated and demonstrably imperfect. Chromatin accessibility, the most
direct test of the regulatory hypothesis, is unavailable for these lines.

---

## Methods (summary)

Full implementations are in `scripts/`, with the SLURM job for each analysis in
`jobs/`, per-figure source data in each figure bundle, and the rationale for each
methodological choice — including the alternatives rejected and what they cost —
in `docs/methodology_rationale.md`.

**Data.** Tahoe-100M (95.6M cells; 50 lines, 379 compounds, three doses, 14
plates) aggregated to 65,918 (line, drug, dose, plate) pseudobulk profiles; LINCS
phase 1 (GSE92742) and phase 2 (GSE70138) Level 4, landmark genes only; PRISM
Repurposing secondary screen (738 lines, ~1,500 compounds, 8 doses, replicate
detection plates); GDSC1 and GDSC2 fitted dose response (978 lines, 542 drugs);
CCLE mutation calls via Cellosaurus, together with CCLE baseline expression
(19,221 genes), gene-level copy number and RPPA protein (214 antibodies) from
the legacy CCLE distribution and the DepMap figshare release; OP3 and sciPlex3
for architecture
replication.

**Decomposition.** Responses are batch-matched deltas against controls of the
same context and batch. The shared response is the **leave-one-context-out** mean
for the same compound and dose. The interaction is the covariance of residuals
between independent replicates of the same context × compound × **dose**, with
same-batch pairs excluded, reported against a **matched cross-context null**.
Perturbation-level indices carry their replicate-pair counts and an `estimable`
flag; with no replicate axis the tool refuses to report an interaction.

**Cross-laboratory analysis.** Per (line, compound) residuals in each dataset;
per-compound Spearman correlation across shared lines for viability, and
residual-profile correlation over shared genes for transcription. Within-lab
ceilings from PRISM replicate plates and GDSC1 vs GDSC2. Identity is verified by
response fingerprint on one random half of the shared compounds and evaluated on
the disjoint half.

**Predictor blocks.** Per compound, 5-fold cross-validated ridge with alpha
selected inside each training fold, solved in the dual when predictors exceed
samples. Blocks are matched so that none wins on capacity or coverage:
synonymous and nonsynonymous variants are restricted to the same 3,435 genes;
baseline expression and gene-level copy number are each reduced to their 2,000
most variable features; and every block is scored on the same cell lines for a
given compound. Copy number, expression and protein are compared on the subset
of lines present in all three.

**Statistics.** Nulls are matched to the confound in each case: cross-context
pairs for the residual offset, size-matched solid panels for the immune
comparison, permuted pairings for identity, burden-conditioned partial
correlation for MSI, and lineage stratification for allele associations.

---

## Figures

1. **Estimating the interaction.** (a) The four failure modes and their measured
   cost. (b) Tahoe's replicate structure: 4.7% of conditions replicated. (c) True
   replicate versus cross-dose pairing against a matched null.
2. **Architecture where it can be measured.** (a) Shared/interaction shares in
   LINCS phase 1 and 2, PRISM, OP3. (b) Residual reproduces across replicates but
   not across compounds or contexts. (c) Dose curve: rise to a peak below
   lethality, then collapse.
3. **What does and does not predict it.** (a) Cross-validated *R*² by predictor
   block — baseline expression, baseline protein, lineage, gene-level copy
   number and nonsynonymous variants — on identical compounds and lines, with
   expression and copy number capped at the same 2,000 most-variable features so
   neither block wins on capacity. Copy number beats mutations (p = 6×10⁻⁴) but
   trails expression by 0.088 (p = 3×10⁻²¹), placing DNA-level features at the
   bottom of the ordering. (b) The synonymous control, plotted as the
   per-compound nonsynonymous-minus-synonymous difference, since both blocks are
   individually negative and the claim is about their gap. (c) Power curve for
   genotype linkage versus context count. (d) Readout decoupling: transcription
   versus viability mechanism rankings.
4. **Cross-laboratory reproducibility.** (a) The ladder from repeat plates to a
   different laboratory. (b) Rank of the identifier-matched line by response
   similarity. (c) Reproducible fraction before and after identity verification,
   on held-out compounds. (d) Transfer versus strength of the line-specific
   component.
5. **Matching strategy, and why the divergence reading was withdrawn.**
   (a) Fingerprint similarity under random, same-tissue, identifier and best-hit
   pairing. (b) Their distributions. (c) Expression rank of the line that
   outranked the identifier match: near chance, and sharing tissue only 7% of
   the time, which is why best-hit failure is read as metric noise rather than
   demonstrated culture divergence.
