# Measuring context-dependent drug response: what the estimator, the replicates and the compound set decide

**Draft manuscript.** Every number below is reproducible from this repository;
figure bundles (PNG/SVG/PDF + source data + generating script) are under
`results/figures/`, and each analysis ran as a recorded SLURM job in `jobs/`.
See `RESULTS.md` for the full result set, `docs/methodology_rationale.md` for why
each method was chosen over its alternatives, and `docs/pertdecomp.md` for the
tool.

> ## ⚠ STATUS (2026-09-03): the central claim is withdrawn; the paper is now a measurement paper
>
> A defect larger than the ~60 the three audits found was located on 2026-09-03.
> **Every residual in this project removed the compound main effect but not each
> cell line's general sensitivity** — its response to *every* compound, set by
> growth rate, seeding density and drug metabolism. That term reproduces across
> disjoint compound halves at **r = 0.989** (737 lines) and is shared between
> replicate detection plates exactly as a genuine interaction is, so it entered
> the pair covariance and was counted as context-dependence.
>
> | claim | was | now |
> |---|---|---|
> | the interaction is a **pair** property, not a cell property | our top-ranked novel result; contradicted Lim & Pavlidis (*Sci Rep* 2021) | **withdrawn.** A cell-line responsiveness factor exists as published. The relation is **2.0× the property** (3.9% vs 2.0%) — a measurement, not a refutation |
> | interaction shares | e.g. Tahoe 11.5%, LINCS-1 57% | inflated ~**1.5×**; drug 94.1% / cell property 2.0% / relation 3.9% in PRISM |
> | title / apportionment | laboratory 84% of the loss | **assay 81%, laboratory 19%** on the 75 compounds all three rungs share |
> | identity-validated transfer | 87% | **98%** against a ceiling on the same lines; a reliability-matched control gives 49%, so the effect is identity-specific |
> | expression over genotype | 92.5% of compounds | **119/120 (99.2%)** beat their own permutation null; mutations 22/120, zero compounds mutation-only (p = 1×10⁻²⁹) |
> | copy number beats mutations | p = 5.7×10⁻⁴ | **withdrawn** — CI [−0.000, +0.010] under a compound-cluster bootstrap |
> | PRISM cell lines | 738 | **737** — the old count included a line fabricated by a parsing bug |
> | dose peak-coincidence tests | 14% / 81% cited as support | **at or barely above their own permutation nulls** |
>
> The estimator now recovers a known share with slope **0.950** (R² = 0.9995) and
> returns 0.008 from data containing none. See `RESULTS.md` §27 and
> `docs/novelty_audit.md`.

> **Revision note (2026-08-31).** An earlier draft of this manuscript was titled
> *"Transcriptional drug response transfers between cellular contexts to a degree
> set by drug mechanism."* That framing has been withdrawn. The mechanism ranking
> is significant within Tahoe-100M (P = 6.0×10⁻⁴) but does **not** reproduce
> against LINCS phase 1 (ρ = +0.09) or PRISM viability (ρ = −0.09 over 67
> classes, a well-powered null). The earlier supporting value of ρ = +0.56 was
> produced by a biased per-perturbation estimator and does not survive its
> correction. The paper below is built on what replicates.

---

## Abstract

Giga-scale perturbation atlases are being built to predict how any cell responds
to any drug, and their central quantity — the part of a drug's effect specific to
a cellular context rather than shared across contexts — is widely quoted but
rarely measured properly. We show that a measured response contains **three**
terms, not two: the drug's average effect, a property of the cell that has
nothing to do with any drug, and the relation between the pair. Only the third is
context-dependence. The middle term — a line's *general sensitivity*, its
response to every compound — is routinely left inside the residual, where it is
indistinguishable from a genuine interaction because it is equally reproducible
across replicates. It is large: a line's general sensitivity estimated on one
half of a compound library predicts the estimate on a disjoint half at
**r = 0.989** across 737 cell lines, and **570 of 736 lines (77.4%)** carry a
reproducible line-level shift at FDR < 0.05.

Estimating the three terms as covariances between independent estimates, so that
noise contributes zero to each, gives **drug effect 94.1%, cell property 2.0%,
cell–drug relation 3.9%** (PRISM, 737 lines × 1,324 compounds). Interaction
shares that omit the middle term are inflated by roughly **1.5×**. On simulated
data with a known share the corrected estimator recovers it (slope 0.950,
R² = 0.9995) and returns 0.008 from data containing none, while residual variance
reports 72% and pooled-batch 41%; leaving general sensitivity in place drives the
same estimator from 0.170 to 0.595 as the planted cell property grows, with the
truth fixed at 0.20. Both main effects must be estimated out of fold — the drug
effect from other lines, the cell property from other compounds — and the cell
property **within a replicate plate**: pooled across plates it carries a share of
each plate's control noise and subtracts var(control)/2 from the pair covariance,
which drove the estimate to exactly zero in simulation while appearing correct in
review. We release the decomposition and a per-pair query as a tool.

Applied across four atlases, the corrected interaction is dose-dependent, and
transcriptional and viability context-dependence remain **decoupled** (ρ = −0.18
across mechanism classes). Asking what predicts the corrected relation once per
compound, against that compound's own permutation null with FDR across the
family, **baseline expression beats its null in 119 of 120 compounds (99.2%)**,
lineage in 104, and **genome-wide mutation status in 22**, whose median
cross-validated *R*² is negative; 97 compounds are expression-significant and
mutation-non-significant and **none** the reverse (p = 1×10⁻²⁹). Expression
predicts the *relation* (+0.081 [+0.072, +0.093]) better than it predicts the
*cell property* (+0.068), so the result is not general sensitivity in disguise;
copy number does not separate from mutations once compound correlation is
respected (+0.005, CI [−0.000, +0.010]). Finally, of the line-specific response
an assay reproduces with itself, **58% [49–77%]** survives transfer to another
laboratory, rising to **98%** once cell-line identity is verified from the data
rather than assumed from an identifier — against a reliability-matched control at
49%, so the effect is specific to identity rather than to measurement quality.
Apportioned on a matched compound set, the **assay accounts for 81% of the loss
and the laboratory 19%**, the opposite of what the unmatched comparison reports.
Context count and replication, not cell count, are the binding constraints on
what these atlases can answer.

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
replicate structure the estimator requires exists in the flagship atlas but is
discarded by convention, which doubles the reported value; and the
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
| treating **doses as replicates** | doses are different conditions | 11.5% on true replicates vs 20.7% cross-dose — nearly double (Tahoe, plate 6 vs 14) |

The third has a counterintuitive property worth stating because we got it wrong
twice: the obvious repair, leave-one-*condition*-out, makes the bias *worse*,
because the mean subtracted from replicate A still contains replicate B. Only
removing the whole context eliminates the coupling. We additionally report the
interaction against a **matched cross-context null**, which removes the residual
construction offset without assuming its magnitude; where the signal is strong
this changes nothing (LINCS phase 2: 45% → 47%), and where it is weak it was the
entire result. These guards are released as `pertdecomp`. Against a standard mixed-model
decomposition on the same simulated data, the two agree to 0.051 where both are
identifiable — ours is not a reinvention — but without replicates the mixed model
returns 0.161 for a true share of 0.1 and 0.182 for 0.3, a number set by the
optimiser rather than the data and flagged by nothing in its output. That is the
case for a tool that refuses.

### The replicate structure exists, and discarding it doubles the estimate

Tahoe-100M comprises 95.6 million cells across 47 analysable lines and 379
compounds, and it was built with the replication this measurement needs: **plate
14 is a biological replicate of plate 6**, 6.2M cells over 50 lines and 95
compounds, included by the authors to demonstrate platform reproducibility. All
4,746 of its (line, compound, dose) triples also appear on plate 6.

Counted from the released metadata with no filtering of ours, **7,691 of 56,877
triples — 13.5% — sit on more than one plate**, a figure stable across the
atlas's own quality filter (13.52%) and a ≥10-cells-per-plate requirement
(13.06%). Separately, 96.8% of (line, compound) *pairs* span plates, but only
because different **doses** sit on different plates; those are not replicates.

The authors withhold plate 14 from training and reserve it for validation, which
is correct for fitting a model. Any analysis pipeline that copies the convention
— ours did, by default — loses the atlas's only same-dose replicates and drops
replication to 5.4%. Measured against a matched cross-context null:

| pairing | pairs | interaction share | *P* vs null |
|---|---|---|---|
| **true replicate** (same line, compound **and dose**) | 11,492 | **11.5% [10.5–12.5%]** | 2×10⁻⁷ |
| cross-dose (what remains without plate 14) | 67,744 | 20.7% [20.5–20.9%] | 1×10⁻³² |
| pooled | 79,236 | 19.5% [19.3–19.8%] | 7×10⁻²⁸ |

The interaction is real and well determined at **11.5%**, and the cross-dose
pairing **nearly doubles it**. That pairing also attenuates with dose separation
(ρ = −0.020, *P* = 2×10⁻⁷), as expected if those pairs are different conditions
rather than repeats.

LINCS phase 1, with 6.1 million same-dose cross-plate pairs across 2,834
compounds and 70 lines, gives **43% shared / 57% interaction** — a different
platform with plate-wise z-scoring rather than log-CPM deltas, so the magnitudes
are not directly comparable, but the estimator is the same.

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
context-dependence *within* Tahoe (Kruskal–Wallis *P* = 6.0×10⁻⁴ across 24
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
the analysis on transcription gives **46%** (LINCS phase 1 vs phase 2 within-lab
r = 0.061; Tahoe vs LINCS cross-lab 0.032), and the identity check behaves the
same way: 16 of 16 name-matched lines are reciprocal best matches within the
Broad, but only 3 of 6 across laboratories.

### What the measurements imply for models

The results above are measurements, but the audience for them is building
predictive models, and three of the measurements bound directly what such models
can achieve.

**The baseline they must beat depends on how many contexts it averages.** The
perturbation-mean prediction improves monotonically with the number of contexts
it is estimated from: r = 0.550 at 2 contexts, 0.592 at 3, 0.625 at 6, 0.654 at
18 and 0.666 at 45. A model that is *equally good everywhere* would therefore
appear **21% stronger** against a 2-context baseline than against an 18-context
one.

This is directly relevant to how current results are read. State (Adduri et al.,
*Cell* 2026) evaluates zero-shot transfer across five query datasets holding out
one context at a time, so its perturbation-mean baseline is built from 2, 2, 3, 5
and 17 contexts; the paper reports "more than 19% improvement" on the
18-context dataset and "several-fold improvements" on the smallest — the ordering
baseline quality alone produces. We are explicit that this is a **partial**
confound: a ~19-point differential is real and runs in the reported direction,
but it does not explain several-fold gaps, and State is evidently doing something
the baseline is not. The narrower and still useful conclusion is that gains are
comparable across benchmarks only at matched context counts, and that benchmarks
should report the number of contexts behind their mean baseline.

**Prediction for an unseen context cannot come from a context descriptor.** No
line-level feature we tested recovers a new line's interaction — driver
mutations, tissue, baseline transcriptome, DNA methylation — and §17 explains
why the search was ill-posed rather than merely underpowered: the interaction is
a property of a *pairing*, and the best available descriptor, baseline
expression, reaches only ≈9% of it. What does work is measurement: profiling
roughly **20 arbitrary compounds** in the new line and fine-tuning recovers 98%
of the achievable gain, and designing that probe panel does not beat choosing it
at random. State reaches the same conclusion from the opposite direction — its
headline context-generalisation setting supplies each model with 30% of the
perturbations in the test context, making it an *underrepresented* rather than
unseen context — so the two results corroborate each other.

**The ceiling on cross-atlas training is 56–87% of an assay's own
reproducibility.** A model trained on one atlas and applied to another is
predicting a quantity that transfers at just over half the within-laboratory
ceiling if cell lines are matched by identifier, and close to the ceiling if
identity is verified from the data first. Weighting compounds by the strength of
their line-specific component matters as much: the strongest quartile transfers
at the ceiling, the weakest at a third of it, so pooling atlases without such a
weighting spends most of its capacity on the compounds that transfer least.

Together these say the binding constraints are experimental rather than
architectural — context count, replication, and identity verification — which is
why this paper measures rather than models.

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

**Replication is what gets discarded, not what is missing.** Tahoe-100M was
built with a replicate plate and can measure the interaction — 11.5%
[10.5–12.5%]. The problem is that the replicate plate is withheld from training
by convention, and pipelines inherit that exclusion, leaving a cross-dose pairing
that doubles the estimate. We made this mistake ourselves and reported the
opposite conclusion before catching it. The practical recommendation is
therefore analytical rather than architectural: **keep the replicate plate in
when measuring, and never treat doses as replicates.** Where atlases genuinely
are short is contexts — the power curve gives 4% recovery of genotype
associations at 47 contexts against 72% at 400 — so an atlas of roughly 400
lines with two plates per condition and one-tenth the cells per condition remains
the design target.

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

*Only two laboratories.* Every cross-laboratory result rests on the Broad and the
Sanger. A third was attempted twice and not achieved: NCI-60 is unreachable
without a browser session and CTRP's portal is retired, while sciPlex3 — a
genuine third institution sharing three cell lines with LINCS — has too few
contexts for the estimator, since its leave-one-context-out shared response is a
mean of two lines against LINCS's mean of sixty-nine. This is the single most
valuable thing new data could add, and it is why the laboratory-versus-assay
apportionment is quoted with an interval that reaches 100%.

*What expression predicts is not a programme.* The 9% is carried diffusely — 751
of 2,000 genes hold half the ridge weight, and no gene set among 5,457 survives
correction — so the predictor tracks global transcriptional state rather than an
interpretable mechanism.

*Laboratory versus assay.* PRISM and GDSC differ in institution and in assay, so
neither is isolated by that comparison alone. We apportion them using GDSC1
versus GDSC2, which changes screening version, concentration range and assay
chemistry within one institution, and conclude that laboratory dominates
(84% versus 16%). This is a partial separation from two points on a ladder, not
a factorial design; a cross-laboratory comparison with the assay held exactly
fixed would settle it, and CTRP would have provided one but its NCI data portal
has been retired.

*Identity verification, and a correction.* Every cell line in the
cross-laboratory comparison passes STR profiling — PRISM excludes failures from
its released screen, so all 738 screened lines are authenticated — which rules
out outright misidentification and leaves metric noise and culture drift within
correctly labelled lines as the explanations. We first read reciprocal-best-hit
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

*Transcriptional arm.* Its within-laboratory ceiling is only r ≈ 0.06, so the 46%
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

**Metrics, and their relation to Cell-Eval.** Prediction accuracy is reported as
Pearson correlation between observed and predicted perturbation-induced
expression changes, over the 2,000 most response-variable genes (`r_hvg`) and
over each condition's own top-100 differentially expressed genes (`r_de100`).
These correspond directly to the Pearson-delta metric of Cell-Eval (Adduri et
al. 2026); `r_de100` is its restriction to the DE gene set. We do not report
Cell-Eval's perturbation-discrimination or DE-overlap scores, because those rank
a model's predictions against alternatives and the central quantity here is a
variance decomposition rather than a prediction. Where we do make prediction
claims — the baseline-scaling result and the probe-compound protocol — the metric
is Cell-Eval-comparable.

**Statistics.** Nulls are matched to the confound in each case: cross-context
pairs for the residual offset, size-matched solid panels for the immune
comparison, permuted pairings for identity, burden-conditioned partial
correlation for MSI, and lineage stratification for allele associations.

---

## Figures

1. **A response is three things, not two.** (a) Variance apportionment into the
   drug effect, the cell property and the cell–drug relation, each measured as a
   covariance between independent estimates so noise contributes zero
   (94.1% / 2.0% / 3.9%; PRISM, 737 lines × 1,324 compounds; log axis, since the
   drug term is 25× the other two combined). (b) A line's general sensitivity
   estimated on one half of the compounds against the estimate on a disjoint
   half, r = 0.989 — it is a property of the culture, not noise. (c) With the
   truth held at 0.20 and only the planted general response growing, the
   corrected estimator returns 0.211 / 0.212 / 0.211 where the uncorrected one
   returns 0.170 / 0.347 / 0.595. (d) The four failure modes on data with a known
   answer, each dumbbell labelled with the estimator that produced it.
   (e) Tahoe's replicate structure from the released metadata: 13.5% of
   (line, drug, dose) triples sit on more than one plate, 5.4% once plate 14 — a
   designed replicate of plate 6 — is dropped under the training convention.
   (f) True replicate versus cross-dose pairing against a matched null.
2. **Architecture where it can be measured.** (a) Shared/interaction shares in
   LINCS phase 1 and 2, PRISM, OP3, under the corrected estimator. (b) The
   residual reproduces across replicates but not across compounds or contexts.
   (c) Dose curve: rise to a peak below lethality, then collapse.
3. **What does and does not predict the relation.** (a) Cross-validated *R*² by
   predictor block — baseline expression, baseline protein, lineage, gene-level
   copy number, nonsynonymous variants — on identical compounds and lines, with
   expression and copy number capped at the same 2,000 most-variable features so
   neither wins on capacity. The copy-number advantage over mutations is
   withdrawn (+0.005, CI [−0.000, +0.010]). (b) The synonymous control, plotted
   as the per-compound nonsynonymous-minus-synonymous difference, since both
   blocks are individually negative and the claim is about their gap.
   (c) Power curve for genotype linkage versus context count. (d) Readout
   decoupling: transcription versus viability mechanism rankings. (e) The same
   question asked once per compound against that compound's own permutation null,
   with FDR across the family — a count that a few strong compounds cannot carry
   the way a median can: expression 119/120, lineage 104/120, mutations 22/120,
   with 97 compounds expression-only and none mutation-only.
4. **Cross-laboratory reproducibility.** (a) The ladder from repeat plates to a
   different laboratory, with the loss apportioned on a matched compound set
   (assay 81%, laboratory 19%). (b) Rank of the identifier-matched line by
   response similarity. (c) Reproducible fraction before and after identity
   verification, each against a ceiling measured on its own lines, beside the
   reliability-matched control that decides whether the selection is picking
   correctly-identified lines or merely well-measured ones. (d) Transfer versus
   strength of the line-specific component. (e) The transcriptional arm, which
   has no cross-laboratory test: after identity validation no (line, compound)
   pair is shared between Tahoe and LINCS.
5. **Matching strategy, and why the divergence reading was withdrawn.**
   (a) Fingerprint similarity under random, same-tissue, identifier and best-hit
   pairing. (b) Their distributions. (c) Expression rank of the line that
   outranked the identifier match: near chance, and sharing tissue only 6% of the
   time, which is why best-hit failure is read as metric noise rather than
   demonstrated culture divergence.
6. **What the measurements imply for how models are read.** (a) Baseline quality
   against the number of contexts it averages, with published benchmark context
   counts marked. (b) Few-shot transfer to an unseen line, against the honest
   floor of the additive prior plus the line's general response estimated from
   the same probe compounds — the gap between the two curves is what the model
   contributes beyond learning that a line responds strongly to everything.
   (c) Apparent model gain as a function of benchmark context count.
