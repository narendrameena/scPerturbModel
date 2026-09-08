# Context-dependence of drug response, measured: dose, transfer, genotype and chromatin

**Draft manuscript — second of two.** Split from the combined draft on 2026-09-07.
The design argument (the detection-limit calculation, the 38-dataset field survey,
the tool and the pre-registration) is now `MANUSCRIPT_DESIGN.md`. This paper holds
the measurements that do not serve that argument.

> **Status: needs its own abstract, introduction and discussion.** The Results
> sections below are carried over verbatim from the combined draft and are
> internally consistent, but the framing material has not been rewritten for this
> narrower scope. The estimator itself is described in `MANUSCRIPT_DESIGN.md`; this
> paper depends on it and should cite it rather than restate it.
>
> Several of these sections are negative or partially withdrawn results
> (§28's generalisation, the mechanism-ranking claim, the copy-number claim, two
> of three published statistics that failed to reproduce). They are kept because
> they are informative, but the framing must be honest that this is in substantial
> part a paper about what could *not* be shown.

---

## Results

### Four further ways the estimate goes wrong, none visible in the output

The interaction is naturally estimated as the covariance of response residuals
between independent replicates of the same context × compound × dose. Four
apparently reasonable shortcuts each corrupt it, and none is visible in the
output. On our data:

| shortcut | why it fails | measured cost |
|---|---|---|
| residual **variance** instead of replicate covariance | noise does not cancel | 66% "interaction" from data containing none; 72% at a true share of 20% |
| pooling **same-batch** comparisons | batch state is shared signal | within-plate pairs agree ~7× better than cross-plate; 41% against a true share of 20%, and 28% from data containing none |
| **in-sample** shared response | residuals sum to zero, forcing E[r_a·r_b] = −σ²/n | per-drug covariance driven negative for 21 of 24 drugs; clamping then yields exact zeros |
| treating **doses as replicates** | doses are different conditions | 0.5% on true replicates vs 9.2% cross-dose — an 18× inflation (Tahoe, plate 6 vs 14) |

The third has a counterintuitive property worth stating because we got it wrong
twice: the obvious repair, leave-one-*condition*-out, makes the bias *worse*,
because the mean subtracted from replicate A still contains replicate B. Only
removing the whole context eliminates the coupling. We additionally report the
interaction against a **matched cross-context null**, which removes the residual
construction offset without assuming its magnitude; where the signal is strong
this changes nothing (LINCS phase 2: 45% → 61%), and where it is weak it was the
entire result. These guards are released as `pertdecomp`. Against a standard mixed-model decomposition on the same simulated data the two
do *not* agree: they differ by 0.179 on average and in opposite directions, ours
over-estimating (mean bias +0.090) and the mixed model under-estimating
(−0.124), so at a true share of 0.3 with replicates it returns 0.045 and misses
a real interaction almost entirely. Neither is unbiased and ours errs upward,
which makes it useful for ordering designs and detecting presence rather than
for quoting a share precisely. Without replicates the mixed model still returns
a number — 0.004 for a true 0.1 and 0.005 for 0.3 — set by the optimiser rather
than the data and flagged by nothing in its output. A near-zero share reads as a
confident negative, which is the case for a tool that refuses.

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

### Transcriptional and viability context-dependence agree weakly, if at all

With every dataset on one corrected estimator, mechanism structures
context-dependence *within* Tahoe (Kruskal–Wallis *P* = 6.0×10⁻⁴ across 24
classes) and within PRISM (*P* = 2×10⁻¹⁶ across 67 classes). Whether the ordering
transfers between a transcriptional and a viability readout depends entirely on
which platform pair is used, and only one of them is adequately powered.

| comparison | classes | ρ | *p* |
|---|---:|---:|---:|
| Tahoe vs LINCS-1 (transcription vs transcription) | 10 | +0.103 | 0.78 |
| Tahoe vs PRISM (transcription vs viability) | 11 | −0.155 | 0.65 |
| **LINCS-1 vs PRISM (transcription vs viability)** | **67** | **+0.279** | **0.023** |

The Tahoe comparisons rest on 10–11 mechanism classes, at which even ρ = 0.6
would not reach significance; they are consistent with anything from strong
agreement to strong disagreement and support no conclusion. **The informative
comparison is LINCS-1 against PRISM**, which shares 67 classes — an order of
magnitude more — and shows *positive* agreement (ρ = +0.279, *p* = 0.023; on all
annotated rather than active compounds, +0.353 over 78 classes, *p* = 0.002). A
consensus transcriptional rank against viability gives +0.190 (*p* = 0.108,
n = 73).

**An earlier version of this section claimed the opposite** — "ρ = −0.09 over 67
shared classes, a well-powered null" — and a later draft restated it as a
demonstrated absence of agreement. Both were wrong, for two different reasons. The
class count was right for this pair but the correlation came from a table
generated before three of its four inputs were re-run on the corrected estimator;
regenerating it flips the sign. Separately, the *Tahoe* comparisons were quoted as
though they carried the power of the 67-class one, which they do not.

What the data now support is weak positive agreement rather than decoupling, with
one caution: six pairwise comparisons were made, Bonferroni at α = 0.05 requires
*p* < 0.0083, and +0.279 does not meet it. **The claim of readout decoupling is
withdrawn.** A transcriptional and a viability context-dependence index should
still not be treated as interchangeable without evidence, but the evidence
available now points weakly toward agreement, not away from it.

Generating table: `docs/source_data/three_platform_mechanism_cdi.csv`, re-run
2026-09-08.

### Molecular state predicts the interaction; genotype does not

At 737 PRISM cell lines the genotype scan recovers the clinical biomarker set de
novo (TP53 with MDM2 inhibitors, BRAF with vemurafenib and dabrafenib, PIK3CA
with alpelisib, KRAS with a MEK inhibitor; 80 associations at FDR < 0.05 over 1.5
million tests), confirming the estimator and the genotype join. But out of
sample, across **120** compounds with 5-fold cross-validated ridge
(`results/tables/expression_architecture.csv`):

| predictor block | median CV *R*² | compounds positive |
|---|---:|---:|
| **baseline expression** (2,000 genes) | **+0.0998** | **100.0%** |
| **baseline protein** (RPPA, 214 antibodies) | **+0.0763** | 98.3% |
| lineage | +0.0251 | 85.0% |
| copy number (2,000 genes) | +0.0075 | 63.3% |
| nonsynonymous variants | +0.0012 | 52.5% |

*Mutational burden and synonymous variants are measurable only in the separate
300-compound genotype partition (`genetic_architecture_summary.csv`), where they
score **−0.0031 / 29.7%** and **−0.0077 / 24.3%** against that analysis's own
lineage baseline of +0.0124 / 68.3%. They are reported separately rather than
appended to the table above, because the two analyses use different compound sets
and different line filters.*

> *Corrected 2026-09-08.* An earlier version of this table gave 150 compounds and
> a "compounds positive" column reading 99.2 / 88.3 / 86.7 / 57.5 / 18.3%. The
> *R*² column was right, but the count was 120 and the percentages came from a
> different analysis — 99.2% is the positive fraction of `r2_cnv_expr` (median
> +0.0801), not of the `r2_expression` printed beside it. The burden and
> synonymous rows were a third source.
>
> **The values above are read from the committed file, not from a fresh run.** A
> re-run of `expression_gap_closure.py` was attempted on 2026-09-08 and was killed
> by a 1-hour timeout at 101 of 120 compounds, so it never rewrote the table; an
> earlier note in this file claiming the re-run "reproduces the file exactly" was
> wrong and is retracted. The file itself dates from 2026-09-03 and postdates the
> estimator correction of that morning, so it is not known to be stale — but it
> has not been independently regenerated, and these numbers may move when it is.
>
> Note that `RESULTS.md` §17 reports this comparison as 4.6× on a restricted line
> set and gives +0.0927 / 92.5% for expression. That analysis is internally
> consistent but is *not* the one tabulated here, and its exact filter has not
> been traced to a committed table. The values above are the ones this paper
> uses, from a single named file.

**Genome-wide mutation status carries no generalisable information; baseline
molecular state does.** Expression beats lineage by 4.0× (+0.0998 against +0.0251) and lineage adds
nothing on top of expression, so lineage was acting as a coarse proxy for
expression state. Copy number appears to beat mutations (+0.0064, p = 5.7×10⁻⁴), as Schlüter &
Schönhuth report, **but that advantage is withdrawn**: the 120 compounds were
treated as independent when their residuals correlate at mean r = +0.23, an
effective n ≈ 11, and a compound-cluster bootstrap returns CI [−0.002, +0.015],
including zero. What survives is that copy number trails expression by 0.092 (p
= 3.3×10⁻²¹) and adds nothing to it — the joint block scores below expression
alone. The
absolute effect remains modest (≈9% of variance) but is highly consistent across
compounds. Using synonymous variants as a
control — silent changes cannot alter a protein but carry identical ancestry,
lineage and germline-contamination structure, and we match both blocks to the
same 3,435 genes — isolates a genuinely mechanistic excess of **+0.0036**
(*P* = 0.025), about 0.3% of the interaction variance.

Subsampling the confirmed PRISM associations quantifies why line-level tests fail
in smaller atlases: recovery is 5% at 47 cell lines, 40% at 250, 78% at 400 and
100% at 600 (the 600-line row is over 9 associations rather than 12, the rest
over 12). The genotype negatives reported in 47-line atlases are a power
ceiling, not a biological absence.

Screening 111,589 allele × compound tests for associations that are *not* known
pharmacogenomics and are invisible to a gene-level indicator yields four
candidates; none survives validation in GDSC, and all three testable ones are
identified as **germline** by a single shared genomic position, no excess
mutational burden, and neighbouring recurrent variants that are synonymous. The
same pipeline recovers 11 of 11 known biomarkers in GDSC (BRAF V600E ×
dabrafenib, *P* = 3×10⁻³¹), so the negative is informative.

### Cross-laboratory transfer is limited by the assay and by cell-line identity

Against within-laboratory ceilings of r = 0.685 (PRISM replicate plates) and
0.365 (GDSC1 vs GDSC2), the line-specific response transfers between institutions
at r = 0.290 — **57.9% [49.0–77.5%] of what the assay achieves with itself**. The
ceiling is itself low: even repeating a measurement in one laboratory, the line's
compound-specific deviation agrees at only r ≈ 0.5.

Cross-atlas comparisons assert that a line in one atlas is the same as a line in
another, from an identifier that is never checked. Verifying it — each line must
be its own best match among all candidates by **response fingerprint**, its
residual across the compounds both atlases share — shows the assertion usually
fails: of 488 COSMIC→DepMap-matched lines, **only 5–12% are their own best
match**, ranking a median 82nd of 971 candidates. Pairing each line with its
single best fingerprint match instead of its identifier raises agreement to
r = 0.419, against 0.235 for identifier matching and 0.048 for a same-tissue
random line — an upper bound the atlas's own design does not let us reach. That is far better than chance
(~486), so identity carries real information; it is simply not unique.

**Verifying identity raises transfer from 57.9% to 98.0% of the ceiling**
[47.0–116.2%], measured against a ceiling computed on the same lines and
compounds; against the all-lines ceiling it reads 110%, which a fraction of a
ceiling cannot be. This is not circular: identity is validated on one random half
of the shared compounds and agreement measured on the disjoint other half. A
reliability-matched control — lines selected for how well they are measured
rather than for identity — reaches only 49%, so the gain is specific to
identity. We initially attributed this to detecting divergent
cultures; baseline expression does not support that (below), so the gain is
reported as a practical filter whose mechanism is unresolved.

A matching ladder shows what each rule buys: random pairing r = −0.002,
same-tissue random pairing 0.048, identifier 0.235, best available partner 0.419.
Identifier matching beats same-tissue pairing decisively (*P* = 1.5×10⁻¹²⁴), so
identity carries information well beyond lineage — but the best available partner
reaches nearly twice what the identifier finds.

Apportioning the loss across the ladder — repeat plates (0.685) → different assay
version within one laboratory (0.365) → different laboratory and assay (0.290) —
**changing assay costs 81% of the total drop and changing laboratory the
remaining 19%** (CI on the laboratory share [−6%, 38%], which includes zero). The
apportionment must be computed on compounds all three rungs share: on the 1,435 /
123 / 187 different compounds each rung happened to cover, it inverts to 84%
laboratory, which is what an earlier draft reported and what its title asserted.

Transfer scales with signal strength (ρ = +0.496, *P* = 1.6×10⁻¹⁰, n = 147): the
strongest quartile of compounds transfers *at* the within-laboratory ceiling, the
weakest at about a third of it. The middle two quartiles are within noise of each
other, so this separates strong from weak rather than ordering all four.

The result is a property of laboratories rather than of a killing assay, though
transcription transfers considerably worse than viability. On the matched set —
114 within-lab pairs sharing a (line, compound) with the 317 cross-lab pairs —
the within-laboratory ceiling is *r* = 0.044 and cross-laboratory agreement
0.012, a **reproducible fraction of 27.3%** against viability's 56%. The identity
check behaves the same way: within the Broad (LINCS phase 1 vs phase 2) **15 of
16** name-matched lines are reciprocal best matches (median rank 1 of 26), but
across laboratories only **2 of 6** (median rank 4 of 62).

The transcriptional ceiling is itself only *r* ≈ 0.034–0.044, so this residual is
barely reproducible even within one laboratory and 27.3% is a ratio of two very
small numbers. Shen (`rs-10846736`, 2026) reaches a compatible conclusion on
Tahoe alone, by measuring per-drug repeat reliability directly (median 0.067).

*These figures were corrected on 2026-09-08. An earlier version reported 46%,
within-lab r = 0.061 and 16 of 16 identity matches; those came from a keying bug
in which an inner loop overwrote the compound key with a cell line, so the two
datasets shared no keys and every cross-laboratory comparison silently returned
zero pairs. Fixed in `383e890` with an AST regression test. Generating tables are
committed to `docs/source_data/`.*

### The same omission inflates another group's statistic, in another modality

A defect in our own code becomes a methodological finding only if the same term
sits inside statistics computed elsewhere. We tested three published
context-specificity statistics on the publishing group's own data, under the rule
that the published number must be reproduced first — a correction applied to a
statistic we cannot reproduce measures our reimplementation, not their result.

Two failed that test outright. Our implementation of the sci-Plex classification
(Srivatsan et al. 2020) gives 74.4% cell-type-dependent against a published 48%,
and our implementation of the CMap panel-conservation statistic (Subramanian et
al. 2017) gives 5.0% against a published 26%. Neither reproduces, so neither says
anything about those papers; and within our own versions the correction moves the
statistic by only 1.01× and 1.15×, the latter in the opposite direction.

The informative case is TRADE (Nadig et al. 2025), whose deposit carries
**lfcSE beside every log2 fold change** for four cell lines given the same 2,052
essential-gene perturbations over 6,642 shared genes. The standard errors make a
noise-corrected version constructible, and four cell lines make the consistent
component estimable as a between-line covariance, which noise cannot inflate
because two lines are measured independently. Both forms are computed on
identical data with the identical correction:

| statistic | cell-line effect left in | removed | shift |
|---|---:|---:|---:|
| raw (noise included) | 61.5% cell-type-dependent | 58.1% | +3.4 points |
| **noise-corrected** (covariance + published SEs) | 66.3% | **54.6%** | **+11.7 points** |

The noise-corrected split moves **3.4× more**, a 1.21-fold inflation of the
cell-type-dependent share. Sampling noise is 38% of the raw per-line variance
here, which is why the raw form is comparatively inert — and why the two
statistics above, computed on raw residuals, barely moved.

**The correction must be a contrast across contexts, and this is easy to get
wrong.** A first version subtracted each line's uncentred mean response across
the other perturbations. Every perturbation in this deposit is an essential-gene
knockdown, and knocking down any essential gene drives a common growth-arrest
programme, so those four vectors correlate **with each other at r = +0.51** and
with the cross-line shared response at **r = +0.80**: one shared programme, not
four line properties. Subtracting it removed signal the statistic is meant to
keep and drove the split the wrong way. Centred across lines, the same vectors
correlate at **r = −0.53** between lines and −0.19 with the shared response, and
behave as context properties should. Any application of this correction should
report that diagnostic.

Our raw replication gives 39% consistent against TRADE's published 56%, so this
does not restate or correct their number. What it establishes is that the effect
is not peculiar to us, to viability, or to chemical perturbation.

### Drug exposure writes a new axis rather than reorganising the programme

The decomposition treats the drug effect and the cell's own state as separable
terms. Whether that is geometrically reasonable is testable: does exposure move
cells along the directions in which untreated cells already differ, or into
directions none of them occupies?

The baseline programme space is the span of principal components of **untreated**
pseudobulk across the Tahoe cell lines. Every component here is an inner product
between deltas measured on **two independent plates**, because measurement noise
is close to isotropic and lands outside any low-dimensional span — a purely noisy
response would otherwise read as pure remodelling. From single observations the
shift share is 8.5%; across plates it is 10.1%.

The reference that makes the number readable is a **held-out untreated cell
line**, a signal that is a shift by construction since it *is* line-to-line
variation. Over 5,345 replicated conditions, with 36 lines fitting the span and
12 held out:

| baseline directions (k) | captures of a drug response | captures of a held-out untreated line |
|---:|---:|---:|
| 1 | 1.7% | 9.7% |
| 5 | 6.2% | 22.1% |
| 10 | 9.0% | 25.0% |
| 20 | 12.2% | 29.6% |
| 35 | 15.7% | 32.8% |

Responses are captured **two to five times less at every k**. Exposure does not
move cells along the axes that separate untreated lines; it adds one.

Three controls fix what kind of change this is. It is **not cytotoxic collapse**:
the off-axis share is 89%, 86% and 90% at 0.05, 0.5 and 5 µM — flat, and already
at full size at the lowest dose, where cells are not dying. **Cell identity is
preserved**: between-line spread, treated against the same lines untreated, has
ratio 1.002 / 1.000 / 1.007 across the three doses (all *P* > 0.5), so lines stay
exactly as far apart as they were. **Co-expression is untouched**: the gene–gene
correlation structure across lines is essentially unchanged by treatment
(ρ = +0.998 over 400 genes; mean absolute change 0.0127 against a
label-permutation null of 0.150 [0.128, 0.173]) — a smaller change than any
random split of the same samples produces.

So exposure remodels expression, but additively: a new drug-specific axis laid on
top of an unchanged programme, rather than reprogramming or convergence onto a
common stress state. This is also the geometric reason the decomposition
separates so cleanly — the drug axis is close to orthogonal to the axes on which
cells differ, which is how a drug main effect can account for 94.1% of the
variance while the cell property and the cell–drug relation remain small and
separable.

### Chromatin, and what a decomposition over all features costs

Two objections apply to everything above. Transcription and viability are both
readouts of one cell state, so their agreement is weaker evidence than it looks;
and every number so far comes from an estimator built for this paper. Spear-ATAC
(Pierce, Greenleaf et al. 2021, via scPerturb) addresses both: **fully crossed** —
all 41 CRISPR perturbations in all three cell lines (K562, GM12878, MCF7), 4–6
independent replicate samples each, 73,344 cells — genetic rather than chemical,
and in a different molecular layer.

The experiment supplies its own positive control. Each perturbation targets a
transcription factor, and most of those factors have a motif in the ChromVar
feature set, so the knockout predicts a fall in that one matched motif. **Six of
72 (line, perturbation) pairs reject at Bonferroni *q* < 0.05 and all six move
down** — GATA1 (Δ = −4.19), NFE2 (−1.89), FOSL1 (−1.52), KLF1 (−0.62) and CEBPB
in K562, NFYB in GM12878, which are the canonical K562 erythroid regulators.

Yet an energy-distance E-test (Peidli et al. 2024) computed over **all** 2,174
motifs rejects for **none** of those 123 pairs, and the same is true of gene
scores over 24,919 genes. One motif moving four standard deviations is invisible
to a statistic that averages over two thousand that did not move. The failure is
dilution, and it applies with equal force to a variance decomposition run over
every feature.

That cost can be measured. Ranking features on one disjoint half of the replicate
samples by their **perturbation main effect**, and computing the decomposition on
the other half — the two quantities are orthogonal by construction and the halves
share no measurement — gives:

| features included | residual (noise) | interaction / reproducible |
|---|---:|---:|
| all 2,174 | 89.1% | **0.3%** |
| 100 responsive, selected out of fold | 64.2% | **91.6%** |

Identical data, identical estimator; only the feature set differs, and the answer
moves from 0.3% to 91.6%. **A context-dependence index computed atlas-wide is
largely reporting the features that carry no signal**, which is a general caution
about how these indices are built, this paper's own included.

On those 100 responsive features the interaction variance is 33.9%, against a
perturbation-label permutation null of **33.0% [24.2–39.5%], *P* = 0.45**. The
result is therefore a null — but a meaningful one, because the main effects are
demonstrably present, the responsive features were identified without
circularity, and the interaction is still at chance. That is consistent with
Tahoe's 0.5% [0.0–1.5%] in transcription. With three contexts the power is
limited and the null interval spans 24–39%, so this is *no evidence of an
interaction* rather than evidence of none.

The methodological claim gains support that does not depend on any of this.
`variancePartition` (Hoffman & Schadt 2016) is the standard variance-component
method in genomics and its model is exactly the split this paper proposes,
`y ~ (1|context) + (1|perturbation) + (1|context:perturbation)`. Its interaction
term is identifiable **only because the design is replicated** — without
replicates it and the residual are one stratum — yet it is routinely fitted to
unreplicated designs.
