# Why each method was chosen

One entry per analysis step: the choice made, the alternatives considered, and
why they were rejected. Where a choice was later found to be wrong, the error is
kept rather than overwritten — the correction is part of the record.

---

## 1. Pseudobulk deltas rather than single-cell modelling for the main analysis

**Chosen.** Mean log1p-CPM profile per (cell line, drug, dose, plate), minus the
plate-matched DMSO profile of the same line.

**Alternatives rejected.**
- *Single-cell likelihood models (scVI/cVAE) as the primary analysis.* Run as
  models 1–2, but not as the measurement instrument: their outputs depend on
  latent dimensionality, likelihood and training seed, so a variance
  decomposition built on them would confound biology with fitting choices.
- *Raw counts without control subtraction.* Cell lines differ enormously at
  baseline; without subtracting the matched control the between-line variance
  swamps the drug effect entirely.

**Why this one.** The question is about the *response*, which is a difference.
Subtracting the same line's own control on the same plate removes both the
baseline identity and the plate offset in one step, with no fitted parameters.

---

## 2. Plate-matched controls, and cross-plate pairs only

**Chosen.** Controls are matched within plate; every replicate comparison uses
two *different* plates.

**Alternatives rejected.**
- *Pooling all replicate pairs.* This is what the first version did, and it was
  wrong: within-plate pairs agree ~7× better than cross-plate pairs (r 0.45 vs
  0.06), so pooling inflated the interaction estimate from 28% to 41%. The
  correction is recorded in RESULTS.md §4.
- *Global (not plate-matched) DMSO reference.* Leaves a plate offset in every
  delta, which then reappears as spurious "context" structure.

**Why this one.** Batch structure survives plate-matched normalisation. Only a
comparison across plates can distinguish reproducible biology from shared
technical state.

---

## 3. Interaction estimated as a covariance between replicates, not a variance

**Chosen.** interaction = mean over independent replicate pairs of
`E[residual_a × residual_b]`, clamped at zero.

**Alternatives rejected.**
- *Variance of the residual.* Contains noise, which does not cancel; it
  overstates interaction by exactly the noise level, and the overstatement grows
  as cells per condition falls.
- *ANOVA / linear mixed model variance components.* Would work, but assumes
  homoscedastic Gaussian residuals across 2,000 genes and thousands of
  conditions, and gives no natural per-perturbation estimate.

**Why this one.** Independent replicates share signal but not noise, so the
cross-replicate covariance estimates the reproducible interaction directly, with
noise cancelling in expectation rather than being modelled. It also degrades
gracefully: with no replicates the estimator is simply refused (the
`pertdecomp` tool raises rather than reporting a number it cannot support).

---

## 4. Leave-one-context-out additive prior

**Chosen.** The shared response for a line is the mean over all *other* lines.

**Alternative rejected.** *Mean over all lines including the target.* The target
line contributes 1/n of its own prediction, which leaks and inflates the
apparent additive share — badly at small n, which is exactly the regime of the
context-count scaling experiment (§2).

**Why this one.** It is the honest version of the same quantity and makes the
5→45-line comparison interpretable.

---

## 5. Context-dependence index as a ratio, with the inert-compound caveat

**Chosen.** CDI = interaction / (additive + interaction).

**Alternatives considered and reported alongside** (`compare_context_metrics.py`).
- *Cross-context profile correlation*, which is what the published literature
  mostly uses. Rejected as the primary metric because it confounds
  context-specificity with the compound's own reproducibility: a weak or noisy
  compound scores "context-specific" for the wrong reason.
- *Disattenuated transfer* (cross-context r ÷ replicate r). Kept as a reported
  cross-check; it is arguably the most defensible published-style metric.

**Known weakness, stated rather than hidden.** A compound that does nothing has
a near-zero numerator *and* denominator and lands at 0 by default. LINCS phase 1
has 35 mechanism classes sitting at exactly 0.000 for this reason. Every
mechanism comparison is therefore repeated on compounds with a detectable shared
effect (§12); the conclusions do not change, but the caveat is real.

---

## 6. Mechanism labels from the Broad Repurposing Hub

**Chosen.** Hub MOA annotations, name-normalised.

**Alternative rejected.** *Tahoe's own drug metadata, matched by name.* Tried
first; it annotated only 23 of 273 LINCS compounds, which is too thin to rank
mechanisms.

**Why this one.** One vocabulary across all four datasets, so rankings are
comparable rather than each dataset being scored against its own labels.

---

## 7. PRISM for the genotype question

**Chosen.** PRISM Repurposing viability, 738 lines.

**Alternatives rejected.**
- *More cells from Tahoe.* Does not help. The genotype question is line-level,
  so power is set by the number of lines (47), not the number of cells (100 M).
- *Additional transcriptional atlases.* None exist with hundreds of genotyped
  contexts; LINCS phase 1 reaches 70 lines, still far short.

**Why this one.** It trades the transcriptional readout for a ~15× increase in
context count, which is the axis that was limiting. The trade is stated
explicitly, because a viability-linked effect need not be transcription-linked.

---

## 8. Cross-replicate rather than cross-dose agreement in PRISM

**Chosen.** Covariance between X1/X2/X3 detection plates at matched compound
*and* dose.

**Alternative rejected — after first being used by mistake.** The first PRISM
pass used agreement *across doses*, which put the interaction at 5% of
reproducible variance. That is wrong twice over: doses are not replicates (a
line's response genuinely differs between doses), and the many no-effect low
doses enter as pure noise. Switching to the replicate axis moved the estimate to
22%, in line with the other datasets.

**Why this one.** It is the same estimator used everywhere else in the project,
so PRISM is comparable to Tahoe and LINCS rather than being a separate method.

---

## 9. Subsampling to a power curve rather than quoting a power calculation

**Chosen.** Take the associations confirmed at 738 lines, resample to
20…600 lines, and count how often each is recovered at the same per-test
threshold.

**Alternative rejected.** *An analytic power calculation.* Requires assuming an
effect-size distribution and a test statistic; the assumption would be doing all
the work, and reviewers would rightly discount it.

**Why this one.** It measures recovery on the real effect sizes, the real
allele frequencies and the real noise, and it answers the exact question asked:
at Tahoe's context count, would we have seen this?

---

## 10. Vectorised Mann–Whitney instead of looping scipy calls

**Chosen.** Rank once per compound, then one matrix multiply against a
(gene × line) membership matrix.

**Alternative rejected.** *`scipy.stats.mannwhitneyu` per test.* Statistically
identical at these sample sizes (scipy uses the same normal approximation), but
1.5 M calls made the subsampling loop unaffordable.

**Why this one.** Same statistic, ~100× faster, which is what made the power
curve possible at all.

---

## 11. Within-compound trend testing for the dose question

**Chosen.** Spearman of CDI on dose *within each compound*, then a Wilcoxon
signed-rank test over compounds.

**Alternative rejected.** *Pooling all (compound, dose) points into one
regression.* Between-compound differences in CDI are large and would dominate;
a compound tested only at high dose would masquerade as a dose effect.

**Why this one.** Each compound is its own control, so the question asked is
strictly "does raising the dose of *this* compound change *its*
context-dependence".

---

## 12. Allele-level testing with an MSI burden control

**Chosen.** Per-protein-change indicators (`BRAF p.V600E`), tested against the
same residuals as the gene-level scan, then a rank-based *partial* correlation
conditioning on each line's total frameshift burden.

**Alternatives rejected.**
- *Gene-level indicators alone.* `BRAF mutant` pools 104 lines of which only 50
  carry V600E; the rest dilute the signal. Eleven compounds show a significant
  V600E association whose gene-level test fails.
- *Quintile matching on burden.* Tried first and discarded: carriers span
  several burden quintiles, so restricting controls to "the carriers'
  quintiles" excluded almost nothing and every hit trivially survived. A control
  that cannot fail is not a control.

**Why this one.** Frameshifts in homopolymer runs (RPL22 K16fs, ACVR2A K437fs)
are microsatellite-instability markers that co-occur in the same lines, so a set
of "independent" frameshift hits can be one MSI association counted many times.
Partial correlation on burden is the test that can actually distinguish them.

---

## 13. Reporting negatives as power statements, not as biology

**Chosen.** Where a line-level predictor fails at n≈47, the claim made is "not
detectable at this context count", never "does not exist".

**Why.** §11 shows the distinction is not academic: effects that are
unambiguously real (BRAF→vemurafenib) are recovered 4% of the time at 47
contexts. Any of the four failed predictors in §10 could be in the same
position.

---

## 14. Matched cross-context null instead of assuming the offset is small

**Chosen.** interaction = (mean covariance of same-context replicate pairs) −
(mean covariance of cross-context pairs, matched on perturbation and batch).

**Alternative rejected.** *Taking the raw covariance and clamping at zero* — what
the code did originally. Residuals are formed against a mean estimated from a
finite number of other contexts, so the mean subtracted from one residual still
contains the other; every pair therefore carries a negative offset of order
−2σ²/(n_ctx−1). Measured: Tahoe's cross-line offset is −0.0022 against a shared
component of 0.059, large enough to drive the raw per-drug covariance negative
for 21 of 24 drugs, which clamping then converted into exact zeros.

**Why it matters more in some datasets than others.** LINCS phase 2: raw
covariance +0.359, offset −0.030 — the correction moves the share 45%→47%, i.e.
nothing. Tahoe: the offset is comparable to the signal. The correction is
therefore not cosmetic where it matters and not disruptive where it does not,
which is the behaviour a bias correction should have.

## 15. Requiring same-dose replicates, and reporting when they do not exist

**Chosen.** Only pairs sharing (context, perturbation, **dose**) on different
plates count as replicates.

**Alternative rejected.** *Treating different doses of the same (line, drug) as
replicates* — the original pairing. Measured on identical (line, drug) sets, the
two disagree at p = 3×10⁻¹²⁰ and in the wrong direction: true replicates covary
at −0.00307, cross-dose pairs at +0.00620. Replicates agreeing *less* than
different doses cannot be reproducible interaction.

**Consequence, reported rather than hidden.** This leaves Tahoe with 6,482
usable pairs from 25 of 379 drugs, which is not enough. The honest output is a
bracket (0% to 21.6%) plus the statement that the atlas cannot resolve it, not a
point estimate from whichever pairing gives a publishable number.

## 16. Screening known pharmacogenomics out before claiming novelty

**Chosen.** An explicit gene→mechanism table of established relationships
(BRAF→RAF/MEK/EGFR, TP53→MDM2, PIK3CA/PTEN→PI3K/AKT, …), applied before any
"novel" claim.

**Alternative rejected.** *Reporting the top of the FDR-sorted list.* It is
dominated by BRAF V600E, which is a positive control.

**And a control that was applied and then distrusted.** Surviving candidates
passed split-sample replication at 94–99%, lineage stratification and MSI-burden
partial correlation — and are still reported as **negative**, because (a) three
sit in MUC16 and HMCN1, ranked 2nd and 21st of 18,739 genes by mutation count,
where recurrent variants are long-gene artefacts, and (b) the split-sample test
is circular, since candidates were selected using all the lines it then splits.
A validation that cannot fail is not a validation — the same lesson as the
quintile matching in §12. Real validation needs PRISM primary, GDSC or CTRP.
