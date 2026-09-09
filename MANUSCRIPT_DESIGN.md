# Replicates, not cells: a design calculation for perturbation atlases

**Draft manuscript.** Every number is reproducible from this repository. Figure
bundles (PNG/SVG/PDF + source data + generating script) are under
`results/figures/`; each analysis ran as a recorded SLURM job in `jobs/`. Full
result set in `RESULTS.md`; method choices justified in
`docs/methodology_rationale.md`; the pre-registration is `docs/PREREGISTRATION.md`
(tag `prereg-2026-09-07`).

> **Scope note.** This paper was split out of a longer draft covering the same
> project's measurements of dose response, cross-laboratory transfer, genotype
> prediction, chromatin and expression remodelling. Those become a separate paper
> (`MANUSCRIPT_MEASUREMENT.md`); nothing is discarded. What is kept here is the
> design argument and exactly the material load-bearing for it — chiefly the
> estimator, since a design calculation is only meaningful if the quantity it
> designs for is estimated correctly.

---

## Abstract

Perturbation atlases profile thousands of compounds across dozens of cellular
contexts to learn how context shapes drug response, and are scaled by cell count.
Tahoe-100M sequenced 95.6 million cells and still cannot resolve the interaction
it was built to measure; **the same budget, spent as six replicates of 301 cells
rather than two of 1,810, would have detected an effect three times smaller than
its own.** Whether an atlas can measure context-dependence is fixed by its design,
not its scale. Because the quantity is a covariance between independent
replicates, its precision is set by replication rather than depth: a condition
measured once contributes nothing however deeply it is sequenced. From context,
perturbation and replicate counts alone, the calculation predicts which of five
published atlases resolved an interaction and which could not, correctly in all
five, with one shared noise constant and no per-atlas tuning. Testing its
precision law on 8,427 LINCS compounds across a 50,000-fold range of sample size
refuted its original form — pairs within a condition share profiles, so the
estimator is a U-statistic whose variance is `a/k + b/k²`, not `b/pairs`; the
corrected form is fitted on half the compounds, validated on the other half, and
leaves all five verdicts intact. The premise holds where tested: 2% of
Tahoe-100M's cells retain 92% of a known interaction, while 10% of its contexts
destroy precision. Applied to all 38 scPerturb datasets, **one can support the
estimate at all**. We release the calculation as a tool and pre-register it.

---

## Introduction

Single-cell and pooled perturbation atlases measure the consequences of thousands
of chemical interventions across dozens to hundreds of cellular contexts,
motivated by predictive "virtual cell" models. The difficulty such models face is
*context transfer*: predicting a compound's effect in a cell state it was not
profiled with. Progress is judged almost entirely by prediction accuracy, and the
picture is contested — benchmarks report that deep models do not beat simple mean
baselines, while method papers report substantial gains.

Underneath that dispute sits a prior question that is not asked. The quantity
these models are trying to learn — the part of a drug's effect specific to a
context rather than shared across contexts — is an interaction term, and
interaction terms have detection limits. Whether a given atlas can estimate one is
a property of how it was built, decidable before any data is collected. That
calculation is not published anywhere we can find, and the designs of the field's
atlases suggest it is not being done. Tahoe-100M spends 95.6 million cells and
replicates 13.5% of its conditions, which leaves it a detection floor of 0.0169
against a true interaction of 0.005 — it misses its own effect by 3.4×. The same
95.6 million cells spread as six replicates of 301 cells give a floor of 0.0018,
which would have resolved it with nearly 3× margin. Nothing about the biology or the
budget changed; only how the cells were spread. Spear-ATAC, at the other extreme,
replicates almost everything across three cell lines and still cannot resolve an
interaction, because three contexts is too few however often each is repeated.

We give the calculation, validate it against five published atlases whose outcomes
we measured independently, test its central premise by subsampling, apply it to
an entire public collection, and release it as a tool. The argument requires one
piece of estimator work, included here because it is load-bearing: a
context-dependence estimate that omits a cell line's general sensitivity is
inflated by roughly 1.5×, so the design calculation would otherwise be
calibrated against the wrong quantity.

---

## Results

### What an atlas must measure, and how to know before building it

The interaction is estimated as a covariance between independent replicates, so
its precision is set by the number of replicate **pairs** —
`n_ctx × n_pert × n_rep(n_rep−1)/2` — and by per-observation noise. **Cell count
does not appear.** A condition measured once contributes no pair however deeply it
is sequenced.

That turns "can this atlas answer the question" into arithmetic. Given only each
atlas's context, perturbation and replicate counts, and told nothing about what
any of them found:

| atlas | replicate pairs | smallest detectable share | observed | predicted | actual |
|---|---:|---:|---:|---|---|
| Tahoe-100M | 4,560 | 0.0169 | 0.005 | **not resolvable** | not resolvable |
| LINCS phase 1 | 177,003 | 0.0027 | 0.57 | resolvable | resolvable |
| OP3 | 2,646 | 0.0222 | 0.331 | resolvable | resolvable |
| sci-Plex 3 | 567 | 0.0479 | 0.302 | resolvable | resolvable |
| Spear-ATAC | 1,230 | 0.0325 | 0.014 | **not resolvable** | not resolvable |

**Five of five**, from three integers each, before any data is examined. A
**single shared noise constant** is used for all five rather than a per-atlas
value: five free parameters fitting five binary outcomes would prove nothing.
The result holds for any shared constant between 0.10 and 0.30, a threefold
range, degrading to 4/5 at 0.35.

The prescription follows directly. **Tahoe needed six replicates per
condition, not two** — at two its floor is 0.0169 against a true 0.005, and even a
third replicate only reaches 0.0097, still short. Spear-ATAC, at three cell lines
and five replicates, sits at 0.0325 against a true 0.014 and would have needed
eleven.

**The premise is testable and holds.** Binomially downsampling Tahoe's counts
while keeping every condition, a biologically specified interaction — MEK
inhibitors suppressing the Pratilas ERK-output signature further in BRAF/RAS-driven
lines than in wild-type ones — retains **92% of its size at 2% of the cells**, a
fiftyfold reduction on a flat, non-decaying curve. Thinning contexts instead
leaves the mean intact but collapses precision: at 10% of contexts the spread
across draws is 0.084 against an effect of 0.10. Tahoe would have obtained the
same answer from about 2 million cells, and the remaining budget spent on
replicates would have moved it across its own detection threshold.

The readout pools across ~120 lines, so per-condition sampling noise averages out
before the contrast is taken; a single-condition estimate would degrade under cell
thinning sooner, and this experiment does not measure how much sooner.

*Calibration.* Against simulation with a known interaction, power above the
predicted threshold is 83–100% and the false-positive rate 0–8% on larger designs;
on the smallest design tested the false-positive rate reaches 17%. The threshold is
also not sharp from below: on some designs an effect *beneath* the predicted floor
is still detected up to 75% of the time, so the floor is conservative rather than
exact. The minimum detectable share is therefore an order-of-magnitude guide, which
suffices here because the five atlases differ by two orders of magnitude in what
they can resolve — but it should not be read as a precise boundary.

### One dataset in thirty-eight can support the estimate at all

Five atlases are few, and they were the five we happened to analyse. We therefore
applied the calculation to **every RNA and protein dataset in scPerturb** (Peidli
et al. 2024; *n* = 38), reading design parameters from deposited metadata without
touching an expression value, using the same shared constant and no per-dataset
tuning. As a check on the reader, our hand-entered row for sci-Plex 3 (3 contexts
× 189 perturbations × 2 replicates) is exactly what independent parsing of the raw
`.h5ad` returns.

| | datasets |
|---|---:|
| scPerturb RNA/protein datasets | 38 |
| more than one context | 4 |
| …and ≥ 1 condition measured twice | 1 |
| …and able to detect a 5% interaction | **1** |

**One of thirty-eight can support a context × perturbation interaction estimate** —
sci-Plex 3, at a floor of 0.048 against its measured 0.302. Twenty-six carry
replicate-like annotation, but in 23 the replicates never cover the same condition
twice, which yields no pair.

A CRISPR screen in one cell line is single-context *by design*, and scoring it as
a failed atlas would be unfair. Splitting by perturbation type leaves the claim
narrower and sharper: of **27 CRISPR datasets none use more than one context**,
while of **ten drug screens three do and one of those replicates a condition**.
The datasets a reanalyst would reach for to ask whether drug response depends on
cell context mostly cannot answer it — not because the effect is absent, but
because the design does not admit the estimate.

**Two caveats bound this.** The deposited metadata cannot distinguish an
independent replicate from a split capture: `batch` denotes a genuine replicate in
some datasets and a 10x capture lane in others. Replogle's ~40 `batch` levels are
gem groups — one transduced pool distributed across captures, sharing
transduction, culture and selection — so that study's union of K562 and RPE1 does
not qualify despite appearing to. Reported replicate counts are therefore an
**upper bound on usable replicates**, and so is the count of qualifying datasets;
the bias runs conservative for the claim. Second, this counts what is recoverable
from deposited metadata, not what was run at the bench. The practical consequence
survives that distinction, because an index computed by a reanalyst — or by the
original authors from the released object — can only use annotation that is
present, and the usual substitute where it is absent is to split one well's cells,
which returns 0.500 from data containing no interaction.

Files are also not studies: scPerturb distributes `ReplogleWeissman2022` and
`TianKampmann2019` as per-context files, so single-context datasets are 34/38
(89%) per file and 17/22 (77%) per study. The qualifying count is one under either
unit.

**Independent convergent evidence.** While this work was in preparation, Shen
(Research Square, 2026-09-01, `rs-10846736`) reported a leakage-safe
leave-one-cell-line-out benchmark on Tahoe-100M and reached a compatible
conclusion by a different route: across 107 replicated drugs the context-specific
residual has a median repeat reliability of 0.067, and model skill tracks that
ceiling (ρ = 0.76), so prediction is "bounded primarily by the reproducibility of
the target signal". That is the same diagnosis as ours for Tahoe specifically, and
it is reassuring that two independent analyses reach it.

The contribution here is what that study does not do. It **measures** reliability
post hoc, on one atlas, from data already collected; we **predict** resolvability
from three design integers before any data exists, validate that prediction across
five atlases and 38 datasets, and release it as a calculator an atlas builder can
run while the design is still changeable. Shen's conclusion — that benchmarks
should report repeatability ceilings — is a diagnostic recommendation; ours is a
prescriptive one, and the two are complementary rather than competing.

**Prospective test.** Predictions about scPerturb's ATAC collection were
registered publicly and pushed *before* the unseen archives were downloaded
(`docs/PREREGISTRATION.md`, tag `prereg-2026-09-07`; the git history shows the
order). The one informative prediction — that `MimitouSmibert2021` would not
support the estimate — held: CD4+ T cells only, 6 perturbations, one replicate,
zero pairs. This is *n* = 1 and is reported as one dataset, not as a confirmation.
A first draft that treated the whole ATAC article as unseen was withdrawn on
discovering three of its files had already been analysed; the withdrawal is
recorded in the registration. A standing registration commits the frozen analysis
script, by SHA-256, to the first three qualifying atlases released after
2026-09-07, with the falsifiable assertion that **an atlas with one replicate per
condition will not resolve a reproducible interaction however many cells it
sequenced.**

### The quantity being designed for, and why the obvious estimator is wrong

A design calculation is only meaningful if the quantity it designs for is
estimated correctly. A response measured in a cellular context contains three
terms, not two:

    y(context, drug, dose) = β(drug, dose) + α(context) + γ(context, drug) + noise

Only γ is context-dependence. α — a context's **general sensitivity**, how it
responds to every compound, set by growth rate, seeding density and drug
metabolism — is almost always left inside the residual, because the standard
construction removes the perturbation mean and nothing else. Whatever makes one
context respond strongly to everything then reads as that context responding
*differently*.

That term is not noise. Estimated on one half of a compound library it predicts
the estimate on a disjoint half at **r = 0.989** across 737 PRISM cell lines, and
**570 of 736 lines (77.4%)** carry a reproducible line-level shift at FDR < 0.05.
Because it is identical across replicate plates, it survives every noise-cancelling
device an estimator applies and lands directly in the replicate covariance that
validates interactions. Ben-David et al. (2018) name the same axis: their most
resistant MCF7 strains are resistant *in general*, through downregulated
drug-metabolism pathways, not through anything drug-specific.

Estimating all three as covariances between independent estimates, so noise
contributes zero to each:

| component | share | variance |
|---|---:|---:|
| drug effect (identical in every context) | 94.1% | 1.6338 |
| **cell property** (general sensitivity) | **2.0%** | 0.0343 |
| **cell–drug relation** (interaction) | **3.9%** | 0.0683 |

PRISM, 737 lines × 1,324 compounds. **Interaction estimates that omit the middle
term are inflated by roughly 1.5×.**

Two rules make the estimates honest. Both main effects must be estimated **out of
fold** — β from other contexts, α from other perturbations — since an α computed
with the perturbation included absorbs the interaction and biases it toward zero.
And α must be estimated **within a replicate batch**: pooled across batches it
carries a share of each batch's control noise and subtracts var(control)/2 from
the pair covariance. In simulation that error drove the estimate to exactly
**0.000** while appearing correct in review. With a planted general response of 0,
0.5 and 1.0 SD the corrected estimator returns 0.211, 0.212 and 0.211 against a
truth of 0.20, where the uncorrected one returns 0.170, 0.347 and 0.595. On
simulated data with a known share the corrected estimator recovers it (slope
0.950, R² = 0.9995) and returns 0.008 from data containing none, while residual
variance reports 72% and pooled-batch 41%.

### The replicate structure exists, and discarding it doubles the estimate

Tahoe-100M's own design illustrates the cost of the calculation not being done.
**13.5% of its (line, drug, dose) triples sit on more than one plate** — the
replicate structure the estimator needs exists in the released data. Under the
training convention that drops plate 14, a designed replicate of plate 6, that
falls to **5.4%**.

Pairing across doses instead of across true replicates — the natural substitute
when replicates appear unavailable — inflates the estimate from **0.5%
[0.0–1.5%], *P* = 0.10** at matched dose to **9.2%**: from a value not
distinguishable from zero to a clearly non-zero one. Two doses of the same
compound in the same line share the compound's dose–response shape, and that
shared shape enters the covariance as if it were context-specific. (Both figures
are from the corrected estimator. Under the uncorrected one the same substitution
reads 11.5% → 20.7%, a doubling; the correction moves both levels but not the
direction.) The matched-dose value being indistinguishable from zero is exactly
what the design calculation predicts from Tahoe's 4,560 pairs.

This is the failure mode the calculation is meant to prevent. The atlas was built
at a scale that made the question look answerable, the replicate structure needed
to answer it was present but thin, and the convenient substitute for it inflates
the answer twofold.

### What this implies for models and benchmarks

Two consequences constrain models rather than experiments.

**A fixed context embedding falls behind as an atlas grows.** The interaction's
effective dimensionality grows at roughly 0.05 directions per context (LINCS,
slope 0.05, *r* = +0.55, *P* = 1×10⁻⁶⁴) and does not saturate, so a
*d*-dimensional context embedding represents a falling fraction as contexts
accumulate — *d* = 5 holds 100% of the interaction at 20 contexts and 47% at 200.
The usable form is a rule: **choose *d* ≥ 0.05 × n_contexts.** This is a property
of the model class rather than of training, so more data does not remove it.

**Benchmarks must report the context count behind their mean baseline.** The
perturbation-mean prediction improves monotonically with the number of contexts it
is estimated from: *r* = 0.550 at 2 contexts, 0.592 at 3, 0.625 at 6, 0.654 at 18,
0.666 at 45. A model equally good everywhere appears **19% stronger** against a
2-context baseline than against an 18-context one, and 21% stronger than against a
45-context one. State (Adduri et al., *Cell*
2026) evaluates across five query datasets whose mean baselines are built from 2,
2, 3, 5 and 17 contexts, and reports its largest gains on the smallest. This is a
**partial** confound — a ~19-point differential runs in the reported direction but
does not explain several-fold gaps — and the useful conclusion is narrower: gains
are comparable across benchmarks only at matched context counts.

Together these say the binding constraints are experimental rather than
architectural: context count, replication, and identity verification.

---

## Discussion

The field is scaling perturbation atlases along the one axis that does not buy the
measurement they are built to make. Cell count enters an interaction estimate only
through per-condition noise, which saturates; replicate pairs enter through a term
that does not. The consequence is concrete and unflattering: Tahoe-100M would have
obtained the same answer from 2% of its cells, and the budget released by that
would have bought the replication its own question needed.

We think the reason the calculation has not been done is that the interaction is
treated as an output of analysis rather than as an object with a detection limit.
Once it is treated as the latter, the question "can this atlas answer it" becomes
arithmetic available before the first cell is sequenced, and the answer for most
of the public corpus is no. That the corpus is dominated by single-context CRISPR
screens is not a criticism of those experiments, which were not built for this
question. It is a caution against computing an index of context-dependence on
them anyway, which is what the absence of a stated detection limit invites.

The most useful form of this work is prescriptive. For an atlas being designed
now, the calculation says how to spend a fixed budget — for Tahoe's own 95.6
million cells over 48 × 1,100 conditions, eight replicates of 226 cells detect
0.0016, eightfold better than what was achieved — and says which of contexts,
perturbations or replicates is binding. We release it as `perturbdesign` with
`plan`, `audit` and `budget` subcommands and a browser calculator, and we
pre-register it against atlases that do not yet exist, because a calculation of
this kind earns its keep prospectively or not at all.

### Limitations

**The mathematics is a standard power calculation.** What is new is its
application, its validation against real atlases with no free parameters, and the
field survey — not the statistics, which are a variance calculation for a
covariance estimator.

**Five atlases, chosen by us.** They were the five this project had analysed. The
38-dataset survey addresses the selection concern for the design parameters but
not for the measured shares, which exist only for the five.

**One informative prospective test.** Part A of the pre-registration is *n* = 1.
The standing registration (Part B) is a commitment, not a result, and cannot be
credited until qualifying atlases appear.

**The snr constant is not measured.** It is fixed at 0.20 and shown to be
insensitive over 0.10–0.30, but it is a stand-in for per-condition noise that a
pilot should estimate directly.

**The budget optimiser's noise curve has two assumed parameters.** Swept over a
fortyfold range the optimum is never one replicate and its direction is stable,
but the exact optimum moves.

**Replicate counts are upper bounds.** Deposited metadata cannot separate an
independent replicate from a split capture.

**No new experiments.** Every measurement here is a reanalysis of public data. The
prescriptions are untested prospectively by construction; that is what the
registration is for.

---

## Methods (summary)

**Decomposition.** `y = β(drug,dose) + α(context) + γ(context,drug) + noise`, each
term estimated as a covariance between independent replicate estimates so noise
contributes zero in expectation. β is estimated from other contexts, α from other
perturbations and **within a replicate plate**. Implementation:
`perturbmodel.celldrug`; 31 unit tests including recovery of a planted share and a
null.

**Design calculation.** Pair count `n_ctx × n_pert × n_rep(n_rep−1)/2`; standard
error `(1 + 1/snr²)/√(pairs × n_feat)` with `n_feat = 2000`; minimum detectable
share `1.96 ×` that. Shared `snr = 0.20` throughout, swept in
`scripts/design_calculator.py`. Implementation: `perturbmodel.design`.

**Field survey.** Design parameters read from `obs` with
`perturbmodel.atlas_meta`, no expression values loaded. Missing markers (`''`,
`None`, `nan`) excluded from level counts. Replicates counted **per condition**,
not per dataset. All 38 scPerturb RNA/protein datasets scored;
`scripts/atlas_design_benchmark.py`.

**Premise test.** Counts binomially downsampled per condition with replicates and
contexts held fixed, and contexts thinned with cells held fixed; readout is the
MEK×MAPK effect of `RESULTS.md` §36. `scripts/cells_vs_replicates.py`.

**Pre-registration.** `docs/PREREGISTRATION.md`, frozen at tag
`prereg-2026-09-07`; frozen artefacts verified by SHA-256 in the test suite;
outcomes appended to `docs/PREREGISTRATION_OUTCOMES.md`.

**Tool.** `perturbdesign` (`perturbmodel.design_cli`); browser calculator
`docs/design_calculator.html`.

---

## Figures

1. **A design calculation for perturbation atlases.** (a) For five published
   atlases, the smallest interaction share each design can detect against what we
   measured in it; the tick or cross is the prediction, made from each atlas's
   context, perturbation and replicate counts alone and correct in all five.
   (b) The sensitivity sweep: 5/5 holds for any shared noise constant in
   0.10–0.30, so the result is not a fitted parameter. (c) The fraction of the
   interaction a fixed *d*-dimensional context embedding can represent as an atlas
   grows; the usable form is *d* ≥ 0.05 × n_contexts.
   *Bundle:* `results/figures/00_manuscript/fig10/`.
2. **Cells or replicates: the premise, tested.** (a) The §36 MEK×MAPK effect under
   binomial downsampling of cells with every condition kept, against the same
   effect under thinning of contexts with cells kept. (b) The two axes relative to
   the full data: 2% of cells retains 92%, while 10% of contexts collapses
   precision to a spread of 0.084 against an effect of 0.10.
   *Bundle:* `results/figures/00_manuscript/cells_vs_replicates/`.
3. **The same calculation across the whole field.** All 38 scPerturb RNA and
   protein datasets, scored from deposited metadata with one shared constant and
   no per-dataset tuning. (a) Contexts against replicates per condition, coloured
   by perturbation type, with the five atlases of Figure 1 as reference points;
   the estimator is defined only in the shaded corner. (b) The funnel: 38
   datasets, 4 multi-context, 1 that also replicates a condition. (c) Drug screens
   alone, since a CRISPR screen in one line is single-context by design.
   *Bundle:* `results/figures/00_manuscript/fig11/`.
4. **The quantity being designed for.** (a) Variance apportionment into drug
   effect, cell property and cell–drug relation, each a covariance between
   independent estimates (94.1% / 2.0% / 3.9%; PRISM, 737 lines × 1,324
   compounds; log axis). (b) A line's general sensitivity on one half of the
   compound library against a disjoint half, *r* = 0.989. (c) With the truth held
   at 0.20 and only the planted general response growing, the corrected estimator
   returns 0.211/0.212/0.211 where the uncorrected returns 0.170/0.347/0.595.
   (d) The four failure modes on data with a known answer.
   *Bundle:* `results/figures/00_manuscript/fig1/` (panels a–d).
5. **Tahoe's replicate structure, and the cost of the substitute.** (a) 13.5% of
   (line, drug, dose) triples sit on more than one plate, falling to 5.4% once
   plate 14 — a designed replicate of plate 6 — is dropped under the training
   convention. (b) True-replicate against cross-dose pairing versus a matched
   null: the substitute doubles the reported interaction, 11.5% → 20.7%.
   *Bundle:* `results/figures/00_manuscript/fig1/` (panels e–f).
