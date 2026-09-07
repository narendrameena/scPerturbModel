# Number audit of `MANUSCRIPT_DESIGN.md`

**2026-09-07.** Every quantitative claim in the design paper checked against its
source — a result table, a rerun of the generating code, or `RESULTS.md`. Prompted
by finding that the five-atlas table's detection floors had been stale for some
time while the verdicts they supported stayed correct, which is exactly the
failure mode a spot-check misses.

**49 sentences carrying numbers.** Method: extract every numeric claim
programmatically, group by source of truth, verify each group against that source
rather than against the surrounding prose.

## Result

| | count |
|---|---:|
| verified against source | 45 |
| **wrong, corrected** | **3** |
| incomplete rather than wrong, expanded | 1 |

No claim was found to be unverifiable.

## Errors found and corrected

### 1. Baseline *r* values misquoted (two of five)

The context-count sweep was quoted as *r* = 0.557 at 2 contexts and 0.660 at 18.
`RESULTS.md`'s table gives **0.550** and **0.654**. The other three (0.592, 0.625,
0.666) were right.

### 2. "21% stronger" applied to the wrong comparison

The paper said a model appears "21% stronger against a 2-context baseline than
against an 18-context one". 21.1% is the handicap of a 2-context baseline against
a **45**-context one; against an 18-context baseline it is **18.9%**. The error was
inherited from `RESULTS.md`'s prose, which states ~21% for the 2-vs-18 comparison
while its own table gives 21.1% for 2-vs-45.

It was also internally inconsistent: the very next sentence described the same
quantity as "a ~19-point differential". Now 19% for 2-vs-18, with the 21% figure
kept and correctly attributed to 2-vs-45.

### 3. Corrected and uncorrected estimators mixed in one paragraph

The cross-dose paragraph said pairing across doses "doubles the reported
interaction, from 11.5% to 20.7%", then two sentences later gave Tahoe's
matched-dose interaction as 0.5%. Both are true, but 11.5% and 20.7% are
**pre-correction** values and 0.5% is **post-correction**; a reader is left unable
to tell whether Tahoe's interaction is 11.5% or 0.5%.

Under the corrected estimator the substitution runs **0.5% → 9.2%**, which is a
stronger result than the doubling: it moves the estimate from not distinguishable
from zero to clearly non-zero. Rewritten to lead with the corrected pair and give
the uncorrected one parenthetically.

## Expanded rather than corrected

**Calibration was reported selectively.** The paper gave power above the threshold
(83–100%) and the false-positive rate (0–8%, 17% on the smallest design), but
omitted that power *below* the threshold reaches **75%** on one design — the floor
is conservative rather than sharp. That is in the generating script's own output
and is unflattering, so it now appears in the text.

## Verified, with source

| claim | value | source |
|---|---|---|
| five-atlas pairs / floors / verdicts | all 15 numbers | `perturbmodel.design`, pinned by `test_five_atlas_table_matches_the_code` |
| snr sweep 5/5 over 0.10–0.30, 4/5 at 0.35 | — | `scripts/design_calculator.py` |
| Tahoe 6 replicates, Spear-ATAC 11 | — | `required_replicates` |
| 2% of cells → 92% of effect | 92% | `results/tables/cells_vs_replicates.csv` |
| 10% of contexts → spread 0.084 vs effect 0.10 | 0.084, 0.1007 | same |
| calibration power / FPR | 83–100%, 0–17% | `design_calculator.py` |
| 38 datasets → 4 multi-context → 1 | 38/4/1/1 | `results/tables/atlas_design_benchmark.csv` |
| 27 CRISPR none multi-context; 10 drug, 3 multi, 1 replicated | — | same |
| sci-Plex 3 floor | 0.048 | same |
| 34/38 (89%) per file, 17/22 (77%) per study | — | `atlas_design_benchmark.py` |
| MimitouSmibert2021: 10,018 cells, 1 context, 6 perts, 1 rep | — | `PREREGISTRATION_OUTCOMES.md` |
| general sensitivity *r* = 0.989, 737 lines | — | `RESULTS.md` |
| 570/736 lines (77.4%) FDR < 0.05 | — | `RESULTS.md:1590` |
| PRISM 94.1% / 2.0% / 3.9%; variances 1.6338 / 0.0343 / 0.0683 | — | `RESULTS.md:1613–1615` |
| corrected estimator slope 0.950, *R*² 0.9995 | recomputed | `estimator_simulation/simulation.csv`, all-seeds fit |
| returns 0.008 at the null | 0.0079 | `recovery.csv` |
| residual variance 72%, pooled-batch 41% at true 0.2 | 0.717, 0.411 | `recovery.csv` |
| planted-response sweep 0.211/0.212/0.211 vs 0.170/0.347/0.595 | — | `RESULTS.md` |
| Tahoe 13.5% of triples replicated | 13.446% recomputed | `RESULTS.md:762`, pseudobulk metadata |
| falls to 5.4% without plate 14 | 5.35% | `RESULTS.md:762` |
| matched-dose 0.5% [0.0–1.5%], *P* = 0.10 | — | `RESULTS.md:1854` |
| dimensionality slope 0.05/context, *r* = +0.55, *P* = 10⁻⁶⁴ | LINCS | `RESULTS.md:2808` |
| *d* = 5 holds 100% at 20 contexts, 47% at 200 | recomputed | `captured_fraction` |
| budget optimum 8 reps × 226 cells → 0.0016 | — | `perturbdesign budget` |
| 95.6M cells, 48 × 1,100 conditions, 13.5% replicated | — | Tahoe-100M dataset card |

## What this changes about the paper's status

None of the three errors touches a headline. The 5/5 validation, the 1-of-38
survey and the premise test are unaffected. Two of the three were in the
model-implications section, which is the paper's least load-bearing part, and the
third was a presentational mixing of two estimators rather than a wrong
measurement.

The pattern worth noting is that **all three were inherited from the combined
draft rather than introduced in the split**, and two trace to `RESULTS.md` prose
disagreeing with its own tables. The tables are right in every case checked. Prose
restating a table is where the errors live, which argues for generating such
sentences from the tables where possible, and for pinning published tables with
tests as the five-atlas table now is.

---

# Number audit of `MANUSCRIPT_MEASUREMENT.md`

**2026-09-07**, same method: extract every numeric claim programmatically, group
by source of truth, verify against that source rather than against surrounding
prose. **52 sentences carrying numbers.**

| | count |
|---|---:|
| verified against source | 46 |
| **wrong, corrected** | **4** |
| **unverifiable, flagged in the text** | **1** |
| internal inconsistency, resolved | 1 (counted in the 4) |

This paper carries a higher error rate than the design paper (4 wrong of 52
versus 3 of 49) and, unlike that one, an error of a more serious kind.

## Errors found and corrected

### 1. A withdrawn claim asserted as fact — the most serious finding

The paper stated: *"Copy number beats mutations (+0.0064, p = 5.7×10⁻⁴), as
Schlüter & Schönhuth report."* That advantage was **withdrawn** by this project's
own audit (`RESULTS.md:1549`): the 150 compounds were treated as independent when
their residuals correlate at mean *r* = +0.23, an effective *n* ≈ 11, and a
compound-cluster bootstrap returns **CI [−0.002, +0.015]**, which includes zero.

The combined draft's status banner recorded the withdrawal correctly. The Results
prose was never updated to match, and the split carried the stale sentence
forward. Now states the withdrawal; what survives — copy number trailing
expression by 0.088 and adding nothing to it — is kept.

### 2. Mechanism-ranking ρ against LINCS phase 1

Stated as **+0.19**; the corrected table (`RESULTS.md:872`) gives **+0.09**
(*p* = 0.80, 10 shared classes). An earlier superseded table in the same file
gives +0.56, produced by a biased per-perturbation estimator and explicitly
retired. Neither value supports a different conclusion — both are null — but the
quoted number was wrong.

### 3. Genotype power sweep wrong on all four points

Stated as 4% / 45% / 72% / 96% at 47 / 250 / 400 / 600 cell lines.
`results/tables/prism_genotype_power.csv` gives **5% / 40% / 78% / 100%**. The
78% at 400 also reconciles with `RESULTS.md` §11's "~400 lines for 80% power",
which the quoted 72% did not. The 600-line row rests on 9 associations rather
than 12; that is now stated.

### 4. The matching ladder given twice, with different numbers

Two sentences in the same section reported the same quantities differently:

| | identifier | best available |
|---|---:|---:|
| first statement | 0.241 | 0.427 |
| second statement | **0.235** | **0.419** |
| `RESULTS.md:1044–1048` | **0.235** | **0.419** |

The second was right; the first is corrected.

## Unverifiable, flagged rather than silently changed

**The predictor-block table** (baseline expression +0.0998 / 99.2%, and six
further rows) is from a **150-compound** run whose output table is not in
`results/tables/`. Not one of its seven *R*² values appears anywhere in
`RESULTS.md`. The traceable analysis there is a **120-compound** version
restricted to lines on which every block is measurable, giving baseline
expression **+0.0927 / 92.5%**.

These are different analyses, so substituting one for the other would be wrong. A
warning block now sits above the table in the paper. Note also that 99.2%
elsewhere in this project is a *different statistic* — 119/120 compounds beating
their own permutation null — so its appearance in a "compounds positive" column
may be a transcription of the wrong quantity.

## Two generations of corrections, and which won

`RESULTS.md` contains an audit table (line ~1545) whose corrections were
themselves later superseded by a rerun. Where they conflict, the rerun is
authoritative and the paper follows it correctly:

| quantity | audit table (superseded) | rerun (current) | paper |
|---|---|---|---|
| assay vs laboratory | 68% / 32% | **81% / 19%** | 81% / 19% ✓ |
| identity-validated transfer | 71% | **98.0%** | 98.0% ✓ |
| reliability-matched control | 59% | **49%** | 49% ✓ |

All three verify against `results/tables/cross_lab_summary.csv`
(0.8089 / 0.1911, 0.9804, 0.4934). The paper also correctly states that the
laboratory share's CI [−6%, 38%] includes zero.

## Verified, with source

Dose (1,443 compounds, ρ = +0.33, 74% positive, *P* = 4×10⁻⁹², +74% vs +36%);
mechanism ρ = −0.18 vs PRISM and the 11–12 class-count correction; 111,589
allele × compound tests; 11 of 11 GDSC biomarkers, BRAF V600E *P* = 3×10⁻³¹;
3,435 matched genes for the synonymous control; 488 COSMIC→DepMap lines;
identifier vs same-tissue *P* = 1.5×10⁻¹²⁴; the full cross-lab ladder
(0.685 / 0.365 / 0.290); transfer vs signal strength ρ = +0.54, *P* = 1.8×10⁻²⁹;
transcription 46% (0.061, 0.032), 16 of 16 lines; sci-Plex 74.4% vs published
48%; correction factors 1.01× and 1.15×; TRADE 2,052 perturbations, 6,642 genes,
61.5% → 58.1%, noise 38%, *r* = +0.51 uncentred and −0.53 centred, 39% vs
published 56%; shift share 8.5% / 10.1%; remodelling 5,345 conditions, 36 lines,
12 held out, off-axis 89/86/90%, spread ratios 1.002/1.000/1.007, co-expression
ρ = +0.998, mean absolute change 0.0127; Spear-ATAC 41 perturbations × 3 lines,
73,344 cells, 6 of 72 at Bonferroni, GATA1 Δ = −4.19, 2,174 motifs, 24,919 genes,
interaction 33.9% against null 33.0% [24.2–39.5%], *P* = 0.45.

## Conclusion across both papers

101 numeric claims audited, **7 wrong, 1 unverifiable**. Every error was
inherited from the combined draft; none was introduced by the split. The
distribution matters: the design paper's three errors were all presentational or
peripheral, while this paper's include a retracted result stated as fact.

The mechanism is consistent across both. **Result tables were right in every case
checked; prose restating them was where errors lived**, and the worst case arose
when a correction was recorded in a status banner but never propagated into the
Results text. Two practical consequences: generate such sentences from tables
where possible, and treat any claim carrying a withdrawal note in one place as
requiring a search for every other place it is stated.
