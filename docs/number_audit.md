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
