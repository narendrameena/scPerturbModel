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

---

# Withdrawal-propagation sweep of `RESULTS.md`

**2026-09-07.** Both manuscript audits found the same mechanism: a correction
recorded in one place, never propagated to the other places the claim appears.
`RESULTS.md` (3,313 lines) is the source both papers draw from and the one file
neither audit checked against anything. This sweep takes every withdrawal or
supersession marker in it and asks whether the withdrawn claim is still asserted
elsewhere in the same file.

**13 withdrawal markers. 9 of the withdrawn claims are still asserted somewhere
in the file, 4 of them in text carrying no supersession banner at all.**

## The serious findings

### 1. A reversed conclusion, stated twice as a considered nuance

§14 says: *"Changing the assay within one laboratory costs 16% of the total loss;
changing laboratory costs the remaining 84%."* The Limitations then reasons from
it: *"the assay contribution is not distinguishable from zero… protocol is not the
problem."*

Both are **inverted**. On a matched compound set it is **assay 81% / laboratory
19%**, and it is the *laboratory* share whose CI [−6%, 38%] includes zero
(`cross_lab_summary.csv`: 0.8089 / 0.1911). Protocol is the larger term, not the
smaller one — the opposite of what the Limitations concludes.

§31 (line ~1896) explicitly records that this "was still inverted in three places
— the Discussion, the Limitations and a Results paragraph — despite having been
corrected in the status banner." **It was still inverted in this file when this
sweep ran.** The fix had been applied to the manuscript and not to its source.

### 2. The copy-number claim asserted as fact, as in the measurement paper

§17 withdraws it; the §17-adjacent text at line ~948 still read *"copy number does
beat mutations (+0.0064, p = 5.7×10⁻⁴) — their claim holds"*. Same error, same
cause, and the direct origin of the measurement paper's worst finding.

### 3. "Not an underpowered comparison" — directly contradicted

§12 argues a null is meaningful because it rests on "67 shared classes, so this is
*not* an underpowered comparison". §31 establishes that only **11–12** classes
clear the minimum class size, that the class count was wrong six-fold, and that
the power claim was unsupportable. §12's sentence stood unmarked.

### 4. The PRISM line count

Corrected to **737** at line ~1541 (one line was a parsing artefact; ten carry
`FAILED_STR`). "738 cell lines" was still asserted in four separate places,
including a dataset-description line.

## One quantity is irreproducible, which the sweep found by accident

The mechanism ranking of Tahoe against LINCS phase 1 is quoted as **+0.09** in
§16's table and **+0.19** in §31's prose. The only saved table,
`results/tables/three_platform_mechanism_cdi.csv`, recomputes to **+0.558**
(n = 10, *p* = 0.093) — the *pre-correction* value §16 explicitly retires.

**The post-correction run's output was never committed, so neither +0.09 nor +0.19
can be reproduced from this repository.** All three candidates are null at n = 10,
so no conclusion turns on it. But this also means the measurement paper's
correction earlier today — +0.19 → +0.09, made on the strength of §16's table —
substituted one unverifiable number for another. That paper now states the null
without a point estimate and records why.

## What was done

- A **superseded-value index** at the top of `RESULTS.md` listing every known
  stale quantity, where it still appears, its current value, and its source of
  truth.
- **Inline markers** at the four places stating a reversed or retracted
  conclusion. Stale sentences are struck through and kept, not deleted: the file
  is a chronological record and destroying it would lose the provenance both
  papers depend on.
- The four "738" assertions corrected to 737.
- The measurement paper's ρ restated as a null without an irreproducible point
  estimate.

## The pattern, across all three audits

| audit | claims | wrong | unverifiable |
|---|---:|---:|---:|
| design paper | 49 | 3 | 0 |
| measurement paper | 52 | 4 | 1 |
| `RESULTS.md` sweep | 13 markers | 9 still asserted | 1 |

Every manuscript error traced to a superseded value copied out of a section that
had not been re-run. The corrections were not missing — they were present, in the
same file, sometimes on a line explicitly complaining that the fix had not
propagated. What was missing was any mechanism forcing a correction to reach every
site of the claim.

Two things would prevent recurrence, and neither is expensive: keep prose that
restates a table generated from that table, and treat "committed the corrected
table" as part of the definition of a correction — the irreproducible ρ exists
because a rerun's output was never saved.

---

# Re-run: the mechanism comparison, 2026-09-08

The sweep found that the Tahoe-vs-LINCS mechanism correlation was quoted as +0.09
(§16) and +0.19 (§31) while the only saved table reproduced +0.558 — the
pre-correction value §16 retires — because the corrected run's output had never
been committed. It has now been re-run and its table committed to
`docs/source_data/`.

**Cause.** `results/tables/three_platform_mechanism_cdi.csv` was generated
2026-08-31. Three of its four inputs were regenerated on 2026-09-01 and
2026-09-03 when the corrected estimator was applied. The table was never
re-generated, so it silently described a superseded state. `results/` is
gitignored, so nothing recorded that it had gone stale.

**Result.** Two of §16's three rows are confirmed; the third changes sign.

| comparison | classes | §16 (2026-09-03 prose) | re-run 2026-09-08 |
|---|---:|---:|---:|
| Tahoe vs LINCS-1 | 10 | +0.09 | **+0.103** (*p* = 0.78) |
| Tahoe vs PRISM | 11 | −0.18 | **−0.155** (*p* = 0.65) |
| LINCS-1 vs PRISM | 67 | −0.09 (n.s.) | **+0.279** (*p* = 0.023) |

§16 was right to one decimal on the Tahoe rows; **§31's +0.19 was wrong** and the
measurement paper's point estimate is restored to +0.103.

## The readout-decoupling claim is withdrawn

The third row is not a rounding difference. LINCS-1 against PRISM is the only
adequately powered comparison in this analysis — 67 shared mechanism classes
against 10–11 for either Tahoe comparison — and it moves from a null to
**positive agreement between a transcriptional and a viability readout**. On all
annotated compounds rather than active ones it is +0.353 over 78 classes
(*p* = 0.002). A consensus transcriptional rank against viability moves from
−0.19 to +0.190.

This is the second time this claim has been wrong, each time for a different
reason. It was first stated as "ρ = −0.09 over 67 shared classes, a well-powered
null"; §31 showed the *Tahoe* comparisons rest on 11–12 classes and restated it as
an absence of evidence. That restatement kept the −0.09, which the re-run now
shows came from the stale table.

**Held back from overstatement:** six pairwise comparisons were made, Bonferroni
at α = 0.05 requires *p* < 0.0083, and +0.279 does not meet it. The defensible
statement is that the evidence for decoupling is withdrawn and what remains is
weak positive agreement — not that agreement is established. The measurement
paper's section is retitled accordingly and now leads with the powered comparison
rather than the underpowered ones.

## The process fix

`docs/source_data/` is a new tracked directory for small tables backing numbers
quoted in either manuscript, with provenance recorded in its README: generating
script, date, each input and its date. `results/` stays gitignored for figures and
large intermediates.

This closes the specific gap the sweep identified — "committed the corrected
table" is now something that can be checked, rather than assumed. The rule it
implements: **a number may appear in a manuscript only if the table that produced
it is committed.**

---

# Full staleness sweep and bug hunt, 2026-09-08

Prompted by the mechanism-table finding: if one table had gone stale relative to
its inputs, others might have too.

## Staleness

Two detectors were run over all 144 tables. **6 stale outputs found; all 6
re-run.**

| detector | found |
|---|---|
| output older than its newest input table | 5 |
| output of a corrected-estimator script predating `celldrug.py` (2026-09-03 11:18) | 1 |
| **remaining after re-run** | **0** |

A third detector — figure bundles older than their own source-data CSV — flagged
98 bundles and was a **false positive**: `save_figure` writes the PNG and the CSV
in one call, and the largest gap across 99 bundles was 37 s. No real figure
staleness.

### What changed on re-run

| table | change |
|---|---|
| `cdi_vs_target_genetics.csv` | **CDI changed for all 254 drugs** (max Δ 0.236); target-annotation columns unchanged |
| `prism_vs_tahoe_cdi.csv` | PRISM CDI max Δ 0.264, Tahoe max Δ 0.036 |
| `sparse_validation*.csv` | `n_flagged` ±1; hypergeometric *p* shifted slightly as a consequence (deterministic, not RNG) |
| `cross_lab_identity.csv` | **not rewritten — see the bug below** |

### A claim changed

`RESULTS.md` said CDI is uncorrelated with target expression level
(ρ = +0.08, *p* = 0.21). Re-run: **ρ = +0.135, *p* = 0.032** — nominally
significant, though not surviving Bonferroni across the three tests
(*p* < 0.017 required). The mutation-frequency and expression-variance
correlations also moved (+0.10 → +0.113; −0.21 → −0.228). Corrected in place.

Separately, `RESULTS.md` quoted "264 drugs" in a summary table while its own prose
and its table both had **254**. Pre-existing, unrelated to the re-run, corrected.

## A confirmed bug: the transcription cross-laboratory arm

`cross_lab_identity.csv` is dated 2026-09-01 — **three script revisions** before
`cross_lab_transcription.py` was updated (`fb4830f`, 2026-09-03) to strip each
line's general response. Re-running the current script:

| quantity | published | re-run |
|---|---|---|
| LINCS p1 vs p2, within-lab median *r* | 0.061 | **0.023** |
| Tahoe vs LINCS, cross-lab pairs | (r = 0.032) | **0 pairs** |
| reproducible fraction | 46% | **not computable** |
| identity check | 16 of 16 | **skipped** (no shared compounds) |

**The zero is a bug, not a null.** 172 (line, compound) pairs demonstrably exist
in the raw inputs — ht29 135, a549 25, hs578t 12, from 136 shared normalised
compound names and 3 shared cell lines — and the pipeline loses all of them
between loading and pairing. `remove_line_effect_profiles` is *not* the cause: its
min-3-compounds filter cannot drop ht29. **Cause not located; not guessed at.**

Both `RESULTS.md` and the measurement paper now carry a do-not-quote banner on
this arm. **The viability arm is unaffected** — it has a committed table
(`cross_lab_summary.csv`) and reproduces exactly.

## Bug hunt

Four failure classes were scanned, chosen because this project has hit three of
them before.

| class | result |
|---|---|
| DataFrame attribute shadowed by a column name (`.shift`, `.ndim` — hit 3× historically) | **0 live instances.** 16 tables do carry colliding column names (`index`, `size`, `median`, `diff`, `rank`, `mode`, `filter`, `sample`, `ndim`), and all 8 dot-access sites resolve to the real attribute or use bracket indexing. |
| unseeded randomness → irreproducible *p*-values | **0.** Four scripts flagged by a first pass were all seeded; the detector required a literal digit and missed `default_rng(seed)` with a loop variable, plus two hits that were in comments. |
| circular table dependencies | **0 real.** Three scripts read a table they also write, all behind a `--replot` flag that redraws a figure without recomputing. |
| stale outputs | 6, all re-run, 0 remaining |

## The process gap this exposes

The mechanism table went stale because `results/` is gitignored and nothing
recorded that its inputs had moved. `docs/source_data/` fixes that for tables
backing manuscript numbers, and now holds three. But the transcription bug shows
the deeper issue: **a script can be corrected without its outputs being
regenerated, and nothing fails.** A `make`-style dependency check — outputs older
than inputs or older than the code that writes them — would have caught all six of
today's cases and the mechanism table in September. It is the obvious next piece
of infrastructure and does not exist yet.

---

# Cross-laboratory transcription arm: fixed and re-run, 2026-09-08

The keying bug is fixed (`383e890`) and the arm re-run. **Every number in it
changed, and the headline nearly halved.**

| quantity | published (buggy) | corrected |
|---|---:|---:|
| reproducible fraction, transcription | 46% | **27.3%** |
| within-lab LINCS p1 vs p2, median *r* | 0.061 | **0.034** |
| within-lab matched ceiling | — | **0.044** |
| cross-lab Tahoe vs LINCS, median *r* | 0.032 | **0.012** |
| cross-lab pairs | 489 | **317** (p1), 172 (p2) |
| identity, within Broad | 16 of 16 | **15 of 16** (94%) |
| identity, across labs | 3 of 6 | **2 of 6** (33%) |

**The fix is confirmed by an independent prediction.** Before fixing, I computed
from the raw inputs that exactly **172** (line, compound) pairs should exist
between Tahoe and LINCS phase 2. The re-run reports **172**. The pre-fix run
reported 0.

**What survives.** The direction of the claim: transcription transfers less well
than viability, and cross-laboratory transfer is worse than within-laboratory.
**What does not:** the magnitude. Transcription now transfers at *half* the
viability arm's rate rather than four-fifths of it, and the within-lab
transcriptional ceiling (*r* ≈ 0.034–0.044) is low enough that 27.3% is a ratio
of two very small numbers.

Corrected in `RESULTS.md`, `MANUSCRIPT_MEASUREMENT.md` and the archived draft,
including three secondary sites that restated the stale values — the propagation
check that the September sweep showed is necessary. Generating tables are now
committed to `docs/source_data/`.

## Tally for the day

| | |
|---|---|
| stale tables found | 6 |
| stale tables re-run | 6 |
| remaining stale | **0** |
| bugs found | **1** (loop-variable shadowing, silent, ~5 days live) |
| bugs fixed | 1, with an AST regression test |
| claims materially changed | **3** — readout decoupling withdrawn; CDI vs target expression level now nominally significant; transcription transfer 46% → 27.3% |
| scooping risks resolved | 2 of 2, one a material partial overlap |

---

# Rerun of the 49 stale-output scripts, 2026-09-08 (in progress)

Snapshot of 142 tables taken first, so every regenerated table is diffed against
what it replaced. Deltas below 1e-4 are classed as float jitter, not results.

## Status

14 of 49 scripts complete, 0 failed. **23 tables re-run: 19 identical, 2 float
noise, 2 genuinely changed.** Stale outputs 124 → 115.

## Changed

**`variance_component_comparison.csv` — §21 reversed.** Reported above; the
mixed-model agreement claim is withdrawn and the tool's failure mode inverts.

**`sparse_benchmark.csv` — one figure.** Permutation jitter only. Every aggregate
claim holds (Higher Criticism rejects 100% at *k* = 4, 50% at *k* = 8, 25% at
*k* = 16, 0% at the null; the pooled test flat at 12% throughout). Precision at
*k* = 4 is **92%**, not the published 98%. Corrected; not quoted in either
manuscript.

## Confirmed unchanged

Nineteen tables reproduce exactly, including several carrying real claims:
`estimator_simulation.csv` (§19's slope 0.950, *R*² 0.9995, 0.008 at the null),
`spectrum_benchmark.csv` (the spectral estimator's rejection),
`two_stage_benchmark.csv` and `two_stage_discovery.csv` (the two-stage method's
rejection), `cell_cycle_log2or.csv` (the paper replication),
`replication_correlations.csv` (plate 6 vs 14), and all dataset-fact tables.

Two tables differ only in the sixth decimal of an energy distance
(`published_methods_etest.csv`, `atac_gene_scores_etest.csv`); *p* and *q* are
identical.

## The predictor-block table, traced

The 2026-09-07 audit flagged the measurement paper's predictor-block table as
unverifiable — its source was not in `results/tables/`. **It is
`expression_architecture.csv`**, written by `expression_gap_closure.py`, and
finding it turns one vague flag into three specific defects:

| | paper | file |
|---|---|---|
| median CV *R*², all five blocks | +0.0998, +0.0763, +0.0251, +0.0075, +0.0012 | **identical** |
| compounds | 150 | **120** |
| expression, % positive | 99.2% | **100.0%** (99.2% is `r2_cnv_expr`, median +0.0801) |
| protein, % positive | 88.3% | **98.3%** |
| lineage, % positive | 86.7% | **85.0%** |
| copy number, % positive | 57.5% | **63.3%** |
| nonsynonymous, % positive | 18.3% | **52.5%** |

So the *R*² column is right and the "compounds positive" column beside it comes
from a different analysis — the two were merged into one table. The last two rows
(mutational burden, synonymous) come from a third source,
`genetic_architecture.py`, whose current values are −0.0031 / 29.7% and
−0.0077 / 24.3% against the paper's −0.0023 / 32.0% and −0.0064 / 22.0%.

The flag in the manuscript now records this precisely. The table will be rebuilt
on one source per column once `expression_gap_closure.py` finishes re-running.

## Note on what this exercise is worth

Two conclusions have now changed from re-running stale tables (§21's mixed-model
agreement, §16's readout decoupling) and one manuscript table has been shown to
merge three sources. The remaining 35 scripts are the heavy ones — Tahoe and
LINCS at full scale — and are where the estimator correction is most likely to
bite, since those are the analyses that use it directly.

## Correction to this file, 2026-09-08

An entry added earlier today said the predictor-block table had been "verified
against a fresh re-run, which reproduces the file exactly". **That was wrong.**
The re-run of `expression_gap_closure.py` was killed by a one-hour timeout at 101
of 120 compounds (`rc=124`) and never rewrote the table; I misread a polling loop
and compared the snapshot against itself. The claim is retracted in the
manuscript and here.

What stands: the table's three defects — n = 120 not 150, a "compounds positive"
column taken from a different analysis, and two rows from a third script — are
real, and were established by reading the committed file, which needs no re-run
to check. What does not stand is any claim that the file has been independently
regenerated. It dates from 2026-09-03, after that morning's estimator correction,
so it is not known to be stale; it is simply unverified. The script has been
restarted without a timeout cap.

The episode is the same failure this whole audit is about — asserting a
verification that did not happen — and it is recorded rather than quietly fixed.

## Cascading staleness, working as intended

Re-running `prism_context_genetics.py` rewrote `prism_decomposition.csv`, which
the freshness check immediately flagged as making two downstream tables stale:
`prism_vs_tahoe_cdi.csv` and `three_platform_mechanism_cdi.csv`. Both were
regenerated.

**The mechanism comparison reproduces exactly on the fresh PRISM inputs**:
Tahoe vs LINCS-1 +0.103 (*p* = 0.777), Tahoe vs PRISM −0.155 (*p* = 0.650),
LINCS-1 vs PRISM **+0.279** (*p* = 0.022). The withdrawal of the
readout-decoupling claim therefore survives an independent regeneration of its
main input, which is the strongest form of confirmation available here.

This is also the first time the dependency check has caught a *cascade* rather
than a single stale file — the failure mode it was built for.

## A filename collision, found 2026-09-08

The freshness audit kept reporting `three_platform_mechanism_cdi.csv` as behaving
oddly — content changing for reasons unrelated to its own inputs. The cause is not
staleness at all.

**Two scripts wrote that filename with different schemas.**
`three_platform_synthesis.py` writes the per-mechanism CDI matrix (179 × 4: one
column per platform). `prism_vs_tahoe.py` wrote a 12 × 6 overlap table
(`m, prism, n_prism, tahoe, n_tahoe, lincs`) to the same path. Whichever ran last
won.

The consequence was silent. `manuscript_figures.py` reads that name and guards on
`{"Tahoe", "PRISM"} <= set(tp.columns)`; under the wrong writer the guard simply
failed, so **a panel disappeared from Figure 3 with no error and no warning**. It
also explains why this table looked mysteriously stale during the audit — the
mechanism comparison's apparent instability was partly this, not just the 08-31
timestamp.

Fixed by renaming `prism_vs_tahoe.py`'s output to `prism_tahoe_lincs_overlap.csv`.
Verified: running `prism_vs_tahoe.py` after `three_platform_synthesis.py` now
leaves the mechanism table intact (ρ = +0.279, *p* = 0.022, n = 67 unchanged).

A test now asserts no two scripts write the same table name, validated both ways —
it names this exact collision on the pre-fix source and passes on the fixed tree.
Writing the detector required care: an f-string like `f"eval{suf}.csv"` parses
into the constants `"eval"` and `".csv"`, and treating the bare extension as a
filename makes twelve unrelated scripts appear to collide.

## A silently empty analysis, found 2026-09-08

`compare_context_metrics.py` crashed on re-run with
`KeyError: np.int64(25625)`. The cause: a condition whose cell line is the only
one in its group has no leave-one-context-out mean, so the residual loop skips it
via `continue` — but the pairing loop below still looks it up in `resid`.

The crash is not the interesting part. **The committed
`context_metric_mechanism_ranks.csv` was empty — 0 rows.** This script exists to
provide the cross-check that `docs/methodology_rationale.md` §5 cites in defence
of using CDI over published-style metrics, and it had been producing nothing at
all. Nothing flagged it: an empty table is not a crash, and no test asserted the
file was non-empty.

With the guard added it yields **4 mechanism classes** — enough to run, far too
few to compare rankings. The correlations at that *n* are ρ = 0.000 (*p* = 1.00)
and ρ = −1.000, neither meaningful. All three metrics do put the same four classes
at the context-specific end (JAK/STAT, adrenoceptor agonist, DNA synthesis/repair,
cyclooxygenase), which is weakly reassuring and nothing more.
`methodology_rationale.md` now describes the cross-check as **attempted and
uninformative** rather than as supporting the choice of metric.

## Cascades cleared, values confirmed

Re-running `drug_context_dependence.py` and `prism_context_genetics.py` made four
critical tables stale in turn. All regenerated, and every number held:

* `drug_context_dependence.csv` re-ran **byte-identical** (367 × 14, no column changed).
* The mechanism comparison reproduces exactly again — LINCS-1 vs PRISM
  **+0.279, *p* = 0.022, n = 67** — now for the second independent regeneration.
* `cdi_vs_target_genetics.csv` reproduces its corrected values (+0.113, −0.228, +0.135).

## Other changes this batch

`tcga_anchoring_survival.csv`: component C2 dropped one cancer type (31 → 30) and
its *z* moved 0.629 → 0.959, *p* 0.529 → 0.338. Non-significant either way, and
this table is not quoted — `RESULTS.md`'s survival numbers come from the
*adjusted* analysis, a different table. No claim moves.

## The predictor table, now genuinely verified

`expression_gap_closure.py` completed without the timeout cap and reproduces all
five blocks **exactly**: +0.0998 / 100.0%, +0.0763 / 98.3%, +0.0251 / 85.0%,
+0.0075 / 63.3%, +0.0012 / 52.5%, n = 120. The verification I claimed prematurely
earlier today has now actually happened, and the manuscript note says so rather
than saying it twice.

## A sign that had been wrong in two documents

`methodology_evidence.csv` re-ran with its D2 row moving from −0.00307 to
+0.00687. Chasing it: `docs/methodology_rationale.md` §15 and `docs/pertdecomp.md`
both argued against cross-dose pairing by quoting "true replicates covary at
**−0.00307**, cross-dose pairs at +0.00620" — replicates agreeing *less* than
different doses.

The negative was an artefact of a superseded prior, and
`methodology_evidence.py` says so in its own source comment: the in-sample and
leave-one-condition-out priors give covariances near −σ²/n, and "the negative
numbers were arithmetic, not absence of signal". Under the leave-one-context-out
prior the value is slightly positive.

Authoritative values, from `tahoe_true_replicates.csv` regenerated today:

| pairing | pairs | raw covariance | share |
|---|---:|---:|---:|
| true replicate (same line, drug, dose) | 11,492 | **+0.00025** | 0.46% |
| cross-dose (previous pairing) | 67,744 | **+0.00583** | 9.2% |

**The argument survives; the number did not.** The ordering that carries it —
replicates agreeing less than cross-dose pairs — holds by a factor of 23. Both
documents corrected, with the reason recorded rather than the value quietly
swapped. Note also that the corrected shares here (0.46% and 9.2%) are exactly
the figures already used in the manuscripts, so those were right while the
supporting covariances beside them were not.

## Confirmed unchanged this batch

* **Methylation null holds**: 0 of 20 features at FDR < 0.10 before and after
  (`epigenome_vs_context.csv`, min *q* 0.578 → 0.643).
* `target_abundance_pgx.csv` shifted on 256 rows but is not quoted in prose.
* `methodology_evidence.csv`'s other rows moved slightly (D3 additive prior
  0.4490 → 0.4422, D6 dose trend 0.3333 → 0.1557) without changing which option
  each decision selects.
