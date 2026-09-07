# Pre-registration: prospective tests of the design calculation

**Registered 2026-09-07.** Frozen at commit `2a27f3b66e64` (filled by the tag
`prereg-2026-09-07`, which is what timestamps this document; GitHub's record of
when that tag was pushed is the evidence that the predictions preceded the data).

RESULTS.md §46–§48 validated a design calculation *retrospectively*: on five
atlases whose outcomes this project had already measured, and on 38 scPerturb
datasets whose designs were read after the calculation existed. Retrospective
validation is worth what its blindness is worth, and a reviewer is entitled to
ask whether the calculation would have said the same thing before the answer was
visible.

This document fixes predictions before the corresponding data is examined. It has
two parts: one that settles within minutes of registration, and one that stays
open for any qualifying atlas released from now on.

Nothing below is edited after the fact. Outcomes are recorded in
`docs/PREREGISTRATION_OUTCOMES.md`, which appends and never rewrites.

---

## Part A — one genuinely unseen ATAC dataset (settles immediately)

**A first draft of this section claimed the whole scPerturb ATAC article was
unseen. That was false and is corrected here rather than quietly deleted.**
Checking before pushing showed that three of the article's six files —
`PierceGreenleaf2021_{GM12878,K562,MCF7}`, i.e. Spear-ATAC — were downloaded on
2026-09-04 and analysed at length in §35 and §39. Spear-ATAC is multi-context and
replicated, which is exactly the case the draft's predictions were counting. I
had already seen part of the answer, so those predictions were not blind and are
withdrawn in full.

What remains genuinely unseen, of article **24160968**:

| file | status at registration |
|---|---|
| `PierceGreenleaf2021_GM12878.zip` | **seen** — Spear-ATAC, §35/§39 |
| `PierceGreenleaf2021_K562.zip` | **seen** — Spear-ATAC, §35/§39 |
| `PierceGreenleaf2021_MCF7.zip` | **seen** — Spear-ATAC, §35/§39 |
| `Liscovitch-BrauerSanjana2021_K562_1.zip` | not downloaded, but the filename states the cell line |
| `Liscovitch-BrauerSanjana2021_K562_2.zip` | not downloaded, but the filename states the cell line |
| `MimitouSmibert2021.zip` | **not downloaded, design unknown** |

Only the last is a clean test, so Part A is a test with *n* = 1 and is reported as
such. A single dataset is weak evidence and will not be described as more.

### Predictions

Reasoning: §48 found 34/38 RNA/protein files single-context, driven by CRISPR
screens, which are single-line by construction. ATAC perturbation experiments are
overwhelmingly CRISPR screens and chromatin assays are lower-throughput per cell,
so I expect the same pattern. Mimitou et al. 2021 is, to my memory, an ASAP-seq
paper, but I do not know its released design and am not confident of the memory.

| # | Prediction | Falsified if |
|---|---|---|
| **A1** | `MimitouSmibert2021` does **not** support a context × perturbation interaction estimate — i.e. it has one context, or fewer than 2 replicates per condition, or both | it has ≥2 contexts *and* ≥2 replicates per condition |
| **A2** | Both `Liscovitch-BrauerSanjana2021` files are single-context | either has ≥2 contexts |

A2 is nearly free, since the filenames say `K562`; it is registered only so the
run is fully specified in advance, and it is **not counted as evidence** whichever
way it falls. A1 is the only prediction in Part A that carries information.

### What a failure would mean

If A1 fails, then a chromatin perturbation atlas exists whose design admits the
estimate and which §48's framing would not have anticipated, and §48's claim about
the field's designs must be softened. With *n* = 1 a success is correspondingly
weak support, and will be reported as one dataset, not as a confirmation.

### Analysis plan, fixed in advance

- Scored by `perturbmodel.atlas_meta.describe`, **unmodified**, at the commit this
  document is tagged with. The column-preference lists `CTX`, `PERT`, `REP` and
  the `MISSING` set are frozen as they stand.
- Detection floors from `perturbmodel.design.min_detectable_share` at the shared
  constant **snr = 0.20**, the value used throughout §46–§48. No per-dataset
  tuning, and no re-tuning of `snr` after seeing the ATAC designs.
- "Supports the estimate" means ≥2 contexts and ≥2 replicates per condition —
  the same rule as §48, not a new one.
- All three unseen files are scored. None is dropped for being inconvenient. A
  file that fails to parse is reported as a parse failure with its error, and A1
  is recorded as **untested** rather than as a success.
- These archives are `.zip`, not `.h5ad`. Whatever unpacking is needed to reach a
  metadata table is permitted, since it precedes the frozen scorer; the scorer
  itself is not adapted to the result.

### What a failure would mean

If A1 or A3 fails, §48's claim that the field's designs almost never admit this
estimate is too strong and must be softened to the RNA/protein collection alone.
If A2 fails, the single-context pattern is specific to RNA CRISPR screens rather
than general, which materially weakens the motivation for the whole design
argument. I will report either outcome in the manuscript.

---

## Part B — a standing prediction for atlases released from 2026-09-07 onward

I could not verify, at registration, the design of any specific unreleased atlas:
this session's search budget was exhausted, the Tahoe-100M dataset card announces
no successor, and I will not assert design parameters I cannot check. Naming a
target I cannot verify would be theatre. So Part B is written to apply to whatever
appears, which also removes any ability to choose a flattering target.

### Eligibility

Any single-cell perturbation atlas whose data and metadata become public on or
after **2026-09-07**, and which reports **≥ 2 contexts**, **≥ 20 perturbations**,
and enough metadata to determine replicate structure. The first three eligible
atlases are the test set. They are taken in order of public release date, not
chosen.

### The prediction

For each eligible atlas, before any interaction is estimated:

1. `perturbdesign audit` is run on the released object to obtain
   (contexts, perturbations, replicates) and a detection floor **f**.
2. The interaction share **s** is then measured with the frozen estimator.
3. **Prediction: `s > f` if and only if the atlas resolves a reproducible
   interaction** — i.e. the sign of `s − f` predicts, atlas by atlas, whether a
   permutation test on the interaction rejects at *P* < 0.05.

**Falsified if the prediction is wrong for 2 or more of the first 3 eligible
atlases.** One miss out of three is consistent with the 83–100% power and 0–17%
false-positive rates measured in §46's calibration and will not be claimed as a
success either; it will be reported as one miss.

### The sharper, riskier prediction

The binary test above is easy to pass when an atlas is far from its threshold.
The informative case is an atlas built near it, so:

4. **For any eligible atlas with one replicate per condition, the measured
   interaction share will not be separable from zero by a permutation test**,
   regardless of how many cells it sequenced. This is the calculation's strongest
   and most falsifiable commitment — a single well-powered counterexample, an
   unreplicated atlas that nonetheless resolves a reproducible interaction, kills
   the premise that pairs and not cells set the precision.

I expect assertion 4 to be the one that is actually tested, because unreplicated
designs are common (§48: 23 of 38).

### Analysis plan, fixed in advance

- The test is executed by `scripts/frozen_prospective_test.py`, **run without
  modification**. Its SHA-256 is recorded below. If the script must change to run
  at all — a format change, an API break — the diff is recorded in the outcomes
  file and the run is flagged as modified.
- Estimator: replicate-covariance decomposition with both main effects taken out
  of fold and the context effect estimated **within a replicate plate** (§46,
  `perturbmodel.celldrug`).
- Null: permutation of perturbation labels within context, 200 draws, the same
  construction as §35 and §39.
- Contexts and perturbations are those the atlas itself annotates. I do not
  re-cluster, re-label, or select a favourable subset.
- No atlas is excluded after its result is known. Exclusions may only be made on
  the eligibility criteria above, which are fixed here.

---

## Standing commitments

- Outcomes are appended to `docs/PREREGISTRATION_OUTCOMES.md`, never edited.
- Failures are reported in the manuscript with the same prominence as successes.
  This project has withdrawn four claims already (§28's generalisation, the
  pair-property claim, the genotype claim, the chromatin corroboration); one more
  is not costly enough to justify hiding it.
- `snr = 0.20` is frozen. If a future result is only obtained by changing it, that
  is a refutation of the calibration and is reported as one.

### Frozen artefacts

| file | SHA-256 |
|---|---|
| `src/perturbmodel/design.py` | `fd2a8814cb39a606baf7d2cb0b493d5473a9551c628779b814a1322fa555799d` |
| `src/perturbmodel/atlas_meta.py` | `d456516ff4fe81424fca7529daaab190a8ef2f58abd89d0b3632b7c45a6419da` |
| `scripts/frozen_prospective_test.py` | `49c33236770a28a88cf4895ba4fd2d67f1a32cafc4886c0837acdbbb7bec9503` |
