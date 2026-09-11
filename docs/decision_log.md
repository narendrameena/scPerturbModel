# Decision log — test before applying

**Opened 2026-09-11.** Every restructuring decision is written here **before** it is
applied, with a test that could fail. The test is run, the result recorded, and only
then is the decision acted on or abandoned.

## Why this exists

The failure that has cost this project the most is not wrong analysis. It is
**judgment applied before it was tested**. The pattern, five times in five days:

| date | judgment applied | tested afterwards | outcome |
|---|---|---|---|
| 09-07 | Replogle's 40 `batch` levels are replicates → headline becomes 2 of 38 | checked `z_gemgroup_UMI` | gem groups, not replicates. Reverted before publishing. |
| 09-09 | the pairs^−1/2 law is refuted → fit a U-statistic correction | checked the bootstrap | the refutation was an identity of the bootstrap. Both withdrawn. |
| 09-09 | claimed a verification while the run was still going | checked the exit code | killed at timeout, never ran. Retracted. |
| 09-10 | contexts and replicates are not interchangeable → invert the prescription | applied the CV control to the right coefficient | *P* = 0.70. Withdrawn same day. |
| 09-10 | "~120 lines" → corrected to "50 lines" | applied the script's own filters | 47 lines. The correction was also wrong. |

In every case the test existed, was cheap, and was run *after* the claim was
committed. The cost was not the error; it was that the error had already been
built on.

## The rule

No decision below is applied until its test has been run and recorded. A decision
whose test cannot fail is not tested — it must be restated until it can.

## Criteria, in priority order

1. **Biological defensibility first.** Where two framings are available, take the
   one whose claim is a statement about biology that could be wrong, over the one
   that is most novel-sounding or most convenient. A result that survives because
   it is vague is not a result.
2. **Survives adversarial check.** Prefer claims that have already been attacked
   and held over claims that are merely unrefuted because nobody looked.
3. **Traceable to a committed table.** No number enters without its generating
   script and output committed.
4. **Scale of evidence.** 737 lines beats 47; 8,427 compounds beats 5 atlases.

---

## Decisions

### D1 — What is the paper's spine?

**Status: PENDING — three agents auditing. Not applied.**

*Options.* (a) Pharmacogenomic power: genotype × drug across 737 PRISM lines, the
finding that molecular state predicts response where genotype does not, and that
~50-line panels are underpowered ~10×. (b) The design/estimator work, retargeted.
(c) Something else an agent proposes.

*The test, fixed before the answer is known.* The spine must clear all four:
1. its central claim is biological and falsifiable — state the experiment that
   would refute it;
2. every number in it recomputes from a committed table;
3. it has been adversarially attacked at least once and survived;
4. it does not depend on any claim withdrawn in §46–§52.

*Refutation condition.* If no option clears all four, the honest output is a
preprint or a specialist journal, and the log says so.

### D2 — Does adding data change anything?

**Status: PENDING. Not applied.**

*The test.* For each proposed dataset, name the specific claim it would change and
the result that would change it. A dataset that would "strengthen" a claim without
a stated way of weakening it does not qualify.

*Refutation condition.* If no addition has a named failure mode, adding data is
volume, not evidence, and is not done.

### D3 — Which genetics claims survive as load-bearing?

**Status: PENDING — agent auditing. Not applied.**

*The test.* Recompute from committed tables; check the multiplicity correction on
1.5M and 111,589 tests; establish whether the "de novo biomarker recovery" is
circular (were the 11 GDSC biomarkers selected before or after seeing which
replicated?).

*Refutation condition.* If the biomarker set was chosen post hoc, "11 of 11
replicated" is selection, not validation, and must be withdrawn.


### D0 — Resolve the §52 provenance flag (retrospective, now closed)

*Decision.* Replace "reported, not reproduced" on three §52 figures with a
confirmed statement, or leave the flag standing.

*Test, fixed before the run finished.* Re-run `cells_vs_replicates.py` at a
different seed with more draws. If 2%-of-cells retention comes back near 92%, the
published figure is real and the audit's "four-draw artefact" claim is wrong. If
it comes back near 100%, the audit is right.

*Result.* Seed 7, 36 draws (12 cells-arm): **102%**. Audit confirmed; the
published 92% is a low-draw excursion. The contexts arm's spread is 0.075 against
0.084 published — large under both, so that half is qualitatively stable.

*Applied.* §52 updated with my own numbers replacing the borrowed ones. This is
the first decision in the log to complete the cycle: stated, tested, then applied.
