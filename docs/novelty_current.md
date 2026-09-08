# What is actually new, as of 2026-09-08

`docs/novelty_audit.md` was written before the estimator correction and the
withdrawals that followed. Its headline conclusion — that the most defensible
novelty is "the pair-property result" — names a claim this project has since
**withdrawn**, and its second-ranked mechanism result does not reproduce across
platforms. This file replaces its verdict; the original is kept because its
literature search is still the record of what was looked for.

**Short answer: 7 claims are novel and survive. 2 are novel but partly
anticipated. 5 were novel and are withdrawn. 4 are informative negatives.**

## Novel and surviving (7)

| # | Claim | Where | Why nothing published covers it |
|---|---|---|---|
| N1 | **A design calculation for perturbation atlases**: the interaction is a covariance between replicates, so precision is set by replicate *pairs* and cell count does not enter. Predicts 5/5 published atlases with one shared constant. | Design §1 | No published power calculation for this quantity. The two scooping candidates are about prediction (Svensson) and a descriptive Tahoe audit (Shen), neither of which computes a detection limit. |
| N2 | **One of 38 scPerturb datasets can support the estimate at all**; of ten drug screens, three use >1 context and one of those replicates. | Design §2 | scPerturb (Peidli 2024) catalogues the datasets but does not score their designs against an estimator's requirements. |
| N3 | **Cells and contexts are not interchangeable**: 2% of Tahoe's cells retain 92% of a known interaction; 10% of contexts collapse precision. | Design §1 | The field's scaling arguments are asserted, not measured. |
| N4 | **Omitting a context's general sensitivity inflates published interaction indices ~1.5×**, demonstrated on another group's deposited data (TRADE) and in simulation. | Design §3, Measurement | The α term is named in the literature (Ben-David 2018; Lim & Pavlidis 2021) but its effect on interaction *estimators* is not quantified anywhere found. |
| N5 | **Tahoe-100M's replicate structure exists but is thin**, and the natural substitute — pairing across doses — moves the estimate from 0.5% (indistinguishable from zero) to 9.2%. | Design §4 | No published Tahoe analysis reports this control. |
| N6 | **A fixed *d*-dimensional context embedding falls behind as an atlas grows**; interaction dimensionality rises ~0.05 directions per context. Usable rule *d* ≥ 0.05 × n_contexts. | Design §5 | Novel as a stated scaling rule. |
| N7 | **Benchmark gains are not comparable across datasets with different context counts.** The perturbation-mean baseline strengthens monotonically (r 0.550 → 0.666, 2 → 45 contexts), so a model equally good everywhere looks 19% stronger against a 2-context baseline. | Design §5 | Explains a standing disagreement between benchmark papers and method papers; not previously stated. |

## Novel but partly anticipated (2)

| # | Claim | Nearest prior work |
|---|---|---|
| A1 | Additive/interaction split of a **chemical** perturbation atlas at scale | TRADE (genetic, 4 lines, 56/44) and sci-Plex (3 lines, 48%) get numerically similar splits. Ours differs in modality, scale (47 contexts) and in validating the split against replicate covariance rather than asserting it from correlations. A reviewer who knows TRADE will ask; the answer is method and scale, not the number. |
| A2 | Molecular state predicts the interaction where genotype does not | Schlüter & Schönhuth (2025) report the copy-number-over-mutations ordering. We reproduce the ordering and then **withdraw** our own version of it under a cluster bootstrap. What survives is expression ≫ genotype, which is consistent with but broader than published work. |

## Novel, then withdrawn (5)

These were the previous audit's leading results. They are listed because a reader
of that file will look for them.

| # | Claim | Status |
|---|---|---|
| W1 | The interaction is a **pair** property, not a cell property — contradicting Lim & Pavlidis | **Withdrawn.** A line-level responsiveness factor exists as published; the relation is 2.0× the property, a measurement rather than a refutation. |
| W2 | Context-dependence is set by drug **mechanism**, reproducibly | **Withdrawn across platforms.** Significant within Tahoe (*P* = 6.0×10⁻⁴) but +0.103 against LINCS-1 at 10 shared classes. |
| W3 | Transcriptional and viability context-dependence are **decoupled** | **Withdrawn 2026-09-08.** The only powered comparison (67 classes) shows +0.279 (*p* = 0.023) — weak positive agreement, not decoupling. |
| W4 | Copy number beats mutations for predicting the interaction | **Withdrawn.** Cluster bootstrap CI [−0.002, +0.015] includes zero. |
| W5 | §28's generalisation of the inflation result to sci-Plex and CMap statistics | **Withdrawn permanently** after three failed replication attempts. |

## Informative negatives (4)

Not novel *claims*, but results a reviewer would count.

- **0 of 3 published statistics reproduced** (sci-Plex 74.4% vs published 48%; CMap 5.0% vs 26%; TRADE 39% vs 56%). Reported as failures to reproduce, with no claim made about those papers.
- **Genotype does not linearly predict the interaction** even at 737 lines; the 47-line negatives in the literature are a power ceiling (80% power needs ~400 lines).
- **Chromatin shows no detectable interaction** on features that demonstrably respond (33.9% against a null of 33.0%, *P* = 0.45) — an absence of evidence at three contexts, not evidence of absence.
- **Two proposed methods rejected on their own benchmarks** (a spectral estimator, a two-stage discovery procedure).

## The honest tally

Counting only what a reviewer would call a new, defended claim: **7 novel and
surviving, 2 novel but anticipated in part**. The design paper carries N1–N3 and
N5–N7 — six of the seven — which is the main argument for it being the stronger
of the two papers.

Against that, five claims that the previous audit ranked as this project's best
work are now withdrawn, and one of them (W3) was withdrawn only after a stale
table was re-run today. That is the correct outcome each time, but it is also the
reason to treat any remaining un-re-run analysis as provisional.

## Unresolved scooping risk

Both items from the original audit remain **unverified by direct fetch**:

1. **Svensson et al.** (bioRxiv, June 2026) — a Tahoe author building the same
   baseline ladder on the same data. The sweep found no residual analysis or
   additive/interaction split, which is the opening, but the preprint has not
   been read in full.
2. **Shen et al.** (Research Square 2026) — reportedly defines a
   "context-robustness" axis on Tahoe-100M. Direct fetch returned HTTP 403.

Neither touches N1 or N2, which is where the design paper's novelty sits. Both
should be read before submission; this session's web-search budget was exhausted
and could not settle them.
