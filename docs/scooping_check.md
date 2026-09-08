# Scooping check, resolved 2026-09-08

`docs/novelty_audit.md` left two risks unverified because direct fetch failed
(rate-limited and HTTP 403). Both are now **verified from published abstracts**
via the Europe PMC REST API, which serves records the publisher sites block.

Abstracts only — full texts were not obtainable — so what follows is a
conservative reading. Where an abstract is silent I record it as *not stated*
rather than as absent.

---

## Risk 1: Svensson et al. — **not a scoop**

**"Back to basics: Observed statistics are sufficient to predict drug responses"**
Svensson V, Khan U, Heydari H, Ubas AA, Thomas N, Merico D, Goodarzi H, Yu J,
Alidoust N, Gandhi S. bioRxiv 2026, doi 10.64898/2026.06.09.731197
(Europe PMC `PPR1251624`).

Introduces **Rhaister**, a predictor operating on screen-level summary
statistics: measure a few perturbations in a new context, predict the rest by
learning how response patterns vary across reference contexts. Also **Rhaister-O**,
zero-shot from baseline expression. Claims to match or exceed virtual-cell models
while training in seconds.

| our claim | touched? |
|---|---|
| N1 design calculation / detection limit | **no** — no power calculation or detection limit is described |
| N2 38-dataset field survey | **no** |
| N3 cells-vs-contexts exchange rate | **no** |
| N5 Tahoe replicate structure | **no** |
| N6, N7 model/benchmark rules | **no** |
| the few-shot "arithmetic floor beats the fine-tuned model" result | **substantially overlapped** |

Rhaister is the same *thesis* as our few-shot finding — that a cheap
summary-statistic method matches expensive models given a few probe compounds.
That result was **already dropped from both papers during the 2026-09-07 split**,
so nothing currently claimed is affected; but it should not be revived without
citing this.

**Verdict: not a scoop of the design paper.** It is a prediction paper; ours is a
measurability paper.

---

## Risk 2: Shen — **two preprints, one a material partial overlap**

Both are single-author (Shen H.) audits of Tahoe-100M on Research Square.

### 2a. `rs-10448056` — the one the old audit flagged

**"Two axes of drug transcriptional response, and a mechanistic correlate that
organizes them: a conservation–divergence audit of Tahoe-100M"**, 2026-07-23,
doi 10.21203/rs.3.rs-10448056/v1 (`PPR1283738`).

Defines **context-robustness**: "the fraction of a drug's response that is a
background-independent conserved 'core' rather than a cell-line-specific
'periphery'". That **is** a per-drug context-specificity score, so the old audit's
suspicion was right. It further reports that "inhibitors of broad signalling hubs
are robust, agents that depend on a specific target or lineage are labile" — the
same shape as our mechanism ranking.

**Impact: low, because we already withdrew that claim.** Our mechanism ranking is
significant within Tahoe but does not transfer across platforms (+0.103 against
LINCS at 10 shared classes), and is listed as withdrawn in
`docs/novelty_current.md`. Shen reaches a similar descriptive result on the same
single atlas; we cannot claim it and no longer do.

### 2b. `rs-10846736` — **the serious one, and it is new**

**"Basal expression predicts context-specific drug responses only where they are
reproducible: a leakage-safe 47-line LOCO benchmark on Tahoe-100M"**,
**2026-09-01**, doi 10.21203/rs.3.rs-10846736/v1 (`PPR1309852`). *Posted one week
before this check; it post-dates the old audit, which could not have found it.*

Its three headline findings, quoted:

1. "on the full benchmark **no model meaningfully beats a drug-mean baseline**"
   (best model +0.004 reliability-weighted residual *R*²).
2. "**repeat measurements explain why**: for the 107 drugs with replicate
   pseudobulk data, the context-specific residual is itself **only weakly
   reproducible** (median repeat reliability 0.067; only 4 drugs exceed 0.2),
   **capping what any context predictor could achieve**."
3. "where the residual *is* reproducible, basal expression *is* predictive: model
   skill tracks the repeatability ceiling (Spearman ρ = 0.76, *p* = 3×10⁻²¹)."

Concluding: "context-specific drug-response prediction from basal expression is
**bounded primarily by the reproducibility of the target signal**, and …
**benchmarks should report repeatability ceilings alongside model skill**."

**This is the closest published work to our argument, and it was not on the
record when the design paper was framed.**

| our claim | status against `rs-10846736` |
|---|---|
| The context-specific component of Tahoe's response is barely reproducible, which bounds what can be learned | **independently established there**, by per-drug repeat reliability rather than by our decomposition |
| N5 Tahoe's replicate structure is thin and load-bearing | **partly anticipated** — he uses the 107 replicated drugs for exactly this purpose |
| N1 a *prospective* design calculation: detection limit from (contexts, perturbations, replicates) before data exists, validated on five atlases | **not touched.** He measures reliability post hoc on one atlas; he does not predict it from a design |
| N2 the 38-dataset field survey | **not touched** — Tahoe only |
| N3 cells-vs-contexts exchange rate (2% of cells → 92% of effect) | **not touched** |
| N6 *d* ≥ 0.05 × n_contexts | **not touched** |
| N7 benchmark gains vs context count | **not touched.** He argues benchmarks should report *repeatability* ceilings; we argue they should report the *context count* behind the mean baseline. Adjacent, not the same |

---

## What this changes

**The design paper's core survives**: N1, N2, N3, N6 and N7 are untouched by any
of the three preprints, and they are the paper's argument. Nothing has to be
withdrawn.

**But its framing has to change.** The motivating observation — that Tahoe's
context-specific signal is too weakly reproducible to support what is being asked
of it — is now published by someone else, five weeks before our submission
window. The design paper must cite `rs-10846736`, state plainly that it reaches a
compatible conclusion on Tahoe by a different route, and position our contribution
as the part he does not do: **turning a post-hoc reliability measurement into a
prospective design calculation, and testing it across atlases rather than within
one.**

That is a defensible and still-novel position. It is a weaker one than "nobody has
noticed this", which is what the paper currently implies.

**Two practical consequences.** The convergence is evidence the question is live
and that others are moving on it, which argues for submitting sooner rather than
adding scope. And the earlier estimate of publication odds was made without this
preprint on the record; it should be revised down, because a reviewer who knows
`rs-10846736` will ask what is left, and the honest answer is "the prospective
calculation and the cross-atlas survey" rather than "the observation itself".

## Method note

Publisher sites blocked both (bioRxiv 429, Research Square 403). The Europe PMC
REST API served all three records:

```
https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:PPR1309852&resultType=core&format=json
```

Worth remembering: when a preprint host blocks automated access, Europe PMC
usually has the record, including the full abstract.
