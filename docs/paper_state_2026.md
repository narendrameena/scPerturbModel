# Adduri et al. 2026, *Cell* — "Predicting cellular responses to perturbation across diverse contexts with State"

Read 2026-09-01. This is the most directly relevant paper to our work: the
flagship model for the problem we have been measuring, from the institute that
produced Tahoe-100M. It deserves a careful reading rather than a citation.

---

## What the paper is

**State** is a two-part architecture from Arc Institute:

- **SE (State Embedding)** — a bidirectional transformer trained on **167 million
  cells** of observational data (scBaseCount 72M, CZ CELLxGENE 59M, Tahoe-100M
  35M) to produce a cell-state embedding. Trained with a loss along two axes:
  across genes within a cell, and across cells for each gene.
- **ST (State Transition)** — operates on *sets* of cells rather than single
  cells, predicting how a population of unperturbed cells shifts under a
  perturbation. Self-attention over the set is load-bearing; removing it degrades
  performance.

The stated motivation is one we agree with entirely: "Current deep learning
methods do not consistently outperform linear models when generalizing
perturbation effects across cellular contexts," which they attribute to two
noise sources — unmodelled biological heterogeneity within a population, and
technical variation across datasets.

They also release **Cell-Eval**, an evaluation suite spanning perturbation
discrimination, Pearson correlation of effects, DE overlap, p-value calibration
and effect-size prediction.

## Headline results

| dataset | perturbation discrimination | Pearson of effects |
|---|---|---|
| Tahoe-100M | +66% | +91% over second best |
| Parse-PBMC | +29% | +48% |
| Replogle-Nadig | +10% | +10% |

Baselines: a linear model, CPA, scVI, GEARS, scGPT, plus two mean baselines —
**"context mean"** (average expression in that context across all perturbations)
and **"perturbation mean"** (average perturbation effect across training contexts
applied to the context's basal expression). The latter is exactly our additive
baseline. Trained on the top 2,000 HVGs, as we do.

They note that for the Pearson metric the runner-up was the *context* mean rather
than the perturbation mean, arguing State is "not trivially like either mean
baseline."

## The three things that matter for us

### 1. Their headline context-generalisation setting is our few-shot setting

The main "held-out context" experiment gives each model **30% of the
perturbations in the test context during training**. That is not an unseen
context; it is an *underrepresented* one. The paper is explicit and reasonable
about this — it reflects the realistic case where a screen covers some contexts
deeply and others sparsely.

This corroborates rather than contradicts our §6: what makes a new context
predictable is **measuring some perturbations in it**, not describing it. We
found ~20 arbitrary compounds recovers 98% of the achievable gain, and that
designing the probe panel does not beat choosing it at random. State's central
result is the same phenomenon at scale with a better model.

### 2. Their reported gains are largest exactly where the baseline is weakest

For the zero-shot arm, five "query datasets" are used, holding out one context at
a time: Srivatsan (**3** cell lines), McFaline-Figueroa (**3** contexts),
Replogle-Nadig (**4** cell lines), Jiang (**6** contexts), Parse-PBMC (**18** cell
types) — 34 contexts total.

The paper states plainly: "For larger datasets like Parse-PBMC, using the State
embedding achieved more than 19% improvement... In smaller genetic-perturbation
datasets (e.g., Jiang et al. and McFaline-Figueroa et al.), ST+SE yielded
**several-fold improvements**."

That ordering is what our §2 predicts for reasons that have nothing to do with
the model. Under a hold-one-context-out protocol the perturbation-mean baseline
is estimated from *m−1* contexts, so it is built from **2 contexts** in the
3-context datasets and **17** in Parse-PBMC. Measuring the baseline's quality as
a function of that number on Tahoe (`additive_vs_n_lines.py`, extended to the
relevant range):

| contexts averaged | baseline r | handicap vs a 45-context baseline |
|---|---|---|
| 2 | 0.550 | **21.1%** |
| 3 | 0.592 | 12.4% |
| 4 | 0.606 | 9.8% |
| 6 | 0.625 | 6.4% |
| 18 | 0.655 | **1.7%** |
| 45 | 0.666 | 0% |

A **fixed** model — no better in one dataset than another — would appear ~21%
stronger against a 2-context baseline than against an 18-context one.

**This is a partial confound, not a refutation.** The differential it creates
(~19 points) is real and runs in exactly the reported direction, but it does not
account for *several-fold* gaps. State is doing something the baseline is not.
The correct reading is that a portion of the spread in reported improvements
across benchmarks is an artefact of how many contexts each benchmark's baseline
was estimated from, and that improvements should be compared at matched context
counts. The paper's own caution supports this: it notes the several-fold gains
"did not consistently extend to other metrics."

### 3. They trained on Tahoe, which makes our plate-14 finding relevant to them

SE is trained on 35M Tahoe cells and ST is pretrained on Tahoe for the transfer
experiments. The Tahoe authors reserve **plate 14 — a designed biological
replicate of plate 6 — for validation rather than training** (§14). Any pipeline
inheriting that split has no same-dose replicates in its training data, which is
the condition under which the context × perturbation term is not identifiable at
all (§21). This does not affect State's *predictions*, which are evaluated
against held-out data directly, but it does bear on any attempt to interpret a
model's learned interaction structure as a measurement.

## What we should change

- **Cite State as the state of the art** and as agreeing with us on the
  diagnosis: deep models have not consistently beaten linear ones at context
  transfer, and heterogeneity is a principal reason.
- **Frame our few-shot result as corroborating theirs**, not competing: their
  30%-of-perturbations setting and our 20-compound probe finding are the same
  claim from different directions.
- **Add the baseline-handicap calibration** to the benchmark-design argument,
  with the explicit statement that it explains part and not all of the spread.
- **Recommend that benchmarks report the number of contexts behind the mean
  baseline**, since that number alone moves apparent model gains by ~20 points.
