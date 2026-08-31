# pertdecomp — what it is and why it exists

A perturbation atlas is usually summarised by how much of the response is shared
between cellular contexts and how much is context-specific. That number is
quoted, compared between studies, and used to argue that context-transfer models
have headroom. `pertdecomp` computes it — and, more importantly, refuses to
compute it when the data cannot support it.

The case for the tool is not that the estimator is novel. It is that **every
naive route to this number is wrong in a way that is invisible in the output**,
and we hit four of them in this project before the tool existed. Each is now a
guard rail with a measured cost.

## The four failure modes, with what they cost us

| # | Naive choice | What it does | Measured on our data |
|---|---|---|---|
| 1 | Estimate interaction as residual **variance** | Absorbs noise, which does not cancel | 82% "interaction" vs a true value near 0 on the same Tahoe conditions |
| 2 | Pool **same-batch** comparisons | Batch state is shared signal | Within-plate pairs agree ~7× better than cross-plate; share 41% vs 28% |
| 3 | Use the **in-sample mean** as the shared response | Residuals sum to zero, forcing E[r_a·r_b] = −σ²/n | Drove per-drug covariance negative for 21 of 24 drugs; clamping then produced exact zeros |
| 4 | Treat **different doses** as replicates | Doses are different conditions, not repeats | Tahoe: true replicates covary −0.0031, cross-dose +0.0062, p = 3×10⁻¹²⁰ |

Failure 3 has a subtlety worth stating because we got it wrong twice: the
obvious fix, leave-one-*condition*-out, makes it **worse**, since the mean
subtracted from replicate A still contains replicate B. Only leaving out the
whole *context* removes the coupling.

## What the tool does about them

- **Replicates are mandatory.** Interaction is the covariance of residuals
  between independent replicates of the same context × perturbation × dose.
  Given no replicate column, it emits a warning and reports the interaction as
  NaN rather than a number.
- **Same-batch pairs are excluded**, and the within/cross-batch inflation is
  reported so the user sees what pooling would have cost.
- **The shared response is leave-one-context-out**, and the interaction is
  reported against a **matched cross-context null**, so the residual-construction
  offset is subtracted rather than assumed negligible. The tool prints the raw
  covariance, the offset, and what the uncorrected answer would have been.
- **Per-perturbation indices carry their pair counts** and an `estimable` flag,
  distinguishing "not context-dependent" from "not measurable here".

## Why this is worth a tool rather than a paragraph in a methods section

Applied across five datasets it changes the conclusion, not just the decimal:

| dataset | contexts | true same-dose replicate pairs | interaction share |
|---|---|---|---|
| LINCS phase 1 | 70 | 6.1 M | **57%** |
| LINCS phase 2 | 30 | 1.4 M | 47% |
| PRISM (viability) | 738 | X1/X2/X3 plates | 22% |
| **Tahoe-100M** | 47 | **6,482, from 25 of 379 drugs** | **not estimable** |
| OP3 (primary immune) | 6 types × 3 donors | 18 per compound | not estimable per compound |

The last two rows are the point. Tahoe-100M is the largest single-cell
perturbation atlas published, and it cannot support the statistic it is most
often used to produce, because only 4.7% of its (line, drug, dose) combinations
appear on more than one plate. A tool that reports "not estimable" is the
difference between that being visible and being quietly papered over with a
cross-dose pairing that yields a publishable 21.6%.

## Use

```python
from perturbmodel.decompose import decompose
res = decompose(adata, context="cell_line", perturbation="drug",
                replicate="plate", dose="dose", control="DMSO")
print(res.report())        # shares, batch inflation, null offset, warnings
res.per_perturbation       # per-compound index with pair counts + estimable flag
```

```bash
python -m perturbmodel.decompose --h5ad atlas.h5ad \
    --context cell_line --perturbation drug --replicate plate \
    --dose dose --control DMSO
```

## Design recommendation this yields

For an atlas intended to study context-transfer, the two axes that bind are the
number of **contexts** (§10–11: 47 gives 4% recovery of effects that 400 recovers
at 72%) and the presence of **true replicates** (§14). Cells per condition binds
neither. An atlas of 400 lines × fewer cells × two plates per condition would
answer what 100 million cells across 47 lines and one plate cannot.
