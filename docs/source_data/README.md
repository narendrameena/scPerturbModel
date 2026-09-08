# Source data for numbers quoted in the manuscripts

`results/` is gitignored — it is large and fully regenerable. That policy is
right for figures and intermediate artefacts and wrong for one narrow case, which
the 2026-09-07 audits found the hard way: a number quoted in a manuscript whose
generating table was never committed cannot be checked by anyone, including us.

The Tahoe-vs-LINCS mechanism correlation was quoted as +0.09 in `RESULTS.md` §16
and +0.19 in §31, while the only saved table reproduced +0.56 — the
*pre-correction* value §16 explicitly retires. The corrected run's output had
never been saved, so neither quoted value was reproducible.

This directory is the fix. **Any small table backing a number quoted in
`MANUSCRIPT_DESIGN.md` or `MANUSCRIPT_MEASUREMENT.md` is committed here**, with
its provenance recorded below. Large outputs stay in `results/`.

## Contents

### `three_platform_mechanism_cdi.csv`, `three_platform_mechanism_rho.csv`

Per-mechanism-class context-dependence index on each platform, and the pairwise
Spearman correlations between platforms.

| | |
|---|---|
| generated | 2026-09-08 |
| script | `scripts/three_platform_synthesis.py` |
| variant saved | active compounds only (compounds with a detectable additive effect) |
| inputs | `results/tables/drug_context_dependence.csv` (Tahoe, matched-null, 2026-09-01)<br>`results/tables/lincs_mechanism_ranking_phase1.csv` (LINCS-1, matched-null, 2026-09-03)<br>`results/tables/lincs_mechanism_ranking_lincs.csv` (LINCS-2, 2026-09-03)<br>`results/tables/prism_decomposition.csv` (PRISM, three-way via `celldrug.general_sensitivity_by_rep`, 2026-09-03) |
| minimum class size | 3 compounds |

All four inputs are on the corrected estimator. The previously committed version
of this table was generated 2026-08-31, **before three of its four inputs were
regenerated**, which is why it reproduced superseded values.

Pairwise Spearman correlations, active compounds only:

| comparison | classes | ρ | *p* |
|---|---:|---:|---:|
| Tahoe vs LINCS-1 | 10 | +0.103 | 0.777 |
| Tahoe vs PRISM | 11 | −0.155 | 0.650 |
| **LINCS-1 vs PRISM** | **67** | **+0.279** | **0.023** |
| LINCS-1 vs LINCS-2 | 19 | +0.511 | 0.026 |
| Tahoe vs LINCS-2 | 6 | −0.200 | 0.704 |
| LINCS-2 vs PRISM | 23 | +0.277 | 0.201 |

Six comparisons were made; Bonferroni at α = 0.05 requires *p* < 0.0083, which
none of them meets. Read the two nominally significant rows accordingly.

### `cross_lab_identity.csv`, `cross_lab_summary.csv`

Cross-laboratory transfer of the line-specific drug response, transcription and
viability arms.

| | |
|---|---|
| generated | 2026-09-08 (transcription arm re-run after the keying fix in `383e890`) |
| script | `scripts/cross_lab_transcription.py`, `scripts/cross_lab_reproducibility.py` |

Transcription arm, corrected:

| comparison | median *r* | pairs |
|---|---:|---:|
| LINCS p1 vs p2 (within-lab) | 0.034 | 5,803 |
| Tahoe vs LINCS p1 (cross-lab) | 0.012 | 317 |
| Tahoe vs LINCS p2 (cross-lab) | 0.022 | 172 |

Matched ceiling 0.044 vs cross-lab 0.012 → **reproducible fraction 27.3%**
(viability: 56%). Identity: 15/16 within the Broad, 2/6 across laboratories.

The pre-fix version of this table reported 46%, within-lab 0.061 and 16/16, from
profiles keyed (line, last-line-seen) instead of (line, compound).
