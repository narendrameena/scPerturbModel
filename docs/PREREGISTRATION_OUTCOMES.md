# Pre-registration outcomes

Append-only. Entries are added below, never edited or removed. The registration
is `docs/PREREGISTRATION.md`, frozen at tag `prereg-2026-09-07`.

---

## 2026-09-07 — Part A settled

**Order of operations, verifiable from git:** the registration was committed and
pushed (`49a1e79`, tag `prereg-2026-09-07`) *before* the three unseen archives
were downloaded. The download command runs in the shell history after the push.
Anyone can check that the predictions predate the data by comparing the tag's
push time to the figshare access.

### Result

| # | Prediction | Outcome |
|---|---|---|
| **A1** | `MimitouSmibert2021` does not support a context × perturbation interaction estimate | **HELD** |
| **A2** | Both `Liscovitch-BrauerSanjana2021` files are single-context | held (registered as *not* evidence) |

Scored by `perturbmodel.atlas_meta` at the frozen commit. The archives ship
`obs.csv` rather than `.h5ad`, so `pandas.read_csv` replaced `read_obs`; the
column-preference lists, the `MISSING` set and every downstream rule were used
unmodified, as the analysis plan permits.

| dataset | cells | context | perturbations | replicates/condition | pairs |
|---|---:|---|---:|---:|---:|
| `MimitouSmibert2021` | 10,018 | CD4+ T cells (1) | 6 | 1 | 0 |
| `Liscovitch-BrauerSanjana2021_K562_1` | 8,723 | K562 (1) | 22 | 1 | 0 |
| `Liscovitch-BrauerSanjana2021_K562_2` | 12,788 | K562 (1) | 84 | 1 | 0 |

All three are single-context with one replicate per condition, so the pair count
is zero and the interaction is not estimable in any of them at any effect size.

### What this is worth

**One informative prediction, one success.** A1 was the only prediction in Part A
carrying information, and *n* = 1 is weak evidence. It is recorded as one
dataset, not as a confirmation of §48. Its value is procedural: the prediction
was fixed, public and falsifiable before the file was fetched, and it could have
failed.

### What went wrong first, and is left on the record

The first draft of Part A claimed the *entire* scPerturb ATAC article was unseen
and registered four predictions (A1–A4) about all six files. That claim was
false: three files — `PierceGreenleaf2021_{GM12878,K562,MCF7}`, i.e. Spear-ATAC —
had been downloaded on 2026-09-04 and analysed in §35 and §39, and Spear-ATAC is
multi-context *and* replicated, which is exactly the case those predictions were
counting. Registering them would have been scoring a test whose answer was
partly known. The error was caught before pushing; the draft predictions are
withdrawn in full and the episode is documented in the registration itself rather
than deleted.

---

## 2026-09-07 — an error found while extending §48, and its correction

Not a registered prediction; recorded here because it changed a published number's
justification and because the first attempt at the fix was itself wrong.

**The observation.** scPerturb distributes some multi-context *studies* as
per-context *files*. Scoring per file therefore undercounts multi-context designs.
Two cases in the 38:

| study | files | contexts per file | union |
|---|---:|---|---:|
| `ReplogleWeissman2022` | 3 | 1, 1, 1 | 2 (K562, RPE1) |
| `TianKampmann2019` | 2 | 1, 1 | 2 (iPSC, iPSC-derived neuron) |

Grouped by study, single-context datasets fall from **34/38 files (89%)** to
**17/22 studies (77%)**.

**The wrong fix.** On seeing that `ReplogleWeissman2022` unions to 2 contexts,
2,056 shared perturbations and ~40 annotated `batch` levels per condition, I was
about to report that §48 had missed a dataset which does support the estimate,
and to raise the headline from 1 to 2.

**Why it was wrong.** Replogle's `batch` is the 10x gem group — the file also
carries `z_gemgroup_UMI`, which names it. A pooled library is transduced once and
the cells are then distributed across gem groups for capture. Two cells with the
same guide in different gem groups share the transduction, the culture, the
selection and the perturbation duration, and differ only in the emulsion. They
are split captures, not independent treatments, and by the criterion stated in
`perturbmodel.atlas_meta` — the same criterion that made §38's inflated indices
diagnosable — they do not constitute usable replicates. `ReplogleWeissman2022`
does **not** support the estimate. `TianKampmann2019` has one replicate
regardless.

**Net effect on §48.** The headline is unchanged: exactly one dataset supports
the estimate. Its justification is now explicit rather than accidental, and the
per-file/per-study distinction is reported as a caveat.

**A limitation this exposes, now stated in §48.** `describe` cannot tell an
independent replicate from a split capture, because the deposited metadata does
not distinguish them — `batch` means both things in different datasets. The
replicate counts it reports are therefore an **upper bound on usable
replicates**, and the count of datasets supporting the estimate is an upper bound
too. The direction of §48's claim is conservative under this bias; the specific
number is not exact.


---

## 2026-09-09 — Amendment to a frozen artefact, declared

`docs/PREREGISTRATION.md` froze `src/perturbmodel/design.py` by SHA-256. **That
file has been modified**, and the test that checks the hash caught it. Declaring
the change rather than re-freezing quietly:

**What changed.** `interaction_se` now uses the U-statistic variance
`u/k + (1−u)/n_pairs` instead of `1/n_pairs`, with `U_FIRST_ORDER = 0.1013`.
See `RESULTS.md` §49: the original form was refuted on 8,427 LINCS compounds
(measured slope −0.370 against a predicted −0.500), and the replacement is fitted
on half those compounds and validated on the other half.

| | SHA-256 |
|---|---|
| frozen at registration | `fd2a8814cb39a606baf7d2cb0b493d5473a9551c628779b814a1322fa555799d` |
| after this amendment | `93632f85ff39afb49290d7f96d4c587e68dbae13c58ea611a206c85dba6d8b01` |

**Why this does not compromise the registration.**

* **Part A already settled**, on 2026-09-08, under the original frozen version.
  Its outcome stands and is not revisited.
* **The change was driven by LINCS calibration data**, which is not the target of
  any registered prediction. No Part B atlas, and no scPerturb dataset named in
  B5/B6, influenced it.
* **The change is a strict generalisation.** `u_first_order=0.0` reproduces the
  original numbers exactly, so any registered prediction can be re-scored under
  either version and both are available.
* **It does not move any registered prediction.** The five-atlas verdicts are
  unchanged (5/5, across the same snr range), and B5/B6 depend on design
  parameters — contexts, replicates — that the variance formula does not touch.

**How Part B is scored from here.** Under the amended version, with the original
reported alongside wherever the two differ. Both hashes are recorded above, so
either can be reconstructed.

**The honest cost.** A registration is worth less once its artefacts move, however
good the reason. The mitigation is that the change is declared, dated, hash-
pinned in both directions, and reversible by a single argument — not that it is
harmless.


---

## 2026-09-10 — Second amendment, and the registered artefact cannot be restored

The amendment declared on 2026-09-09 was made on the strength of `RESULTS.md`
§49, which is **now withdrawn**: the refutation it rested on was an artefact of an
i.i.d. pair bootstrap, and the correction fitted to it deflated the variance at
`n_rep = 2` where the corrected term must be zero. `U_FIRST_ORDER` is reverted
to 0.

`src/perturbmodel/design.py` is therefore on its **third** hash:

| | SHA-256 |
|---|---|
| at registration (2026-09-07) | `fd2a8814cb39a606baf7d2cb0b493d5473a9551c628779b814a1322fa555799d` |
| first amendment (2026-09-09) | `93632f85ff39afb49290d7f96d4c587e68dbae13c58ea611a206c85dba6d8b01` |
| **after reverting (2026-09-10)** | `bbf3b6875a2720379b758c69090e16fd6a23dfe95d4716d670134a55826965c6` |

**The registration-time artefact cannot be restored by reverting the constant.**
The added `n_profiles` helper and the commentary explaining the episode remain in
the file, deliberately — deleting them would erase the record of what happened.
Numerically the module is now equivalent to its registered state
(`min_detectable_share` reproduces `0.0169 / 0.0027 / 0.0222 / 0.0479 / 0.0325`),
but it is not byte-identical and cannot be made so without destroying the audit
trail.

**The honest cost, restated and now larger.** A frozen artefact has been modified
twice in three days, the first time on the strength of a result that was not real.
The registration's remaining value is that Part A settled before any of this, and
that every change is declared with hashes in both directions. That is worth
something. It is worth less than a registration whose artefacts never moved.


## 2026-09-10 (second entry) — `atlas_meta.py` amended, declared

`describe()` gained an exact pair count (`n_pairs_exact`) alongside the existing
median-based one. The median implies zero pairs for an unbalanced atlas — most
conditions measured once, a few many times — and nine of the 38 scPerturb
datasets are in that state, carrying up to 18,007 real pairs while the tool told
each of them "UNRESOLVABLE AT ANY EFFECT SIZE".

**This does not disturb any registered prediction.** Part A settled on 2026-09-08
under the original file, and all nine affected datasets are single-context, so
their verdicts are unchanged and the "1 of 38" headline is unaffected — I checked
every multi-context dataset and the exact and median counts agree. B5/B6 depend on
context and replicate counts, which the change does not touch. The median field is
retained, so the registered scoring rule still evaluates identically.

| | SHA-256 |
|---|---|
| at registration | `d456516ff4fe81424fca7529daaab190a8ef2f58abd89d0b3632b7c45a6419da` |
| after this amendment | `c34819f0829b85ac11b7653c2edeffd9e65bcd5d59ae1e745cca87e9dcdaf4ae` |
