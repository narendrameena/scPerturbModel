# Novelty audit: what is already published, claim by claim

Assessed against the literature sweep in `docs/related_work_perturbation_models.md`
(which included Europe PMC full-text phrase searches) plus full readings of MAP,
XPert and the GLP1R GWAS. Verdicts are graded, and the two items I could not
verify directly are flagged as such.

---

## The closest published precedent: TRADE

**Nadig, Replogle, Pogson et al., "Transcriptome-wide analysis of differential
expression in perturbation atlases", Nature Genetics 57:1228 (2025).**

TRADE splits perturbations into **56% cross-cell-type-consistent and 44%
cell-type-dependent**. That is uncomfortably close to our **59% / 41%** — the
single most important precedent to cite and distinguish, and it appeared in the
journal we were considering.

What differs, and why our claim is not contained in theirs:

| | TRADE | this work |
|---|---|---|
| perturbations | **genetic** (CRISPR) | **chemical**, 379 compounds |
| contexts | 4 cell lines | 47 lines + 6 immune types + 3 lines |
| method | correlation of perturbation profiles between lines | variance decomposition into additive + interaction |
| dose | none | 3 doses, used as the replication control |
| residual | not analysed | the object of study |
| replicate validation | none | cross-plate / cross-donor covariance |
| reproducibility of the split | asserted from correlations | established against permuted nulls |

So a similar *number* exists, obtained a different way, in a different
perturbation modality, with 4 contexts rather than 47. Our contribution is the
decomposition and the demonstration that the non-shared part is a **pair**
property — which TRADE does not test and cannot, at 4 lines.

**Also close, in weaker form:** sci-Plex (Srivatsan et al., *Science* 2020)
reported that of 4,308 DE genes, **48% responded cell-type-dependently** and 22%
identically across its three lines — again a comparable fraction, from three
contexts, with no decomposition and no replicate-based validation.

---

## Claim-by-claim verdict

| # | Claim | Verdict | Nearest prior work |
|---|---|---|---|
| 1 | Additive/interaction split of a **chemical** perturbation atlas | **Novel**, but numerically anticipated | TRADE (genetic, 4 lines, 56/44); sci-Plex (3 lines, 48%) |
| 2 | The interaction is a **pair** property — reproduces across replicates but not across drugs within a line | **Novel**, and contradicts a published assertion | Lim & Pavlidis (*Sci Rep* 2021) assert a line-level "cell-line responsiveness" factor; Subramanian (*Cell* 2017) report only 26% of compounds give similar signatures panel-wide |
| 3 | Dose replication as an internal positive control for residual reproducibility | **Novel** | none found; the factorial vocabulary is absent from this field |
| 4 | Plate structure inflates reproducibility ~7× in Tahoe-100M despite plate-matched normalisation | **Novel** | no published Tahoe analysis reports this control |
| 5 | Rank-1 viability model falsified; ≥6 reproducible programmes | **Novel** | McFarland MIX-Seq (*Nat Commun* 2020) uses a rank-1 structure implicitly but never states or tests it |
| 6 | Context-dependence set by **mechanism**, nuclear receptors highest | **Partly anticipated** — see risk below | Niepel (*Nat Commun* 2017): RTK/MAPK/PI3K inhibitors cell-type-specific, chaperone/cell-cycle shared, but only 6 lines and no index; Subramanian 2017 (26% conserved) |
| 7 | Uniform population shift, not pre-existing subpopulation, at 24 h | **Novel for transcription**; do not overstate | Persister literature (Shaffer, Goyal, Iyer 2025: 70–90% of survive/die set by pre-existing state) concerns **survival over days**, a different question |
| 8 | Additive baseline strengthens with the number of contexts averaged | **Novel** | none found; it explains the standing disagreement between Nature Methods 2025 / DrEval and individual method papers |
| 9 | Post-hoc adaptation fails, fine-tuning succeeds; ~20 arbitrary compounds | **Novel** | von Kügelgen (arXiv 2504.18522) formalises the additive-latent-shift assumption and calls it "a modelling assumption requiring validation" — we validate it empirically |
| 10 | Context count, not cell count, is the binding constraint | **Novel as stated** | Jones, Tsherniak & McFarland (bioRxiv 2020) note limited power "to detect subtle cell-line-specific transcriptional responses" and call for a larger panel |

---

## Two live risks I could not fully verify

1. **Svensson et al., "Back to basics: Observed statistics are sufficient to
   predict drug responses"** (bioRxiv, June 2026; doi 10.64898/2026.06.09.731197).
   A Tahoe-100M author builds our exact baseline ladder — global mean, context
   mean, perturbation mean, additive — on this dataset. The sweep found no
   residual analysis or additive/interaction split in it, which is our opening,
   but this is the same team on the same data and the preprint should be re-read
   in full before submission. *Direct fetch was rate-limited; verdict rests on
   the sweep.*

2. **Shen et al. (Research Square 2026, rs-10448056)** reportedly defines a
   "context-robustness" axis on Tahoe-100M. If that is a per-drug
   context-specificity score, it overlaps claim 6 — our mechanism ranking —
   though probably not the decomposition. *Direct fetch returned HTTP 403; this
   must be checked manually before submission.* This is the most consequential
   unresolved item in this audit.

---

## Honest summary

Nothing we report is a duplicate of a published result, but two of our headline
numbers have close published cousins. The 59/41 split is numerically very near
TRADE's 56/44, and a reviewer who knows that paper will ask what is new; the
answer is the modality (chemical, 379 compounds), the scale (47 contexts versus
4), the method (decomposition with replicate-based validation rather than
profile correlation), and above all claim 2 — that the non-shared part is a
*pair* property — which TRADE cannot test at four lines and which contradicts the
line-level "responsiveness" factor asserted elsewhere.

The most defensible novelty, in order: the pair-property result (2), the
benchmark-design result (8), the adaptation dichotomy (9), the plate confound
(4), and the rank-1 falsification (5). The mechanism ranking (6) is the most
attractive result but carries the largest unresolved scooping risk and rests on
a single atlas.
