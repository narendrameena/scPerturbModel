# Methodological positioning: what we do differently, and what it buys

An honest audit of our methods against the published literature, written to be
argued with. Section 3 lists where we are *not* stronger, because that is where
review will concentrate.

---

## 1. Methodological choices that are genuinely distinctive

### 1.1 Dose replication as an internal positive control
Each (line, drug) pair appears at three doses. We use agreement *between doses
of the same pair* to separate reproducible line-specific signal from noise:
independent noise does not covary across doses, so the covariance estimates the
real interaction.

*Why it matters.* Every other treatment of this problem reports how large a
residual is; none establishes that the residual is **real**. Without this, "the
residual is small" and "the residual is noise" are indistinguishable — and they
lead to opposite conclusions. Neither MAP, XPert, CPA, chemCPA, PRnet, STATE nor
the Tahoe-100M paper itself uses replicate structure this way.

### 1.2 Cross-plate stratification (and the confound it exposed)
All reproducibility estimates are restricted to comparisons made on *different
plates*. We found this is not a nicety: within-plate agreement is **+0.45–0.48**
versus **+0.06–0.07** across plates — a sevenfold inflation that survives
plate-matched DMSO normalisation, the standard correction.

*Why it matters.* Any analysis of this atlas that pools within- and cross-plate
comparisons overstates reproducibility by ~7×. We know of no published analysis
of Tahoe-100M that reports this control. It also changed our own headline
(interaction share 41% → 28% at dev47 scale), which is why it is stated as a
correction in `RESULTS.md` rather than silently applied.

### 1.3 Split-half validation of latent components
Components are fitted on one half of a line×drug pair's plates and scored on the
other half. Only components that reproduce are interpreted.

*Why it matters.* Papers reporting programs from NMF/PCA/factor models
(expiMap, VEGA, OntoVAE, cNMF applications) typically report variance explained,
which is guaranteed to be positive and says nothing about reproducibility. Our
procedure demonstrably discriminates: components 1, 2, 5, 6, 7, 10 reproduce
(r = 0.33–0.48) while 3, 4 and 8 come back **negative**. Reporting the failures
is the point — it shows the filter has teeth.

### 1.4 A model architecturally nested inside its own baseline
`prediction = additive_prior + residual`, with the residual **zero-initialised**.
The model therefore *starts* at the additive baseline and can only depart from it
where data supports that. The no-context ablation lands on the baseline exactly.

*Why it matters.* The standard practice — train a deep model, compare to a
baseline afterwards — cannot distinguish "the model learned context" from "the
model learned the drug main effect slightly differently". Given that
[Nature Methods 2025](https://www.nature.com/articles/s41592-025-02772-6) and
DrEval both find deep models failing to beat mean baselines, nesting the
baseline inside the architecture converts an empirical comparison into a
structural guarantee.

### 1.5 Isolating architectural from informational failure
For an unseen line we ran the *same probe data* through two adaptation regimes:
fitting only the line embedding with the residual head frozen, versus fine-tuning
the whole model. Frozen: **zero** gain at every k, including an oracle using all
~30 available probe drugs (n = 8,420 conditions). Fine-tuned: **+0.034**,
recovering 98% of the ceiling.

*Why it matters.* The field reports that unseen contexts are hard (XPert: 121%
performance drop in cold-cell; scPerturBench; Systema). None separates *the
information is absent* from *the architecture cannot use it*. Holding data
constant and varying only whether weights may adapt answers that directly, and
the answer — architectural, not informational — is actionable.

### 1.6 Varying one design parameter of the benchmark itself
We recomputed the additive baseline while varying **only** the number of lines
averaged: r 0.649 (5 lines) → 0.718 (45), **+10.6%**.

*Why it matters.* This reconciles a real contradiction. MAP reports +12.3% for
its model over the best baseline using six Tahoe lines; Nature Methods 2025,
DrEval, XPert's own mean baselines and this work (40+ contexts) find mean
baselines nearly unbeatable. The handicap a sparse benchmark imposes on the
baseline is the same order as the gain attributed to the model. We are not aware
of anyone testing baseline strength as a function of benchmark design.

### 1.7 Nulls matched to the confound, not to nothing
Expression-matched random gene sets for target-abundance tests; permuted-pair
nulls for residual correlations; within-cancer-type z-scoring before pooling
TCGA; sequencing depth as an explicit control in the epigenome test; the
responsive-gene universe (not the genome) as enrichment background.

*Why it matters.* Each of these, omitted, produces a false positive we would
have reported. The genome-background version of the enrichment, for instance,
recovers the responsive-gene signature and calls it a property of the component.

### 1.8 Portability enforced by an abstract replicate axis
The decomposition is written against three roles — *context*, *perturbation*,
*independent replicate* — rather than against Tahoe's schema, so it runs
unchanged on other atlases by naming which column plays which role. Applied to
OP3 (6 immune cell types, donors as replicates) and sciPlex3 (3 cancer lines,
experimental replicates), the structure reproduces: interaction share 41% / 33%
/ 30%, residual across replicates +0.062 / +0.135 / +0.048, residual across
perturbations +0.019 / +0.007 / +0.026.

*Why it matters.* The replicate axis is the load-bearing requirement, and making
it explicit is what allows the same estimator to be applied to plates, donors
and experimental repeats without reinterpretation. It also produced a free
internal check: our §2 claim that a prior estimated from fewer contexts is
noisier predicts that the additive share should rise as contexts fall, and the
three datasets — 47, 6 and 3 contexts — land in exactly that order
(59% → 67% → 70%).

---

## 2. What better biology this actually buys

Ranked by how defensible each claim is.

**Strong — well powered, replicated across subsets.**
1. **Context-dependence is a property of the drug's mechanism** (Kruskal–Wallis
   p = 2×10⁻³, 367 drugs): nuclear-receptor agonists most context-specific
   (0.363), glucose-transporter and RAF inhibitors most conserved (0.087–0.119).
   This is prospectively useful — it says how broad a screening panel a compound
   needs, from its mechanism. The top class is mechanistically expected
   (nuclear receptors read out each cell's enhancer landscape), which makes it a
   positive control and a result simultaneously.
2. **The architecture of drug response**: ~59% of reproducible response is the
   drug's average effect, ~41% is line×drug interaction; the interaction
   reproduces across doses (+0.062) but not across drugs (+0.019) or lines
   (−0.002). It is therefore an *interaction*, not a line property — which
   explains, in one statement, why line-level descriptors cannot work.
   **Replicated in two independent atlases** — OP3 (primary human immune cells,
   donors as replicates) and sciPlex3 (Trapnell lab, sci-Plex) — across three
   platforms, three labs and three kinds of replicate. This is now the most
   robust claim we have.
3. **A working protocol for a new cell model**: ~20 *arbitrary* compounds plus
   fine-tuning recovers 98% of the achievable gain; probe-panel design does not
   help (random ties or beats chemically-diverse, MOA-diverse, most-discriminative
   and strongest-effect panels).

**Moderate — real but bounded.**
4. **It is not a viability artefact.** The rank-1 model implicit in MIX-Seq is
   falsified: component 1 carries only 34% of reproducible residual variance and
   is uncorrelated with cell-cycle arrest (r = −0.05). At least six components
   reproduce. Answers the first deflationary question a reviewer asks.
5. **At the cell level it is a uniform population shift, not a pre-existing
   subpopulation** (20.5M cells; median shape excess 0.91× the noise floor).
   Corollary: pseudobulk is adequate for this question — worth saying in a field
   investing heavily in single-cell perturbation atlases.
6. **Named response programs** (epithelial-vs-neuronal identity, ERBB4/PI3K,
   immune-vs-ECM, EMT, ribosome-vs-antimicrobial, chromatin/senescence), though
   the enrichment q-values are modest (10⁻²–10⁻¹).

**Framing — a design insight rather than a discovery.**
7. **The binding constraint is context count, not cell count.** 100M cells, ~47
   usable contexts; four independent line-level predictors fail, plausibly one
   power ceiling hit four times. This reframes what these atlases should optimise.

---

## 3. Where we are *not* stronger — expect review here

- **No new data and no wet-lab validation.** Every Nature Genetics computational
  paper surveyed ships a prospective experimental test (15 designed enhancers,
  5,850 STARR-seq elements, a validated EGFR variant). We ship none. This is the
  single largest gap and it is not closable computationally.
- **We lose on prediction accuracy, and should not compete there.** MAP and
  XPert are better predictors than anything we built. Our model's advantage is
  interpretability of attribution, not performance.
- **Our positives are outnumbered by our negatives.** Four failed line-level
  predictors, a failed survival association, no coherence of the programs in
  tumours. Defensible individually; cumulatively they make a "discovery" framing
  hard to sustain.
- **n = 47 cuts both ways.** It excuses our negatives but also caps any positive
  line-level claim we might later want to make.
- **The MOA annotations are GPT-4o-derived and demonstrably imperfect** — we
  independently caught Clobetasol propionate ("unclear" → corticosteroid) and
  Capmatinib (MET inhibitor labelled EGFR/ERBB). Our mechanism-level result
  inherits that noise, which biases toward the null but weakens precision.
- **Partly single-dataset.** The *decomposition* now replicates in OP3 and
  sciPlex3, so the architectural claim is multi-atlas. But everything downstream
  of it remains Tahoe-only: the drug-mechanism CDI ranking, the six named
  programs, the few-shot protocol, the baseline-scaling curve, the TCGA
  anchoring and the identity audit. Those are the results a reviewer would most
  reasonably ask to see repeated, and OP3's 147 compounds with MOA annotation
  would support at least the mechanism ranking.

---

## 4. The honest one-line summary

We are not a better predictor and we have no new experiments. What we have is a
**measurement of the problem**: what fraction of drug response is transferable,
which drugs transfer, what the non-transferable part is made of, why no cell-line
descriptor recovers it, what it costs to obtain it for a new context, and how
much of the field's disagreement about baselines is benchmark design. Those are
answers the prediction papers assume rather than establish.
