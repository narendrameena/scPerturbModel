# Related work: single-cell perturbation-response models (surveyed 2026-08-27)

Web survey of the field, with our own dev-subset results positioned against it.

## 1. Latent-space arithmetic (first generation)
- **scGen** (Lotfollahi 2019): VAE + additive perturbation vector in latent space.
- Our finding: the same idea via scVI latent shifts scored r_de100 0.62 vs 0.90
  for the plain additive baseline — latent arithmetic through a bottleneck loses
  drug-specific signal (consistent with why the field moved on).

## 2. Conditional / disentanglement models
- **CPA** (Lotfollahi 2023): compositional perturbation autoencoder — additive
  latent composition of drug/dose/covariate embeddings; adversarial disentangling.
- **chemCPA** (Hetzel 2022): drug-structure encoder → unseen-compound prediction.
- **biolord** (Piran, [Nat Biotech 2024](https://www.nature.com/articles/s41587-023-02079-x)):
  disentanglement of known/unknown attributes; reported to beat chemCPA
  (r² 0.76 vs 0.51 on unseen-drug tasks).
- **PRnet** (Qi, [Nat Comms 2024](https://www.nature.com/articles/s41467-024-53457-1)):
  conditional encoder–decoder for novel chemical perturbations; evaluated on
  Drug_unseen / Cell_line_unseen / Both_unseen splits (same split philosophy as
  our hard splits); reported to outperform chemCPA.
- Our model 1 (context-residual) is in this family but *anchored to the additive
  baseline by construction* (zero-init residual), which guarantees >= baseline.

## 3. Graph / GRN-informed
- **GEARS** (Roohani 2023): GNN over GO/co-expression graphs; combinatorial
  *genetic* perturbations, unseen-gene extrapolation.
- **CellOracle**, **scLAMBDA**: GRN-based TF-perturbation simulation.

## 4. Optimal transport, diffusion, flow matching (2023–2026 wave)
- **CellOT / CondOT**: learn control→perturbed transport maps.
- **Multi-conditional diffusion transformer** (Hu,
  [Quant Biol 2026](https://doi.org/10.1002/qub2.70016)): beats PRnet on all
  three unseen splits.
- **scDFM** ([arXiv 2026](https://arxiv.org/pdf/2602.07103)) distributional flow
  matching; **PerturbDiff** ([arXiv 2026](https://arxiv.org/pdf/2602.19685));
  **PRESCRIBE** ([arXiv 2025](https://arxiv.org/pdf/2510.07964)) Bayesian
  uncertainty; RAG-style retrieval prediction
  ([arXiv 2026](https://arxiv.org/pdf/2603.07233)).

## 5. Foundation models / virtual cell models — including models trained ON Tahoe-100M
- **STATE** (Arc Institute,
  [announcement](https://arcinstitute.org/news/virtual-cell-model-state)):
  multi-scale set transformer over cell *populations*; trained on >100M cells
  including **Tahoe-100M**, Parse-PBMC, Replogle-Nadig. Arc's first virtual
  cell model ([Virtual Cell Initiative](https://arcinstitute.org/virtual-cell-initiative)).
- **Tahoe-x1** (Tahoe Therapeutics,
  [bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.10.23.683759v1),
  [code](https://github.com/tahoebio/tahoe-x1)): perturbation-trained single-cell
  foundation models up to **3B parameters**, trained on Tahoe-100M; emphasizes
  gene essentiality and pathway-informative gene representations.
- Earlier general FMs repurposed for perturbation: **scGPT**, **Geneformer**,
  **scFoundation**, **scBERT**, **UCE**; simple scalable virtual cell:
  **OCOO-T** ([arXiv 2026](https://arxiv.org/pdf/2606.12838)).
- Ecosystem: Tahoe + Arc + Biohub announced a joint effort to generate the
  largest perturbation dataset for virtual cell models
  ([Arc news](https://arcinstitute.org/news/tahoe-arc-biohub)).

## 6. Benchmarks and critiques — the field's reality check
- **"Deep-learning-based gene perturbation effect prediction does not yet
  outperform simple linear baselines"**
  ([Nature Methods 2025](https://www.nature.com/articles/s41592-025-02772-6)):
  GEARS/scGPT/foundation models vs additive & mean-response baselines — deep
  models fail to beat them on genetic-perturbation benchmarks.
- **Benchmarking generalizable perturbation prediction**
  ([Nature Methods 2025/26](https://www.nature.com/articles/s41592-025-02980-0)):
  27 methods × 29 datasets across test scenarios.
- **PerturBench** (NeurIPS): standardized perturbation-model benchmark.
  Model library: [Benchmarking-Single-Cell-Perturbation](https://github.com/xianglin226/Benchmarking-Single-Cell-Perturbation).
- **scPerturb** (Peidli 2024): harmonized datasets + E-distance metrics
  (we use both).

## Where our results sit
1. Our **additive-baseline dominance** (r_de100 0.78–0.90 seen-condition) and
   **scVI-latent-shift failure** independently replicate the field's central
   critique — on *chemical* perturbations in Tahoe-100M, at single-lab scale.
2. Our **context-residual model** is a direct answer to that critique:
   architecture cannot fall below the baseline, and its gain (88% of held-out
   conditions) is attributable to context by ablation.
3. Our **unseen-line triple negative** (learned-embedding n/a, mutation/organ =
   additive, control-expression = additive at 46 training lines) matches why
   STATE/Tahoe-x1 push context breadth: cross-context transfer is the field's
   open problem, and static or baseline-expression descriptors do not solve it
   at ~50-line scale.
4. Obvious next reads/comparisons for us: STATE's population-level architecture
   (our CVAE evaluates population match with E-distance already), Tahoe-x1's
   evaluation protocol on the very same dataset, and PRnet/biolord's
   unseen-drug splits vs our ECFP prior-dropout model.
