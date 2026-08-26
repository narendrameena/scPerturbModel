# Project roadmap — perturbation response model on Tahoe-100M

Purpose: build a model, based on single-cell perturbation data, that predicts how a
cell's transcriptomic state changes under a (drug, dose) perturbation in a given
cellular context — and evaluate it honestly against strong baselines.

## Phase 0 — Environment & data access
- [ ] `conda env create -f envs/environment.yml`
- [ ] `scripts/download_tahoe100m.py --mode list` then `--mode metadata`
- [ ] Explore metadata tables in `notebooks/01_metadata_exploration.ipynb`
      (drugs, doses, MOAs, targets, cell lines, driver mutations, cell counts/condition)
- [ ] Decide storage format: streamed parquet vs. converted AnnData/zarr shards on scratch

## Phase 1 — Subset-first development loop
Working at 100M-cell scale from day one is a mistake. Define a small, fixed dev subset:
- [ ] e.g. 5–10 cell lines × ~50 drugs (spanning MOAs + E-distance range) × all doses,
      plus matched DMSO controls, capped ~2–5M cells → `data/interim/dev_subset/`
- [ ] QC mirroring the paper (≥700 UMIs, <20% mito, ≥250 genes, singlets, ≥50 cells/condition)
- [ ] Reproduce paper sanity checks on the subset (E-distance dose trend, Dabrafenib/BRAF specificity)

## Phase 2 — Baselines (must-beat)
- [ ] Non-parametric: predicted perturbed state = control cells + mean pseudobulk delta
      of that (drug, cell line) — the "additive shift" baseline
- [ ] Linear: ridge regression on pseudobulk deltas with drug + cell line + dose factors
- [ ] scVI latent-space shift model (train scVI, model perturbation as latent arithmetic)
- [ ] Published methods for comparison where feasible: CPA, scGen/chemCPA, GEARS-style
      (drug-structure-aware), optionally a pretrained foundation model (scGPT / STATE / scFoundation)

## Phase 3 — The model
Candidate directions (pick after Phase 2 results):
1. **Conditional generative model**: scVI-style decoder conditioned on (cell line, drug
   embedding, dose); drug embedding from chemical structure (ECFP/ChemBERTa) so unseen
   drugs are addressable.
2. **Latent transport**: learn a dose-conditioned map control→perturbed in latent space
   (optimal transport / flow matching), context-conditioned.
3. **Transformer over gene tokens** fine-tuned on perturbation prediction (foundation-model
   route) — expensive; justify with Phase 2 evidence.
Key design requirements:
- dose–response must be explicit (monotone-ish embedding or continuous conditioning)
- context (cell line / mutation status) as input, so context-dependence is testable
- distribution-level prediction (not just mean shift) — the atlas exists because
  heterogeneity matters (cell cycle substates, resistant subpopulations)

## Phase 4 — Evaluation protocol (fix BEFORE training the big model)
Splits, hardest→easiest generalization:
- [ ] **Held-out condition**: seen drug + seen line, unseen combination
- [ ] **Held-out cell line** (context generalization)
- [ ] **Held-out drug** (chemistry generalization; needs structure-aware drug encoder)
- [ ] **Held-out plate 14** (replication/technical robustness — paper does this)
Metrics:
- [ ] E-distance realized vs predicted; energy/MMD distance between predicted and true
      cell populations; pseudobulk delta correlation (all genes + top-k DE genes)
- [ ] DE-gene overlap (precision/recall of top DE genes vs plate-matched DMSO)
- [ ] Biological probes: BRAF/KRAS context specificity, CDK inhibitor cell-cycle shifts
- [ ] Report against baselines; a model that doesn't beat the additive-shift baseline
      on unseen combinations isn't learning context-dependence

## Phase 5 — Scale up & applications
- [ ] Scale training to full 95.6M-cell filtered atlas (multi-GPU; jobs/ SLURM scripts)
- [ ] Transfer tests on external data (`data/external/`): sciPlex3, Replogle Perturb-seq
      via scPerturb; DepMap sensitivity correlation
- [ ] Possible aging angle (project context): map drug-induced transcriptomic shifts
      onto aging signatures; look for perturbations that move cells along/against them
- [ ] Write-up: figures in results/figures, tables in results/tables

## Practical notes
- Keep every experiment driven by a YAML in `configs/` + git commit hash; log to
  `results/` and `logs/` with the config name.
- Heavy jobs through `jobs/*.sbatch`; never run multi-hour training on login nodes.
- Track data lineage: raw → interim → processed transformations only via scripts in
  `scripts/` (no untracked notebook-side effects).
