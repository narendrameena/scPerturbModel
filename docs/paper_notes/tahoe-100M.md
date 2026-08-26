# Paper notes — Tahoe-100M (Zhang et al., bioRxiv 2025.02.20.639398v3)

**Title:** Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent
Gene Function and Cellular Modeling
**Groups:** Tahoe Therapeutics + Parse Biosciences. Corresponding: Johnny Yu, Hani Goodarzi.
**Data:** Hugging Face (`tahoebio/Tahoe-100M`).

## One-paragraph summary

The largest public single-cell perturbation dataset to date: ~100 million scRNA-seq
profiles measuring how ~1,100 small-molecule drug–dose treatments (379 distinct drugs,
mostly approved oncology agents, 3 doses each) reshape the transcriptomes of 50 cancer
cell lines after 24 h. The "Mosaic" platform grows all 50 lines together as mixed
spheroid "cell villages" in each well, so every drug is applied to all lines
simultaneously; cells are assigned back to their line of origin by SNP demultiplexing.
This design massively reduces batch effects and yields 52,886 cell line × drug × dose
conditions with a median of ~1,287 cells each — explicitly built as training data for
AI/foundation models of the cell.

## Experimental design (Mosaic platform)

1. 50 commercially available cancer cell lines pooled at equal ratios into 3D spheroids
   ("Mosaic tumors") with Cultrex, seeded in 14 deep 96-well plates.
2. Each well = one drug treatment (compounds from MedChemExpress, ≤1% DMSO); wells
   H11/H12 on every plate = DMSO vehicle controls. 24 h exposure.
3. Dissociation → fixation (Parse Evercode v3) → Parse **GigaLab** combinatorial
   barcoding (3 rounds) → ~100k-cell sublibraries (1,786 total) → Ultima UG-100
   sequencing (~1.4 trillion reads).
4. Cell line identity recovered by **demuxlet** against a curated SNP reference
   (dbSNP, exonic/UTR, gnomAD AF > 0.1; >98% accuracy, >90% singlets).
   Treatment identity from the first barcode (well position).

## Scale / QC numbers

- 153M cells sequenced → **100.6M pass minimal filters** → **95.6M pass full filters**
  (used downstream).
- Full filters: ≥700 UMIs, <20% mito, |z| ≤ 3 on UMI and mito%, ≥250 genes,
  demuxlet singlet. Condition-level: ≥50 cells per cell line × drug condition.
- 3 low-abundance lines dropped (NCI-H661, NCI-H596, NCI-H2122) → 47 lines analyzed.
- Mean ~2,288 transcripts/cell (median 1,890). Genes: Ensembl 109 / GRCh38; 62,710 genes.
- 47 lines from 13 organs (mostly lung, bowel, pancreas, skin); TP53/KRAS/CDKN2A
  altered in ~half. 180/379 drugs classified into 25 MOAs; 325 target genes.
- Scale comparison: 31× more drug-perturbed conditions than sciPlex3; 29× more
  observations per condition than Replogle 2022 Perturb-seq.
- Plate 14 replicates plate 6: matched pseudobulk Pearson r median 0.975 vs 0.915 unmatched.

## Computational approach (their pipeline — worth mirroring)

- **scVI** (scvi-tools), 10-d latent, negative binomial likelihood, trained on plates
  1–13 (89.4M cells; 90/10 train/val, 10 epochs ≈ >800M examples). Plate 14 held out
  for validation/criticism. Model + AnnData "minified" for interactive use (41 GB).
- **LISI** to rank drivers of variation: cell line identity ≫ cell cycle ≫ drug/dose.
  → drug effects need differential analysis, not raw clustering.
- **E-distance** (Peidli/scPerturb) in scVI latent space, per cell line vs plate-matched
  DMSO_TF controls; median across plates then across lines. Dose-dependent; largest for
  proteasome/HDAC/PI3K-AKT inhibitors and outliers (homo)harringtonine, dinaciclib.
- **Differential gene-set scores**: pseudobulk per (cell line, drug, dose); Vision scores
  over 7,390 MSigDB gene sets on scVI-normalized expression; 250 Monte Carlo samples;
  ~62k comparisons → PCA → tSNE (cuML) and neural-network MDS. E-distance gradient
  center→periphery; MOAs partially separate.
- **Cell cycle**: Scanpy score_genes_cell_cycle (Regev lab signature), per-plate
  log2 odds ratio vs control.

## Key biological validations (sanity checks for any model)

- **Dabrafenib** (BRAF-V600E inhibitor) suppresses a RAS/RAF activity signature
  (KRAS.600.LUNG.BREAST.V1) specifically in BRAF-V600E lines, not KRAS-mutant lines.
- **RMC-6236** (pan-RAS inhibitor) hits KRAS-mutant lines most; separates from
  **Adagrasib** (KRAS-G12C-specific) which only affects G12C lines.
- **Palbociclib** (CDK4/6) → G1 arrest; **Dinaciclib** (CDK1/2/5/9) → G2/M arrest.
- HDAC inhibitors Belinostat/Panobinostat → G2/M arrest; Tucidinostat → G0/G1
  (negative G2/M enrichment); Carbamazepine (weak HDAC activity) → no effect.
- Microtubule inhibitors → G2/M arrest (except Paclitaxel, Tubulin Inhibitor 6 minimal).

These are ideal **held-out qualitative benchmarks**: a good perturbation model should
reproduce mutation-context-dependent responses and MOA-specific cell cycle shifts.

## Why it matters for modeling (their framing)

1. Diverse, balanced training data across contexts/MOAs → context-dependent gene function.
2. Multiple drugs per pathway/target → fine-grained MOA representations.
3. Cell-village design unifies contexts in one batch structure → genetic/epigenetic
   context modulation is learnable rather than confounded.
Supports self-/semi-supervised and transfer learning; explicit "virtual cell" ambitions
(cites Bunne et al. 2024).

## Caveats / things to watch

- Cancer cell lines only, 24 h single timepoint, small-molecule perturbations only.
- Drug/MOA/target annotations were generated with **GPT-4o + MedChemExpress scraping**
  (only 97/380 manually inspected) — treat MOA labels as noisy; consider re-curating
  (e.g., against ChEMBL/DrugBank) before using as supervision.
- Doses vary per drug (typically 0.05 / 0.5 / 5 µM); dose is a covariate, not uniform.
- Demuxlet doublet/ambiguity filtering removes low-coverage cells non-randomly.
- Preprint (not peer-reviewed at this version); "no reuse" copyright note on the PDF,
  but the dataset itself is publicly released on Hugging Face.
