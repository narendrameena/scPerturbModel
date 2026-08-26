#!/usr/bin/env python3
"""Phase 3: chemical features for every drug from canonical SMILES.

ECFP4 (Morgan radius-2) 1024-bit fingerprints -> data/processed/drug_ecfp.npz
(keys: drugs [str array], fp [n_drugs x 1024 float32]). DMSO_TF and unparseable
SMILES are excluded (controls are not modeled as perturbations).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "processed"

dr = pd.read_parquet(ROOT / "data/metadata/metadata/drug_metadata.parquet")
gen = AllChem.GetMorganGenerator(radius=2, fpSize=1024)

drugs, fps, failed = [], [], []
for r in dr.itertuples():
    mol = Chem.MolFromSmiles(r.canonical_smiles) if r.canonical_smiles else None
    if mol is None:
        failed.append(r.drug)
        continue
    fp = gen.GetFingerprintAsNumPy(mol).astype(np.float32)
    drugs.append(r.drug)
    fps.append(fp)

np.savez_compressed(OUT / "drug_ecfp.npz",
                    drugs=np.array(drugs), fp=np.stack(fps))
print(f"{len(drugs)} drugs fingerprinted -> {OUT/'drug_ecfp.npz'}")
if failed:
    print(f"failed SMILES ({len(failed)}): {failed[:10]}")
