#!/usr/bin/env python3
"""Download Tabula Sapiens per-tissue h5ads (CELLxGENE) for the cell-type audit.

Independent reference atlas: normal human tissues, different lab/technology
from Tahoe-100M. Only tissues matching Tahoe organ annotations are fetched
(~15 GB total) into data/external/tabula_sapiens/.

Run on a login node (compute nodes may lack internet). Resumable: existing
complete files are skipped.
"""
import json
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "external" / "tabula_sapiens"
OUT.mkdir(parents=True, exist_ok=True)
CID = "e5f58829-1a66-40b5-a624-9046778e74f5"  # Tabula Sapiens collection

# TS tissue dataset title suffix -> Tahoe Organ label(s) it audits
TISSUES = {
    "Lung": ["Lung"],
    "Large_Intestine": ["Bowel"],
    "Pancreas": ["Pancreas"],
    "Skin": ["Skin"],
    "Mammary": ["Breast"],
    "Stomach": ["Esophagus/Stomach"],
    "Uterus": ["Uterus", "Cervix"],
    "Liver": ["Liver"],
    "Bladder": ["Bladder/Urinary Tract"],
    "Kidney": ["Kidney"],
    "Ovary": ["Ovary/Fallopian Tube"],
    "Neural": ["CNS/Brain", "Peripheral Nervous System"],
}


def main():
    url = f"https://api.cellxgene.cziscience.com/curation/v1/collections/{CID}"
    with urllib.request.urlopen(url, timeout=120) as r:
        col = json.load(r)

    mapping = {}
    for d in col["datasets"]:
        title = d.get("title", "")
        for suffix in TISSUES:
            if title == f"Tabula Sapiens - {suffix}":
                h5 = [a for a in d["assets"] if a["filetype"] == "H5AD"]
                if h5:
                    mapping[suffix] = (h5[0]["url"], h5[0].get("filesize", 0))

    manifest = {}
    for suffix, (u, size) in sorted(mapping.items()):
        dest = OUT / f"{suffix}.h5ad"
        manifest[suffix] = {"file": dest.name, "organs": TISSUES[suffix],
                            "bytes": size}
        if dest.exists() and dest.stat().st_size == size:
            print(f"skip {suffix} (complete)")
            continue
        print(f"downloading {suffix} ({size/1e9:.2f} GB) ...", flush=True)
        tmp = dest.with_suffix(".part")
        with urllib.request.urlopen(u, timeout=300) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh, length=1 << 22)
        tmp.rename(dest)
    with open(OUT / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"done: {len(mapping)} tissues -> {OUT}")


if __name__ == "__main__":
    main()
