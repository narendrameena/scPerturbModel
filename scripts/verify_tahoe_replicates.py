#!/usr/bin/env python3
"""Verify the 4.7% replicate figure directly from Tahoe's released metadata.

The claim that Tahoe-100M replicates only 4.7% of its (cell line, drug, dose)
combinations is the most consequential and most contestable statement in this
work, and until now it was computed from OUR pseudobulk build. If that build
dropped conditions -- through cell-count thresholds, plate exclusions or the
delta construction -- the figure would be an artefact of our processing rather
than a property of the atlas, and it would be dismissed on exactly that ground.

This script recomputes it from `obs_metadata.parquet` as released: 100.6 million
cell records with their plate, cell line and drug-and-concentration annotation,
with no filtering of ours applied. It reports the count three ways, because the
honest number depends on which cells one is willing to use:

  all cells               every record in the released metadata
  pass_filter == 'full'   the atlas's own recommended QC filter
  >= 10 cells per plate   additionally requiring enough cells for a usable
                          pseudobulk profile on each plate, which is the
                          condition under which a replicate is actually usable

If the three agree, the figure is a property of the experimental design. If they
diverge, the difference is itself the finding and belongs in the paper.

Outputs: results/tables/tahoe_replicate_structure.csv
"""
import argparse
import ast
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata" / "obs_metadata.parquet"
TAB = ROOT / "results" / "tables"


def parse_conc(s):
    """drugname_drugconc is a stringified list of (name, conc, unit) tuples."""
    try:
        v = ast.literal_eval(s) if isinstance(s, str) else s
        if isinstance(v, (list, tuple)) and len(v):
            return float(v[0][1])
    except Exception:
        pass
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cells", type=int, default=10)
    args = ap.parse_args()

    cols = ["plate", "drug", "cell_line", "drugname_drugconc", "pass_filter"]
    print(f"reading {META.name} ...", flush=True)
    d = pd.read_parquet(META, columns=cols)
    print(f"  {len(d):,} cell records", flush=True)

    # concentration is inside the annotation string; parse it once on the
    # unique values rather than per cell
    uniq = d.drugname_drugconc.astype(str).drop_duplicates()
    cmap = {s: parse_conc(s) for s in uniq}
    d["conc"] = d.drugname_drugconc.astype(str).map(cmap)
    print(f"  {d.conc.notna().mean():.1%} of records have a parsed "
          f"concentration; {d.conc.nunique()} distinct values", flush=True)

    rows = []
    for label, sub in (("all cells", d),
                       ("pass_filter == full",
                        d[d.pass_filter.astype(str) == "full"])):
        g = sub.groupby(["cell_line", "drug", "conc"], observed=True)
        n_plates = g.plate.nunique()
        n_tri = len(n_plates)
        n_rep = int((n_plates > 1).sum())
        pair = sub.groupby(["cell_line", "drug"], observed=True).plate.nunique()
        rows.append({"filter": label, "n_triples": n_tri,
                     "n_replicated": n_rep, "frac_replicated": n_rep / n_tri,
                     "n_pairs": len(pair),
                     "frac_pairs_multiplate": float((pair > 1).mean())})
        print(f"  {label:24s} {n_tri:,} (line, drug, dose) triples, "
              f"{n_rep:,} on >1 plate = {n_rep/n_tri:.2%}", flush=True)

    # a replicate is only usable if each plate carries enough cells to build a
    # profile; this is the number that matters for the estimator
    sub = d[d.pass_filter.astype(str) == "full"]
    cnt = (sub.groupby(["cell_line", "drug", "conc", "plate"], observed=True)
           .size().rename("n_cells").reset_index())
    ok = cnt[cnt.n_cells >= args.min_cells]
    per = ok.groupby(["cell_line", "drug", "conc"], observed=True).plate.nunique()
    rows.append({"filter": f">= {args.min_cells} cells per plate",
                 "n_triples": len(per), "n_replicated": int((per > 1).sum()),
                 "frac_replicated": float((per > 1).mean()),
                 "n_pairs": np.nan, "frac_pairs_multiplate": np.nan})
    print(f"  >= {args.min_cells} cells/plate      {len(per):,} triples, "
          f"{int((per > 1).sum()):,} on >1 plate = {(per > 1).mean():.2%}")

    R = pd.DataFrame(rows)
    R.to_csv(TAB / "tahoe_replicate_structure.csv", index=False)
    print("\n" + R.round(4).to_string(index=False))

    f = R.frac_replicated
    print(f"\nrange across filters: {f.min():.2%} to {f.max():.2%}")
    if f.max() - f.min() < 0.02:
        print("  The three agree, so the figure is a property of the "
              "experimental design\n  and not of any processing choice of "
              "ours.")
    else:
        print("  The filters disagree; the spread is itself the result and "
              "should be reported\n  as a range rather than a point estimate.")
    print(f"\nfor contrast, (line, drug) pairs spanning >1 plate: "
          f"{R.frac_pairs_multiplate.dropna().iloc[0]:.1%} — high only because "
          f"different\nDOSES sit on different plates, which is the confusion "
          f"the paper is about.")


if __name__ == "__main__":
    main()
