#!/usr/bin/env python3
"""Phase 1: materialize the dev subset from the raw Tahoe-100M shards.

Scans data/raw/tahoe-100m/data/train-*.parquet, keeps cells whose
(cell_line_id, drug) match data/interim/dev_subset/selection.yaml (drugs +
DMSO_TF controls, all doses/plates), joins the drug concentration from
sample_metadata, and writes matching cells to
data/interim/dev_subset/data/train-XXXXX.parquet (same schema + conc columns,
so TahoeStreamDataset(local_dir=data/interim/dev_subset) just works).

Resumable: an input shard whose output part (or empty-marker) exists is skipped.

Usage:
  python scripts/build_dev_subset.py [--workers 32] [--limit N]
"""
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "tahoe-100m" / "data"
META = ROOT / "data" / "metadata" / "metadata"
OUT = ROOT / "data" / "interim" / "dev_subset"   # overridden by --subset-dir
OUT_DATA = OUT / "data"
OUT_MARK = OUT / ".done"  # empty-output markers for resumability


def load_selection():
    sel = yaml.safe_load(open(OUT / "selection.yaml"))
    lines = [c["cellosaurus"] for c in sel["cell_lines"]]
    drugs = list(sel["drugs"]) + list(sel.get("controls", ["DMSO_TF"]))
    return lines, drugs


def dose_lookup() -> pd.DataFrame:
    sm = pd.read_parquet(META / "sample_metadata.parquet")

    def parse(s):
        try:
            entries = ast.literal_eval(s)
            name, conc, unit = entries[0]
            return float(conc), str(unit)
        except Exception:
            return float("nan"), ""

    parsed = sm["drugname_drugconc"].map(parse)
    sm["conc"] = [p[0] for p in parsed]
    sm["conc_unit"] = [p[1] for p in parsed]
    return sm.set_index("sample")[["conc", "conc_unit"]]


def process_shard(args):
    shard, out_path, mark_path, lines, drugs, dose_records = args
    t = pq.read_table(shard)
    mask = pc.and_(pc.is_in(t["cell_line_id"], value_set=pa.array(lines)),
                   pc.is_in(t["drug"], value_set=pa.array(drugs)))
    t = t.filter(mask)
    if t.num_rows == 0:
        mark_path.touch()
        return shard.name, 0
    samples = t["sample"].to_pylist()
    nan = (float("nan"), "")
    conc = [dose_records.get(s, nan)[0] for s in samples]
    unit = [dose_records.get(s, nan)[1] for s in samples]
    t = t.append_column("conc", pa.array(conc, pa.float32()))
    t = t.append_column("conc_unit", pa.array(unit, pa.string()))
    tmp = out_path.with_suffix(".tmp")
    pq.write_table(t, tmp, compression="zstd")
    tmp.rename(out_path)
    return shard.name, t.num_rows


def main():
    global OUT, OUT_DATA, OUT_MARK
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=None, help="debug: only N shards")
    ap.add_argument("--subset-dir", default="data/interim/dev_subset",
                    help="dir containing selection.yaml; outputs go under it")
    args = ap.parse_args()
    OUT = ROOT / args.subset_dir
    OUT_DATA, OUT_MARK = OUT / "data", OUT / ".done"

    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_MARK.mkdir(parents=True, exist_ok=True)
    lines, drugs = load_selection()
    dl = dose_lookup()
    dose_records = {s: (float(r.conc), str(r.conc_unit)) for s, r in dl.iterrows()}

    shards = sorted(RAW.glob("train-*.parquet"))[: args.limit]
    todo = []
    for s in shards:
        out_path = OUT_DATA / s.name
        mark = OUT_MARK / s.name
        if out_path.exists() or mark.exists():
            continue
        todo.append((s, out_path, mark, lines, drugs, dose_records))
    print(f"{len(shards)} shards found, {len(todo)} to process "
          f"({len(shards) - len(todo)} already done)")

    total = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_shard, t) for t in todo]
        for i, f in enumerate(as_completed(futs), 1):
            name, n = f.result()
            total += n
            if i % 100 == 0 or i == len(futs):
                print(f"[{i}/{len(futs)}] +{total:,} cells so far", flush=True)

    # condition-level summary over ALL output shards
    counts = {}
    for p in sorted(OUT_DATA.glob("train-*.parquet")):
        t = pq.read_table(p, columns=["cell_line_id", "drug", "conc", "plate"])
        df = t.to_pandas().groupby(["cell_line_id", "drug", "conc", "plate"],
                                   observed=True).size()
        for k, v in df.items():
            counts[k] = counts.get(k, 0) + v
    summary = pd.Series(counts).rename("n_cells").reset_index()
    summary.columns = ["cell_line_id", "drug", "conc", "plate", "n_cells"]
    summary.to_csv(OUT / "condition_counts.csv", index=False)
    print(f"\nDev subset: {summary.n_cells.sum():,} cells, "
          f"{len(summary):,} (line,drug,conc,plate) groups")
    print(f"-> {OUT_DATA}  |  summary: {OUT / 'condition_counts.csv'}")


if __name__ == "__main__":
    main()
