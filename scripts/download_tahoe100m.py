#!/usr/bin/env python3
"""Access the Tahoe-100M dataset on Hugging Face (tahoebio/Tahoe-100M).

Modes:
  list      List all files in the HF dataset repo (see what exists before downloading).
  metadata  Download only small annotation/metadata files into data/metadata/.
  preview   Stream the first N expression records without downloading the dataset.
  full      Snapshot the ENTIRE dataset into data/raw/ (hundreds of GB - run via jobs/).

Examples:
  python scripts/download_tahoe100m.py --mode list
  python scripts/download_tahoe100m.py --mode metadata
  python scripts/download_tahoe100m.py --mode preview -n 5
  python scripts/download_tahoe100m.py --mode full
"""
import argparse
import sys
from pathlib import Path

REPO_ID = "tahoebio/Tahoe-100M"
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Core annotation tables + docs (~2.3 GB, dominated by obs_metadata.parquet).
# Deliberately excludes metadata/pseudobulk_differential_expression/ (~89 GB, 1026
# shards) and data/ (~338 GB, 3388 shards) — those come with --mode full.
def is_core_metadata(f: str) -> bool:
    if f.startswith("metadata/") and "pseudobulk" not in f:
        return True
    return f in ("README.md", "LICENSE.md") or f.startswith("tutorials/")


def repo_files():
    from huggingface_hub import list_repo_files
    return sorted(list_repo_files(REPO_ID, repo_type="dataset"))


def mode_list(_args):
    for f in repo_files():
        print(f)


def mode_metadata(_args):
    from huggingface_hub import hf_hub_download
    out_dir = PROJECT_ROOT / "data" / "metadata"
    out_dir.mkdir(parents=True, exist_ok=True)
    picked = [f for f in repo_files() if is_core_metadata(f)]
    if not picked:
        print("No metadata-like files matched; run --mode list and adjust is_core_metadata().")
        return
    for f in picked:
        print(f"-> {f}")
        hf_hub_download(REPO_ID, f, repo_type="dataset", local_dir=out_dir)
    print(f"\nDone: {len(picked)} files in {out_dir}")


def mode_preview(args):
    from datasets import load_dataset
    ds = load_dataset(REPO_ID, streaming=True, split="train")
    for i, record in enumerate(ds):
        if i >= args.n:
            break
        print(f"--- record {i} ---")
        for k, v in record.items():
            s = str(v)
            print(f"  {k}: {s[:120]}{'...' if len(s) > 120 else ''}")


def mode_full(_args):
    from huggingface_hub import snapshot_download
    out_dir = PROJECT_ROOT / "data" / "raw" / "tahoe-100m"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Snapshotting {REPO_ID} -> {out_dir} (~429 GB total; resumable).")
    snapshot_download(REPO_ID, repo_type="dataset", local_dir=out_dir, max_workers=8)
    print("Done.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["list", "metadata", "preview", "full"],
                    default="list")
    ap.add_argument("-n", type=int, default=3, help="records to show in preview mode")
    args = ap.parse_args()
    {"list": mode_list, "metadata": mode_metadata,
     "preview": mode_preview, "full": mode_full}[args.mode](args)


if __name__ == "__main__":
    sys.exit(main())
