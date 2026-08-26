"""PyTorch datasets for Tahoe-100M.

Verified record schema (HF `tahoebio/Tahoe-100M`, config `expression_data`,
95,624,334 records = the paper's full-filter cells):

  genes              sequence<int64>   token IDs of non-zero genes (map via gene_metadata)
  expressions        sequence<float32> raw counts aligned with `genes`
  drug               string            "DMSO_TF" = vehicle control
  sample             string            join key -> sample_metadata (has drug concentration)
  BARCODE_SUB_LIB_ID string            unique cell ID, index into obs_metadata
  cell_line_id       string            Cellosaurus ID (join -> cell_line_metadata)
  moa-fine           string            fine MOA label (GPT-curated; noisy)
  canonical_smiles   string            drug structure
  pubchem_cid        string            "" for DMSO controls
  plate              string            "1".."14"

Quirk: the FIRST entry of `genes`/`expressions` is a marker/CLS token
(expressions[0] < 0) and must be stripped before use — see strip_marker().

Patterns:
  - development: materialized dev subset in data/interim/ -> map-style Dataset
  - full scale:  local parquet shards in data/raw/tahoe-100m/data/ -> IterableDataset
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

import torch
from torch.utils.data import IterableDataset, get_worker_info

META_FIELDS = (
    "drug", "sample", "BARCODE_SUB_LIB_ID", "cell_line_id",
    "moa-fine", "canonical_smiles", "pubchem_cid", "plate",
)


def strip_marker(genes: list[int], expressions: list[float]) -> tuple[list[int], list[float]]:
    """Drop the leading CLS/marker token (flagged by a negative first count)."""
    if expressions and expressions[0] < 0:
        return genes[1:], expressions[1:]
    return genes, expressions


class TahoeStreamDataset(IterableDataset):
    """Stream cells from local parquet shards (preferred) or the HF hub.

    Yields dicts: {'genes': LongTensor (g,), 'expressions': FloatTensor (g,),
    <metadata str fields>}. Densification to a fixed gene space belongs in the
    collate_fn (gene vocabulary = 62,710 genes from gene_metadata.parquet).
    """

    def __init__(
        self,
        local_dir: str | Path | None = "data/raw/tahoe-100m",
        repo_id: str = "tahoebio/Tahoe-100M",
        split: str = "train",
    ):
        super().__init__()
        self.local_dir = Path(local_dir) if local_dir else None
        self.repo_id = repo_id
        self.split = split

    def _stream(self):
        from datasets import load_dataset

        if self.local_dir is not None and (self.local_dir / "data").exists():
            shards = sorted(str(p) for p in (self.local_dir / "data").glob("train-*.parquet"))
            ds = load_dataset("parquet", data_files={"train": shards},
                              streaming=True, split="train")
        else:
            ds = load_dataset(self.repo_id, streaming=True, split=self.split)
        worker = get_worker_info()
        if worker is not None:  # shard the stream across DataLoader workers
            ds = ds.shard(num_shards=worker.num_workers, index=worker.id)
        return ds

    def __iter__(self) -> Iterator[dict]:
        for record in self._stream():
            genes, expressions = strip_marker(record["genes"], record["expressions"])
            out: dict = {
                "genes": torch.as_tensor(genes, dtype=torch.long),
                "expressions": torch.as_tensor(expressions, dtype=torch.float32),
            }
            for f in META_FIELDS:
                out[f] = record.get(f)
            yield out
