"""Push the PETai activity label release to the HuggingFace Hub.

Publishes the four tables of `combined/` as four named configs of a single
dataset repo, uploads `sequences.fasta` and `summary.json` as plain files, and
writes the dataset card from `DATASET_CARD.md`.

The tables are small (13,092 label rows, 1,219 sequences), so everything is
loaded in memory and pushed in one pass. Expect well under a minute.

Run `--dry-run` first. It builds every config and prints the schema and row
counts without contacting the Hub.

Usage:
    python push_to_hf.py <combined_dir> --repo petadex/petai-activity-labels --dry-run
    python push_to_hf.py <combined_dir> --repo petadex/petai-activity-labels
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Columns that must stay typed rather than being inferred per-shard. `value`
# spans several units and is read as float; the condition columns are sparse
# and would otherwise flip between float and object across configs.
LABEL_DTYPES = {
    "seq_md5": "string",
    "source": "string",
    "task": "string",
    "variable": "string",
    "unit": "string",
    "substrate": "string",
    "value": "float64",
    "is_negative": "boolean",
    "ph": "float64",
    "temperature_c": "float64",
    "time_h": "float64",
    "ref": "string",
    "curator": "string",
}

CONFIGS = {
    "labels": "labels.csv",
    "sequences": "sequences.csv",
    "split": "split.csv",
    "clusters": "clusters.csv",
}

EXTRA_FILES = ["sequences.fasta", "summary.json"]


def load_config(combined: Path, name: str, filename: str) -> pd.DataFrame:
    df = pd.read_csv(combined / filename, dtype="string", keep_default_na=False)
    df = df.replace({"": None})
    if name == "labels":
        for col, dtype in LABEL_DTYPES.items():
            if col not in df.columns:
                continue
            if dtype == "boolean":
                df[col] = df[col].str.strip().str.lower().isin({"true", "1", "yes", "t"})
            elif dtype == "float64":
                df[col] = pd.to_numeric(df[col], errors="coerce")
    if name == "sequences":
        for col in ("length", "n_labels", "n_negative", "n_refs"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        if "is_fragment" in df.columns:
            df["is_fragment"] = (
                df["is_fragment"].str.strip().str.lower().isin({"true", "1", "yes", "t"})
            )
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("combined", type=Path, help="the release directory")
    ap.add_argument("--repo", required=True, help="e.g. petadex/petai-activity-labels")
    ap.add_argument("--private", action="store_true", help="create the repo private")
    ap.add_argument("--dry-run", action="store_true", help="build but do not push")
    ap.add_argument(
        "--card",
        type=Path,
        default=Path(__file__).parent / "DATASET_CARD.md",
        help="dataset card to upload as README.md",
    )
    args = ap.parse_args()

    frames = {n: load_config(args.combined, n, f) for n, f in CONFIGS.items()}
    for name, df in frames.items():
        print(f"[{name}] {len(df):,} rows x {len(df.columns)} cols")
        print("   " + ", ".join(f"{c}:{df[c].dtype}" for c in df.columns))

    if args.dry_run:
        print("\ndry run, nothing pushed")
        return 0

    try:
        from datasets import Dataset
        from huggingface_hub import HfApi
    except ImportError as exc:  # pragma: no cover
        print(f"missing dependency: {exc}", file=sys.stderr)
        print("needs `datasets` and `huggingface_hub`", file=sys.stderr)
        return 1

    for name, df in frames.items():
        print(f"pushing config {name} ...")
        Dataset.from_pandas(df, preserve_index=False).push_to_hub(
            args.repo, config_name=name, private=args.private
        )

    api = HfApi()
    for filename in EXTRA_FILES:
        src = args.combined / filename
        if not src.exists():
            continue
        print(f"uploading {filename} ...")
        api.upload_file(
            path_or_fileobj=str(src),
            path_in_repo=f"raw/{filename}",
            repo_id=args.repo,
            repo_type="dataset",
        )

    if args.card.exists():
        print("uploading dataset card ...")
        api.upload_file(
            path_or_fileobj=str(args.card),
            path_in_repo="README.md",
            repo_id=args.repo,
            repo_type="dataset",
        )

    print(f"\ndone: https://huggingface.co/datasets/{args.repo}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
