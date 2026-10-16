"""Push the PETai v2 label release to the HuggingFace Hub.

Issue 117. Adapted from the issue 116 push script. Differences from v1:

- five configs rather than four, adding `normalised`
- `labels` carries `orig_seq_md5` and `remapped` from the fragment remap
- `sequences` and `split` carry `component`, which is the fold assignment unit
- an extra raw file, `targets_summary.json`

The tables are small, so everything is loaded in memory and pushed in one pass.

Run `--dry-run` first. It builds every config and prints the schema and row
counts without contacting the Hub.

Nothing in this script has been run against the Hub. The v1 repo at
petadex/petai-activity-labels is public and its card overstates the homology
guarantee, which `DATASET_CARD_v2.md` corrects. Pushing replaces that card, so
read it before running this.

Usage:
    python push_v2_to_hf.py <v2_dir> --repo petadex/petai-activity-labels --dry-run
    python push_v2_to_hf.py <v2_dir> --repo petadex/petai-activity-labels
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

BOOL_TRUE = {"true", "1", "yes", "t"}

# Columns that must stay typed rather than being inferred per-shard. `value`
# spans several units and is read as float; the condition columns are sparse and
# would otherwise flip between float and object across configs.
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
    "orig_seq_md5": "string",
    "remapped": "boolean",
}

# The normalised table is the labels table plus the step 3 columns.
NORMALISED_DTYPES = dict(LABEL_DTYPES, **{
    "unit_key": "string",
    "target": "string",
    "value_canonical": "float64",
    "unit_canonical": "string",
    "flag": "string",
    "split": "string",
})

CONFIGS = {
    "labels": "labels.csv",
    "sequences": "sequences.csv",
    "split": "split.csv",
    "clusters": "clusters.csv",
    "normalised": "normalised/labels_normalised.csv",
}

EXTRA_FILES = {
    "sequences.fasta": "sequences.fasta",
    "summary.json": "summary.json",
    "targets_summary.json": "normalised/targets_summary.json",
}


def coerce(df: pd.DataFrame, dtypes: dict[str, str]) -> pd.DataFrame:
    for col, dtype in dtypes.items():
        if col not in df.columns:
            continue
        if dtype == "boolean":
            df[col] = df[col].str.strip().str.lower().isin(BOOL_TRUE)
        elif dtype == "float64":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_config(root: Path, name: str, filename: str) -> pd.DataFrame:
    df = pd.read_csv(root / filename, dtype="string", keep_default_na=False)
    df = df.replace({"": None})
    if name == "labels":
        df = coerce(df, LABEL_DTYPES)
    elif name == "normalised":
        df = coerce(df, NORMALISED_DTYPES)
    elif name == "sequences":
        for col in ("length", "n_labels", "n_negative", "n_refs",
                    "n_absorbed_fragments"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        if "is_fragment" in df.columns:
            df["is_fragment"] = (
                df["is_fragment"].str.strip().str.lower().isin(BOOL_TRUE)
            )
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("v2", type=Path, help="the v2 release directory")
    ap.add_argument("--repo", required=True, help="e.g. petadex/petai-activity-labels")
    ap.add_argument("--private", action="store_true", help="create the repo private")
    ap.add_argument("--dry-run", action="store_true", help="build but do not push")
    ap.add_argument(
        "--card",
        type=Path,
        default=Path(__file__).parent / "DATASET_CARD_v2.md",
        help="dataset card to upload as README.md",
    )
    args = ap.parse_args()

    missing = [f for f in CONFIGS.values() if not (args.v2 / f).exists()]
    if missing:
        print(f"missing input files: {missing}", file=sys.stderr)
        print("run remap_fragments.py, split_components.py and "
              "normalise_units.py first", file=sys.stderr)
        return 1

    frames = {n: load_config(args.v2, n, f) for n, f in CONFIGS.items()}
    for name, df in frames.items():
        print(f"[{name}] {len(df):,} rows x {len(df.columns)} cols")
        print("   " + ", ".join(f"{c}:{df[c].dtype}" for c in df.columns))

    # Guard against pushing a split that leaks. This is the defect v2 exists to
    # fix, so it is checked here rather than trusted.
    seqs = frames["sequences"]
    for level in ("cluster_30", "cluster_50", "cluster_90", "component"):
        if level not in seqs.columns:
            continue
        spanning = seqs.groupby(level)["split"].nunique()
        n = int((spanning > 1).sum())
        print(f"   leakage check {level}: {n} groups span folds")
        if n:
            print(f"refusing to push: {n} {level} groups span folds", file=sys.stderr)
            return 1

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
    for target, source in EXTRA_FILES.items():
        src = args.v2 / source
        if not src.exists():
            continue
        print(f"uploading {target} ...")
        api.upload_file(
            path_or_fileobj=str(src),
            path_in_repo=f"raw/{target}",
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
