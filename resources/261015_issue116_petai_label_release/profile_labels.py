"""Profile and validate the PETai activity label release.

Reads the six files in the `combined/` release directory and reproduces every
number quoted in the Issue 116 notebook: release counts, the per-fold negative
deficit that makes the shipped split unusable for evaluation, the fragment
remap recovery, and the unit collisions in the value column.

Nothing is written. This is a read-only audit so the published figures can be
re-derived from the release alone.

Usage:
    python profile_labels.py <combined_dir> [--json out.json]
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

TRUE = {"true", "1", "yes", "t"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def count_fasta(path: Path) -> int:
    with path.open(encoding="utf-8") as fh:
        return sum(1 for line in fh if line.startswith(">"))


def audit(combined: Path) -> dict:
    labels = read_csv(combined / "labels.csv")
    sequences = read_csv(combined / "sequences.csv")
    split = read_csv(combined / "split.csv")
    clusters = read_csv(combined / "clusters.csv")
    summary = json.loads((combined / "summary.json").read_text())
    by_md5 = {s["seq_md5"]: s for s in sequences}

    out: dict = {}

    # --- Release shape, checked against the shipped summary.json ---------
    out["counts"] = {
        "labels": len(labels),
        "sequences": len(sequences),
        "split_rows": len(split),
        "cluster_rows": len(clusters),
        "fasta_records": count_fasta(combined / "sequences.fasta"),
        "refs": len({l["ref"] for l in labels if l["ref"]}),
    }
    out["summary_agrees"] = {
        "labels": summary["labels_total"] == len(labels),
        "sequences": summary["sequences_total"] == len(sequences),
        "clean": summary["sequences_labelled_clean"] == len(split),
    }

    # --- Blocker 1: the split cannot evaluate ----------------------------
    # Negatives are the scarce class. Count them per fold, and count the halo
    # screen separately since it supplies most of them.
    neg_labels: collections.Counter = collections.Counter()
    neg_seqs: dict[str, set] = collections.defaultdict(set)
    halo: collections.Counter = collections.Counter()
    fold_seqs: dict[str, set] = collections.defaultdict(set)

    for row in labels:
        seq = by_md5.get(row["seq_md5"])
        if seq is None:
            continue
        fold = seq["split"] or "DROPPED"
        fold_seqs[fold].add(row["seq_md5"])
        if row["is_negative"].strip().lower() in TRUE:
            neg_labels[fold] += 1
            neg_seqs[fold].add(row["seq_md5"])
        if row["source"] == "halo":
            halo[fold] += 1

    out["negatives"] = {
        "labels_by_fold": dict(neg_labels),
        "sequences_by_fold": {k: len(v) for k, v in neg_seqs.items()},
        "halo_labels_by_fold": dict(halo),
        "fold_sizes": {k: len(v) for k, v in fold_seqs.items()},
    }

    # --- Step 1: fragment remap recovery ---------------------------------
    # A sequence with no split is a fragment or a contained duplicate. Most
    # name a `contained_in` parent that is itself in a fold, so their labels
    # can be moved onto the parent rather than discarded.
    dropped = [m for m, s in by_md5.items() if not s["split"]]
    labels_per_seq = collections.Counter(l["seq_md5"] for l in labels)
    remappable = [
        m for m in dropped if by_md5.get(by_md5[m]["contained_in"], {}).get("split")
    ]
    out["fragments"] = {
        "dropped_sequences": len(dropped),
        "remappable_sequences": len(remappable),
        "orphaned_labels": sum(labels_per_seq[m] for m in dropped),
        "recoverable_labels": sum(labels_per_seq[m] for m in remappable),
    }
    out["fragments"]["lost_labels"] = (
        out["fragments"]["orphaned_labels"] - out["fragments"]["recoverable_labels"]
    )

    # --- Step 3: the unit collisions -------------------------------------
    # `value` is a single float column spanning several incompatible units, so
    # it is not one regression target.
    out["units"] = dict(collections.Counter(l["unit"] for l in labels).most_common(20))
    out["substrates"] = dict(
        collections.Counter(l["substrate"] for l in labels).most_common(12)
    )
    out["tasks"] = dict(collections.Counter(l["task"] for l in labels))
    out["condition_coverage"] = {
        col: sum(1 for l in labels if l[col] not in ("", "NA", "nan"))
        for col in ("ph", "temperature_c", "time_h")
    }

    # --- Cluster-size asymmetry, the reason CIs go over clusters ---------
    per_fold_clusters: dict[str, collections.Counter] = collections.defaultdict(
        collections.Counter
    )
    for seq in sequences:
        if seq["split"]:
            per_fold_clusters[seq["split"]][seq["cluster_30"]] += 1
    out["clusters_30"] = {}
    for fold, counter in per_fold_clusters.items():
        sizes = sorted(counter.values(), reverse=True)
        out["clusters_30"][fold] = {
            "clusters": len(sizes),
            "sequences": sum(sizes),
            "mean_size": round(sum(sizes) / len(sizes), 2),
            "max_size": sizes[0],
            "top3_share_pct": round(100 * sum(sizes[:3]) / sum(sizes), 2),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("combined", type=Path, help="the release directory")
    ap.add_argument("--json", type=Path, help="also write the audit as JSON")
    args = ap.parse_args()

    report = audit(args.combined)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if args.json:
        args.json.write_text(json.dumps(report, indent=2, ensure_ascii=False))
        print(f"\nwrote {args.json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
