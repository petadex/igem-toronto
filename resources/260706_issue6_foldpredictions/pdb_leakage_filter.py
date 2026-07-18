"""
Ground-truth Note 2, Method 1 (date-based leakage filter).

Of the cores that have a PDB crystal structure (from pdb_to_orfid_map.py), how
many are usable as *leakage-free* ground truth for the ESMFold2 base-vs-finetune
comparison?

Only ESMFold2's structure-training cutoff matters here: ESMFold2 is the only
component trained on 3D coordinates, so it is the only place a crystal structure
could be memorized. ESMC is a sequence-only LM and cannot leak a structure, so
its cutoff is irrelevant to structural ground truth (see experiment notes).

ESMFold2 (per the paper / HF card):
    training cutoff   2021-09-30   -> released after = never in training weights
    eval cutoff       2023-01-13   -> released after = not even in the val set
The 2021-09-30 .. 2023-01-13 window is the validation set (mild leakage only).

We pull each PDB's initial_release_date from RCSB (batched GraphQL), flag it
against both cutoffs, and propagate to ORFids: an ORFid is date-clean if it has
at least one post-cutoff PDB. Split by provenance so 'direct' (high-confidence)
and 'via_uniprot' (weaker link) clean cores are counted separately.

Outputs (this directory):
  pdb_dates.tsv             pdb_id, initial_release_date, deposit_date   (cached)
  pdb_orfid_map_dated.tsv   the map + release date + the two post-cutoff flags

Method 2 (40 pid near-neighbor to training chains) is NOT done here -- for a
single enzyme family a literal 40 pid cut removes nearly everything; that needs a
separate design decision.
"""

import argparse
import json
import time
import urllib.request
from pathlib import Path

import polars as pl

RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
TRAIN_CUTOFF = "2021-09-30"   # ESMFold2 training cutoff
EVAL_CUTOFF = "2023-01-13"    # ESMFold2 evaluation cutoff
GQL = ("query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id "
       "rcsb_accession_info{initial_release_date deposit_date}}}")


def fetch_dates(pdb_ids, chunk=250):
    """Batch-fetch release/deposit dates from RCSB GraphQL. Returns
    {pdb_id: (initial_release_date, deposit_date)}; missing/obsolete -> absent."""
    out = {}
    ids = sorted(pdb_ids)
    for i in range(0, len(ids), chunk):
        batch = ids[i:i + chunk]
        body = json.dumps({"query": GQL, "variables": {"ids": batch}}).encode()
        req = urllib.request.Request(
            RCSB_GRAPHQL, data=body,
            headers={"Content-Type": "application/json"})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    payload = json.load(r)
                break
            except Exception as e:
                if attempt == 2:
                    raise
                print(f"    retry {i}-{i+len(batch)} after {e}")
                time.sleep(3)
        for e in (payload.get("data", {}).get("entries") or []):
            info = e.get("rcsb_accession_info") or {}
            out[e["rcsb_id"].upper()] = (
                (info.get("initial_release_date") or "")[:10],
                (info.get("deposit_date") or "")[:10],
            )
        print(f"    fetched {min(i + chunk, len(ids)):,}/{len(ids):,}")
    return out


def load_or_fetch_dates(pdb_ids, dates_path, force):
    if dates_path.exists() and not force:
        print(f"[*] Reusing cached {dates_path.name} (pass --force to refetch)")
        return pl.read_csv(dates_path, separator='\t',
                           schema_overrides={c: pl.Utf8 for c in
                                             ("pdb_id", "initial_release_date", "deposit_date")})
    print(f"[*] Fetching release dates for {len(pdb_ids):,} PDBs from RCSB...")
    start = time.time()
    d = fetch_dates(pdb_ids)
    df = pl.DataFrame(
        {"pdb_id": list(pdb_ids),
         "initial_release_date": [d.get(p, ("", ""))[0] for p in pdb_ids],
         "deposit_date": [d.get(p, ("", ""))[1] for p in pdb_ids]}
    ).sort("pdb_id")
    df.write_csv(dates_path, separator='\t')
    missing = df.filter(pl.col("initial_release_date") == "").height
    print(f"[+] Dates fetched in {time.time() - start:.0f}s "
          f"({missing:,} PDBs returned no date -- obsolete/withdrawn)")
    return df


def main():
    script_dir = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path, default=script_dir / "pdb_orfid_map.tsv")
    ap.add_argument("--out-dir", type=Path, default=script_dir)
    ap.add_argument("--force", action="store_true", help="refetch dates from RCSB")
    args = ap.parse_args()

    dates_path = args.out_dir / "pdb_dates.tsv"
    dated_map_path = args.out_dir / "pdb_orfid_map_dated.tsv"

    m = pl.read_csv(args.map, separator='\t', schema_overrides={"orfid": pl.Utf8})
    pdb_ids = m["pdb_id"].unique().to_list()

    dates = load_or_fetch_dates(pdb_ids, dates_path, args.force)

    # Attach dates + post-cutoff flags to every (core, pdb) row -------------
    dated = m.join(dates, on="pdb_id", how="left").with_columns([
        (pl.col("initial_release_date") > TRAIN_CUTOFF).alias("post_train_2021_09_30"),
        (pl.col("initial_release_date") > EVAL_CUTOFF).alias("post_eval_2023_01_13"),
    ])
    dated.write_csv(dated_map_path, separator='\t')

    # PDB-level counts -------------------------------------------------------
    pdb_lvl = dated.select(
        ["pdb_id", "initial_release_date",
         "post_train_2021_09_30", "post_eval_2023_01_13"]).unique()
    n_pdb = pdb_lvl.height
    n_pdb_dated = pdb_lvl.filter(pl.col("initial_release_date") != "").height

    def orfid_clean(flag_col, provenance=None):
        d = dated.filter(pl.col(flag_col))
        if provenance:
            d = d.filter(pl.col("provenance") == provenance)
        return set(d["orfid"].to_list())

    print("=" * 64)
    print(f"distinct PDBs: {n_pdb:,} (with a release date: {n_pdb_dated:,})")
    print(f"total ground-truth ORFids (any PDB): {dated['orfid'].n_unique():,}")
    print("-" * 64)
    for label, flag in (("TRAIN cutoff  > 2021-09-30 (not in training weights)", "post_train_2021_09_30"),
                        ("EVAL  cutoff  > 2023-01-13 (never seen at all)", "post_eval_2023_01_13")):
        n_pdb_post = pdb_lvl.filter(pl.col(flag)).height
        any_clean = orfid_clean(flag)
        direct_clean = orfid_clean(flag, "direct")
        via_clean = orfid_clean(flag, "via_uniprot")
        via_only = via_clean - direct_clean
        print(f"{label}")
        print(f"    clean PDBs:               {n_pdb_post:,}")
        print(f"    clean ORFids (any):       {len(any_clean):,}")
        print(f"      via a direct PDB:       {len(direct_clean):,}")
        print(f"      only via_uniprot PDB:   {len(via_only):,}")
        print("-" * 64)

    # Year histogram ---------------------------------------------------------
    print("release-year histogram (distinct PDBs):")
    hist = (pdb_lvl.filter(pl.col("initial_release_date") != "")
            .with_columns(pl.col("initial_release_date").str.slice(0, 4).alias("yr"))
            .group_by("yr").len().sort("yr"))
    for row in hist.iter_rows(named=True):
        print(f"    {row['yr']}: {row['len']:,}")
    print("=" * 64)
    print(f"[+] wrote {dates_path.name} and {dated_map_path.name}")


if __name__ == "__main__":
    main()
