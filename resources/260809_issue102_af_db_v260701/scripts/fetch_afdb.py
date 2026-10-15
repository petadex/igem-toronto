#!/usr/bin/env python3
"""Download AlphaFold DB structures in resumable batches and stage them to S3.

Run 1 downloaded the whole set to local disk first and needed ~169 GB of it. This script
downloads a batch, syncs it to S3, deletes it, and moves on, so peak disk is one batch
(~8 GB at the default 20,000 files) rather than the whole set. That is what lets the job run
on a modest scratch volume attached to an existing instance.

Resumability follows the ledger pattern already used elsewhere in this repo (see
resources/260729_issue6d_finalmetrics/): a completed batch is a row in batches.csv, and a
restart skips every batch already recorded there. The ledger is mirrored to S3 after each
batch so an instance failure loses at most one batch of progress.

Three ledger files are written under --ledger-dir:
    batches.csv   batch_idx, n_requested, n_downloaded, bytes, seconds, status, finished_at
    hits.txt      one UniProt AC per successfully downloaded structure
    misses.txt    one UniProt AC per requested-but-absent structure

hits.txt is the input to build_unfolded_fasta.py -- it is the authoritative record of what
AFDB actually had, which is not knowable up front unless build_manifest.py was given an index.

Prerequisites: an authenticated `gcloud` (the public-datasets AlphaFold bucket denies anonymous
reads -- verified 403 on both the bucket index and a known-good object) and AWS credentials with
write access to the target prefix.

Usage:
    python fetch_afdb.py --manifest manifests/afdb_manifest.txt \
        --workdir /mnt/afdb/staging \
        --s3 s3://petadex-protein-structures/af_db/v260701/structures/ \
        --ledger-dir /mnt/afdb/ledger \
        [--batch-size 20000] [--jobs 32] [--max-batches 1]   # --max-batches 1 = pilot
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

AC_FROM_NAME_RE = re.compile(r"^AF-([A-Z0-9]+)-F\d+-model_v\d+\.cif$")

# Resolved once at startup by main(). Full paths rather than bare names so the calls do not
# depend on the parent shell's PATHEXT resolution.
GCLOUD = "gcloud"
AWS = "aws"


def ac_from_uri(uri):
    return AC_FROM_NAME_RE.match(os.path.basename(uri)).group(1)


def read_manifest(path):
    with open(path, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def completed_batches(ledger_csv):
    """Batch indices already recorded as done, so a restart can skip them."""
    done = set()
    if not os.path.exists(ledger_csv):
        return done
    with open(ledger_csv, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            if row.get("status") == "ok":
                done.add(int(row["batch_idx"]))
    return done


def run(cmd, **kwargs):
    return subprocess.run(cmd, **kwargs)


def download_batch(uris, dest, jobs):
    """Fetch one batch with gcloud storage. -c continues past missing objects."""
    os.makedirs(dest, exist_ok=True)
    cmd = [GCLOUD, "storage", "cp", "-c", "-I", dest]
    env = dict(os.environ)
    # gcloud reads parallelism from config; setting it per-invocation keeps the shared VM's
    # global gcloud config untouched.
    env["CLOUDSDK_STORAGE_PROCESS_COUNT"] = str(jobs)
    env["CLOUDSDK_STORAGE_THREAD_COUNT"] = str(jobs)
    proc = run(cmd, input="\n".join(uris) + "\n", text=True, env=env,
               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # A non-zero exit with -c just means some objects were absent, which is expected.
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()[-3:]
        print(f"    gcloud exit {proc.returncode} (expected when objects are missing): {tail}")
    return proc.returncode


def sync_batch(src, s3_prefix):
    proc = run([AWS, "s3", "sync", src, s3_prefix, "--only-show-errors"])
    if proc.returncode != 0:
        raise RuntimeError(f"aws s3 sync failed with exit {proc.returncode}")


def mirror_ledger(ledger_dir, s3_ledger):
    if not s3_ledger:
        return
    run([AWS, "s3", "sync", ledger_dir, s3_ledger, "--only-show-errors"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--workdir", required=True, help="scratch dir for the current batch")
    ap.add_argument("--s3", required=True, help="destination S3 prefix for structures/")
    ap.add_argument("--ledger-dir", required=True)
    ap.add_argument("--s3-ledger", default=None, help="optional S3 prefix to mirror the ledger to")
    ap.add_argument("--batch-size", type=int, default=20000)
    ap.add_argument("--jobs", type=int, default=32)
    ap.add_argument("--max-batches", type=int, default=None, help="stop after N batches (pilot runs)")
    ap.add_argument("--keep-local", action="store_true", help="do not delete a batch after syncing")
    args = ap.parse_args()

    global GCLOUD, AWS
    resolved = {}
    for tool in ("gcloud", "aws"):
        path = shutil.which(tool)
        if path is None:
            sys.exit(f"[-] '{tool}' not found on PATH.")
        resolved[tool] = path
    GCLOUD, AWS = resolved["gcloud"], resolved["aws"]

    os.makedirs(args.ledger_dir, exist_ok=True)
    ledger_csv = os.path.join(args.ledger_dir, "batches.csv")
    hits_path = os.path.join(args.ledger_dir, "hits.txt")
    misses_path = os.path.join(args.ledger_dir, "misses.txt")

    uris = read_manifest(args.manifest)
    batches = [uris[i:i + args.batch_size] for i in range(0, len(uris), args.batch_size)]
    done = completed_batches(ledger_csv)
    todo = [i for i in range(len(batches)) if i not in done]
    if args.max_batches is not None:
        todo = todo[: args.max_batches]

    print(f"[*] {len(uris):,} URIs in {len(batches):,} batches of {args.batch_size:,}")
    print(f"[*] {len(done):,} already complete; running {len(todo):,} now.")

    if not os.path.exists(ledger_csv):
        with open(ledger_csv, "w", encoding="utf-8", newline="") as f:
            csv.writer(f).writerow(
                ["batch_idx", "n_requested", "n_downloaded", "bytes", "seconds", "status", "finished_at"]
            )

    total_hits = total_miss = total_bytes = 0
    run_start = time.time()

    for n, idx in enumerate(todo, 1):
        batch = batches[idx]
        dest = os.path.join(args.workdir, f"batch_{idx:06d}")
        shutil.rmtree(dest, ignore_errors=True)
        t0 = time.time()
        print(f"[{n}/{len(todo)}] batch {idx}: {len(batch):,} objects -> {dest}")

        try:
            download_batch(batch, dest, args.jobs)

            got = {}
            for name in os.listdir(dest):
                m = AC_FROM_NAME_RE.match(name)
                if m:
                    got[m.group(1)] = os.path.getsize(os.path.join(dest, name))
            requested = {ac_from_uri(u) for u in batch}
            missing = sorted(requested - set(got))
            nbytes = sum(got.values())

            sync_batch(dest, args.s3)

            with open(hits_path, "a", encoding="utf-8", newline="\n") as f:
                for ac in sorted(got):
                    f.write(ac + "\n")
            with open(misses_path, "a", encoding="utf-8", newline="\n") as f:
                for ac in missing:
                    f.write(ac + "\n")

            elapsed = time.time() - t0
            with open(ledger_csv, "a", encoding="utf-8", newline="") as f:
                csv.writer(f).writerow(
                    [idx, len(batch), len(got), nbytes, round(elapsed, 1), "ok",
                     datetime.now(timezone.utc).isoformat(timespec="seconds")]
                )

            total_hits += len(got)
            total_miss += len(missing)
            total_bytes += nbytes
            rate = len(batch) / elapsed if elapsed else 0
            print(f"    {len(got):,} hits / {len(missing):,} misses, "
                  f"{nbytes / 1e9:.2f} GB, {elapsed:,.0f}s ({rate:,.0f} obj/s)")
        finally:
            if not args.keep_local:
                shutil.rmtree(dest, ignore_errors=True)

        mirror_ledger(args.ledger_dir, args.s3_ledger)

    total_elapsed = time.time() - run_start
    print("[+] Run complete.")
    print(f"    batches run : {len(todo):,}")
    print(f"    hits        : {total_hits:,}")
    print(f"    misses      : {total_miss:,}")
    print(f"    downloaded  : {total_bytes / 1e9:.2f} GB in {total_elapsed / 60:,.1f} min")
    remaining = len(batches) - len(done) - len(todo)
    if remaining > 0 and total_elapsed > 0:
        per_batch = total_elapsed / len(todo)
        print(f"    {remaining:,} batches remain -> ~{remaining * per_batch / 3600:,.1f} h at this rate")


if __name__ == "__main__":
    main()
