"""
foldseek_novelty.py -- figure Db: structural novelty of a fold run against a reference structure DB.

Runs ON an EC2 VM in the bucket's region (see the notebook's Db VM section). Pulls a run's predicted
structures down to local disk, searches all of them against a Foldseek target database in ONE pass,
and writes a resumable per-ORF ledger of the best structural hit plus a mergeable aggregate:

    python foldseek_novelty.py \
        --run s3://petadex-protein-structures/esmfold2_paramsweep/s68_l10/ \
        --target ~/db/pdb --out nov_s68_l10 --threads 16

`--run` is an `s3://bucket/prefix/` or a local dir laid out the way esmfold2_local_predictor.py
writes one (`<root>/structures/<name>.cif`). `--target` is a database built by
`foldseek databases PDB ~/db/pdb ~/tmp`.

WHY THIS IS A DIFFERENT SHAPE OF JOB FROM conf_metrics.py (Da)
  Da reads ~4 KB per ORF and never touches a structure, so it streams from S3 with a big thread pool
  and never lands anything on disk. Db needs the structures THEMSELVES: Foldseek wants every query in
  one database it can mmap, because the entire point of the tool is that the target index is loaded
  once and amortised over the whole query set. Per-structure `easy-search` would reload the target DB
  630k times and throw away the only optimisation that matters.

  So the pipeline is staged, and each stage is separately resumable:

      fetch     S3 -> local structures/        (the slow, restartable part: ~100 GB at full scale)
      createdb  structures/ -> queryDB          (foldseek createdb)
      search    queryDB x targetDB -> alnDB     (foldseek search; the actual compute)
      reduce    alignments -> novelty.tsv + aggregates.json

  --stage runs a subset; the default runs all four and skips any whose output already exists.

READ THIS BEFORE INTERPRETING THE OUTPUT -- three ways to get a wrong novelty number
  1. "NO HIT" IS NOT PROOF OF NOVELTY. Foldseek's speed comes from a k-mer prefilter that discards
     most of the target DB before any alignment happens. It is a heuristic: a query with no reported
     hit means "the prefilter found nothing", not "nothing similar exists". That error points in
     exactly the direction that flatters a novelty claim, so it cannot be left unexamined. Mitigations
     built in here: --evalue defaults to a PERMISSIVE 10 (not Foldseek's 1e-3) so weak hits are still
     reported, and --exhaustive re-runs a subset with the prefilter disabled. Validate the no-hit set
     with --exhaustive before any novelty number goes in the writeup.
  2. TM-SCORE NORMALISATION IS NOT A DETAIL. TM-score is asymmetric. `qtmscore` is normalised by the
     QUERY length, `ttmscore` by the target. For "does any solved structure resemble MY protein",
     query-normalised is the honest one -- otherwise a small PDB domain matching a fragment of a long
     ORF looks like a strong hit and the ORF is wrongly called non-novel. `alntmscore` additionally
     has an open report of returning the same value as `ttmscore`
     (github.com/steineggerlab/foldseek/issues/312), so all three are recorded and the CHOICE is made
     in the notebook, on evidence, not here.
  3. BEST-BY-E-VALUE != BEST-BY-TM. Foldseek ranks by E-value; the novelty claim is about TM. They
     disagree often enough to matter (a short, near-perfect domain match can outrank a long, moderate
     whole-chain one). Both best hits are recorded per ORF, with the field that decides "the" score
     chosen by --score-field.

WHAT THE SCORE DOES AND DOESN'T LICENSE
  This measures RESEMBLANCE TO A REFERENCE SET, not correctness and not function. A structure with no
  PDB neighbour may be a genuinely new fold, or may be a bad prediction -- the two are indistinguish-
  able from this number alone. That is exactly why Db is plotted against confidence: the interesting
  quadrant is CONFIDENTLY novel (high pLDDT, low best-hit TM), and neither axis means much by itself.
  Pass --conf-tsv to join Da's per-ORF confidence in and get that 2D aggregate for free.

  Note also that the target DB choice IS the claim. Novelty vs PDB ("no solved structure looks like
  this") and vs AFDB ("nothing predicted looks like this") are different statements, and ESMFold2 was
  trained on AFDB, so a low AFDB hit is entangled with the training set in a way a low PDB hit is not.
  This script is agnostic; the notebook has to say which DB it ran and why.

BUILT FOR THE FULL CENTROID SET, same conventions as conf_metrics.py / pair_tm.py:
  * `--shard K/N` fans one command across a fleet; each shard owns a disjoint ORFid band, its own
    local structure dir, its own query DB and its own ledger, so workers never clobber each other;
  * ledgers are append-only and re-read on start, so an interrupted run RESUMES;
  * `aggregates*.json` (fixed-bin 1D + 2D histograms, Welford moments, EXACT threshold counts,
    reservoir subsample) is O(1) in N and mergeable with `--merge`, so the figure costs the same at
    121 ORFs or 2 million;
  * `--s3-out` mirrors ledger + aggregates back to S3 so a replaced spot instance resumes from its own.

Deps on the VM: foldseek on PATH (static binary is fine -- there is no build step, unlike MolProbity),
plus `pip install boto3`. numpy is NOT required.
"""
import argparse
import csv
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor


# ---------------------------------------------------------------------------
# the columns we ask Foldseek for
# ---------------------------------------------------------------------------
# Order matters: convertalis emits them positionally, so this tuple IS the TSV schema. All three
# TM normalisations are requested deliberately -- see note 2 in the module docstring.
FORMAT_FIELDS = (
    "query", "target", "fident", "alnlen", "qlen", "tlen",
    "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "prob", "lddt",
    "alntmscore", "qtmscore", "ttmscore",
)

# Fields eligible for --score-field, i.e. the number that becomes "the" novelty score.
TM_FIELDS = ("qtmscore", "ttmscore", "alntmscore")

STRUCT_EXTS = (".cif", ".cif.gz", ".pdb", ".pdb.gz", ".mmcif", ".mmcif.gz")


# ---------------------------------------------------------------------------
# ORFid <-> filename
# ---------------------------------------------------------------------------
# Two naming conventions exist in the bucket and both must round-trip:
#   * the original local predictor wrote  structures/orf<ORFid>.cif
#   * the --s3-only fleet path writes     structures/<ORFid>.cif      (bare id)
# On top of that, `foldseek createdb` derives its entry name from the file and may or may not keep
# the extension, and appends the chain label -- so the SAME structure can appear as any of
# `orf12345.cif`, `orf12345.cif_A`, `orf12345_A`, `12345_A`, `12345`. All of them have to collapse to
# the same join key, because ORFid is what joins Db to Da's conf.tsv.
#
# The extension is therefore OPTIONAL here. That costs a little safety on the chain suffix, so the
# suffix is only stripped when it is unambiguous (see below) rather than whenever it matches.
_EXT = r"\.(?:cif|mmcif|pdb)"
_CHAIN = r"[A-Za-z0-9]{1,4}"
# Order matters: the ext+chain form must be tried BEFORE the bare-ext form, or `x.cif_A` fails the
# `$`-anchored ext match, falls through to chain-stripping, and leaves a stray `.cif` on the id.
_EXT_CHAIN_RX = re.compile(rf"^(?P<stem>.+){_EXT}_(?P<chain>{_CHAIN})$", re.IGNORECASE)
_EXT_RX = re.compile(rf"^(?P<stem>.+){_EXT}$", re.IGNORECASE)
_CHAIN_RX = re.compile(rf"^(?P<stem>.+)_(?P<chain>{_CHAIN})$")
# A residual alphabetic extension means this was never one of our structures (`junk.txt`). Digits are
# allowed through so a versioned accession like WP_012345.1 survives.
_FOREIGN_EXT_RX = re.compile(r"\.[A-Za-z]{1,6}$")


def orfid_from_name(name):
    """ORFid from a structure filename or a foldseek query entry name. None if unrecognisable.

    Handles every shape the two bucket conventions and foldseek's own entry naming can produce:
    `orf12345.cif`, `orf12345.cif_A`, `orf12345_A`, `12345_A`, `12345`, plus `.gz`.

    Anything unrecognisable returns None rather than a guess: a silently mangled id would join to the
    WRONG ORF in conf.tsv, which is far worse than a visible gap.
    """
    s = os.path.basename(str(name)).strip()
    if not s:
        return None
    if s.lower().endswith(".gz"):
        s = s[:-3]

    m = _EXT_CHAIN_RX.match(s)
    if m:
        s = m.group("stem")
    else:
        m = _EXT_RX.match(s)
        if m:
            s = m.group("stem")
        else:
            # No extension to disambiguate, so only strip a chain suffix when the stem holds no
            # further underscore -- otherwise an ORFid that legitimately contains one is truncated.
            m = _CHAIN_RX.match(s)
            if m and "_" not in m.group("stem"):
                s = m.group("stem")

    if _FOREIGN_EXT_RX.search(s):
        return None
    if s[:3].lower() == "orf" and len(s) > 3:
        s = s[3:]
    return s or None


def is_structure(name):
    return name.lower().endswith(STRUCT_EXTS)


def shard_of(orfid, n):
    """Stable shard assignment. Uses a content hash, NOT python's salted hash(), so the same id lands
    in the same shard on every process and every machine."""
    import hashlib
    return int(hashlib.md5(str(orfid).encode()).hexdigest()[:8], 16) % n


# ---------------------------------------------------------------------------
# fetching structures out of S3 (or a local dir)
# ---------------------------------------------------------------------------
class StructSource:
    """Lists and downloads <root>/structures/*. `root` is s3://bucket/prefix/ or a local dir.

    One boto3 client PER THREAD -- clients are only nominally thread-safe and sharing one across a
    pool serialises on its connection pool, which is the whole point of the pool.
    """

    def __init__(self, root, retries=4):
        self.root, self.retries = root, retries
        self.is_s3 = str(root).startswith("s3://")
        self._tl = threading.local()
        if self.is_s3:
            self.bucket, _, self.prefix = root[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")

    @property
    def s3(self):
        if getattr(self._tl, "c", None) is None:
            import boto3
            self._tl.c = boto3.client("s3")
        return self._tl.c

    def list_structures(self):
        """[(orfid, key_or_path, size)] for every structure under <root>/structures/."""
        out = []
        if not self.is_s3:
            d = os.path.join(self.root, "structures")
            if not os.path.isdir(d):
                d = self.root
            for n in sorted(os.listdir(d)):
                if not is_structure(n):
                    continue
                oid = orfid_from_name(n)
                if oid:
                    p = os.path.join(d, n)
                    out.append((oid, p, os.path.getsize(p)))
            return out

        pfx = "/".join(p for p in (self.prefix, "structures/") if p)
        tok = None
        while True:
            kw = {"Bucket": self.bucket, "Prefix": pfx}
            if tok:
                kw["ContinuationToken"] = tok
            r = self.s3.list_objects_v2(**kw)
            for o in r.get("Contents", []):
                name = o["Key"].rsplit("/", 1)[-1]
                if not is_structure(name):
                    continue
                oid = orfid_from_name(name)
                if oid:
                    out.append((oid, o["Key"], o["Size"]))
            if not r.get("IsTruncated"):
                break
            tok = r["NextContinuationToken"]
        return out

    def download(self, key, dest):
        """Download to `dest` atomically -- a half-written CIF left by a killed process would be fed
        to createdb on the next run and poison the query DB."""
        if not self.is_s3:
            shutil.copyfile(key, dest)
            return
        tmp = dest + ".part"
        for attempt in range(self.retries):
            try:
                self.s3.download_file(self.bucket, key, tmp)
                os.replace(tmp, dest)
                return
            except Exception:                                     # noqa: BLE001
                if os.path.exists(tmp):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass
                if attempt == self.retries - 1:
                    raise
                time.sleep(0.4 * 2 ** attempt)


def stage_fetch(src, dest_dir, shard, workers, limit=None, log=print):
    """Mirror the run's structures into `dest_dir`. Resumable: an existing non-empty file is skipped.

    Returns the list of ORFids this worker owns -- which is also the DENOMINATOR for the novelty
    figure, and the only way a no-hit ORF can be reported at all.
    """
    os.makedirs(dest_dir, exist_ok=True)
    listing = src.list_structures()
    log(f"[fetch] {len(listing)} structures under {src.root}")

    if shard:
        k, n = shard
        listing = [t for t in listing if shard_of(t[0], n) == k]
        log(f"[fetch] shard {k}/{n} owns {len(listing)}")
    if limit:
        listing = listing[:limit]
        log(f"[fetch] --limit {limit}")

    # Duplicate ORFids would silently become one query. Keep the first, report the collision.
    seen, uniq = {}, []
    for oid, key, size in listing:
        if oid in seen:
            continue
        seen[oid] = True
        uniq.append((oid, key, size))
    if len(uniq) != len(listing):
        log(f"[fetch] WARNING: {len(listing) - len(uniq)} duplicate ORFids collapsed")

    todo = []
    for oid, key, size in uniq:
        ext = ".cif.gz" if str(key).endswith(".gz") else os.path.splitext(str(key))[1] or ".cif"
        dest = os.path.join(dest_dir, f"{oid}{ext}")
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            continue
        todo.append((key, dest))

    log(f"[fetch] {len(todo)} to download, {len(uniq) - len(todo)} already local")
    if todo:
        t0, done, errs = time.time(), [0], []
        lock = threading.Lock()

        def one(job):
            key, dest = job
            try:
                src.download(key, dest)
            except Exception as e:                                # noqa: BLE001
                with lock:
                    errs.append((key, repr(e)))
            with lock:
                done[0] += 1
                if done[0] % 500 == 0:
                    el = time.time() - t0
                    log(f"[fetch] {done[0]}/{len(todo)}  {done[0]/max(el,1e-9):.1f}/s")

        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(one, todo))
        log(f"[fetch] done in {time.time()-t0:.1f}s, {len(errs)} errors")
        for k, e in errs[:10]:
            log(f"[fetch]   ERR {k}: {e}")

    return [oid for oid, _, _ in uniq]


# ---------------------------------------------------------------------------
# shelling out to foldseek
# ---------------------------------------------------------------------------
def sh(cmd, log=print, check=True):
    """Run a command, streaming nothing but reporting the tail on failure."""
    log(f"[cmd] {' '.join(str(c) for c in cmd)}")
    t0 = time.time()
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    dt = time.time() - t0
    if p.returncode != 0:
        log(f"[cmd] FAILED rc={p.returncode} after {dt:.1f}s")
        log((p.stdout or "")[-2000:])
        log((p.stderr or "")[-2000:])
        if check:
            raise SystemExit(f"foldseek command failed: {' '.join(str(c) for c in cmd)}")
    else:
        log(f"[cmd] ok in {dt:.1f}s")
    return p


def foldseek_version(binary="foldseek"):
    try:
        p = subprocess.run([binary, "version"], capture_output=True, text=True)
        return (p.stdout or p.stderr or "").strip().splitlines()[0] if p.returncode == 0 else None
    except FileNotFoundError:
        return None


def stage_createdb(struct_dir, query_db, binary, log=print):
    # createdb writes several sidecar files; the base file existing is foldseek's own resume signal.
    if os.path.exists(query_db) and os.path.getsize(query_db) > 0:
        log(f"[createdb] {query_db} exists, skipping")
        return
    os.makedirs(os.path.dirname(os.path.abspath(query_db)), exist_ok=True)
    sh([binary, "createdb", struct_dir, query_db], log=log)


def stage_search(query_db, target_db, aln_db, tmp_dir, binary, threads, evalue, max_seqs,
                 alignment_type, exhaustive, extra, log=print):
    if os.path.exists(aln_db) or os.path.exists(aln_db + ".dbtype"):
        log(f"[search] {aln_db} exists, skipping")
        return
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = [binary, "search", query_db, target_db, aln_db, tmp_dir,
           "--threads", threads,
           "-e", evalue,
           "--max-seqs", max_seqs,
           "-a", "1",                       # keep backtrace: required for the TM/LDDT columns
           "--alignment-type", alignment_type]
    if exhaustive:
        # Disables the k-mer prefilter -> every query is aligned against every target. Orders of
        # magnitude slower; only ever for validating a no-hit subset (see note 1).
        cmd += ["--exhaustive-search", "1"]
    cmd += list(extra)
    sh(cmd, log=log)


def stage_convert(query_db, target_db, aln_db, out_tsv, binary, threads, log=print):
    if os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0:
        log(f"[convertalis] {out_tsv} exists, skipping")
        return
    sh([binary, "convertalis", query_db, target_db, aln_db, out_tsv,
        "--threads", threads,
        "--format-output", ",".join(FORMAT_FIELDS)], log=log)


# ---------------------------------------------------------------------------
# reduce: alignments -> one row per ORF
# ---------------------------------------------------------------------------
def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def best_hits(aln_tsv, score_field, log=print):
    """Fold the alignment TSV down to the best hit per ORF, streaming.

    The TSV is one row per (query, target) pair and is the only big intermediate -- at full scale it
    can be tens of GB, so it is never held in memory. Two winners are tracked per ORF because
    Foldseek ranks by E-value while the novelty claim is about TM (note 3).
    """
    idx = {name: i for i, name in enumerate(FORMAT_FIELDS)}
    best = {}
    n_rows, n_bad = 0, 0
    bad_names = []          # kept so an id-parsing failure names itself instead of looking like
                            # "the search found nothing" -- these two are indistinguishable in the
                            # totals, and one is a bug while the other is the headline result.

    with open(aln_tsv, newline="") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            n_rows += 1
            if len(parts) < len(FORMAT_FIELDS):
                n_bad += 1
                if len(bad_names) < 5:
                    bad_names.append(f"<{len(parts)} cols> {line[:80]}")
                continue
            oid = orfid_from_name(parts[idx["query"]])
            if oid is None:
                n_bad += 1
                if len(bad_names) < 5:
                    bad_names.append(parts[idx["query"]])
                continue

            ev = _f(parts[idx["evalue"]])
            tm = _f(parts[idx[score_field]])
            rec = best.get(oid)
            if rec is None:
                rec = best[oid] = {"n_hits": 0, "by_evalue": None, "by_tm": None}
            rec["n_hits"] += 1

            if ev is not None and (rec["by_evalue"] is None
                                   or ev < rec["by_evalue"]["_ev"]
                                   # tie on E-value -> prefer the higher bitscore, then higher TM
                                   or (ev == rec["by_evalue"]["_ev"]
                                       and (_f(parts[idx["bits"]]) or -1) > rec["by_evalue"]["_bits"])):
                rec["by_evalue"] = {"_ev": ev,
                                    "_bits": _f(parts[idx["bits"]]) or -1.0,
                                    **{k: parts[idx[k]] for k in FORMAT_FIELDS}}
            if tm is not None and (rec["by_tm"] is None or tm > rec["by_tm"]["_tm"]):
                rec["by_tm"] = {"_tm": tm, **{k: parts[idx[k]] for k in FORMAT_FIELDS}}

    log(f"[reduce] {n_rows} alignment rows, {len(best)} ORFs with >=1 hit, {n_bad} unparsable")
    if n_bad:
        log(f"[reduce] WARNING: {n_bad}/{n_rows} rows had an unparsable query name. Samples:")
        for b in bad_names:
            log(f"[reduce]   {b!r}")
        log("[reduce]   ^ these rows found REAL hits that are being thrown away -- fix "
            "orfid_from_name() before trusting any novelty number from this run.")
    return best


LEDGER_FIELDS = [
    "orfid", "status", "n_hits",
    # winner by the TM field named in --score-field: this is the novelty score
    "tm_target", "tm_score", "qtmscore", "ttmscore", "alntmscore",
    "tm_evalue", "tm_prob", "tm_lddt", "tm_fident", "tm_alnlen", "tm_qlen", "tm_tlen",
    # winner by foldseek's own ranking, kept because it disagrees with the above often enough
    "ev_target", "ev_evalue", "ev_bits", "ev_qtmscore",
    # joined from Da, if --conf-tsv was given
    "mean_plddt", "ptm",
]


def write_ledger(path, owned_ids, best, score_field, conf, log=print):
    """One row per OWNED ORF -- including the ones with no hit.

    This is the correctness property the whole figure rests on: a query Foldseek returned nothing for
    is the MOST novel structure in the set, and if it is dropped here the novelty distribution is
    silently truncated at exactly the end that matters.
    """
    n_hit = n_nohit = 0
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_FIELDS, delimiter="\t")
        w.writeheader()
        for oid in owned_ids:
            c = conf.get(oid, {})
            row = {"orfid": oid, "mean_plddt": c.get("mean_plddt"), "ptm": c.get("ptm")}
            rec = best.get(oid)
            if rec is None:
                n_nohit += 1
                row.update(status="no_hit", n_hits=0)
            else:
                n_hit += 1
                t, e = rec["by_tm"], rec["by_evalue"]
                row.update(status="ok", n_hits=rec["n_hits"])
                if t:
                    row.update(tm_target=t["target"], tm_score=t[score_field],
                               qtmscore=t["qtmscore"], ttmscore=t["ttmscore"],
                               alntmscore=t["alntmscore"],
                               tm_evalue=t["evalue"], tm_prob=t["prob"], tm_lddt=t["lddt"],
                               tm_fident=t["fident"], tm_alnlen=t["alnlen"],
                               tm_qlen=t["qlen"], tm_tlen=t["tlen"])
                if e:
                    row.update(ev_target=e["target"], ev_evalue=e["evalue"],
                               ev_bits=e["bits"], ev_qtmscore=e["qtmscore"])
            w.writerow(row)
    log(f"[reduce] ledger {path}: {n_hit} with a hit, {n_nohit} no_hit "
        f"({100.0*n_nohit/max(len(owned_ids),1):.1f}% -- validate these with --exhaustive)")


def read_conf_tsv(path, log=print):
    """Da's conf.tsv -> {orfid: {mean_plddt, ptm}}. The 'free join' that gives Db its second axis."""
    if not path:
        return {}
    out = {}
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            oid = row.get("orfid") or row.get("orf_id") or row.get("id")
            if not oid:
                continue
            out[oid] = {"mean_plddt": _f(row.get("mean_plddt")), "ptm": _f(row.get("ptm"))}
    log(f"[reduce] joined confidence for {len(out)} ORFs from {path}")
    return out


# ---------------------------------------------------------------------------
# aggregates: O(1) in N, mergeable across shards
# ---------------------------------------------------------------------------
TM_BINS = 50          # 0..1
PLDDT_BINS = 50       # 0..100
TM_THRESHOLDS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
PLDDT_THRESHOLDS = (70.0, 80.0, 90.0)
RESERVOIR_N = 5000


def _bin(v, lo, hi, n):
    if v is None:
        return None
    i = int((v - lo) / (hi - lo) * n)
    return min(max(i, 0), n - 1)


def build_aggregates(ledger_path, score_field, target_db, run_root, log=print):
    agg = {
        "score_field": score_field,
        "target_db": str(target_db),
        "run": str(run_root),
        "n_total": 0, "n_no_hit": 0, "n_scored": 0, "n_with_conf": 0,
        "tm_hist": [0] * TM_BINS,
        "tm_bins": [0.0, 1.0, TM_BINS],
        "plddt_hist": [0] * PLDDT_BINS,
        "plddt_bins": [0.0, 100.0, PLDDT_BINS],
        # joint histogram, row-major [tm][plddt] -- the confidently-novel quadrant lives here
        "joint_hist": [[0] * PLDDT_BINS for _ in range(TM_BINS)],
        # EXACT counts, accumulated rather than read off the bins
        "tm_at_or_below": {str(t): 0 for t in TM_THRESHOLDS},
        "tm_at_or_above": {str(t): 0 for t in TM_THRESHOLDS},
        # "confidently novel": pLDDT >= P and TM <= T (no_hit counts as TM = 0)
        "novel_and_confident": {f"tm<={t}|plddt>={p}": 0
                                for t in TM_THRESHOLDS for p in PLDDT_THRESHOLDS},
        # Welford
        "tm_n": 0, "tm_mean": 0.0, "tm_m2": 0.0,
        "reservoir": [], "_res_seen": 0,
    }
    rng = random.Random(0)

    with open(ledger_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            agg["n_total"] += 1
            plddt = _f(row.get("mean_plddt"))
            if plddt is not None:
                agg["n_with_conf"] += 1
                b = _bin(plddt, 0.0, 100.0, PLDDT_BINS)
                if b is not None:
                    agg["plddt_hist"][b] += 1

            if row.get("status") == "no_hit":
                agg["n_no_hit"] += 1
                # A no-hit ORF is not missing data -- it is the strongest novelty observation there
                # is, so it enters the distribution at TM = 0 rather than being skipped.
                tm = 0.0
            else:
                tm = _f(row.get("tm_score"))
                if tm is None:
                    continue
            agg["n_scored"] += 1

            tb = _bin(tm, 0.0, 1.0, TM_BINS)
            if tb is not None:
                agg["tm_hist"][tb] += 1
                pb = _bin(plddt, 0.0, 100.0, PLDDT_BINS) if plddt is not None else None
                if pb is not None:
                    agg["joint_hist"][tb][pb] += 1

            for t in TM_THRESHOLDS:
                if tm <= t:
                    agg["tm_at_or_below"][str(t)] += 1
                if tm >= t:
                    agg["tm_at_or_above"][str(t)] += 1
                if plddt is not None:
                    for p in PLDDT_THRESHOLDS:
                        if tm <= t and plddt >= p:
                            agg["novel_and_confident"][f"tm<={t}|plddt>={p}"] += 1

            agg["tm_n"] += 1
            d = tm - agg["tm_mean"]
            agg["tm_mean"] += d / agg["tm_n"]
            agg["tm_m2"] += d * (tm - agg["tm_mean"])

            agg["_res_seen"] += 1
            pt = [row.get("orfid"), tm, plddt]
            if len(agg["reservoir"]) < RESERVOIR_N:
                agg["reservoir"].append(pt)
            else:
                j = rng.randrange(agg["_res_seen"])
                if j < RESERVOIR_N:
                    agg["reservoir"][j] = pt

    agg["tm_sd"] = math.sqrt(agg["tm_m2"] / (agg["tm_n"] - 1)) if agg["tm_n"] > 1 else 0.0
    log(f"[agg] {agg['n_scored']} scored, mean {score_field} = {agg['tm_mean']:.4f} "
        f"+/- {agg['tm_sd']:.4f}, {agg['n_no_hit']} no_hit")
    return agg


def merge_aggregates(paths, log=print):
    """Bin-for-bin union of shard aggregates. Histograms and exact counts add; Welford combines by
    the parallel form; reservoirs are concatenated and resampled down."""
    out = None
    for p in paths:
        with open(p) as fh:
            a = json.load(fh)
        if out is None:
            out = a
            continue
        for k in ("n_total", "n_no_hit", "n_scored", "n_with_conf"):
            out[k] += a[k]
        for k in ("tm_hist", "plddt_hist"):
            out[k] = [x + y for x, y in zip(out[k], a[k])]
        out["joint_hist"] = [[x + y for x, y in zip(r1, r2)]
                             for r1, r2 in zip(out["joint_hist"], a["joint_hist"])]
        for k in ("tm_at_or_below", "tm_at_or_above", "novel_and_confident"):
            for kk in out[k]:
                out[k][kk] += a[k].get(kk, 0)
        na, nb = out["tm_n"], a["tm_n"]
        if nb:
            delta = a["tm_mean"] - out["tm_mean"]
            tot = na + nb
            out["tm_m2"] = out["tm_m2"] + a["tm_m2"] + delta * delta * na * nb / tot
            out["tm_mean"] = (out["tm_mean"] * na + a["tm_mean"] * nb) / tot
            out["tm_n"] = tot
        out["reservoir"] = out["reservoir"] + a["reservoir"]
        out["_res_seen"] = out.get("_res_seen", 0) + a.get("_res_seen", 0)

    if out is None:
        raise SystemExit("nothing to merge")
    if len(out["reservoir"]) > RESERVOIR_N:
        out["reservoir"] = random.Random(0).sample(out["reservoir"], RESERVOIR_N)
    out["tm_sd"] = math.sqrt(out["tm_m2"] / (out["tm_n"] - 1)) if out["tm_n"] > 1 else 0.0
    log(f"[merge] {len(paths)} shards -> {out['n_scored']} scored, mean {out['tm_mean']:.4f}")
    return out


# ---------------------------------------------------------------------------
# mirroring results back to S3
# ---------------------------------------------------------------------------
def s3_put(local, s3_uri, log=print):
    import boto3
    bucket, _, key = s3_uri[len("s3://"):].partition("/")
    key = key.rstrip("/") + "/" + os.path.basename(local)
    boto3.client("s3").upload_file(local, bucket, key)
    log(f"[s3] {local} -> s3://{bucket}/{key}")


# ---------------------------------------------------------------------------
# --check-env: settle it by measurement before committing to a batch
# ---------------------------------------------------------------------------
def check_env(args, log=print):
    """Verify the binary, the target DB, and a real end-to-end search on ONE structure.

    Same role --check-env played for MolProbity's engine question: the cheap experiment that answers
    'will the batch work, and does the output mean what I think' before hours are spent.
    """
    ok = True

    ver = foldseek_version(args.foldseek)
    log(f"[check] foldseek: {ver or 'NOT FOUND ON PATH'}")
    ok &= ver is not None

    tgt = os.path.expanduser(args.target)
    have_tgt = os.path.exists(tgt) or os.path.exists(tgt + ".dbtype")
    log(f"[check] target db {tgt}: {'ok' if have_tgt else 'MISSING -- run: foldseek databases PDB %s tmp' % tgt}")
    ok &= have_tgt

    try:
        import boto3  # noqa: F401
        log("[check] boto3: ok")
    except ImportError:
        log("[check] boto3: MISSING (pip install boto3) -- only needed for s3:// runs")
        ok &= not str(args.run).startswith("s3://")

    if not ok:
        log("[check] FAILED -- fix the above before running a batch")
        return 1

    # one structure, all four stages, in a throwaway dir
    probe = os.path.join(args.out, "_check")
    shutil.rmtree(probe, ignore_errors=True)
    os.makedirs(probe, exist_ok=True)
    src = StructSource(args.run)
    ids = stage_fetch(src, os.path.join(probe, "structures"), None, 4, limit=1, log=log)
    if not ids:
        log("[check] FAILED -- no structures found under the run prefix")
        return 1
    log(f"[check] probe ORFid = {ids[0]}")

    stage_createdb(os.path.join(probe, "structures"), os.path.join(probe, "qdb"), args.foldseek, log=log)
    stage_search(os.path.join(probe, "qdb"), tgt, os.path.join(probe, "aln"),
                 os.path.join(probe, "tmp"), args.foldseek, args.threads, args.evalue,
                 args.max_seqs, args.alignment_type, False, args.foldseek_arg, log=log)
    aln_tsv = os.path.join(probe, "aln.tsv")
    stage_convert(os.path.join(probe, "qdb"), tgt, os.path.join(probe, "aln"), aln_tsv,
                  args.foldseek, args.threads, log=log)

    best = best_hits(aln_tsv, args.score_field, log=log)
    if not best:
        log("[check] search returned NO HITS for the probe structure.")
        log("[check]   Not necessarily broken -- but on the s68_l10 stand-in it probably IS, since")
        log("[check]   those ORFs were selected for HAVING a PDB entry and should self-hit at TM~1.")
        return 1
    rec = list(best.values())[0]
    t = rec["by_tm"]
    log(f"[check] best hit: {t['target']}  qtm={t['qtmscore']}  ttm={t['ttmscore']}  "
        f"alntm={t['alntmscore']}  e={t['evalue']}  lddt={t['lddt']}")
    log("[check] ^ compare the three TM columns: if alntmscore == ttmscore, use qtmscore "
        "(see issue #312 note in the docstring).")
    log("[check] OK")
    return 0


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Figure Db: Foldseek best-hit novelty for a fold run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--run", help="s3://bucket/prefix/ or local dir written by the fold predictor")
    ap.add_argument("--target", default="~/db/pdb",
                    help="foldseek target DB (build with: foldseek databases PDB ~/db/pdb ~/tmp)")
    ap.add_argument("--out", default="nov_out", help="output directory")
    ap.add_argument("--conf-tsv", help="Da's conf.tsv, joined in for the confidence axis")

    ap.add_argument("--stage", default="all",
                    choices=["all", "fetch", "createdb", "search", "reduce"],
                    help="run one stage only; every stage skips work that already exists")
    ap.add_argument("--shard", help="K/N -- this worker owns ORFids hashing to shard K of N")
    ap.add_argument("--limit", type=int, help="only the first N structures (smoke tests)")

    ap.add_argument("--foldseek", default="foldseek", help="path to the foldseek binary")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 8)
    ap.add_argument("--dl-workers", type=int, default=32, help="S3 download concurrency")
    ap.add_argument("--evalue", default="10",
                    help="PERMISSIVE by default (foldseek's own default is 1e-3): for a novelty "
                         "claim a missed weak hit is a false novelty, which is the costly direction")
    ap.add_argument("--max-seqs", default="1000", help="prefilter depth per query")
    ap.add_argument("--alignment-type", default="2",
                    help="2 = 3Di+AA (default, fast, still yields real TM columns); "
                         "1 = TMalign (global, slow)")
    ap.add_argument("--exhaustive", action="store_true",
                    help="disable the k-mer prefilter. Orders of magnitude slower -- use it to "
                         "VALIDATE a no-hit subset, never for the full set")
    ap.add_argument("--score-field", default="qtmscore", choices=TM_FIELDS,
                    help="which TM normalisation is 'the' novelty score; query-normalised is the "
                         "defensible one for 'does anything resemble MY protein'")
    ap.add_argument("--foldseek-arg", action="append", default=[],
                    help="extra raw argument passed through to foldseek search (repeatable)")

    ap.add_argument("--s3-out", help="s3://bucket/prefix/ to mirror ledger + aggregates to")
    ap.add_argument("--check-env", action="store_true", help="verify the VM on ONE structure, then exit")
    ap.add_argument("--merge", nargs="+", help="merge these aggregates*.json into aggregates.merged.json")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    logf = open(os.path.join(args.out, "run.log"), "a")

    def log(*a):
        msg = " ".join(str(x) for x in a)
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    if args.merge:
        merged = merge_aggregates(args.merge, log=log)
        p = os.path.join(args.out, "aggregates.merged.json")
        with open(p, "w") as fh:
            json.dump(merged, fh)
        log(f"[merge] wrote {p}")
        if args.s3_out:
            s3_put(p, args.s3_out, log=log)
        return 0

    if not args.run:
        ap.error("--run is required (except with --merge)")
    if args.check_env:
        return check_env(args, log=log)

    shard = None
    tag = ""
    if args.shard:
        k, n = args.shard.split("/")
        shard = (int(k), int(n))
        tag = f".shard{k}of{n}"
        log(f"[main] shard {k}/{n}")

    tgt = os.path.expanduser(args.target)
    struct_dir = os.path.join(args.out, f"structures{tag}")
    query_db = os.path.join(args.out, f"qdb{tag}")
    aln_db = os.path.join(args.out, f"aln{tag}")
    aln_tsv = os.path.join(args.out, f"aln{tag}.tsv")
    tmp_dir = os.path.join(args.out, f"tmp{tag}")
    ledger = os.path.join(args.out, f"novelty{tag}.tsv")
    agg_path = os.path.join(args.out, f"aggregates{tag}.json")
    ids_path = os.path.join(args.out, f"owned_ids{tag}.txt")

    want = ("fetch", "createdb", "search", "reduce") if args.stage == "all" else (args.stage,)
    t_start = time.time()

    if "fetch" in want:
        ids = stage_fetch(StructSource(args.run), struct_dir, shard, args.dl_workers,
                          limit=args.limit, log=log)
        with open(ids_path, "w") as fh:
            fh.write("\n".join(ids))
    if "createdb" in want:
        stage_createdb(struct_dir, query_db, args.foldseek, log=log)
    if "search" in want:
        stage_search(query_db, tgt, aln_db, tmp_dir, args.foldseek, args.threads,
                     args.evalue, args.max_seqs, args.alignment_type, args.exhaustive,
                     args.foldseek_arg, log=log)
        stage_convert(query_db, tgt, aln_db, aln_tsv, args.foldseek, args.threads, log=log)
    if "reduce" in want:
        # The owned-id list is the denominator and comes from the FETCH listing, not from the
        # alignments -- reading it back off the alignments would define away every no-hit ORF.
        if os.path.exists(ids_path):
            with open(ids_path) as fh:
                owned = [l.strip() for l in fh if l.strip()]
        else:
            owned = sorted({orfid_from_name(n) for n in os.listdir(struct_dir) if is_structure(n)} - {None})
        conf = read_conf_tsv(args.conf_tsv, log=log)
        best = best_hits(aln_tsv, args.score_field, log=log)
        write_ledger(ledger, owned, best, args.score_field, conf, log=log)
        agg = build_aggregates(ledger, args.score_field, tgt, args.run, log=log)
        agg["wall_s"] = round(time.time() - t_start, 1)
        agg["n_owned"] = len(owned)
        with open(agg_path, "w") as fh:
            json.dump(agg, fh)
        log(f"[main] wrote {ledger} and {agg_path}")
        if args.s3_out:
            s3_put(ledger, args.s3_out, log=log)
            s3_put(agg_path, args.s3_out, log=log)

    log(f"[main] done in {time.time()-t_start:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
