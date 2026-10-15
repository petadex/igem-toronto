"""
fetch_clusters.py -- Stage 0, step 1:  90pid cluster id(s)  ->  per-cluster ORF FASTA.

Issue #84.  Self-contained; imports nothing from this repo.

WHY THIS IS ONE SCRIPT AND NOT TWO
----------------------------------
Getting the ORF ids is a cheap SQL query.  Getting the SEQUENCES means a single
streaming pass over s3://petadex/logan/petadex.catalytic_orfs.v1.1.fa, which is
~100 GB.  That pass is the dominant cost of the whole pipeline and it costs the
same for 1 cluster as for 50 -- so the ids for EVERY cluster you want must be
collected BEFORE the pass, and the output partitioned AFTER it.  Splitting this
into "one run per cluster" would multiply the only expensive step by N.

    phase 1  schema   introspect real column names (they are inconsistent -- see below)
    phase 2  ids      SQL -> ids.txt (union over all clusters) + members.tsv
    phase 3  extract  ONE seqkit pass over the ORF corpus -> all.fa
    phase 4  anchors  ONE awk pass over the CORES corpus -> anchors.tsv
    phase 5  split    all.fa -> clusters/<cid>/cluster.fa + centroid.fa + orf_meta.tsv

Phases 3 and 4 are separate passes over two different S3 objects: the ORF corpus has the
sequences, the cores corpus has the HMM hit regions.  Both are streamed, never stored.

Run on the folding VM, not locally: RDS, the ORF corpus and the structure bucket
are all in us-east-1, and the 100 GB pass over the public internet is painful.

A NOTE ON IDENTIFIERS
---------------------
No underscore after "90", anywhere:  table `90pid_enzyme_clusters`, column
`90pid_enzyme_id` (the latter in both that table and petadex_clustering).  Both
START WITH A DIGIT, which is illegal in Postgres unless they were created
double-quoted -- in which case they must stay quoted forever.  So every
identifier here goes through quote_ident().  Names are still resolved against
information_schema rather than hardcoded, because a mismatch then reports the
real names instead of surfacing as a bare 42703/42P01 mid-query; the `schema`
phase prints what it found.  Run that first.

USAGE
-----
    python fetch_clusters.py schema
    python fetch_clusters.py ids 13 27 41 --out-dir out/
    python fetch_clusters.py extract --out-dir out/ --run
    python fetch_clusters.py split --out-dir out/
"""

import argparse
import os
import re
import subprocess
import sys

# --------------------------------------------------------------------------- #
# config
# --------------------------------------------------------------------------- #

DB = {
    "host":     "petadex.ccz9y6yshbls.us-east-1.rds.amazonaws.com",
    "port":     5432,
    "database": "petadex",
    "user":     "readonly_user",
    "password": "petadex",
}

ORF_CORPUS = "s3://petadex/logan/petadex.catalytic_orfs.v1.1.fa"
# Only feeds the pv progress bar.  Measured 2026-08-23 via head_object:
# 116,555,393,240 bytes = 116.56 GB = 108.55 GiB.  The "100g" in the issue-82
# notebook was a guess and makes pv run past 100% with a meaningless ETA.
CORPUS_SIZE_HINT = "117g"

# The CORES file, not the ORF file.  Its headers carry the PETadex HMM hit region as a
# `/start-end` suffix -- `>1|WP_054022242.1|||||/65-261` -- which is the anchor
# domain_from_pae.py needs to decide WHICH parsed domain is the enzyme.  The ORF headers
# carry contig coordinates instead, so this is a separate pass.  13.91 GB measured
# 2026-08-23, about 1/8 of the ORF corpus.
CORE_CORPUS = "s3://petadex/logan/petadex.complete_catalytic_cores.v1.1.fa.zst"
CORE_SIZE_HINT = "14g"

MEMBER_TABLE = "petadex_clustering"

# Confirmed names first, then a loose fallback.  Still resolved against
# information_schema rather than hardcoded, so a mismatch reports the real names
# instead of surfacing as a bare 42703/42P01 mid-query.
#   table  90pid_enzyme_clusters  (90pid_enzyme_id, centroid_orf_id, date_clustered)
#   table  petadex_clustering     (90pid_enzyme_id, orf_id, ...)
CENTROID_TABLE_PATTERNS = [r"^90pid_enzyme_clusters$", r"90.*pid.*enzyme.*cluster"]
CLUSTER_COL_PATTERNS    = [r"^90pid_enzyme_id$", r"90.*pid.*enzyme.*id"]
ORFID_COL_PATTERNS      = [r"^orf_id$", r"orf.*id"]
CENTROID_COL_PATTERNS   = [r"^centroid_orf_id$", r"centroid.*orf.*id"]


# --------------------------------------------------------------------------- #
# db helpers
# --------------------------------------------------------------------------- #

def connect():
    try:
        import psycopg2
    except ImportError:
        sys.exit("psycopg2 not installed.  pip install psycopg2-binary")
    return psycopg2.connect(**DB)


def quote_ident(name):
    """Postgres-safe identifier.  Mandatory here: a name starting with a digit
    is only legal double-quoted."""
    return '"' + name.replace('"', '""') + '"'


def columns_of(cur, table):
    cur.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name = %s ORDER BY ordinal_position", (table,))
    return [r[0] for r in cur.fetchall()]


def resolve(cols, patterns, table, what):
    """Pick the column matching the earliest pattern that hits.  Fails loudly with
    the real column list rather than letting a wrong guess raise 42703 later."""
    for pat in patterns:
        for c in cols:
            if re.search(pat, c, re.I):
                return c
    sys.exit("could not find the %s column in %s.\n"
             "  columns present: %s\n"
             "  patterns tried : %s\n"
             "  pass it explicitly with the matching --*-col flag."
             % (what, table, cols, patterns))


def resolve_table(cur, patterns, override=None):
    """Find the real table name.  Returns None if nothing matches -- the centroid
    table is optional (only Plan 2 needs it), so callers decide how to react."""
    if override:
        return override
    cur.execute("SELECT table_name FROM information_schema.tables "
                "WHERE table_schema NOT IN ('pg_catalog', 'information_schema')")
    names = [r[0] for r in cur.fetchall()]
    for pat in patterns:
        for t in names:
            if re.search(pat, t, re.I):
                return t
    return None


# --------------------------------------------------------------------------- #
# phase 1 -- schema
# --------------------------------------------------------------------------- #

def phase_schema(args):
    with connect() as conn, conn.cursor() as cur:
        ctable = resolve_table(cur, CENTROID_TABLE_PATTERNS, args.centroid_table)
        print("resolved centroid table:", ctable or "NOT FOUND")
        print()
        for t in [MEMBER_TABLE] + ([ctable] if ctable else []):
            cols = columns_of(cur, t)
            print("%-28s %s" % (t, cols if cols else "NOT FOUND"))
        print()
        mcols = columns_of(cur, MEMBER_TABLE)
        if mcols:
            print("resolved cluster col :",
                  resolve(mcols, CLUSTER_COL_PATTERNS, MEMBER_TABLE, "cluster-id"))
            print("resolved orf_id col  :",
                  resolve(mcols, ORFID_COL_PATTERNS, MEMBER_TABLE, "orf-id"))
        if ctable:
            ccols = columns_of(cur, ctable)
            print("resolved centroid col:",
                  resolve(ccols, CENTROID_COL_PATTERNS, ctable, "centroid"))


# --------------------------------------------------------------------------- #
# phase 2 -- ids
# --------------------------------------------------------------------------- #

def phase_ids(args):
    cids = [int(c) for c in args.cluster_ids]
    os.makedirs(args.out_dir, exist_ok=True)

    with connect() as conn, conn.cursor() as cur:
        mcols = columns_of(cur, MEMBER_TABLE)
        cluster_col = args.cluster_col or resolve(
            mcols, CLUSTER_COL_PATTERNS, MEMBER_TABLE, "cluster-id")
        orf_col = args.orf_col or resolve(
            mcols, ORFID_COL_PATTERNS, MEMBER_TABLE, "orf-id")

        cur.execute(
            "SELECT %s, %s FROM %s WHERE %s = ANY(%%s)"
            % (quote_ident(cluster_col), quote_ident(orf_col),
               quote_ident(MEMBER_TABLE), quote_ident(cluster_col)),
            (cids,))
        members = cur.fetchall()

        centroids = {}
        ctable = resolve_table(cur, CENTROID_TABLE_PATTERNS, args.centroid_table)
        if ctable:
            ccols = columns_of(cur, ctable)
            ccluster = args.cluster_col_centroid or resolve(
                ccols, CLUSTER_COL_PATTERNS, ctable, "cluster-id")
            ccent = args.centroid_col or resolve(
                ccols, CENTROID_COL_PATTERNS, ctable, "centroid")
            cur.execute(
                "SELECT %s, %s FROM %s WHERE %s = ANY(%%s)"
                % (quote_ident(ccluster), quote_ident(ccent),
                   quote_ident(ctable), quote_ident(ccluster)),
                (cids,))
            centroids = {int(a): int(b) for a, b in cur.fetchall() if b is not None}
        else:
            print("WARNING: centroid table not found (patterns: %s) -- no centroids "
                  "written; preprocessing Plan 2 needs them.  Pass --centroid-table."
                  % CENTROID_TABLE_PATTERNS, file=sys.stderr)

    if not members:
        sys.exit("no members found for cluster ids %s" % cids)

    mpath = os.path.join(args.out_dir, "members.tsv")
    with open(mpath, "w") as fh:
        fh.write("cluster_id\torf_id\tis_centroid\n")
        for cid, orf in sorted(members):
            is_cent = int(centroids.get(int(cid)) == int(orf))
            fh.write("%s\t%s\t%d\n" % (cid, orf, is_cent))

    ids = sorted({int(o) for _, o in members})
    ipath = os.path.join(args.out_dir, "ids.txt")
    with open(ipath, "w") as fh:
        fh.write("\n".join(str(i) for i in ids) + "\n")

    per = {}
    for cid, _ in members:
        per[int(cid)] = per.get(int(cid), 0) + 1
    print("%d cluster(s), %d unique ORF ids" % (len(cids), len(ids)))
    for cid in sorted(per):
        mark = ("  centroid orf%d" % centroids[cid]) if cid in centroids else "  (no centroid)"
        print("  cluster %-8s %6d ORFs%s" % (cid, per[cid], mark))
    missing = [c for c in cids if c not in per]
    if missing:
        print("  WARNING: no rows for cluster id(s) %s" % missing, file=sys.stderr)
    print("\nwrote %s\nwrote %s" % (ipath, mpath))


# --------------------------------------------------------------------------- #
# phase 3 -- extract  (the one expensive pass)
# --------------------------------------------------------------------------- #

def phase_extract(args):
    ids_path = os.path.join(args.out_dir, "ids.txt")
    out_fa = os.path.join(args.out_dir, "all.fa")
    if not os.path.exists(ids_path):
        sys.exit("%s not found -- run the `ids` phase first" % ids_path)

    # pv is only a progress bar -- drop it rather than fail if it isn't installed.
    import shutil
    pv = " | pv -s %s" % CORPUS_SIZE_HINT if shutil.which("pv") else ""
    if not pv:
        print("# note: pv not found -- running without a progress bar", file=sys.stderr)
    if not shutil.which("seqkit"):
        sys.exit("seqkit not found on PATH -- it is required.  Install the static "
                 "binary from github.com/shenwei356/seqkit/releases (no apt needed).")

    cmd = ("aws s3 cp %s - --no-sign-request"
           "%s"
           " | seqkit grep --id-regexp '^([^|]+)' -f %s -o %s"
           % (ORF_CORPUS, pv, ids_path, out_fa))

    n = sum(1 for _ in open(ids_path))
    print("# one streaming pass over the corpus for all %d ORF ids\n%s\n" % (n, cmd))
    if not args.run:
        print("(dry run -- pass --run to execute)")
        return
    rc = subprocess.call(["bash", "-o", "pipefail", "-c", cmd])
    if rc != 0:
        sys.exit("extract failed (exit %d)" % rc)
    print("wrote %s" % out_fa)


# --------------------------------------------------------------------------- #
# phase 3b -- anchors  (the HMM hit region, from the CORES corpus)
# --------------------------------------------------------------------------- #

def phase_anchors(args):
    """Stream the cores corpus and keep only the header of each wanted ORF, parsing the
    `/start-end` suffix into an anchors TSV.

    Headers only -- no seqkit, no sequences, output is a few KB.  awk does the whole job:
    ids.txt is read first (NR==FNR) into a set, then every '>' line whose field 1 is in
    that set contributes one row."""
    import shutil
    ids_path = os.path.join(args.out_dir, "ids.txt")
    out_tsv = os.path.join(args.out_dir, "anchors.tsv")
    if not os.path.exists(ids_path):
        sys.exit("%s not found -- run the `ids` phase first" % ids_path)
    if not shutil.which("zstd"):
        sys.exit("zstd not found on PATH -- the cores corpus is zstd-compressed and "
                 "seqkit/awk cannot read it directly.  `sudo apt-get install -y zstd`.")

    pv = " | pv -s %s" % CORE_SIZE_HINT if shutil.which("pv") else ""
    awk = (r"""awk -F'|' '"""
           r"""NR==FNR { want[$1]=1; next } """
           r"""/^>/ { id=substr($1,2); """
           r"""       if (id in want) { r=$NF; sub(/^.*\//,"",r); """
           r"""                         if (r ~ /^[0-9]+-[0-9]+$/) print id "\t" r } }' """
           + ids_path + " -")
    cmd = ("{ printf 'orfid\\tcore_range\\n'; "
           "aws s3 cp %s - --no-sign-request%s | zstd -dc | %s ; } > %s"
           % (CORE_CORPUS, pv, awk, out_tsv))

    n = sum(1 for _ in open(ids_path))
    print("# one streaming pass over the CORES corpus for all %d ORF ids\n%s\n"
          % (n, cmd))
    if not args.run:
        print("(dry run -- pass --run to execute)")
        return
    rc = subprocess.call(["bash", "-o", "pipefail", "-c", cmd])
    if rc != 0:
        sys.exit("anchors failed (exit %d)" % rc)

    got = max(0, sum(1 for _ in open(out_tsv)) - 1)
    print("wrote %s  (%d/%d ORFs have an anchor)" % (out_tsv, got, n))
    if got < n:
        print("  NOTE: %d ORF(s) had no core record; domain_from_pae.py will report "
              "them as no_anchor and skip them." % (n - got), file=sys.stderr)


# --------------------------------------------------------------------------- #
# phase 4 -- split + header metadata
# --------------------------------------------------------------------------- #

# >{orf_id}|{accession}|{run_accession}|{contig_id}|{start_nt}|{end_nt}|{strand}
# Verified on ispetase_family.fa: (end_nt - start_nt) / 3 == len(aa), exactly.
# start_nt == 0 means the ORF begins at the contig's first base, i.e. it may be
# 5'-truncated by the assembly.  Carried through here because recovering it later
# would cost another 100 GB pass; whether it becomes the truncation filter is a
# decision for assemble_cores.py.
HEADER_FIELDS = ["orf_id", "accession", "run_accession",
                 "contig_id", "start_nt", "end_nt", "strand"]


def parse_header(h):
    f = h.split("|")
    f += [""] * (len(HEADER_FIELDS) - len(f))
    return dict(zip(HEADER_FIELDS, f[:len(HEADER_FIELDS)]))


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def phase_split(args):
    all_fa = os.path.join(args.out_dir, "all.fa")
    mpath = os.path.join(args.out_dir, "members.tsv")
    for p in (all_fa, mpath):
        if not os.path.exists(p):
            sys.exit("%s not found -- run the earlier phases first" % p)

    orf2cluster, centroid_of = {}, {}
    with open(mpath) as fh:
        next(fh)
        for line in fh:
            cid, orf, is_cent = line.rstrip("\n").split("\t")
            orf2cluster[orf] = cid
            if is_cent == "1":
                centroid_of[cid] = orf

    handles, counts, meta_rows, seen = {}, {}, [], set()
    cdir = os.path.join(args.out_dir, "clusters")

    for header, seq in read_fasta(all_fa):
        rec = parse_header(header)
        orf = rec["orf_id"]
        seen.add(orf)
        start = rec["start_nt"]
        meta_rows.append([rec[k] for k in HEADER_FIELDS] +
                         [str(len(seq)),
                          "1" if start.isdigit() and int(start) == 0 else "0"])

        cid = orf2cluster.get(orf)
        if cid is None:
            continue
        if cid not in handles:
            os.makedirs(os.path.join(cdir, cid), exist_ok=True)
            handles[cid] = (open(os.path.join(cdir, cid, "cluster.fa"), "w"),
                            open(os.path.join(cdir, cid, "cluster.fold.fa"), "w"))
            counts[cid] = 0
        # cluster.fa keeps the full header (provenance); cluster.fold.fa uses the BARE
        # ORF id, because the folder names its outputs after the header -- a full header
        # becomes `1219585044__SRR6391592_110893_1_709_1.json`, which no longer matches
        # anything keyed by orf_id downstream.
        handles[cid][0].write(">%s\n%s\n" % (header, seq))
        handles[cid][1].write(">orf%s\n%s\n" % (orf, seq))
        counts[cid] += 1
        if centroid_of.get(cid) == orf:
            with open(os.path.join(cdir, cid, "centroid.fa"), "w") as fh:
                fh.write(">%s\n%s\n" % (header, seq))
            with open(os.path.join(cdir, cid, "centroid.fold.fa"), "w") as fh:
                fh.write(">orf%s\n%s\n" % (orf, seq))

    for pair in handles.values():
        for fh in pair:
            fh.close()

    meta_path = os.path.join(args.out_dir, "orf_meta.tsv")
    with open(meta_path, "w") as fh:
        fh.write("\t".join(HEADER_FIELDS + ["aa_len", "contig_edge_5p"]) + "\n")
        for row in meta_rows:
            fh.write("\t".join(row) + "\n")

    for cid in sorted(counts):
        cent = "yes" if os.path.exists(os.path.join(cdir, cid, "centroid.fa")) else "NO"
        print("  cluster %-8s %6d seqs   centroid: %s" % (cid, counts[cid], cent))
    print("\nwrote %s" % meta_path)

    missing = set(orf2cluster) - seen
    if missing:
        print("WARNING: %d requested ORF ids absent from the corpus (first few: %s)"
              % (len(missing), sorted(missing)[:5]), file=sys.stderr)


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="phase", required=True)

    s = sub.add_parser("schema", help="print real table/column names and exit")
    s.add_argument("--centroid-table")
    s.set_defaults(func=phase_schema)

    s = sub.add_parser("ids", help="SQL -> ids.txt + members.tsv")
    s.add_argument("cluster_ids", nargs="+",
                   help="90pid cluster id(s) -- pass ALL you want; one pass serves all")
    s.add_argument("--out-dir", default="out")
    s.add_argument("--cluster-col")
    s.add_argument("--orf-col")
    s.add_argument("--centroid-table")
    s.add_argument("--cluster-col-centroid")
    s.add_argument("--centroid-col")
    s.set_defaults(func=phase_ids)

    s = sub.add_parser("extract", help="the single ~100 GB streaming pass")
    s.add_argument("--out-dir", default="out")
    s.add_argument("--run", action="store_true",
                   help="actually run it (default: print the command)")
    s.set_defaults(func=phase_extract)

    s = sub.add_parser("anchors", help="cores corpus -> anchors.tsv (HMM hit regions)")
    s.add_argument("--out-dir", default="out")
    s.add_argument("--run", action="store_true",
                   help="actually run it (default: print the command)")
    s.set_defaults(func=phase_anchors)

    s = sub.add_parser("split", help="all.fa -> per-cluster FASTA + orf_meta.tsv")
    s.add_argument("--out-dir", default="out")
    s.set_defaults(func=phase_split)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
