"""
Reverse-map PDB ground-truth structures to the ORFids they cover, using the
catalytic-cores file.

Context (Issue #6, ground-truth Note 1): everything we fold comes from the ORFs
file, but a subset of cores carry an accession that resolves to a solved PDB
structure. Field 1 of both the cores file and the ORFs file is the *same* ORFid,
so once we know which PDB a core corresponds to we have PDB -> ORFid directly.

Two paths produce a PDB for a core (see --scope in the plan):

  direct       cores field 2 is literally a PDB accession (e.g. 5YNS_A). These
               are the "other" bucket from split_accession_types.py. Highest
               confidence: the PDB was annotated to that exact sequence. Chain
               is preserved.

  via_uniprot  cores field 2 is GenBank/RefSeq/UniProt; we map it to a UniProt
               accession (reusing the existing *_to_uniprot.tsv tables), then map
               UniProt -> PDB via the idmapping_selected.tab.gz PDB cross-ref
               column (col 5), whose entries are 'PDBID:CHAIN' so chain is kept.
               No residue ranges come from idmapping; the sequence-identity
               caveat lives in the `provenance` column.

Deliberately out of scope: leakage / 2023-cutoff / FoldSeek filtering, and any
join against the ORFs file or the folded predictions. Those run downstream on the
output of this script.

Outputs (into this script's directory):
  pdb_orfid_map.tsv   long form, one row per (core, pdb)
  pdb_to_orfids.json  grouped pdb_id -> [orfids]
  + a console summary with a reconciliation check against other_accessions.csv

Rerunnable: staging files and the extracted uniprot_to_pdb.tsv are cached and
skipped unless --force is passed.
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path

import polars as pl
import zstandard as zstd

# ---------------------------------------------------------------------------
# Classification regexes -- copied verbatim from
# resources/260523_issue72_af_db_download/scripts/split_accession_types.py so the
# PDB ("other") bucket here matches exactly what was validated before.
# ---------------------------------------------------------------------------
UNIPROT_PAT = re.compile(r'^([A-N_R-Z][0-9][A-Z0-9]{3}[0-9]|[O,P,Q][0-9][A-Z0-9]{4})(?:[A-Z0-9]{4})?$')
REFSEQ_PAT = re.compile(r'^(WP|NP|XP|YP)_\d+$')
GENBANK_PAT = re.compile(r'^[A-Z]{3}\d{5,7}$')

# idmapping_selected.tab.gz column indices (0-based), per local_accession_mapping.py
IDMAP_UNIPROT_COL = 0
IDMAP_PDB_COL = 5


def classify(acc_raw):
    """Return (acc_type, acc_clean). acc_clean is version-stripped; for PDB the
    raw value (with chain) is returned so the chain can be recovered later."""
    clean = acc_raw.split('.')[0].strip()
    if UNIPROT_PAT.match(clean):
        return 'uniprot', clean
    if REFSEQ_PAT.match(clean):
        return 'refseq', clean
    if GENBANK_PAT.match(clean):
        return 'genbank', clean
    return 'pdb', acc_raw


def parse_header(header):
    """header includes the leading '>'. Cores headers look like
    >ORFid|accession|||||/65-261 -- return (orfid, acc_raw, core_range)."""
    parts = header[1:].split('|')
    orfid = parts[0].strip()
    acc_raw = parts[1].strip() if len(parts) > 1 else ''
    core_range = ''
    last = parts[-1].strip()
    if last.startswith('/'):
        core_range = last[1:]
    return orfid, acc_raw, core_range


def iter_core_headers(cores_path):
    """Stream only the '>' header lines out of the .zst cores file. Works on raw
    bytes and decodes header lines only, so the ~300M sequence lines are skipped
    cheaply."""
    dctx = zstd.ZstdDecompressor()
    with open(cores_path, 'rb') as fh:
        reader = dctx.stream_reader(fh)
        leftover = b''
        while True:
            chunk = reader.read(1 << 24)  # 16 MB
            if not chunk:
                break
            data = leftover + chunk
            lines = data.split(b'\n')
            leftover = lines.pop()  # possibly-incomplete final line
            for ln in lines:
                if ln[:1] == b'>':
                    yield ln.decode('utf-8', 'replace')
        if leftover[:1] == b'>':
            yield leftover.decode('utf-8', 'replace')


def stream_cores(cores_path, direct_path, nonpdb_path):
    """Pass 1: split accessioned cores into a direct-PDB staging file and a
    non-PDB (expansion candidate) staging file."""
    print(f"[*] Pass 1: streaming cores headers from {cores_path}")
    start = time.time()
    scanned = 0          # headers seen
    accessioned = 0      # headers with a non-empty accession
    n_direct = 0
    n_nonpdb = 0

    with open(direct_path, 'w', encoding='utf-8') as f_direct, \
         open(nonpdb_path, 'w', encoding='utf-8') as f_nonpdb:
        f_direct.write("orfid\tpdb_id\tchain\tpdb_raw\tcore_range\n")
        f_nonpdb.write("orfid\tacc_clean\tacc_type\tcore_range\n")

        for header in iter_core_headers(cores_path):
            scanned += 1
            if scanned % 20_000_000 == 0:
                print(f"    -> {scanned:,} headers scanned "
                      f"({accessioned:,} accessioned)...")

            orfid, acc_raw, core_range = parse_header(header)
            if not acc_raw:
                continue  # bare Logan/ESMAtlas ORF -- no ground truth possible
            accessioned += 1

            acc_type, acc_clean = classify(acc_raw)
            if acc_type == 'pdb':
                if '_' in acc_raw:
                    pdb_id, chain = acc_raw.split('_', 1)
                else:
                    pdb_id, chain = acc_raw, ''
                pdb_id = pdb_id.upper()
                f_direct.write(f"{orfid}\t{pdb_id}\t{chain}\t{acc_raw}\t{core_range}\n")
                n_direct += 1
            else:
                f_nonpdb.write(f"{orfid}\t{acc_clean}\t{acc_type}\t{core_range}\n")
                n_nonpdb += 1

    elapsed = time.time() - start
    print(f"[+] Pass 1 done in {elapsed:.0f}s: {scanned:,} headers, "
          f"{accessioned:,} accessioned "
          f"({n_direct:,} direct-PDB, {n_nonpdb:,} non-PDB candidates)")


def build_acc_to_uniprot(mappings_dir):
    """Combine the existing GenBank/RefSeq -> UniProt tables into one acc->uniprot
    frame (identity for native-UniProt accessions is handled at join time)."""
    frames = []
    for name in ("genbank_to_uniprot.tsv", "refseq_to_uniprot.tsv"):
        path = mappings_dir / name
        if not path.exists():
            print(f"[!] warning: {path} not found, skipping")
            continue
        df = pl.read_csv(path, separator='\t')
        # tables use columns From / Entry
        df = df.rename({df.columns[0]: "acc_clean", df.columns[1]: "uniprot"})
        frames.append(df.select(["acc_clean", "uniprot"]))
    if not frames:
        return pl.DataFrame({"acc_clean": [], "uniprot": []},
                            schema={"acc_clean": pl.Utf8, "uniprot": pl.Utf8})
    return pl.concat(frames).unique(subset=["acc_clean"])


def extract_uniprot_to_pdb(idmapping_path, out_path, uniprot_set, force):
    """One streaming pass over idmapping_selected.tab.gz: UniProt (col 0) -> PDB
    (col 5), filtered to the UniProt IDs our cores actually reach. Cached."""
    if out_path.exists() and not force:
        print(f"[*] Reusing cached {out_path.name} (pass --force to rebuild)")
        return pl.read_csv(out_path, separator='\t')

    print(f"[*] Extracting UniProt -> PDB from {idmapping_path} "
          f"(filtering to {len(uniprot_set):,} UniProt IDs)...")
    start = time.time()

    def namer(cols):
        return ['uniprot' if i == IDMAP_UNIPROT_COL
                else 'pdb' if i == IDMAP_PDB_COL
                else f'c{i}' for i in range(len(cols))]

    q = (
        pl.scan_csv(
            idmapping_path,
            separator='\t',
            has_header=False,
            with_column_names=namer,
            infer_schema_length=0,
        )
        .select(['uniprot', 'pdb'])
        .filter(pl.col('pdb').is_not_null() & (pl.col('pdb') != ''))
        .filter(pl.col('uniprot').is_in(list(uniprot_set)))
        .with_columns(pl.col('pdb').str.split(';'))
        .explode('pdb')
        # idmapping PDB entries are 'PDBID:CHAIN' (e.g. 1A0J:A); older forms use
        # 'PDBID_CHAIN'. Split the chain off the 4-char entry id either way.
        .with_columns(pl.col('pdb').str.strip_chars().str.replace(r'\.\d+$', ''))
        .with_columns([
            pl.col('pdb').str.replace(r'[:_].*$', '').str.to_uppercase().alias('pdb_id'),
            pl.col('pdb').str.extract(r'[:_]([A-Za-z0-9]+)$', 1).alias('chain'),
        ])
        .filter(pl.col('pdb_id') != '')
        .select(['uniprot', 'pdb_id', 'chain'])
        .unique()
    )
    u2p = q.collect(engine='streaming')
    u2p.write_csv(out_path, separator='\t')
    print(f"[+] UniProt -> PDB: {u2p.height:,} pairs in {time.time() - start:.0f}s "
          f"-> {out_path.name}")
    return u2p


COLUMNS = ["orfid", "pdb_id", "chain", "pdb_raw", "provenance",
           "uniprot", "source_acc", "source_acc_type", "core_range"]


def build_direct_df(direct_path):
    df = pl.read_csv(direct_path, separator='\t',
                     schema_overrides={"orfid": pl.Utf8, "chain": pl.Utf8})
    return df.select([
        pl.col("orfid"),
        pl.col("pdb_id"),
        pl.col("chain"),
        pl.col("pdb_raw"),
        pl.lit("direct").alias("provenance"),
        pl.lit(None, dtype=pl.Utf8).alias("uniprot"),
        pl.col("pdb_raw").alias("source_acc"),
        pl.lit("pdb").alias("source_acc_type"),
        pl.col("core_range"),
    ])


def build_via_uniprot_df(nonpdb_path, acc2uni, u2p):
    df = pl.read_csv(nonpdb_path, separator='\t', schema_overrides={"orfid": pl.Utf8})
    # acc_clean -> uniprot (identity for native uniprot accessions)
    df = df.join(acc2uni, on="acc_clean", how="left")
    df = df.with_columns(
        pl.when(pl.col("acc_type") == "uniprot")
        .then(pl.col("acc_clean"))
        .otherwise(pl.col("uniprot"))
        .alias("uniprot")
    ).filter(pl.col("uniprot").is_not_null())
    # uniprot -> pdb
    df = df.join(u2p, on="uniprot", how="inner")  # u2p: uniprot, pdb_id, chain
    return df.select([
        pl.col("orfid"),
        pl.col("pdb_id"),
        pl.col("chain"),
        pl.when(pl.col("chain").is_not_null())
          .then(pl.concat_str([pl.col("pdb_id"), pl.lit(":"), pl.col("chain")]))
          .otherwise(pl.col("pdb_id")).alias("pdb_raw"),
        pl.lit("via_uniprot").alias("provenance"),
        pl.col("uniprot"),
        pl.col("acc_clean").alias("source_acc"),
        pl.col("acc_type").alias("source_acc_type"),
        pl.col("core_range"),
    ])


def reconcile_direct(direct_df, other_accessions_path):
    """Sanity check: distinct direct PDB IDs vs the 'other' bucket that was
    validated to be all-PDB previously."""
    direct_ids = set(direct_df["pdb_id"].unique().to_list())
    if not other_accessions_path.exists():
        print(f"[!] {other_accessions_path} not found; skipping reconciliation")
        return
    other_ids = set()
    with open(other_accessions_path, 'r', encoding='utf-8') as f:
        for line in f:
            tok = line.strip()
            if not tok or tok.lower() in ("accession", "id", "entry", "from"):
                continue
            tok = tok.split('.')[0].split('_')[0].upper()
            if tok:
                other_ids.add(tok)
    print(f"    reconciliation: direct distinct PDBs={len(direct_ids):,}, "
          f"other_accessions.csv distinct PDBs={len(other_ids):,}, "
          f"overlap={len(direct_ids & other_ids):,}")
    only_other = other_ids - direct_ids
    if only_other:
        print(f"      note: {len(only_other):,} PDBs in other_accessions.csv not "
              f"seen as direct hits (expected if that file was built from a "
              f"different cores version)")


def main():
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]

    ap = argparse.ArgumentParser(description="Reverse-map PDBs to ORFids via the cores file.")
    ap.add_argument("--cores", type=Path,
                    default=repo_root / "petadex.complete_catalytic_cores.v1.1.fa.zst")
    ap.add_argument("--idmapping", type=Path,
                    default=repo_root / "idmapping_selected.tab.gz")
    ap.add_argument("--mappings-dir", type=Path, default=repo_root / "accessions" / "mappings")
    ap.add_argument("--accessions-dir", type=Path, default=repo_root / "accessions")
    ap.add_argument("--out-dir", type=Path, default=script_dir)
    ap.add_argument("--force", action="store_true",
                    help="rebuild cached staging files and uniprot_to_pdb.tsv")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    direct_path = args.out_dir / "direct_pdb_orfid.tsv"
    nonpdb_path = args.out_dir / "nonpdb_orfid_acc.tsv"
    u2p_path = args.out_dir / "uniprot_to_pdb.tsv"
    long_path = args.out_dir / "pdb_orfid_map.tsv"
    json_path = args.out_dir / "pdb_to_orfids.json"

    for required in (args.cores, args.idmapping):
        if not required.exists():
            sys.exit(f"[-] required input not found: {required}")

    # Pass 1 (cached) --------------------------------------------------------
    if direct_path.exists() and nonpdb_path.exists() and not args.force:
        print(f"[*] Reusing cached staging files (pass --force to rebuild)")
    else:
        stream_cores(args.cores, direct_path, nonpdb_path)

    # acc -> uniprot, then the set of UniProt IDs our cores reach -------------
    acc2uni = build_acc_to_uniprot(args.mappings_dir)
    nonpdb = pl.read_csv(nonpdb_path, separator='\t', schema_overrides={"orfid": pl.Utf8})
    reached = nonpdb.join(acc2uni, on="acc_clean", how="left").with_columns(
        pl.when(pl.col("acc_type") == "uniprot")
        .then(pl.col("acc_clean"))
        .otherwise(pl.col("uniprot"))
        .alias("uniprot")
    ).filter(pl.col("uniprot").is_not_null())
    uniprot_set = set(reached["uniprot"].unique().to_list())

    # uniprot -> pdb (cached) ------------------------------------------------
    u2p = extract_uniprot_to_pdb(args.idmapping, u2p_path, uniprot_set, args.force)

    # Assemble long form -----------------------------------------------------
    print("[*] Joining and assembling the PDB -> ORFid map...")
    direct_df = build_direct_df(direct_path)
    via_df = build_via_uniprot_df(nonpdb_path, acc2uni, u2p)
    full = pl.concat([direct_df.select(COLUMNS), via_df.select(COLUMNS)]).unique()
    full = full.sort(["pdb_id", "provenance", "orfid"])
    full.write_csv(long_path, separator='\t')

    # Grouped JSON -----------------------------------------------------------
    grouped = (
        full.group_by("pdb_id")
        .agg(pl.col("orfid").unique().sort().alias("orfids"))
        .sort("pdb_id")
    )
    mapping = {row["pdb_id"]: row["orfids"] for row in grouped.iter_rows(named=True)}
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(mapping, f, indent=2)

    # Summary ----------------------------------------------------------------
    n_direct = full.filter(pl.col("provenance") == "direct").height
    n_via = full.filter(pl.col("provenance") == "via_uniprot").height
    print("-" * 60)
    print(f"rows (core,pdb pairs):     {full.height:,}")
    print(f"  direct:                  {n_direct:,}")
    print(f"  via_uniprot:             {n_via:,}")
    print(f"distinct PDB ids:          {full['pdb_id'].n_unique():,}")
    print(f"distinct ORFids:           {full['orfid'].n_unique():,}")
    reconcile_direct(direct_df, args.accessions_dir / "other_accessions.csv")
    print("-" * 60)
    print(f"[+] wrote {long_path.name} and {json_path.name} to {args.out_dir}")


if __name__ == "__main__":
    main()
