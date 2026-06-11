"""
Download already-folded structures from the ESM Atlas (Biohub) for sequences in a FASTA file.

Verified facts about s3://esm-protein-atlas (public, --no-sign-request):
  - Folds live in two Lance datasets under v1/folds/:
        folds_atlas.lance :     6,600,755 structures (representative/dedup set)
        folds_1B.lance    : 1,095,530,880 structures (the full ~1.1B set)
  - Each row schema: header, protein_hash, sequence, ptm, mean_plddt,
                     per_residue_plddt, pae, structure_blob
  - protein_hash == md5(sequence).hexdigest()   (exact sequence, as written in FASTA)
  - Both datasets have a BTree scalar index on protein_hash -> indexed IN() lookups are fast.
  - structure_blob == brotli( msgpack-numpy( {sequence, atom37_positions, atom37_mask,
                       residue_index, confidence(=pLDDT), ...} ) )

Install:  pip install pylance pyarrow msgpack numpy brotli zstandard boto3 tqdm

Output model (content-addressed):
  - PDBs are keyed by protein_hash:  <outdir>/<hash>.pdb  (or s3://bucket/prefix/<hash>.pdb)
    One object per UNIQUE sequence, so the input's heavy ORF redundancy collapses.
  - --orf-map writes the full orf_id -> protein_hash index (every ORF record, streamed).
  - found_hashes.tsv lists which hashes actually have a structure (drives resume + the join).
  - "which ORFs have a structure?"  ==  orf_map rows whose protein_hash is in found_hashes.
"""
import os, io, time, hashlib
import numpy as np
import lance, brotli, msgpack
import zstandard as zstd
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from itertools import islice
try:
    from tqdm import tqdm
except ImportError:                                   # graceful no-op if tqdm is absent
    class tqdm:
        def __init__(self, iterable=None, **kw): self.iterable = iterable
        def __iter__(self): return iter(self.iterable or [])
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def update(self, n=1): pass
        def set_postfix(self, **kw): pass
        @staticmethod
        def write(msg, *a, **k): print(msg)

STORAGE = {"aws_skip_signature": "true"}          # anonymous public-bucket access
DATASETS = {                                       # short name -> Lance dataset URI
    "atlas": "s3://esm-protein-atlas/v1/folds/folds_atlas.lance",   #       6.6M representative set
    "1b":    "s3://esm-protein-atlas/v1/folds/folds_1B.lance",      # 1,095M full set
}
DATASET_ORDER = ["atlas", "1b"]                    # default query order (cheaper set first)
OUTDIR = "pdbs"
BATCH = 4000                                       # hashes per indexed query

# --- atom37 -> PDB ----------------------------------------------------------
ATOM37 = ['N','CA','C','CB','O','CG','CG1','CG2','OG','OG1','SG','CD','CD1','CD2',
'ND1','ND2','OD1','OD2','SD','CE','CE1','CE2','CE3','NE','NE1','NE2','OE1','OE2',
'CH2','NH1','NH2','OH','CZ','CZ2','CZ3','NZ','OXT']
AA3 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
'H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
'T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'UNK'}

def _arr(d):
    return np.frombuffer(bytes(d[b'data']), dtype=np.dtype(d[b'type'])).reshape(tuple(d[b'shape']))

def blob_to_pdb(blob):
    o = msgpack.unpackb(brotli.decompress(blob), raw=False, strict_map_key=False)
    seq   = o['sequence']
    pos_c = _arr(o['atom37_positions']).astype(np.float32)     # [n_present, 3]
    mask  = _arr(o['atom37_mask']).astype(bool)                # [n_res, 37]
    resid = _arr(o['residue_index'])                           # 1-based already
    conf  = _arr(o['confidence']).astype(np.float32)           # per-residue pLDDT
    nres  = mask.shape[0]
    full = np.zeros((nres, 37, 3), np.float32); full[mask] = pos_c
    bf = conf*100 if conf.size and conf.max() <= 1.0 else conf
    out, serial = [], 1
    for i in range(nres):
        rn = AA3.get(seq[i], 'UNK')
        for j in range(37):
            if not mask[i, j]:
                continue
            x, y, z = full[i, j]; name = ATOM37[j]
            an = (" " + name) if len(name) < 4 else name
            out.append(f"ATOM  {serial:>5} {an:<4} {rn:>3} A{int(resid[i]):>4}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{bf[i]:6.2f}          {name[0]:>2}")
            serial += 1
    out += ["TER", "END"]
    return "\n".join(out) + "\n"

_CIF_HEADER = ("loop_\n"
    "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
    "_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
    "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
    "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
    "_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n")

def blob_to_cif(blob, name="structure"):
    # Same atom37 -> coordinates decode as blob_to_pdb, emitted as an mmCIF atom_site
    # loop. mmCIF has no 99,999-atom / numbering limits, so it's the safer container.
    o = msgpack.unpackb(brotli.decompress(blob), raw=False, strict_map_key=False)
    seq   = o['sequence']
    pos_c = _arr(o['atom37_positions']).astype(np.float32)
    mask  = _arr(o['atom37_mask']).astype(bool)
    resid = _arr(o['residue_index'])
    conf  = _arr(o['confidence']).astype(np.float32)
    nres  = mask.shape[0]
    full = np.zeros((nres, 37, 3), np.float32); full[mask] = pos_c
    bf = conf*100 if conf.size and conf.max() <= 1.0 else conf
    out = [f"data_{name}", "#", _CIF_HEADER.rstrip("\n")]
    serial = 1
    for i in range(nres):
        rn = AA3.get(seq[i], 'UNK'); ri = int(resid[i])
        for j in range(37):
            if not mask[i, j]:
                continue
            x, y, z = full[i, j]; an = ATOM37[j]
            out.append(f"ATOM {serial} {an[0]} {an} {rn} A 1 {ri} "
                       f"{x:.3f} {y:.3f} {z:.3f} 1.00 {bf[i]:.2f} {ri} A")
            serial += 1
    out.append("#")
    return "\n".join(out) + "\n"

# --- FASTA -> set of unique protein_hashes (+ optional full orf_id->hash index) ---
def parse_fasta(path, orf_map_path=None):
    # Streams the FASTA once. Returns the SET of unique protein_hashes to look up.
    # If orf_map_path is given, also streams every 'orf_id\tprotein_hash' row there
    # (the complete index over all ORF records) -- memory-flat, so it scales to 100+ GB.
    hashes = set()
    hid, seq = None, []
    mf = None
    if orf_map_path:
        mf = (zstd.open(orf_map_path, "wt") if orf_map_path.endswith(".zst")
              else open(orf_map_path, "w"))
        mf.write("orf_id\tprotein_hash\n")

    def flush():
        if hid and seq:
            # complete ORFs carry a trailing stop codon '*'; ESMFold never folds the
            # stop, so the Atlas protein_hash is md5 over the sequence WITHOUT it.
            s = "".join(seq).upper().rstrip("*")
            h = hashlib.md5(s.encode()).hexdigest()
            hashes.add(h)
            if mf:
                mf.write(f"{hid}\t{h}\n")

    try:
        with tqdm(total=_input_size(path), unit="B", unit_scale=True,
                  unit_divisor=1024, desc="parse") as pbar:
            for line in _iter_lines(path):
                pbar.update(len(line) + 1)            # ~bytes read (progress only, approximate)
                line = line.strip()
                if line.startswith(">"):
                    flush()
                    hid = line[1:].split()[0].split("|")[0]   # orf_id = first header field
                    seq = []
                elif line:
                    seq.append(line)
            flush()
    finally:
        if mf:
            mf.close()
    return hashes

def _parse_s3(uri):
    bucket, _, key = uri[5:].partition("/")           # uri == s3://bucket/prefix...
    return bucket, key.strip("/")

def _s3_read_client():
    # Anonymous client for PUBLIC input buckets (e.g. petadex) -- no credentials needed,
    # matching the anonymous Atlas access. (Signed boto3 is only used for WRITING output.)
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    return boto3.client("s3", config=Config(signature_version=UNSIGNED))

def _input_size(path):
    # Total bytes of the input, for the parse progress bar; None if unknown (bar still
    # runs, just without a percentage / ETA).
    try:
        if path.startswith("s3://"):
            bucket, key = _parse_s3(path)
            return _s3_read_client().head_object(Bucket=bucket, Key=key)["ContentLength"]
        return os.path.getsize(path)
    except Exception:
        return None

def _iter_lines(path):
    # Yields text lines from a local path OR an s3:// URI (so the 100+ GB input never
    # has to land on the VM's disk); handles .zst transparently in both cases.
    if path.startswith("s3://"):
        bucket, key = _parse_s3(path)
        body = _s3_read_client().get_object(Bucket=bucket, Key=key)["Body"]
        if path.endswith(".zst"):
            reader = zstd.ZstdDecompressor().stream_reader(body)
            yield from io.TextIOWrapper(reader, encoding="utf-8")
        else:
            for raw in body.iter_lines(chunk_size=1 << 20):
                yield raw.decode("utf-8", "replace")
    else:
        opener = (zstd.open if path.endswith(".zst") else open)
        with opener(path, "rt") as f:
            yield from f

def batched(it, n):
    it = list(it)
    for i in range(0, len(it), n):
        yield it[i:i+n]

def hashes_from_map(path):
    # Load the SET of unique protein_hashes from a previously written orf_map
    # (orf_id<TAB>protein_hash, optionally .zst). Lets the download skip re-parsing the
    # full FASTA: reads ~2-3 GB instead of 116 GB, and does no md5 work.
    hashes = set()
    it = _iter_lines(path)
    next(it, None)                                    # skip the header row
    for line in tqdm(it, desc="load orf_map", unit="row", unit_scale=True):
        line = line.rstrip("\n")
        if line:
            hashes.add(line.rsplit("\t", 1)[-1])      # protein_hash is the last column
    return hashes

def _scan_retry(ds, columns, filt, attempts=4):
    # Lance reads the dataset in byte ranges over S3; cross-region range fetches
    # occasionally time out. Retry the batch with backoff, then give up (return None)
    # so one transient blip can't kill the run -- skipped hashes are simply picked up
    # on the next resume (download) or slightly undercounted (count).
    for i in range(attempts):
        try:
            return ds.scanner(columns=columns, filter=filt).to_table().to_pylist()
        except Exception as e:
            if i == attempts - 1:
                tqdm.write(f"  WARN batch failed after {attempts} tries, skipping: {e}")
                return None
            time.sleep(2 ** i)                        # 1s, 2s, 4s backoff

def _ibatch(seq, n):                                  # lazy fixed-size slices of an indexable
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def _bounded(ex, fn, items, max_inflight):
    # Run fn over items with at most max_inflight futures alive at once, yielding each
    # result as it completes. Never materializes all futures/chunks, so memory stays flat
    # no matter how many items there are (the 124M-hash run would OOM otherwise).
    it = iter(items)
    inflight = set()
    for x in islice(it, max_inflight):
        inflight.add(ex.submit(fn, x))
    while inflight:
        done, inflight = wait(inflight, return_when=FIRST_COMPLETED)
        for f in done:
            yield f.result()
            nxt = next(it, None)
            if nxt is not None:
                inflight.add(ex.submit(fn, nxt))

def count_hits(hashes, dsets, workers=16):
    """Dry-run: query protein_hash only (no structure_blob, no files) and report how
    many sequences are present in the Atlas, broken down by dataset. Batches are queried
    concurrently; the `remaining`/`hits` state is updated only on the main thread."""
    n_total = len(hashes)
    remaining = hashes                                # operate in place -- no duplicate set
    per_ds = {}

    def query_chunk(ds, chunk):                       # hash-only: no blob over the wire
        q = ",".join(f"'{h}'" for h in chunk)
        tbl = _scan_retry(ds, ["protein_hash"], f"protein_hash IN ({q})")
        return [r["protein_hash"] for r in tbl] if tbl else []

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for u, ds in dsets:                           # representative set first, then 1B
            if not remaining:
                break
            name = u.rstrip("/").split("/")[-1]
            hits = 0
            pending = list(remaining)                 # one snapshot (~1 GB at 124M); safe vs mutation
            n_batches = (len(pending) + BATCH - 1) // BATCH
            with tqdm(total=n_batches, desc=f"count {name}", unit="batch") as pbar:
                for found in _bounded(ex, lambda c, ds=ds: query_chunk(ds, c),
                                      _ibatch(pending, BATCH), workers * 2):
                    for h in found:
                        if h in remaining:
                            remaining.discard(h); hits += 1
                    pbar.update(1)
                    pbar.set_postfix(hits=hits, unmatched=len(remaining))
            per_ds[name] = hits
    total = n_total - len(remaining)
    print(f"\nHIT RATE: {total}/{n_total} = {100*total/max(n_total,1):.2f}% present in Atlas")
    for name, n in per_ds.items():
        print(f"  {n:>9} from {name}")
    print(f"  {len(remaining):>9} not found in any dataset")

def download(hashes, dsets, outdir=OUTDIR, s3_dest=None, fmt="cif", workers=16):
    # Writes one structure per unique sequence, keyed by protein_hash, to
    # <outdir>/<hash>.<fmt> or s3://bucket/prefix/<hash>.<fmt>. Atlas queries + decode
    # + upload run concurrently across `workers` threads (the work is I/O-bound).
    # found_hashes.tsv is the resume ledger + ORF-lookup join key; it is written only
    # from the main thread as batches complete, so no file-write locking is needed.
    serialize = (lambda blob, h: blob_to_cif(blob, h)) if fmt == "cif" \
                else (lambda blob, h: blob_to_pdb(blob))
    ctype = {"cif": "chemical/x-mmcif", "pdb": "chemical/x-pdb"}[fmt]
    os.makedirs(outdir, exist_ok=True)
    found_path = os.path.join(outdir, "found_hashes.tsv")
    done = set()
    if os.path.exists(found_path):
        with open(found_path) as fh:
            next(fh, None)                            # skip header
            done = {line.split("\t", 1)[0] for line in fh if line.strip()}
    n0 = len(hashes)
    hashes.difference_update(done)                    # in place -- avoid a 124M-set copy
    remaining = hashes
    print(f"{n0} unique sequences | {len(done)} already done | "
          f"{len(remaining)} to fetch | {workers} workers")

    s3 = bucket = prefix = None
    if s3_dest:
        import boto3
        # one shared low-level client (boto3 clients are thread-safe for put_object)
        s3 = boto3.client("s3")
        bucket, prefix = _parse_s3(s3_dest)

    def put(h, text):
        if s3:
            key = f"{prefix}/{h}.{fmt}" if prefix else f"{h}.{fmt}"
            s3.put_object(Bucket=bucket, Key=key, Body=text.encode(), ContentType=ctype)
        else:
            with open(os.path.join(outdir, f"{h}.{fmt}"), "w") as fh:
                fh.write(text)

    def process_chunk(ds, chunk):
        # one worker: query a batch, then decode + upload every hit; return found hashes
        q = ",".join(f"'{h}'" for h in chunk)
        tbl = _scan_retry(ds, ["protein_hash", "structure_blob"], f"protein_hash IN ({q})")
        if not tbl:
            return []
        got = []
        for r in tbl:
            h = r["protein_hash"]
            try:
                put(h, serialize(r["structure_blob"], h))
                got.append(h)
            except Exception as e:                    # a bad structure must not kill the batch
                tqdm.write(f"  WARN {h}: {e}")
        return got

    found = 0
    write_header = not os.path.exists(found_path)
    with open(found_path, "a") as ff, ThreadPoolExecutor(max_workers=workers) as ex:
        if write_header:
            ff.write("protein_hash\tsource_dataset\n")
        for u, ds in dsets:                           # try representative set first, then 1B
            if not remaining:
                break
            src = u.rstrip("/").split("/")[-1]
            pending = list(remaining)                 # one snapshot; bounded submission below
            n_batches = (len(pending) + BATCH - 1) // BATCH
            with tqdm(total=n_batches, desc=f"download {src}", unit="batch") as pbar:
                for got in _bounded(ex, lambda c, ds=ds: process_chunk(ds, c),
                                    _ibatch(pending, BATCH), workers * 2):
                    for h in got:
                        if h in remaining:
                            ff.write(f"{h}\t{src}\n")
                            remaining.discard(h); found += 1
                    ff.flush()                        # crash-safe: persist ledger per batch
                    pbar.update(1)
                    pbar.set_postfix(found=found, unmatched=len(remaining))
    dest = s3_dest if s3_dest else outdir
    print(f"\nDone: {found} structures written to {dest} ; "
          f"{len(remaining)} sequences not in the atlas.")
    print(f"Found-set ledger: {found_path}")

def main(fasta, count_only=False, sample=None, seed=0,
         orf_map=None, s3_dest=None, outdir=OUTDIR, fmt="cif", workers=16,
         hashes_from=None, datasets=None):
    if hashes_from:                                   # reuse a banked orf_map -> no re-parse
        hashes = hashes_from_map(hashes_from)
        print(f"{len(hashes)} unique hashes loaded from {hashes_from}")
    else:
        hashes = parse_fasta(fasta, orf_map_path=orf_map)
        print(f"{len(hashes)} unique sequences parsed from {fasta}")
        if orf_map:
            print(f"full orf_id -> protein_hash index written to {orf_map}")
    if sample is not None and sample < len(hashes):
        import random
        hashes = set(random.Random(seed).sample(list(hashes), sample))
        print(f"sampled down to {len(hashes)} sequences (seed={seed})")
    names = [n.strip() for n in datasets.split(",")] if datasets else DATASET_ORDER
    bad = [n for n in names if n not in DATASETS]
    if bad:
        raise SystemExit(f"unknown --datasets {bad}; choose from {list(DATASETS)}")
    print(f"querying datasets: {', '.join(names)}")
    dsets = [(DATASETS[n], lance.dataset(DATASETS[n], storage_options=STORAGE)) for n in names]
    if count_only:
        count_hits(hashes, dsets, workers=workers)
    else:
        download(hashes, dsets, outdir=outdir, s3_dest=s3_dest, fmt=fmt, workers=workers)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Fetch ESM Atlas structures (keyed by protein_hash) for a FASTA.")
    ap.add_argument("fasta", nargs="?", default="proteins.fasta")
    ap.add_argument("--count", action="store_true",
                    help="dry-run: query protein_hash only and report Atlas hit rate "
                         "(no blob download, no PDB writes)")
    ap.add_argument("--sample", type=int, default=None, metavar="N",
                    help="randomly subsample N sequences before querying (seeded)")
    ap.add_argument("--seed", type=int, default=0, help="random seed for --sample")
    ap.add_argument("--orf-map", default=None, metavar="PATH",
                    help="write the full orf_id->protein_hash index here (.tsv or .tsv.zst); "
                         "join against found_hashes.tsv to see which ORFs have a structure")
    ap.add_argument("--s3", dest="s3_dest", default=None, metavar="s3://bucket/prefix",
                    help="upload <hash>.<format> objects here instead of writing locally")
    ap.add_argument("--outdir", default=OUTDIR,
                    help=f"local output dir for structures + found_hashes.tsv (default {OUTDIR})")
    ap.add_argument("--format", dest="fmt", choices=["cif", "pdb"], default="cif",
                    help="structure file format (default cif)")
    ap.add_argument("--workers", type=int, default=16,
                    help="concurrent Atlas-query/upload threads (default 16)")
    ap.add_argument("--hashes-from", dest="hashes_from", default=None, metavar="ORF_MAP",
                    help="load unique hashes from a saved orf_map (.tsv/.tsv.zst) instead of "
                         "re-parsing the FASTA -- much faster for the download run")
    ap.add_argument("--datasets", default=None, metavar="LIST",
                    help="comma list of Atlas datasets to query: atlas,1b "
                         "(default both, atlas first; e.g. --datasets 1b for the 1B set only)")
    args = ap.parse_args()
    main(args.fasta, count_only=args.count, sample=args.sample, seed=args.seed,
         orf_map=args.orf_map, s3_dest=args.s3_dest, outdir=args.outdir,
         fmt=args.fmt, workers=args.workers,
         hashes_from=args.hashes_from, datasets=args.datasets)
