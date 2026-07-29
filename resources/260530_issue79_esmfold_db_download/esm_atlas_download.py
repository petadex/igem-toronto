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

Output model (content-addressed; one object per UNIQUE sequence keyed by protein_hash, so
the input's heavy ORF redundancy collapses). Byte-for-byte the SAME LAYOUT that
esmfold2_local_predictor.py writes (issue #6), so Atlas and freshly-folded structures are
interchangeable downstream:
  - <outdir>/structures/<hash>.<fmt>   all-atom coords (atom37 -> mmCIF/PDB; pLDDT in B-factor)
  - <outdir>/metrics/<hash>.json       run block + scalars + per_residue_plddt + pae, one file
                                       per structure (scalars FIRST, so a Range-GET of the
                                       first few KB reads them without the arrays)
  - <outdir>/tracker.csv               the scalar LEDGER, columns identical to the predictor's
                                       TRACKER_COLS (resume table + ORF-join key; mirrored to
                                       <s3-prefix>/tracker.csv every ~5 min + at exit and
                                       seeded back from S3 on a fresh box)
  - --orf-map writes the full orf_id -> protein_hash index (every ORF record, streamed).
  - "which ORFs have a structure?"  ==  orf_map rows whose protein_hash is in tracker.csv.

The id namespace is the protein_hash: the Atlas is content-addressed and has no notion of a
PETadex ORFid, so `id` == `protein_hash` on every row. The predictor keys `id` by ORFid but
carries the same `protein_hash` column, computed the same way (md5 of the uppercased sequence
with the trailing '*' stripped) -- so joining the two sources on `protein_hash` works uniformly.

Fields the Atlas cannot supply are null/empty rather than invented: `esmc` (no backbone id),
`timing.*` (these were folded by Biohub, not by us), `iptm` (single chains). `run_label` carries
which Lance dataset the row came from. There is no `has_pae` column any more (the predictor has
none either) -- test `confidence.pae is not None` in the JSON.

pTM (global) and PAE ([L,L]) are read from the Atlas 'ptm'/'pae' columns. Earlier versions of
this script kept only pLDDT (in the CIF B-factor) and discarded pTM/PAE. To backfill them for a
pre-existing download, delete tracker.csv so the affected hashes are re-fetched.
"""
import os, io, csv, json, time, hashlib
import numpy as np
import zstandard as zstd
# brotli, msgpack, lance and boto3 are imported lazily (inside the functions that use
# them) so this module imports -- and --help / the pure helpers run -- without the heavy
# S3/Lance stack installed. They are only needed for an actual Atlas query.
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
S3_DEST = "s3://petadex-protein-structures/esmatlas/"   # canonical sink (issue #79)

# --- atom37 -> PDB ----------------------------------------------------------
ATOM37 = ['N','CA','C','CB','O','CG','CG1','CG2','OG','OG1','SG','CD','CD1','CD2',
'ND1','ND2','OD1','OD2','SD','CE','CE1','CE2','CE3','NE','NE1','NE2','OE1','OE2',
'CH2','NH1','NH2','OH','CZ','CZ2','CZ3','NZ','OXT']
AA3 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
'H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
'T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'UNK'}

def _arr(d):
    return np.frombuffer(bytes(d[b'data']), dtype=np.dtype(d[b'type'])).reshape(tuple(d[b'shape']))

_CIF_HEADER = ("loop_\n"
    "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
    "_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
    "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
    "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
    "_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n")

def _decode_blob(blob):
    # brotli + msgpack-numpy -> the atom37 structure, expanded to full [L,37,3] coords.
    # pLDDT is normalised to the 0-100 scale HERE so every downstream writer is consistent.
    import brotli, msgpack
    o = msgpack.unpackb(brotli.decompress(blob), raw=False, strict_map_key=False)
    seq   = o['sequence']
    pos_c = _arr(o['atom37_positions']).astype(np.float32)     # [n_present, 3]
    mask  = _arr(o['atom37_mask']).astype(bool)                # [n_res, 37]
    resid = _arr(o['residue_index'])                           # 1-based already
    conf  = _arr(o['confidence']).astype(np.float32)           # per-residue pLDDT
    nres  = mask.shape[0]
    full  = np.zeros((nres, 37, 3), np.float32); full[mask] = pos_c
    plddt = conf * 100 if conf.size and conf.max() <= 1.0 else conf   # -> 0-100
    return {"seq": seq, "pos": full, "mask": mask, "resid": resid, "plddt": plddt, "nres": nres}

def _cif_from_decoded(d, name="structure"):
    # atom37 -> mmCIF atom_site loop (pLDDT in B_iso_or_equiv). mmCIF has no 99,999-atom /
    # numbering limits, so it's the safer container.
    seq, pos, mask, resid, bf = d["seq"], d["pos"], d["mask"], d["resid"], d["plddt"]
    out = [f"data_{name}", "#", _CIF_HEADER.rstrip("\n")]
    serial = 1
    for i in range(mask.shape[0]):
        rn = AA3.get(seq[i], 'UNK'); ri = int(resid[i])
        for j in range(37):
            if not mask[i, j]:
                continue
            x, y, z = pos[i, j]; an = ATOM37[j]
            out.append(f"ATOM {serial} {an[0]} {an} {rn} A 1 {ri} "
                       f"{x:.3f} {y:.3f} {z:.3f} 1.00 {bf[i]:.2f} {ri} A")
            serial += 1
    out.append("#")
    return "\n".join(out) + "\n"

def _pdb_from_decoded(d):
    seq, pos, mask, resid, bf = d["seq"], d["pos"], d["mask"], d["resid"], d["plddt"]
    out, serial = [], 1
    for i in range(mask.shape[0]):
        rn = AA3.get(seq[i], 'UNK')
        for j in range(37):
            if not mask[i, j]:
                continue
            x, y, z = pos[i, j]; nm = ATOM37[j]
            an = (" " + nm) if len(nm) < 4 else nm
            out.append(f"ATOM  {serial:>5} {an:<4} {rn:>3} A{int(resid[i]):>4}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{bf[i]:6.2f}          {nm[0]:>2}")
            serial += 1
    out += ["TER", "END"]
    return "\n".join(out) + "\n"

# Back-compat one-shot wrappers (decode + format), used anywhere that still holds a raw blob.
def blob_to_cif(blob, name="structure"):
    return _cif_from_decoded(_decode_blob(blob), name)

def blob_to_pdb(blob):
    return _pdb_from_decoded(_decode_blob(blob))

def _as_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def _coerce_pae(pae_val, nres):
    # Normalise the Atlas 'pae' column to an [L,L] float32 (or None). VERIFIED on a live
    # row: the column is a *serialized numpy blob* -- np.savez bytes (a zip whose first
    # member is the PAE array), NOT a plain list. So bytes are np.load()-ed; lists/ndarrays
    # are still handled for robustness. Never raises -> a bad PAE can't drop the structure.
    if pae_val is None:
        return None
    try:
        if isinstance(pae_val, (bytes, bytearray, memoryview)):
            obj = np.load(io.BytesIO(bytes(pae_val)), allow_pickle=False)
            if hasattr(obj, "files"):                 # NpzFile -> take its first array
                obj = obj[obj.files[0]]
            arr = np.asarray(obj, dtype=np.float32)
        else:
            arr = np.asarray(pae_val, dtype=np.float32)
    except Exception:
        return None
    if arr.ndim == 1:
        if nres and arr.size == nres * nres:
            arr = arr.reshape(nres, nres)
        else:
            return None
    return arr if arr.ndim == 2 else None

# tracker.csv columns -- IDENTICAL to TRACKER_COLS in esmfold2_local_predictor.py (that file is
# the source of truth; duplicated here to keep this script single-file and import-free). The rich
# per-structure record (per-residue pLDDT + PAE) lives in metrics/<hash>.json, NOT here.
TRACKER_COLS = ["id", "run_label", "esmfold2", "esmc", "seq_len", "protein_hash",
                "truncated", "wall_s_median", "residues_per_s", "mean_plddt", "ptm", "iptm",
                "status", "error"]

def build_metrics(h, src, d, pae, mean_plddt, ptm):
    # The per-structure metrics/<hash>.json, shaped exactly like the predictor's build_metrics():
    # scalars first, then the big per-residue/PAE arrays, so a Range-GET of the head reads the
    # scalars alone. Arrays are rounded to 2 dp to match the cif B-factor precision and keep the
    # JSON small. Atlas-unavailable fields are null rather than fabricated.
    return {
        "id": h,
        # Same key names as the predictor's run_config, so run.label / run.esmfold2 /
        # run.num_sampling_steps read uniformly across both sources. The Atlas does not publish
        # the parameters these structures were folded with, so those are null rather than guessed.
        "run": {"label": src, "esmfold2": "esm-atlas-v1", "esmc": None,
                "num_loops": None, "num_sampling_steps": None, "num_diffusion_samples": None,
                "seed": None, "dtype": None,
                "source": "esm-atlas", "dataset": src,
                "note": "precomputed by Biohub; fold parameters not published"},
        "seq_len": d["nres"],
        "protein_hash": h,
        "truncated": False,
        "timing": {"wall_s_median": None, "residues_per_s": None},
        "confidence": {
            "mean_plddt": mean_plddt,
            "ptm": ptm,
            "iptm": None,
            "per_residue_plddt": [round(float(x), 2) for x in np.asarray(d["plddt"]).tolist()],
            "pae": ([[round(float(x), 2) for x in r] for r in pae.tolist()]
                    if pae is not None else None),
        },
        "status": "ok",
        "error": "",
    }

def build_tracker_row(h, src, d, mean_plddt, ptm):
    # One scalar ledger row, same columns/order as the predictor's Row.as_dict().
    return {"id": h, "run_label": src, "esmfold2": "esm-atlas-v1", "esmc": "",
            "seq_len": d["nres"], "protein_hash": h, "truncated": 0,
            "wall_s_median": "", "residues_per_s": "",
            "mean_plddt": mean_plddt, "ptm": (ptm if ptm is not None else ""), "iptm": "",
            "status": "ok", "error": ""}

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

def _s3_put_file(s3, bucket, key, path, ctype="text/csv"):
    # Mirror a local file to s3://bucket/key. Used to keep tracker.csv (the resume ledger +
    # aggregate confidence table) durable off an ephemeral VM. Re-uploads the whole file, so
    # it is called on a time interval (not per batch) -- cheap relative to the structure data.
    with open(path, "rb") as fh:
        s3.put_object(Bucket=bucket, Key=key, Body=fh.read(), ContentType=ctype)

def _s3_get_file(s3, bucket, key, path):
    # Pull s3://bucket/key down to a local path; return True if it existed (seeds resume
    # from a previous run on a fresh instance). Signed client: reads our own private bucket.
    try:
        obj = s3.get_object(Bucket=bucket, Key=key)
    except Exception:
        return False
    with open(path, "wb") as fh:
        fh.write(obj["Body"].read())
    return True

def download(hashes, dsets, outdir=OUTDIR, s3_dest=None, fmt="cif", workers=16):
    # Writes, per unique sequence (keyed by protein_hash), the SAME {structures/, metrics/,
    # tracker.csv} set as esmfold2_local_predictor.py:
    #   <outdir>/structures/<hash>.<fmt>   all-atom coords (pLDDT in the B-factor)
    #   <outdir>/metrics/<hash>.json       run block + scalars + per_residue_plddt + pae
    #   <outdir>/tracker.csv               TRACKER_COLS scalar ledger (also the resume table
    #                                      + ORF join key, via its protein_hash column)
    # Atlas query + decode + upload run concurrently across `workers` threads (I/O-bound).
    # tracker.csv is appended only from the main thread as batches complete -> no locking.
    ctype = {"cif": "chemical/x-mmcif", "pdb": "chemical/x-pdb"}[fmt]
    tracker_path = os.path.join(outdir, "tracker.csv")
    os.makedirs(outdir, exist_ok=True)

    s3 = bucket = prefix = tracker_key = None
    if s3_dest:
        import boto3
        s3 = boto3.client("s3")                       # thread-safe put_object
        bucket, prefix = _parse_s3(s3_dest)
        tracker_key = "/".join(p for p in (prefix, "tracker.csv") if p)
        if not os.path.exists(tracker_path) and _s3_get_file(s3, bucket, tracker_key, tracker_path):
            print(f"resume ledger seeded from s3://{bucket}/{tracker_key}")
    else:
        os.makedirs(os.path.join(outdir, "structures"), exist_ok=True)
        os.makedirs(os.path.join(outdir, "metrics"), exist_ok=True)

    done = set()
    if os.path.exists(tracker_path):                  # resume: id (== protein_hash) is column 0
        with open(tracker_path) as fh:
            next(fh, None)                            # skip header
            done = {line.split(",", 1)[0] for line in fh if line.strip()}
    n0 = len(hashes)
    hashes.difference_update(done)                    # in place -- avoid a 124M-set copy
    remaining = hashes
    print(f"{n0} unique sequences | {len(done)} already done | "
          f"{len(remaining)} to fetch | {workers} workers")

    ctypes = {"structures": ctype, "metrics": "application/json"}

    def write_out(subdir, name, data):
        if s3:
            key = "/".join(p for p in (prefix, subdir, name) if p)
            body = data.encode() if isinstance(data, str) else data
            s3.put_object(Bucket=bucket, Key=key, Body=body,
                          ContentType=ctypes.get(subdir, "application/octet-stream"))
        else:
            with open(os.path.join(outdir, subdir, name),
                      "w" if isinstance(data, str) else "wb") as fh:
                fh.write(data)

    def process_chunk(ds, chunk, src):
        # one worker: query a batch, decode, write structure + metrics json, return ledger rows.
        # 'ptm' and 'pae' are pulled alongside the structure_blob (previously dropped).
        q = ",".join(f"'{h}'" for h in chunk)
        tbl = _scan_retry(ds, ["protein_hash", "structure_blob", "ptm", "pae"],
                          f"protein_hash IN ({q})")
        if not tbl:
            return []
        rows = []
        for r in tbl:
            h = r["protein_hash"]
            try:
                d = _decode_blob(r["structure_blob"])
                text = _cif_from_decoded(d, h) if fmt == "cif" else _pdb_from_decoded(d)
                pae = _coerce_pae(r.get("pae"), d["nres"])
                ptm = _as_float(r.get("ptm"))
                ptm = round(ptm, 4) if ptm is not None else None
                mean_plddt = round(float(np.mean(d["plddt"])), 3)
                write_out("structures", f"{h}.{fmt}", text)
                write_out("metrics", f"{h}.json",
                          json.dumps(build_metrics(h, src, d, pae, mean_plddt, ptm),
                                     separators=(",", ":")))
                rows.append(build_tracker_row(h, src, d, mean_plddt, ptm))
            except Exception as e:                    # a bad structure must not kill the batch
                tqdm.write(f"  WARN {h}: {e}")
        return rows

    found = 0
    last_sync = time.time()
    write_header = not os.path.exists(tracker_path) or os.path.getsize(tracker_path) == 0
    with open(tracker_path, "a", newline="") as mf, ThreadPoolExecutor(max_workers=workers) as ex:
        writer = csv.DictWriter(mf, fieldnames=TRACKER_COLS)
        if write_header:
            writer.writeheader()
        for u, ds in dsets:                           # try representative set first, then 1B
            if not remaining:
                break
            src = u.rstrip("/").split("/")[-1]
            pending = list(remaining)                 # one snapshot; bounded submission below
            n_batches = (len(pending) + BATCH - 1) // BATCH
            with tqdm(total=n_batches, desc=f"download {src}", unit="batch") as pbar:
                for rows in _bounded(ex, lambda c, ds=ds, src=src: process_chunk(ds, c, src),
                                     _ibatch(pending, BATCH), workers * 2):
                    for row in rows:
                        h = row["id"]
                        if h in remaining:
                            writer.writerow(row)
                            remaining.discard(h); found += 1
                    mf.flush()                        # crash-safe: persist ledger per batch
                    if tracker_key and time.time() - last_sync > 300:   # mirror to S3 ~every 5 min
                        _s3_put_file(s3, bucket, tracker_key, tracker_path)
                        last_sync = time.time()
                    pbar.update(1)
                    pbar.set_postfix(found=found, unmatched=len(remaining))
    if tracker_key:                                   # final mirror so the ledger is durable off-box
        _s3_put_file(s3, bucket, tracker_key, tracker_path)
    dest = s3_dest if s3_dest else outdir
    print(f"\nDone: {found} structures (+ metrics JSONs) written to {dest} ; "
          f"{len(remaining)} sequences not in the atlas.")
    where = f"s3://{bucket}/{tracker_key}" if tracker_key else tracker_path
    print(f"Tracker + resume ledger: {tracker_path}" + (f" (mirrored to {where})" if tracker_key else ""))

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
    import lance
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
                         "join it against tracker.csv's protein_hash column to see which ORFs "
                         "have a structure")
    ap.add_argument("--s3", dest="s3_dest", default=S3_DEST, metavar="s3://bucket/prefix/",
                    help=f"upload structures/ + metrics/ here (default {S3_DEST}); "
                         f"pass --s3 '' to write locally instead (tracker.csv is always written "
                         f"locally and, when uploading, mirrored to <prefix>/tracker.csv)")
    ap.add_argument("--outdir", default=OUTDIR,
                    help=f"local output dir for structures/ + metrics/ + tracker.csv "
                         f"(default {OUTDIR})")
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
