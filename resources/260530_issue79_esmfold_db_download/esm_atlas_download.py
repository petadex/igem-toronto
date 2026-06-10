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

Install:  pip install pylance pyarrow msgpack numpy brotli zstandard boto3

Output model (content-addressed):
  - PDBs are keyed by protein_hash:  <outdir>/<hash>.pdb  (or s3://bucket/prefix/<hash>.pdb)
    One object per UNIQUE sequence, so the input's heavy ORF redundancy collapses.
  - --orf-map writes the full orf_id -> protein_hash index (every ORF record, streamed).
  - found_hashes.tsv lists which hashes actually have a structure (drives resume + the join).
  - "which ORFs have a structure?"  ==  orf_map rows whose protein_hash is in found_hashes.
"""
import os, hashlib, time
import numpy as np
import lance, brotli, msgpack
import zstandard as zstd

STORAGE = {"aws_skip_signature": "true"}          # anonymous public-bucket access
DATASETS = [
    "s3://esm-protein-atlas/v1/folds/folds_atlas.lance",
    "s3://esm-protein-atlas/v1/folds/folds_1B.lance",
]
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

    opener = (zstd.open if path.endswith(".zst") else open)
    try:
        with opener(path, "rt") as f:
            for line in f:
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

def batched(it, n):
    it = list(it)
    for i in range(0, len(it), n):
        yield it[i:i+n]

def count_hits(hashes, dsets):
    """Dry-run: query protein_hash only (no structure_blob, no PDB writes) and
    report how many sequences are present in the Atlas, broken down by dataset."""
    remaining = set(hashes)
    n_total = len(hashes)
    per_ds = {}
    t0 = time.time()
    for u, ds in dsets:                               # representative set first, then 1B
        if not remaining:
            break
        name = u.rstrip("/").split("/")[-1]
        hits = 0
        for chunk in batched(remaining, BATCH):
            q = ",".join(f"'{h}'" for h in chunk)
            tbl = ds.scanner(columns=["protein_hash"],   # hash-only: no blob over the wire
                             filter=f"protein_hash IN ({q})").to_table().to_pylist()
            for r in tbl:
                h = r["protein_hash"]
                if h in remaining:
                    remaining.discard(h); hits += 1
            checked = n_total - len(remaining)
            print(f"  [{name}] hits {hits} | "
                  f"{checked/max(time.time()-t0,1e-9):.0f} resolved/s | "
                  f"{len(remaining)} unmatched", end="\r")
        per_ds[name] = hits
        print()
    total = n_total - len(remaining)
    print(f"\nHIT RATE: {total}/{n_total} = {100*total/max(n_total,1):.2f}% present in Atlas")
    for name, n in per_ds.items():
        print(f"  {n:>9} from {name}")
    print(f"  {len(remaining):>9} not found in any dataset")

def download(hashes, dsets, outdir=OUTDIR, s3_dest=None, fmt="cif"):
    # Writes one structure per unique sequence, keyed by protein_hash, to
    # <outdir>/<hash>.<fmt> or s3://bucket/prefix/<hash>.<fmt>. Appends each resolved
    # hash to found_hashes.tsv, the resume ledger + join key for ORF lookup.
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
    remaining = set(hashes) - done
    print(f"{len(hashes)} unique sequences | {len(done)} already done | {len(remaining)} to fetch")

    s3 = bucket = prefix = None
    if s3_dest:
        import boto3
        s3 = boto3.client("s3")
        bucket, prefix = _parse_s3(s3_dest)

    def put(h, text):
        if s3:
            key = f"{prefix}/{h}.{fmt}" if prefix else f"{h}.{fmt}"
            s3.put_object(Bucket=bucket, Key=key, Body=text.encode(), ContentType=ctype)
        else:
            with open(os.path.join(outdir, f"{h}.{fmt}"), "w") as fh:
                fh.write(text)

    found = 0; t0 = time.time()
    write_header = not os.path.exists(found_path)
    with open(found_path, "a") as ff:
        if write_header:
            ff.write("protein_hash\tsource_dataset\n")
        for u, ds in dsets:                           # try representative set first, then 1B
            if not remaining:
                break
            src = u.rstrip("/").split("/")[-1]
            for chunk in batched(remaining, BATCH):
                q = ",".join(f"'{h}'" for h in chunk)
                tbl = ds.scanner(columns=["protein_hash", "structure_blob"],
                                 filter=f"protein_hash IN ({q})").to_table().to_pylist()
                for r in tbl:
                    h = r["protein_hash"]
                    if h not in remaining:
                        continue
                    put(h, serialize(r["structure_blob"], h))
                    ff.write(f"{h}\t{src}\n")
                    remaining.discard(h); found += 1
                ff.flush()                            # crash-safe: persist ledger per batch
                print(f"  found {found} | {found/max(time.time()-t0,1e-9):.0f}/s | "
                      f"{len(remaining)} unmatched", end="\r")
    dest = s3_dest if s3_dest else outdir
    print(f"\nDone: {found} structures written to {dest} ; "
          f"{len(remaining)} sequences not in the atlas.")
    print(f"Found-set ledger: {found_path}")

def main(fasta, count_only=False, sample=None, seed=0,
         orf_map=None, s3_dest=None, outdir=OUTDIR, fmt="cif"):
    hashes = parse_fasta(fasta, orf_map_path=orf_map)
    print(f"{len(hashes)} unique sequences parsed from {fasta}")
    if orf_map:
        print(f"full orf_id -> protein_hash index written to {orf_map}")
    if sample is not None and sample < len(hashes):
        import random
        hashes = set(random.Random(seed).sample(list(hashes), sample))
        print(f"sampled down to {len(hashes)} sequences (seed={seed})")
    dsets = [(u, lance.dataset(u, storage_options=STORAGE)) for u in DATASETS]
    if count_only:
        count_hits(hashes, dsets)
    else:
        download(hashes, dsets, outdir=outdir, s3_dest=s3_dest, fmt=fmt)

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
                    help="upload <hash>.pdb objects here instead of writing them locally")
    ap.add_argument("--outdir", default=OUTDIR,
                    help=f"local output dir for structures + found_hashes.tsv (default {OUTDIR})")
    ap.add_argument("--format", dest="fmt", choices=["cif", "pdb"], default="cif",
                    help="structure file format (default cif)")
    args = ap.parse_args()
    main(args.fasta, count_only=args.count, sample=args.sample, seed=args.seed,
         orf_map=args.orf_map, s3_dest=args.s3_dest, outdir=args.outdir, fmt=args.fmt)
