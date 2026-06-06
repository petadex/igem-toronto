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

Install:  pip install pylance pyarrow msgpack numpy brotli
"""
import os, sys, hashlib, time
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

# --- FASTA -> {md5: id} -----------------------------------------------------
def read_fasta(path):
    h2id, hid, seq = {}, None, []
    def flush():
        if hid and seq:
            s = "".join(seq).upper()
            h2id.setdefault(hashlib.md5(s.encode()).hexdigest(), hid)

    opener = (zstd.open if path.endswith(".zst") else open)
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                flush(); 
                hid = line[1:].split()[0]; 
                seq = []
            elif line:
                seq.append(line)
    flush()
    return h2id

def batched(it, n):
    it = list(it)
    for i in range(0, len(it), n):
        yield it[i:i+n]

def main(fasta):
    os.makedirs(OUTDIR, exist_ok=True)
    h2id = read_fasta(fasta)
    done = {fn[:-4] for fn in os.listdir(OUTDIR) if fn.endswith(".pdb")}
    # checks which protein accessions are still missing from the Atlas, and which are already present in OUTDIR/ (from previous runs)
    todo = [h for h, pid in h2id.items() if pid not in done]
    print(f"{len(h2id)} unique sequences | {len(todo)} still to fetch")
    dsets = [lance.dataset(u, storage_options=STORAGE) for u in DATASETS]
    found = 0; t0 = time.time()
    remaining = set(todo)
    for ds in dsets:                                  # try representative set first, then 1B
        if not remaining:
            break
        for chunk in batched(remaining, BATCH):
            q = ",".join(f"'{h}'" for h in chunk)
            tbl = ds.scanner(columns=["protein_hash", "structure_blob"],
                             filter=f"protein_hash IN ({q})").to_table().to_pylist()
            for r in tbl:
                h = r["protein_hash"]
                pid = h2id.get(h)
                if not pid:
                    continue
                pid = pid.replace("/", "_")
                with open(os.path.join(OUTDIR, f"{pid}.pdb"), "w") as fh:
                    fh.write(blob_to_pdb(r["structure_blob"]))
                remaining.discard(h); found += 1
            print(f"  found {found} | {found/max(time.time()-t0,1e-9):.0f}/s | "
                  f"{len(remaining)} unmatched", end="\r")
    print(f"\nDone: {found} structures written to {OUTDIR}/ ; "
          f"{len(remaining)} sequences not in the atlas.")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "proteins.fasta")
