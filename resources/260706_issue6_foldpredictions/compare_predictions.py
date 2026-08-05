"""
compare_predictions.py -- score ESMFold2 predictions against experimental PDB ground truth.

For each ORF in the subset map, fetch its predicted structure (STREAMED from S3, or from a local
dir for testing), load the matching experimental PDB (local, indexed by accession), and measure
structural similarity with TM-align. TM-align is a sequence-INDEPENDENT structural superposition, so
it finds the best matching sub-structure automatically -- which transparently handles the fact that
we fold the FULL ORF but the PDB only covers the catalytic core region (no manual core_range slicing
needed). We report TM-score normalised by the experimental structure (= how well the prediction
reproduces the real fold), plus RMSD and aligned/total lengths, and clean-vs-control summaries.

Layout assumed:
  * predictions in S3:  <prefix>/structures/orf<ORFid>.cif      (bare-orfid keyed, from --s3-only runs)
  * ground-truth PDBs:  <pdb-dir>/<PDB_ACCESSION>.cif           (local)
  * mapping + role:     paramsweep_subset_map.tsv               (orfid, role, pdb_id, core_range, ...)

Run once per hyperparameter config (one --s3-pred prefix -> one --out TSV); compare the per-role
median TM across configs to choose parameters. The `control` (leakage) rows are the metric's ceiling
-- if they don't come out high, the SCORING is wrong, not the model. Safe to run mid-fold: ORFs not
yet uploaded are marked `pred_missing` and skipped, so just re-run as folds finish.

Deps:   pip install tmtools biopython boto3 numpy      (S3 read needs signed AWS creds)
Example (run from this resources dir):
  python compare_predictions.py \
      --map paramsweep_subset_map.tsv \
      --pdb-dir ../../notebooks/pdb_structures \
      --s3-pred s3://petadex-protein-structures/esmfold2-centroids/paramsweep/s100_l20/ \
      --label s100_l20 --out cmp_s100_l20.tsv

NEXT (active-site metric, not wired yet): active_site_triad_rmsd() is a documented seam -- fill it
once catalytic Ser-His-Asp residues are identified per PDB (SITE/annotation parse or geometric
detection) and mapped onto the prediction via the TM-align residue alignment.
"""
import argparse
import csv
import os
import statistics
import sys
import tempfile
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# structure loading + TM-align
# ---------------------------------------------------------------------------
def load_model_chains(path):
    """{chain_id: (CA_coords Nx3, one_letter_seq)} for the first model's protein chains.

    We extract residues ourselves rather than tmtools.get_residue_data because the latter maps
    3-letter -> 1-letter through a strict dict and raises KeyError on non-standard residues -- most
    importantly 'UNK', which ESMFold2 emits wherever the ORF had an ambiguous 'X'. That KeyError
    would drop the whole (otherwise perfectly good) chain -> pred_parse_fail. Bio's seq1 maps any
    unknown residue to 'X' (undef_code) instead, and TM-align is coordinate-based, so an occasional
    'X' is harmless for the superposition. Also picks up modified residues (MSE, etc.) as 'X'."""
    from tmtools.io import get_structure
    from Bio.SeqUtils import seq1
    s = get_structure(str(path), format="mmcif")           # our cif + PDB cif are both mmCIF
    model = next(s.get_models())
    chains = {}
    for chain in model.get_chains():
        coords, seq = [], []
        for res in chain.get_residues():
            if "CA" in res:                                # protein residue (incl. UNK/MSE)
                coords.append(res["CA"].coord)
                seq.append(seq1(res.resname, undef_code="X") or "X")
        if len(seq) >= 20:                                 # ignore peptides/ions/ligand chains
            chains[chain.id] = (np.asarray(coords, dtype=float), "".join(seq))
    return chains


def best_tm(pred, ref_chains):
    """tm_align the prediction against each reference chain; keep the best (by TM on the reference).
    Returns (ref_chain_id, tm_ref, tm_pred, rmsd, ref_len) or None."""
    from tmtools import tm_align
    pc, ps = pred
    best = None
    for cid, (rc, rs) in ref_chains.items():
        try:
            r = tm_align(pc, rc, ps, rs)
        except Exception:                                  # noqa: BLE001
            continue
        tm_ref = float(r.tm_norm_chain2)                   # normalised by ref (chain2) length
        if best is None or tm_ref > best[1]:
            best = (cid, tm_ref, float(r.tm_norm_chain1), float(r.rmsd), len(rs))
    return best


def active_site_triad_rmsd(pred_path, ref_path, ref_chain):
    """SEAM (not wired): RMSD of the catalytic Ser-His-Asp triad between prediction and PDB.
    Needs the triad residue ids per PDB (parse SITE/catalytic annotations or detect geometrically:
    Ser-OG near His-NE2, His-ND1 near Asp/Glu-OD/OE), then map them onto the prediction through the
    TM-align residue alignment before superposing on the triad. Returns None until implemented."""
    return None


# ---------------------------------------------------------------------------
# prediction source: S3 stream (default) or local dir (testing)
# ---------------------------------------------------------------------------
class PredSource:
    def __init__(self, s3_pred=None, pred_dir=None):
        self.s3 = self.bucket = self.prefix = None
        self.pred_dir = Path(pred_dir) if pred_dir else None
        if s3_pred:
            import boto3
            self.s3 = boto3.client("s3")                   # signed; default cred chain (aws cli/env)
            self.bucket, _, self.prefix = s3_pred[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")

    def get_cif_text(self, oid):
        """CIF text for orf<oid>, or None if not present yet."""
        name = f"orf{oid}.cif"
        if self.pred_dir is not None:
            for p in (self.pred_dir / "structures" / name, self.pred_dir / name):
                if p.exists():
                    return p.read_text()
            return None
        key = "/".join(x for x in (self.prefix, "structures", name) if x)
        try:
            return self.s3.get_object(Bucket=self.bucket, Key=key)["Body"].read().decode()
        except Exception as e:                             # noqa: BLE001
            code = getattr(e, "response", {}).get("Error", {}).get("Code")
            if code in ("NoSuchKey", "404", "NoSuchBucket") or e.__class__.__name__ == "NoSuchKey":
                return None
            raise


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
COLS = ["orfid", "role", "pdb_id", "ref_chain", "pred_len", "ref_len",
        "tm_ref", "tm_pred", "rmsd", "status", "error"]


def score_one(d, src, pdbdir):
    oid, role, pdbid = d["orfid"], d.get("role", ""), d.get("pdb_id", "")
    rec = {c: "" for c in COLS}
    rec.update(orfid=oid, role=role, pdb_id=pdbid, status="error")

    cif = src.get_cif_text(oid)
    if cif is None:
        rec["status"] = "pred_missing"
        return rec
    refp = pdbdir / f"{pdbid}.cif"
    if not refp.exists():
        rec["status"] = "ref_missing"
        return rec

    tf = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        tf.write(cif)
        tf.close()
        pchains = load_model_chains(tf.name)
    finally:
        os.unlink(tf.name)
    if not pchains:
        rec["status"] = "pred_parse_fail"
        return rec
    pred = max(pchains.values(), key=lambda cs: len(cs[1]))   # predictions are single-chain

    rchains = load_model_chains(refp)
    if not rchains:
        rec["status"] = "ref_parse_fail"
        return rec

    b = best_tm(pred, rchains)
    if b is None:
        rec["status"] = "align_fail"
        return rec
    cid, tm_ref, tm_pred, rmsd, rlen = b
    rec.update(ref_chain=cid, pred_len=len(pred[1]), ref_len=rlen,
               tm_ref=round(tm_ref, 4), tm_pred=round(tm_pred, 4), rmsd=round(rmsd, 3), status="ok")
    return rec


def print_summary(rows, label):
    ok = [r for r in rows if r["status"] == "ok"]
    print("\n" + "=" * 62)
    print(f"config '{label}'  scored {len(ok)}/{len(rows)} ok")
    for role in ("clean", "control"):
        v = [r["tm_ref"] for r in ok if r["role"] == role and isinstance(r["tm_ref"], float)]
        if v:
            print(f"  {role:<8} n={len(v):<4} TM_ref median={statistics.median(v):.3f} "
                  f"mean={sum(v)/len(v):.3f} min={min(v):.3f} max={max(v):.3f}")
    from collections import Counter
    bad = Counter(r["status"] for r in rows if r["status"] != "ok")
    if bad:
        print("  non-ok:", ", ".join(f"{k}={n}" for k, n in bad.items()))
    if bad.get("pred_missing"):
        print("  (some predictions not in S3 yet -- rerun after the fold finishes)")
    print("=" * 62)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", required=True, help="paramsweep_subset_map.tsv (orfid, role, pdb_id, ...)")
    ap.add_argument("--pdb-dir", required=True, help="local dir of ground-truth <PDB>.cif files")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--s3-pred", metavar="s3://bucket/prefix/",
                   help="stream predictions from <prefix>/structures/orf<id>.cif")
    g.add_argument("--pred-dir", help="read predictions from a local dir instead (testing)")
    ap.add_argument("--out", required=True, help="per-ORF results TSV")
    ap.add_argument("--label", default="", help="config label for the summary")
    ap.add_argument("--limit", type=int, help="only score the first N rows")
    a = ap.parse_args()

    for m in ("tmtools", "Bio", "numpy") + (("boto3",) if a.s3_pred else ()):
        try:
            __import__(m)
        except ImportError:
            sys.exit(f"missing dependency '{m}'. Install:  pip install tmtools biopython boto3 numpy")

    src = PredSource(a.s3_pred, a.pred_dir)
    pdbdir = Path(a.pdb_dir)
    if not pdbdir.is_dir():
        sys.exit(f"--pdb-dir not found: {pdbdir}")
    rows_in = list(csv.DictReader(open(a.map), delimiter="\t"))
    if a.limit:
        rows_in = rows_in[:a.limit]
    print(f"scoring {len(rows_in)} ORFs from {a.map}  (label={a.label or '-'})")

    out = []
    for d in rows_in:
        rec = score_one(d, src, pdbdir)
        out.append(rec)
        tm = rec["tm_ref"] if rec["tm_ref"] != "" else "-"
        print(f"  orf{rec['orfid']:<9} {rec['role']:<7} {rec['pdb_id']:<6} "
              f"{rec['status']:<14} TM_ref={tm}")

    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")   # .tsv -> tab, matches the reader
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {a.out}")
    print_summary(out, a.label)


if __name__ == "__main__":
    main()
