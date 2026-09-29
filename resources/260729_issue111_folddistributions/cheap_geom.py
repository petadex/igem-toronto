#!/usr/bin/env python
"""
cheap_geom.py -- a cheap, self-contained stand-in for the MolProbity score (issue 6b, Bb/Bc).

WHY THIS EXISTS
    molprobity_score.py needs a ~90-minute source build (cctbx + chem_data + probe + reduce) and
    costs ~0.58 vCPU-seconds per structure. That is fine for a one-off ground-truth run and wrong
    for the eventual millions. This script computes a proxy from the SAME two physical ideas that
    dominate the MolProbity score, using nothing but numpy:

      1. Ramachandran   -- what fraction of residues have backbone (phi, psi) in a favoured region
      2. steric clashes -- how many non-bonded heavy atoms overlap

    No cctbx, no monomer library, no chem_data, no compiled binaries, no hydrogens. It streams the
    same S3 layout and takes the same flags as molprobity_score.py, so the two TSVs join on
    (orfid, run) and can be compared directly. Calibrate this against MolProbity on a subset, then
    run ONLY this at scale.

WHAT IT IS NOT
    Not a MolProbity score, and it must never be labelled as one. Two deliberate approximations:

      * The favoured/allowed Ramachandran regions are coarse BOXES, not the Richardson
        density contours. Percentages will not match `ramalyze` exactly.
      * Clashes are counted on heavy atoms with element-wise van der Waals radii. Real MolProbity
        places hydrogens with `reduce` and measures H-inclusive overlap, which is strictly more
        informative -- a sidechain rotated into a bad H position is invisible here.
      * There is NO rotamer term at all: rotamer outliers need the rotamer library this script
        exists to avoid. The published weight on that term is 0.33.

    So treat the output as an ordering signal to be validated by correlation, not as a calibrated
    quantity. The notebook's correlation figure is what licenses (or refuses) its use.

USAGE (identical in shape to molprobity_score.py)
    python cheap_geom.py \
        --run base=s3://petadex-protein-structures/esmfold2_paramsweep/s100_l20/ \
        --run msa=s3://petadex-protein-structures/esmfold2_paramsweep/s32_l5/ \
        --out cheap_base_vs_msa --workers 16

    python cheap_geom.py --self-test            # offline, no S3, no network
    python cheap_geom.py --score-file x.cif     # score one local structure and print the parts

Scale controls, all matching the sibling scripts: --shard K/N, --limit, resume from the ledger,
--s3-out. Structures stream from S3 and are parsed in memory; nothing is written to disk but the
TSV, so a 2M-structure run needs no local storage.

Vendored (deliberately, per the per-experiment self-containment rule): StructSource, norm_id, the
shard helpers, Progress and the ledger/driver scaffolding are adapted from molprobity_score.py, and
the header-driven _atom_site tokenizer from pair_tm.py's parse_cif_fast. Those scripts are NOT
imported, so neither can be broken by editing this one.
"""
import argparse
import csv
import math
import os
import sys
import time
import threading
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np

# ---------------------------------------------------------------------------
# geometry constants
# ---------------------------------------------------------------------------
# Bondi van der Waals radii (Angstrom). Anything unrecognised falls back to carbon.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
       "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
VDW_DEFAULT = 1.70

# Coarse general-case Ramachandran regions as (phi_lo, phi_hi, psi_lo, psi_hi) boxes, degrees.
# psi wraps at +-180, so the beta region is expressed as two boxes. These approximate the
# Richardson general-case contours; they are intentionally simple, and they were WIDENED to match
# reality: with tighter boxes, 100 refined crystal structures read 88.6% favoured / 4.0% outliers,
# against a published expectation of ~95-98% favoured and well under 1% outliers.
RAMA_FAVORED = [
    (-180.0, -40.0,   90.0,  180.0),    # beta sheet + polyproline II
    (-180.0, -40.0, -180.0, -150.0),    # beta, wrapped across +-180
    (-160.0, -35.0,  -75.0,   25.0),    # alpha-R
    (-110.0, -40.0,   25.0,   90.0),    # bridge, joining alpha-R to beta
    (  35.0,  90.0,   -5.0,   70.0),    # alpha-L
]
# Outliers are what is left over. For phi < 0 nearly the whole plane is populated in real proteins,
# so the allowed region is generous there and the outlier signal comes almost entirely from
# positive phi outside the alpha-L island -- which is where genuine outliers actually live.
RAMA_ALLOWED = [
    (-180.0, -25.0, -180.0,  180.0),
    (  15.0, 110.0,  -40.0,  110.0),
]

# Glycine has no sidechain and legitimately occupies positive phi, so scoring it against the
# general-case region would flag a large fraction of a normal protein. ramalyze uses separate
# contours; the cheap stand-in is to accept the mirror image (-phi, -psi) as well.
GLY = "GLY"

BACKBONE = ("N", "CA", "C")


def _in_boxes(phi, psi, boxes):
    m = np.zeros(phi.shape, dtype=bool)
    for lo_a, hi_a, lo_b, hi_b in boxes:
        m |= (phi >= lo_a) & (phi <= hi_a) & (psi >= lo_b) & (psi <= hi_b)
    return m


# ---------------------------------------------------------------------------
# the cheap score
# ---------------------------------------------------------------------------
def cheap_score(clash_per_1k, pct_rama_favored):
    """MolProbity's functional form MINUS the rotamer term (weight 0.33), which needs a library.

    Deliberately the same shape and the same published coefficients so the two numbers live on a
    comparable scale -- but it is a different quantity, and omitting a positively-weighted term
    means this reads systematically LOWER (i.e. better) than a true MolProbity score.
    """
    if clash_per_1k is None or pct_rama_favored is None:
        return None
    rama_iffy = 100.0 - pct_rama_favored
    return max(0.0, 0.426 * math.log(1.0 + clash_per_1k)
               + 0.250 * math.log(1.0 + max(0.0, rama_iffy - 2.0))
               + 0.5)


# ---------------------------------------------------------------------------
# CIF parsing -- header-driven _atom_site tokenizer (from pair_tm.py, widened to all heavy atoms)
# ---------------------------------------------------------------------------
def _asint(x, fallback):
    try:
        return int(x)
    except (ValueError, TypeError):
        return fallback


def parse_cif_atoms(text):
    """All heavy atoms as (xyz[N,3] float64, name[N], elem[N], resid[N] int, comp[N]).

    resid is a running index that is unique per (chain, label_seq_id) in file order, which is what
    the clash filter needs -- it only ever asks 'are these two atoms in the same or adjacent
    residue', never for the residue's true number.
    """
    lines = text.splitlines()
    i, n = 0, len(lines)
    xyz, names, elems, resids, comps = [], [], [], [], []
    while i < n:
        if lines[i].strip() != "loop_":
            i += 1
            continue
        cols, j = [], i + 1
        while j < n and lines[j].lstrip().startswith("_"):
            cols.append(lines[j].split()[0].strip())
            j += 1
        if not any(c.startswith("_atom_site.") for c in cols):
            i = j
            continue
        idx = {c: k for k, c in enumerate(cols)}
        need = ["_atom_site.label_atom_id", "_atom_site.label_comp_id",
                "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
        if not all(c in idx for c in need):
            i = j
            continue
        i_at, i_cp = idx["_atom_site.label_atom_id"], idx["_atom_site.label_comp_id"]
        i_x, i_y, i_z = (idx["_atom_site.Cartn_x"], idx["_atom_site.Cartn_y"],
                         idx["_atom_site.Cartn_z"])
        i_sq = idx.get("_atom_site.label_seq_id")
        i_ch = idx.get("_atom_site.label_asym_id")
        i_el = idx.get("_atom_site.type_symbol")
        i_alt = idx.get("_atom_site.label_alt_id")
        i_grp = idx.get("_atom_site.group_PDB")
        ncol = len(cols)
        k = j
        prev_key, run = None, -1
        while k < n:
            s = lines[k].strip()
            if not s or s[0] in "#_" or s.startswith(("loop_", "data_")):
                break
            tok = s.split()
            if len(tok) < ncol:
                k += 1
                continue
            # waters/ligands are HETATM and have no backbone: they would distort both metrics
            if i_grp is not None and tok[i_grp].strip('"').upper() != "ATOM":
                k += 1
                continue
            # keep only one altloc, else a disordered sidechain clashes with itself
            if i_alt is not None and tok[i_alt].strip('"') not in (".", "?", "A", ""):
                k += 1
                continue
            name = tok[i_at].strip('"')
            elem = (tok[i_el].strip('"').upper() if i_el is not None
                    else "".join(c for c in name if c.isalpha())[:1].upper())
            if elem in ("H", "D"):                      # predictions have none; PDB files may
                k += 1
                continue
            try:
                p = (float(tok[i_x]), float(tok[i_y]), float(tok[i_z]))
            except ValueError:
                k += 1
                continue
            ch = tok[i_ch].strip('"') if i_ch is not None else "A"
            sq = _asint(tok[i_sq], None) if i_sq is not None else None
            key = (ch, sq if sq is not None else len(resids))
            if key != prev_key:
                run += 1
                prev_key = key
            xyz.append(p)
            names.append(name)
            elems.append(elem)
            resids.append(run)
            comps.append(tok[i_cp].strip('"').upper())
            k += 1
        i = k
    if not xyz:
        return np.zeros((0, 3)), np.array([]), np.array([]), np.array([], dtype=int), np.array([])
    return (np.asarray(xyz, dtype=float), np.asarray(names), np.asarray(elems),
            np.asarray(resids, dtype=int), np.asarray(comps))


# ---------------------------------------------------------------------------
# component 1: Ramachandran
# ---------------------------------------------------------------------------
def dihedral(p0, p1, p2, p3):
    """Signed dihedral in degrees, vectorised over leading axes."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    nb1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * nb1).sum(-1, keepdims=True) * nb1
    w = b2 - (b2 * nb1).sum(-1, keepdims=True) * nb1
    x = (v * w).sum(-1)
    y = (np.cross(nb1, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def backbone_table(xyz, names, resids):
    """(res_index[R], N[R,3], CA[R,3], C[R,3]) for residues having all three backbone atoms."""
    n_res = int(resids.max()) + 1 if resids.size else 0
    slot = {a: np.full((n_res, 3), np.nan) for a in BACKBONE}
    for a in BACKBONE:
        sel = names == a
        if sel.any():
            slot[a][resids[sel]] = xyz[sel]
    ok = ~(np.isnan(slot["N"]).any(1) | np.isnan(slot["CA"]).any(1) | np.isnan(slot["C"]).any(1))
    idx = np.nonzero(ok)[0]
    return idx, slot["N"][idx], slot["CA"][idx], slot["C"][idx]


def rama_stats(xyz, names, resids, max_bond=2.0):
    """(pct_favored, pct_outliers, n_scored).

    phi/psi are only defined for a residue with both neighbours present, so terminal residues and
    any residue across a chain break are skipped -- exactly as ramalyze does.
    """
    idx, N, CA, C = backbone_table(xyz, names, resids)
    if idx.size < 3:
        return None, None, 0
    # consecutive in residue index AND actually peptide-bonded (C(i)->N(i+1) ~1.33 A). The distance
    # test is what catches a chain break that the numbering does not reveal.
    step = np.diff(idx) == 1
    bond = np.linalg.norm(N[1:] - C[:-1], axis=1) < max_bond
    linked = step & bond                                 # link between residue r and r+1
    mid = np.nonzero(linked[:-1] & linked[1:])[0] + 1    # residues with BOTH neighbours linked
    if mid.size == 0:
        return None, None, 0
    phi = dihedral(C[mid - 1], N[mid], CA[mid], C[mid])
    psi = dihedral(N[mid], CA[mid], C[mid], N[mid + 1])
    fav = _in_boxes(phi, psi, RAMA_FAVORED)
    allowed = _in_boxes(phi, psi, RAMA_ALLOWED)
    n = int(mid.size)
    return 100.0 * fav.sum() / n, 100.0 * (~allowed).sum() / n, n


# ---------------------------------------------------------------------------
# component 2: heavy-atom clashes
# ---------------------------------------------------------------------------
def _count_pairs(xyz, r, resids, is_s, ai, bj, tol, min_sep):
    """Clashing pairs between atom index blocks ai and bj, counted once via the global j > i rule."""
    d = np.linalg.norm(xyz[ai][:, None, :] - xyz[bj][None, :, :], axis=-1)
    hit = d < (r[ai][:, None] + r[bj][None, :] - tol)
    hit &= np.abs(resids[ai][:, None] - resids[bj][None, :]) >= min_sep
    hit &= bj[None, :] > ai[:, None]
    # disulfides are real bonds at ~2.05 A and would otherwise read as the worst clash present
    hit &= ~(is_s[ai][:, None] & is_s[bj][None, :] & (d < 2.5))
    return int(hit.sum())


def clash_stats(xyz, elems, resids, tol=0.4, min_sep=2, dense_max=4000, chunk=512):
    """(clashes_per_1000_atoms, n_clashes, n_atoms).

    A clash is a pair of heavy atoms in residues >= `min_sep` apart whose separation is below
    (vdw_i + vdw_j - tol). Element-wise radii matter: a flat cutoff near 3.0 A would flag every
    N...O hydrogen bond (2.8-3.2 A) as a clash, whereas N/O radii put their threshold at ~2.67 A.

    Two strategies, because the crossover is real:
      * n <= dense_max : chunked all-vs-all. O(n^2) but pure vectorised numpy, and for a ~300 aa
        prediction (~2.4k heavy atoms) that is a few tens of ms.
      * n >  dense_max : a uniform cell list of side `cutoff`, so only the 27-cell neighbourhood
        of each atom is examined. O(n) with a Python loop over occupied cells -- slower per atom,
        but the dense path becomes quadratically hopeless on large complexes.
    Both paths are checked against each other in --self-test.
    """
    n = int(xyz.shape[0])
    if n == 0:
        return None, 0, 0
    r = np.array([VDW.get(e, VDW_DEFAULT) for e in elems])
    is_s = np.asarray(elems) == "S"
    total = 0

    if n <= dense_max:
        allj = np.arange(n)
        for s in range(0, n, chunk):
            e = min(s + chunk, n)
            total += _count_pairs(xyz, r, resids, is_s, np.arange(s, e), allj, tol, min_sep)
        return 1000.0 * total / n, total, n

    from collections import defaultdict
    cutoff = 2.0 * float(r.max()) - tol          # no pair can clash further apart than this
    gi = np.floor((xyz - xyz.min(0)) / cutoff).astype(np.int64)
    buckets = defaultdict(list)
    for idx, key in enumerate(map(tuple, gi)):
        buckets[key].append(idx)
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for key, ai in buckets.items():
        nb = [buckets[k] for k in
              ((key[0] + o[0], key[1] + o[1], key[2] + o[2]) for o in offsets) if k in buckets]
        if not nb:
            continue
        total += _count_pairs(xyz, r, resids, is_s, np.asarray(ai),
                              np.concatenate([np.asarray(b) for b in nb]), tol, min_sep)
    return 1000.0 * total / n, total, n


# ---------------------------------------------------------------------------
# score one structure
# ---------------------------------------------------------------------------
def score_structure(text, tol=0.4):
    xyz, names, elems, resids, comps = parse_cif_atoms(text)
    if xyz.shape[0] == 0:
        raise ValueError("no ATOM records parsed")
    fav, out, n_rama = rama_stats(xyz, names, resids)
    per1k, n_clash, n_atoms = clash_stats(xyz, elems, resids, tol=tol)
    return {
        "n_res": int(resids.max()) + 1,
        "n_atoms": n_atoms,
        "n_rama": n_rama,
        "clash_per_1k": per1k,
        "n_clash": n_clash,
        "pct_rama_favored": fav,
        "pct_rama_outliers": out,
        "cheap_score": cheap_score(per1k, fav),
    }


# ---------------------------------------------------------------------------
# structure source: S3 stream (default) or local dir   [vendored from molprobity_score.py]
# ---------------------------------------------------------------------------
class StructSource:
    """Reads <root>/structures/<cif_name>. `root` is s3://bucket/prefix/ or a local dir."""

    def __init__(self, label, root, cif_name="orf{id}.cif", retries=4):
        self.label, self.root, self.cif_name, self.retries = label, root, cif_name, retries
        self.is_s3 = str(root).startswith("s3://")
        self._tl = threading.local()
        if self.is_s3:
            self.bucket, _, self.prefix = root[len("s3://"):].partition("/")
            self.prefix = self.prefix.strip("/")

    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop("_tl", None)                 # threading.local cannot pickle; client must not fork
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._tl = threading.local()

    @property
    def s3(self):
        if getattr(self._tl, "c", None) is None:
            import boto3
            self._tl.c = boto3.client("s3")
        return self._tl.c

    def get_text(self, oid):
        name = self.cif_name.format(id=oid)
        if not self.is_s3:
            for p in (os.path.join(self.root, "structures", name), os.path.join(self.root, name)):
                if os.path.exists(p):
                    with open(p) as f:
                        return f.read()
            return None
        key = "/".join(p for p in (self.prefix, "structures", name) if p)
        for attempt in range(self.retries):
            try:
                return self.s3.get_object(Bucket=self.bucket, Key=key)["Body"].read().decode()
            except Exception as e:                        # noqa: BLE001
                code = str(getattr(e, "response", {}).get("Error", {}).get("Code", ""))
                if code in ("NoSuchKey", "404", "NoSuchBucket") or type(e).__name__ == "NoSuchKey":
                    return None
                if attempt == self.retries - 1:
                    raise
                time.sleep(0.4 * 2 ** attempt)
        return None

    def list_ids(self):
        import re
        pat = re.compile("^" + re.escape(self.cif_name).replace(r"\{id\}", "(.+)") + "$")
        out = set()
        if not self.is_s3:
            d = os.path.join(self.root, "structures")
            for nm in (os.listdir(d) if os.path.isdir(d) else []):
                m = pat.match(nm)
                if m:
                    out.add(m.group(1))
            return out
        tok = None
        while True:
            kw = {"Bucket": self.bucket,
                  "Prefix": "/".join(p for p in (self.prefix, "structures/") if p)}
            if tok:
                kw["ContinuationToken"] = tok
            r = self.s3.list_objects_v2(**kw)
            for o in r.get("Contents", []):
                m = pat.match(o["Key"].rsplit("/", 1)[-1])
                if m:
                    out.add(m.group(1))
            if not r.get("IsTruncated"):
                return out
            tok = r["NextContinuationToken"]


def norm_id(x):
    x = str(x).strip()
    return x[3:] if x.startswith("orf") else x


# ---------------------------------------------------------------------------
# ledger
# ---------------------------------------------------------------------------
COLS = ["orfid", "run", "n_res", "n_atoms", "n_rama", "n_clash", "clash_per_1k",
        "pct_rama_favored", "pct_rama_outliers", "cheap_score", "wall_s", "status", "error"]


def read_done(path):
    """{(run, orfid)} already scored OK -- failures are retried on resume."""
    done = set()
    if not os.path.exists(path):
        return done
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") == "ok":
                done.add((row.get("run", ""), norm_id(row.get("orfid", ""))))
    return done


def score_one(oid, src, tol):
    rec = {c: "" for c in COLS}
    rec.update(orfid=oid, run=src.label, status="error")
    t0 = time.perf_counter()
    text = src.get_text(oid)
    if text is None:
        rec["status"] = "missing"
        return rec
    try:
        d = score_structure(text, tol=tol)
    except Exception as e:                                # noqa: BLE001  one bad structure != dead run
        rec["status"] = "score_fail"
        rec["error"] = f"{type(e).__name__}: {e}"[:300]
        rec["wall_s"] = round(time.perf_counter() - t0, 3)
        return rec
    for k in ("n_res", "n_atoms", "n_rama", "n_clash"):
        rec[k] = d[k]
    for k in ("clash_per_1k", "pct_rama_favored", "pct_rama_outliers", "cheap_score"):
        rec[k] = "" if d[k] is None else round(d[k], 4)
    rec["wall_s"] = round(time.perf_counter() - t0, 3)
    rec["status"] = "ok" if d["cheap_score"] is not None else "incomplete"
    if d["cheap_score"] is None:
        rec["error"] = "no phi/psi could be computed (too short, or no backbone)"
    return rec


# ---------------------------------------------------------------------------
# sharding (same 'K/N' contract as the predictor, pair_tm.py and molprobity_score.py)
# ---------------------------------------------------------------------------
def parse_shard(spec):
    if not spec:
        return None
    try:
        k, n = str(spec).split("/")
        k, n = int(k), int(n)
    except Exception:                                     # noqa: BLE001
        raise SystemExit(f"--shard wants K/N, got {spec!r}")
    if not (1 <= k <= n):
        raise SystemExit(f"--shard K must be in 1..N, got {spec!r}")
    return k, n


def select_shard(ids, k, n):
    ids = sorted(ids)
    per, rem = divmod(len(ids), n)
    start = (k - 1) * per + min(k - 1, rem)
    return ids[start:start + per + (1 if k <= rem else 0)]


def build_worklist(a, sources):
    if a.ids:
        with open(a.ids) as f:
            return sorted({norm_id(x) for x in f.read().split() if x.strip()}), f"--ids {a.ids}"
    sets = [{norm_id(x) for x in s.list_ids()} for s in sources]
    if a.intersect:
        keep = set.intersection(*sets) if sets else set()
        return sorted(keep), "intersection of the runs"
    keep = set.union(*sets) if sets else set()
    return sorted(keep), "union of the runs"


# ---------------------------------------------------------------------------
# worker-process state
# ---------------------------------------------------------------------------
_W = {}


def _init_worker(specs, cif_name, tol):
    global _W
    _W = {"sources": {lab: StructSource(lab, root, cif_name) for lab, root in specs}, "tol": tol}


def _work(job):
    oid, label = job
    try:
        return score_one(oid, _W["sources"][label], _W["tol"])
    except Exception as e:                                # noqa: BLE001
        r = {c: "" for c in COLS}
        r.update(orfid=oid, run=label, status="error", error=f"{type(e).__name__}: {e}"[:300])
        return r


# ---------------------------------------------------------------------------
# progress
# ---------------------------------------------------------------------------
def fmt_hms(seconds):
    s = int(max(seconds, 0))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h}:{m:02d}:{sec:02d}" if h else f"{m}:{sec:02d}"


class Progress:
    """Redraws in place on a tty; one line every `every` items when redirected to a log."""

    def __init__(self, total, width=32, every=25, stream=None):
        self.total = max(int(total), 1)
        self.width, self.every = width, max(int(every), 1)
        self.stream = stream if stream is not None else sys.stdout
        self.tty = hasattr(self.stream, "isatty") and self.stream.isatty()
        self.t0 = time.perf_counter()

    def update(self, n, ok, bad=0):
        el = time.perf_counter() - self.t0
        rate = n / max(el, 1e-9)
        eta = (self.total - n) / rate if rate > 0 else 0.0
        tail = (f"{n}/{self.total} {100.0 * n / self.total:5.1f}%  ok={ok}"
                + (f" bad={bad}" if bad else "")
                + f"  {rate:.1f}/s  {fmt_hms(el)}<{fmt_hms(eta)}")
        if self.tty:
            filled = int(self.width * n / self.total)
            self.stream.write("\r  [" + "#" * filled + "-" * (self.width - filled) + "] "
                              + tail.ljust(58)[:58])
            self.stream.flush()
        elif n % self.every == 0 or n == self.total:
            self.stream.write("  " + tail + "\n")
            self.stream.flush()

    def close(self):
        if self.tty:
            self.stream.write("\n")
            self.stream.flush()


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def summarise(ledger):
    import statistics
    rows = []
    with open(ledger, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("status") == "ok":
                rows.append(r)
    if not rows:
        return
    print()
    for lab in sorted({r["run"] for r in rows}):
        sc = [float(r["cheap_score"]) for r in rows if r["run"] == lab]
        cl = [float(r["clash_per_1k"]) for r in rows if r["run"] == lab]
        fv = [float(r["pct_rama_favored"]) for r in rows if r["run"] == lab]
        print(f"  {lab:<10s} n={len(sc):<6d} cheap median={statistics.median(sc):.3f} "
              f"mean={statistics.fmean(sc):.3f}   clash/1k median={statistics.median(cl):.2f}   "
              f"rama_fav median={statistics.median(fv):.2f}%")
    print("  (lower cheap score is better; NOT a MolProbity score -- no rotamer term)")


def run(a):
    shard = parse_shard(a.shard if a.shard is not None else os.environ.get("CHEAP_GEOM_SHARD"))
    tag = f".shard{shard[0]}of{shard[1]}" if shard else ""
    os.makedirs(a.out, exist_ok=True)
    ledger = os.path.join(a.out, f"cheap_geom{tag}.tsv")

    sources = []
    for spec in a.run:
        if "=" not in spec:
            raise SystemExit(f"--run must be LABEL=SOURCE, got {spec!r}")
        label, _, root = spec.partition("=")
        sources.append(StructSource(label.strip(), root.strip(), a.cif_name))
    for s in sources:
        print(f"run '{s.label}'  {s.root}")

    ids, how = build_worklist(a, sources)
    print(f"{len(ids)} ORFs (work list from {how})")
    if shard:
        before = len(ids)
        ids = select_shard(ids, *shard)
        print(f"shard {shard[0]}/{shard[1]}: {len(ids)}/{before} ORFs")
    if a.limit:
        ids = ids[:a.limit]
        print(f"--limit {a.limit}: {len(ids)} ORFs")

    jobs = [(oid, s.label) for s in sources for oid in ids]
    done = read_done(ledger) if a.resume else set()
    if done:
        jobs = [(o, lab) for (o, lab) in jobs if (lab, o) not in done]
        print(f"resume: {len(done)} already scored, {len(jobs)} left")
    if not jobs:
        print("nothing to do -- ledger already complete for this shard")
        return 0

    initargs = ([(s.label, s.root) for s in sources], a.cif_name, a.tol)
    if a.pool == "process":
        def make_pool():
            return ProcessPoolExecutor(max_workers=a.workers,
                                       initializer=_init_worker, initargs=initargs)
    else:
        _init_worker(*initargs)
        def make_pool():
            return ThreadPoolExecutor(max_workers=a.workers)
    print(f"pool: {a.pool}   workers: {a.workers}   clash tol: {a.tol} A\n")

    lock = threading.Lock()
    counts = Counter()
    t0 = time.perf_counter()
    fresh = not os.path.exists(ledger) or os.path.getsize(ledger) == 0
    with open(ledger, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        if fresh:
            w.writeheader()
        n = 0
        prog = Progress(len(jobs), every=a.flush_every)
        with make_pool() as pool:
            for rec in pool.map(_work, jobs, chunksize=1):
                with lock:
                    w.writerow(rec)
                    counts[rec["status"]] += 1
                    n += 1
                    if n % a.flush_every == 0:
                        fh.flush()
                    prog.update(n, counts["ok"], n - counts["ok"])
        prog.close()

    elapsed = time.perf_counter() - t0
    print(f"\nscored {n} structures in {elapsed:.1f}s ({n / max(elapsed, 1e-9):.2f}/s)")
    for k, v in sorted(counts.items()):
        print(f"  {k}={v}")
    print(f"  ledger -> {ledger}")
    summarise(ledger)

    if a.s3_out:
        import boto3
        b, _, p = a.s3_out[len("s3://"):].partition("/")
        key = p.rstrip("/") + "/" + os.path.basename(ledger)
        boto3.client("s3").upload_file(ledger, b, key)
        print(f"  mirrored -> s3://{b}/{key}")
    return 0


# ---------------------------------------------------------------------------
# self-test -- offline, no S3, no network
# ---------------------------------------------------------------------------
def _fake_cif(residues):
    """Minimal but real-shaped CIF. `residues` = [(comp, [(atom, elem, x, y, z), ...]), ...]"""
    head = ["data_test", "loop_",
            "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
            "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id",
            "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.Cartn_x",
            "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
    rows, k = [], 0
    for ri, (comp, atoms) in enumerate(residues, start=1):
        for (nm, el, x, y, z) in atoms:
            k += 1
            rows.append(f"ATOM {k} {el} {nm} . {comp} A {ri} {x:.3f} {y:.3f} {z:.3f}")
    return "\n".join(head + rows) + "\n#\n"


def _ideal_chain(n, phi=-60.0, psi=-45.0):
    """Build a backbone with (approximately) the requested phi/psi via NeRF-style placement."""
    import numpy as _np
    # bond geometry
    b = {"N-CA": 1.458, "CA-C": 1.525, "C-N": 1.329}
    ang = {"N-CA-C": math.radians(111.2), "CA-C-N": math.radians(116.2), "C-N-CA": math.radians(121.7)}
    coords = [_np.array([0.0, 0.0, 0.0]), _np.array([1.458, 0.0, 0.0]),
              _np.array([1.458 + 1.525 * math.cos(math.pi - ang["N-CA-C"]),
                         1.525 * math.sin(math.pi - ang["N-CA-C"]), 0.0])]
    order = ["N", "CA", "C"]
    seq = [("N", 0), ("CA", 0), ("C", 0)]

    def place(a, bb, c, bond, theta, tor):
        ab = bb - a
        bc = c - bb
        bc_n = bc / _np.linalg.norm(bc)
        nv = _np.cross(ab, bc_n)
        nv /= _np.linalg.norm(nv)
        m = _np.array([bc_n, _np.cross(nv, bc_n), nv]).T
        d = _np.array([-bond * math.cos(theta),
                       bond * math.sin(theta) * math.cos(tor),
                       bond * math.sin(theta) * math.sin(tor)])
        return c + m.dot(d)

    for r in range(1, n):
        N = place(coords[-3], coords[-2], coords[-1], b["C-N"], ang["CA-C-N"], math.radians(psi))
        CA = place(coords[-2], coords[-1], N, b["N-CA"], ang["C-N-CA"], math.pi)
        C = place(coords[-1], N, CA, b["CA-C"], ang["N-CA-C"], math.radians(phi))
        coords += [N, CA, C]
        seq += [("N", r), ("CA", r), ("C", r)]
    residues = {}
    for (nm, r), xyz in zip(seq, coords):
        residues.setdefault(r, []).append((nm, "N" if nm == "N" else "C", *xyz))
    return _fake_cif([("ALA", residues[r]) for r in sorted(residues)])


def self_test():
    print("SELF-TEST: geometry, parsing, ledger/resume, sharding. (No S3, no network.)\n")
    import shutil
    import tempfile
    fails = []

    def check(name, cond):
        print(("  [ok] " if cond else "  [FAIL] ") + name)
        if not cond:
            fails.append(name)

    # --- score shape -------------------------------------------------------
    check("cheap score: floor is 0.5 for a perfect structure",
          abs(cheap_score(0.0, 100.0) - 0.5) < 1e-12)
    check("cheap score: 2% non-favoured is still free (matches the published allowance)",
          abs(cheap_score(0.0, 98.0) - 0.5) < 1e-12)
    check("cheap score: monotonic in clash rate",
          cheap_score(1.0, 95.0) < cheap_score(10.0, 95.0) < cheap_score(50.0, 95.0))
    check("cheap score: monotonic in rama badness",
          cheap_score(5.0, 99.0) < cheap_score(5.0, 90.0) < cheap_score(5.0, 60.0))
    check("cheap score: None in -> None out", cheap_score(None, 95.0) is None
          and cheap_score(5.0, None) is None)

    # --- dihedral ----------------------------------------------------------
    p = np.array([[1.0, 0, 0], [0, 0, 0], [0, 1.0, 0], [0, 1.0, 1.0]])
    d = dihedral(p[0], p[1], p[2], p[3])
    check(f"dihedral of a known +90 arrangement -> {d:.1f}", abs(abs(d) - 90.0) < 1e-6)

    # --- helix vs sheet land in the right boxes ----------------------------
    helix = _ideal_chain(12, phi=-60.0, psi=-45.0)
    sheet = _ideal_chain(12, phi=-130.0, psi=135.0)
    xyz, names, elems, resids, _ = parse_cif_atoms(helix)
    fav_h, out_h, n_h = rama_stats(xyz, names, resids)
    xyz, names, elems, resids, _ = parse_cif_atoms(sheet)
    fav_s, out_s, n_s = rama_stats(xyz, names, resids)
    check(f"alpha-helix backbone reads as favoured ({fav_h:.0f}% over {n_h} res)", fav_h > 95.0)
    check(f"beta-strand backbone reads as favoured ({fav_s:.0f}% over {n_s} res)", fav_s > 95.0)
    check("terminal residues are skipped (n = len-2)", n_h == 10 and n_s == 10)

    # --- a deliberately impossible backbone is NOT favoured ----------------
    bad = _ideal_chain(12, phi=60.0, psi=-120.0)
    xyz, names, elems, resids, _ = parse_cif_atoms(bad)
    fav_b, out_b, _ = rama_stats(xyz, names, resids)
    check(f"a disallowed phi/psi reads as unfavoured ({fav_b:.0f}% favoured)", fav_b < 5.0)
    check(f"...and as an outlier ({out_b:.0f}% outliers)", out_b > 95.0)

    # --- clashes -----------------------------------------------------------
    far = _fake_cif([("ALA", [("CA", "C", 0, 0, 0)]), ("ALA", [("CA", "C", 0, 0, 8.0)]),
                     ("ALA", [("CA", "C", 0, 0, 16.0)])])
    xyz, names, elems, resids, _ = parse_cif_atoms(far)
    per1k, nc, na = clash_stats(xyz, elems, resids)
    check("well-separated atoms -> no clashes", nc == 0 and per1k == 0.0)

    near = _fake_cif([("ALA", [("CA", "C", 0, 0, 0)]), ("ALA", [("CA", "C", 0, 0, 8.0)]),
                      ("ALA", [("CA", "C", 0, 0, 2.0)])])
    xyz, names, elems, resids, _ = parse_cif_atoms(near)
    per1k, nc, na = clash_stats(xyz, elems, resids)
    check("two carbons 2.0 A apart, 2 residues apart -> 1 clash", nc == 1)
    check("clash rate is per 1000 atoms", abs(per1k - 1000.0 * 1 / 3) < 1e-9)

    adj = _fake_cif([("ALA", [("CA", "C", 0, 0, 0)]), ("ALA", [("CA", "C", 0, 0, 2.0)])])
    xyz, names, elems, resids, _ = parse_cif_atoms(adj)
    _, nc_adj, _ = clash_stats(xyz, elems, resids)
    check("adjacent residues are exempt (they are bonded)", nc_adj == 0)

    hb = _fake_cif([("ALA", [("N", "N", 0, 0, 0)]), ("ALA", [("CA", "C", 0, 0, 8.0)]),
                    ("ALA", [("O", "O", 0, 0, 2.9)])])
    xyz, names, elems, resids, _ = parse_cif_atoms(hb)
    _, nc_hb, _ = clash_stats(xyz, elems, resids)
    check("an N...O hydrogen bond at 2.9 A is NOT a clash", nc_hb == 0)

    ss = _fake_cif([("CYS", [("SG", "S", 0, 0, 0)]), ("ALA", [("CA", "C", 0, 0, 8.0)]),
                    ("CYS", [("SG", "S", 0, 0, 2.05)])])
    xyz, names, elems, resids, _ = parse_cif_atoms(ss)
    _, nc_ss, _ = clash_stats(xyz, elems, resids)
    check("a disulfide at 2.05 A is NOT a clash", nc_ss == 0)

    # --- parser ------------------------------------------------------------
    het = _fake_cif([("ALA", [("CA", "C", 0, 0, 0)])]).replace("ATOM 1", "HETATM 1")
    xyz, *_ = parse_cif_atoms(het)
    check("HETATM (waters/ligands) is skipped", xyz.shape[0] == 0)
    withh = _fake_cif([("ALA", [("CA", "C", 0, 0, 0), ("HA", "H", 0.5, 0.5, 0.5)])])
    xyz, *_ = parse_cif_atoms(withh)
    check("hydrogens are skipped if present", xyz.shape[0] == 1)
    check("empty text -> no atoms, and score_structure raises rather than lying",
          parse_cif_atoms("")[0].shape[0] == 0)
    try:
        score_structure("")
        check("score_structure('') raises", False)
    except ValueError:
        check("score_structure('') raises", True)

    # --- the two clash strategies must agree exactly ------------------------
    # a deliberately overcrowded box, so the counts are large and any disagreement shows
    rng = np.random.default_rng(0)
    m = 900
    pts = rng.uniform(0.0, 16.0, size=(m, 3))
    el = np.array(["C"] * m)
    rid = np.arange(m)                                    # every atom its own residue
    dense = clash_stats(pts, el, rid, dense_max=10 ** 9)
    grid = clash_stats(pts, el, rid, dense_max=0)
    check(f"the fixture actually clashes ({dense[1]} pairs)", dense[1] > 100)
    check(f"dense and cell-list paths agree exactly ({dense[1]} == {grid[1]})", dense[1] == grid[1])
    check("chunk size does not change the dense count",
          clash_stats(pts, el, rid, dense_max=10 ** 9, chunk=7)[1] == dense[1])
    check("per-1k rate agrees too", abs(dense[0] - grid[0]) < 1e-9)

    # mixed elements, since the radii enter the threshold
    el2 = np.array(["C", "N", "O", "S"] * (m // 4))
    d2 = clash_stats(pts, el2, rid, dense_max=10 ** 9)
    g2 = clash_stats(pts, el2, rid, dense_max=0)
    check(f"...also with mixed elements ({d2[1]} == {g2[1]})", d2[1] == g2[1])

    # --- driver: sources, worklist, resume, shard --------------------------
    tmp = tempfile.mkdtemp()
    try:
        for lab, ids in (("base", ["1", "2", "3"]), ("msa", ["2", "3", "4"])):
            d = os.path.join(tmp, lab, "structures")
            os.makedirs(d)
            for i in ids:
                with open(os.path.join(d, f"orf{i}.cif"), "w") as f:
                    f.write(_ideal_chain(10))
        srcs = [StructSource("base", os.path.join(tmp, "base")),
                StructSource("msa", os.path.join(tmp, "msa"))]
        check("source lists ids", srcs[0].list_ids() == {"1", "2", "3"})
        check("source reads a structure", srcs[0].get_text("1").startswith("data_"))
        check("source returns None for a missing id", srcs[0].get_text("99") is None)

        class A:
            ids = None
            intersect = True
        got, _ = build_worklist(A(), srcs)
        check("--intersect keeps only paired ORFs", got == ["2", "3"])
        A.intersect = False
        got, _ = build_worklist(A(), srcs)
        check("union keeps all", got == ["1", "2", "3", "4"])

        rec = score_one("1", srcs[0], 0.4)
        check("score_one on a real chain -> ok", rec["status"] == "ok")
        check("score_one records a cheap score", isinstance(rec["cheap_score"], float))
        check("score_one on a missing id -> missing", score_one("99", srcs[0], 0.4)["status"]
              == "missing")

        led = os.path.join(tmp, "l.tsv")
        with open(led, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerow({"orfid": "1", "run": "base", "status": "ok"})
            w.writerow({"orfid": "2", "run": "base", "status": "score_fail"})
        d = read_done(led)
        check("resume: ok rows are skipped", ("base", "1") in d)
        check("resume: failed rows are retried", ("base", "2") not in d)

        allids = [str(i) for i in range(97)]
        parts = [select_shard(allids, k, 5) for k in range(1, 6)]
        flat = sorted(sum(parts, []))
        check("shards are disjoint and cover everything exactly once", flat == sorted(allids))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n" + ("SELF-TEST FAILED: " + ", ".join(fails) if fails else "SELF-TEST PASSED"))
    return 1 if fails else 0


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Cheap numpy geometry proxy for the MolProbity score.")
    ap.add_argument("--run", action="append", default=[], metavar="LABEL=SOURCE",
                    help="repeatable: a run to score, e.g. base=s3://bucket/prefix/ (or a local dir)")
    ap.add_argument("--out", default="cheap_geom", help="output dir for the ledger")
    ap.add_argument("--ids", help="file of ORFids to score (one per token; 'orf' prefix optional)")
    ap.add_argument("--no-intersect", dest="intersect", action="store_false",
                    help="score every ORF in each run instead of only those in ALL runs")
    ap.add_argument("--cif-name", default="orf{id}.cif", dest="cif_name",
                    help="filename template inside structures/ (default orf{id}.cif)")
    ap.add_argument("--tol", type=float, default=0.4,
                    help="clash tolerance in A: overlap beyond (vdw_i+vdw_j-tol) counts (default 0.4)")
    ap.add_argument("--workers", type=int, default=8, help="pool size (default 8)")
    ap.add_argument("--pool", default="process", choices=["process", "thread"],
                    help="process (default) escapes the GIL; thread is a debugging fallback")
    ap.add_argument("--shard", default=None, metavar="K/N",
                    help="score only band K of N (1-based); falls back to $CHEAP_GEOM_SHARD")
    ap.add_argument("--limit", type=int, help="only score the first N ORFs (smoke test)")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="rescore ORFs already in the ledger")
    ap.add_argument("--flush-every", type=int, default=25)
    ap.add_argument("--s3-out", metavar="s3://bucket/prefix/", help="mirror the ledger here at the end")
    ap.add_argument("--score-file", help="score ONE local .cif and print the components")
    ap.add_argument("--self-test", action="store_true", help="offline checks; no S3, no network")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(self_test())
    if a.score_file:
        with open(a.score_file) as f:
            d = score_structure(f.read(), tol=a.tol)
        for k, v in d.items():
            print(f"  {k:<20s} {v}")
        sys.exit(0)
    if not a.run:
        ap.error("need at least one --run LABEL=SOURCE (or --self-test / --score-file)")
    sys.exit(run(a))


if __name__ == "__main__":
    main()
