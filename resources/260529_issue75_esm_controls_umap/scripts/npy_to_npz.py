#!/usr/bin/env python3
"""Bundle an embed_shard.py output (float16 .npy + parallel ids .txt) into a
single .npz, matching the embed_fasta.py schema (accession_ids, embeddings) and
adding a parallel `variants` array recovered from the id prefixes.

The control ids look like  <variant>_<localid>__<globalidx>  e.g.
    rand_95th_family_0__0        -> variant "rand_95th_family"
    rand_fragment_21633__194188  -> variant "rand_fragment"
    rand_uniprot_Q5JDK7__388380  -> variant "rand_uniprot"   (localid is an accession)
    shuffled_64918__517839       -> variant "shuffled"
so the variant is the id with the trailing __<globalidx> and then the trailing
_<localid> token removed.

Usage:
    python npy_to_npz.py --emb controls.npy --ids controls_ids.txt --out controls.npz [--compress]
"""
import argparse
import re

import numpy as np

_GLOBAL = re.compile(r"__\d+$")   # trailing global index
_LOCAL = re.compile(r"_[^_]+$")   # trailing local id / accession token


def variant_of(acc):
    return _LOCAL.sub("", _GLOBAL.sub("", acc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True, help="Input float16 .npy (n, dim)")
    ap.add_argument("--ids", required=True, help="Input accession ids .txt (one per row)")
    ap.add_argument("--out", required=True, help="Output .npz")
    ap.add_argument("--compress", action="store_true",
                    help="Use np.savez_compressed (slower; float16 embeddings compress poorly)")
    args = ap.parse_args()

    emb = np.load(args.emb, mmap_mode="r")
    with open(args.ids) as fh:
        ids = [line.rstrip("\n") for line in fh]

    if len(ids) != emb.shape[0]:
        raise SystemExit(f"row mismatch: {emb.shape[0]} embeddings vs {len(ids)} ids")

    accession_ids = np.array(ids, dtype=object)
    variants = np.array([variant_of(a) for a in ids], dtype=object)

    names, first, counts = np.unique(variants, return_index=True, return_counts=True)
    print(f"{emb.shape[0]} rows, dim {emb.shape[1]}, dtype {emb.dtype}")
    for i in np.argsort(first):  # report in on-disk order
        print(f"  {names[i]:<24} {counts[i]:>7}  rows {first[i]}..{first[i] + counts[i] - 1}")

    saver = np.savez_compressed if args.compress else np.savez
    print(f"Writing {args.out} ({'compressed' if args.compress else 'uncompressed'}) ...")
    saver(
        args.out,
        accession_ids=accession_ids,
        embeddings=np.asarray(emb),  # materialize the memmap for the writer
        variants=variants,
    )
    print("done")


if __name__ == "__main__":
    main()
