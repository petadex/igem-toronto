# Embed and Plot the UMAP Pilot Controls

Downstream half of the ESM landscape control experiment (issue #75). Takes the six
control CSVs generated in issue #34, embeds them with the issue #32 ESM-2 pipeline,
bundles the embedding shard into a labelled `.npz`, and co-UMAPs all 517,840 controls
with the 64,730 real family centroids in a single coordinate system.

---

## Files

| File | Purpose |
|---|---|
| `scripts/npy_to_npz.py` | Bundle the embedding shard (`controls.npy` + `controls_ids.txt`) into a single `.npz`, recovering a parallel `variants` label array from the id prefixes |
| `scripts/run_umap.py` | Fit a 2D UMAP over one or more `.npz` embedding files; write coordinates CSV + plot |
| `plots/umap.png` | UMAP projection of 64,730 real centroids + 517,840 controls |

---

## Pipeline

```bash
# 1. (issue #32 pipeline) embed the control CSVs from #34
#    esm2_t30_150M_UR50D, mean-pooled, 640-dim, float16
#    -> controls.npy (517840, 640)  +  controls_ids.txt (parallel ids, same order)

# 2. Bundle the shard into a labelled .npz.
#    Variants are recovered from id prefixes, e.g.
#      rand_95th_family_0__0        -> 'rand_95th_family'
#      rand_fragment_21633__194188  -> 'rand_fragment'
#      rand_uniprot_Q5JDK7__388380  -> 'rand_uniprot'
#      shuffled_64918__517839       -> 'shuffled'
python3 scripts/npy_to_npz.py --emb controls.npy --ids controls_ids.txt --out controls.npz

# 3. Co-UMAP controls with the real family centroids (one shared fit).
#    family_embeddings.npz has no `variants` column, so it takes the --labels value 'real';
#    controls.npz is labelled by its own `variants` array.
python3 scripts/run_umap.py family_embeddings.npz controls.npz \
    --labels real --output umap_output.csv --plot plots/umap.png
```

UMAP parameters (in `run_umap.py`): `n_neighbors=15`, `min_dist=0.1`, `metric=euclidean`,
`random_state=42`.

---

## Inputs / outputs (S3)

Embeddings are large and live in S3, not the repo:

- `controls.npy` — 632 MB (517840 × 640 float16) — raw shard
- `controls.npz` — 587 MB — bundled + labelled
- `family_embeddings.npz` — 147 MB (64730 × 640 float32) — real centroids
- `umap_output.csv` — 33 MB — final coordinates (`id, label, x, y`)

```bash
aws s3 sync ./ s3://petadex/esm_embeddings/esm_controls/ \
    --exclude "*" --include "*.npz" --include "*.csv"
```

---

## Procedural pitfalls

These bit (or nearly bit) us during the run — read before re-running:

1. **Id ↔ embedding row order must be preserved.** `npy_to_npz.py` pairs row *i* of
   `controls.npy` with line *i* of `controls_ids.txt` by position — there is no key
   join. If the embedding job completes batches out of order, or silently drops
   sequences (e.g. over a max length), the two desync and every point is mislabelled.
   The script hard-fails on a count mismatch (`row mismatch: N vs M`), but an
   equal-count *reorder* would pass silently. Keep embedding output in input order.

2. **Variant labels are recovered by string surgery, not stored.** `variant_of()`
   strips the trailing `__<globalidx>` then the trailing `_<localid>` token via regex.
   This assumes the id is `<variant>_<localid>__<globalidx>` and that neither the
   variant name nor the local id contains an unexpected underscore. A UniProt accession
   with an underscore, or a `family_id` with one, would be mis-split into the wrong
   variant. Always eyeball the printed per-variant counts (they should be 64,730 each,
   194,190 for fragments) before trusting the bundle.

3. **dtype mismatch between the two `.npz` files.** `controls.npz` is float16 (saved
   from the shard) while `family_embeddings.npz` is float32. `run_umap.py` casts both
   to float32 before stacking, so this is handled — but if you stack them yourself,
   `np.vstack` on mixed dtypes will upcast/copy and can blow memory unexpectedly.

4. **`np.savez` materialises the memmap.** `npy_to_npz.py` opens the `.npy` with
   `mmap_mode='r'` but `np.savez` calls `np.asarray(emb)`, pulling the full 632 MB into
   RAM to write it. The `--compress` flag is available but near-useless here: float16
   embeddings are high-entropy and barely compress, while compression roughly doubles
   write time. Default (uncompressed) is the right call.

5. **The UMAP is a *single joint fit*, and controls outnumber real ~8:1.** All 582,570
   points are fit together so they share one coordinate system — but the 517,840
   controls dominate the manifold, so the layout of the real centroids is shaped largely
   by the controls around them. This is intended (we want one comparable space), but it
   means you cannot read absolute distances as if the real points were embedded alone.

6. **`random_state=42` forces single-threaded UMAP.** Setting a seed disables UMAP's
   parallelism, so fitting 582k × 640 is slow and memory-heavy. Worth it for
   reproducibility, but budget the time; do not assume it hung.

7. **Plot overdraw hides minority classes.** `run_umap.py` scatters groups in
   `groupby` (alphabetical) order with `alpha=0.4`. Whatever is drawn last sits on top,
   so dense sets can occlude sparser ones. If a control looks "missing" from the plot,
   check the coordinates CSV before concluding it overlaps another set.

8. **Random-UniProt under-delivery is silent.** `generate_rand_uniprot.py` skips failed
   REST batches with a warning and filters entries lacking a sequence, so the fetched
   count can fall short of `n` (obsolete/deleted accessions return nothing). Verify the
   `rand_uniprot` row count after generation, not after embedding.

---

## Prerequisites

- `numpy`, `pandas`, `umap-learn`, `matplotlib`
- The issue #32 ESM-2 embedding pipeline (step 1)
- The control CSVs from issue #34
