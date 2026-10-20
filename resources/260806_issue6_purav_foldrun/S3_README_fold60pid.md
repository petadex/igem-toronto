# PETadex v1.1: structures of the 60% identity cluster centroids

Structures for the centroids of the 1,814,296 60%-identity clusters in PETadex v1.1: ESMFold2 all-atom
predictions, plus the AlphaFold DB model for centroids that already have one. 1,813,242 of the 1,814,296
have a structure; only sequences over 2,500 aa do not. iGEM Toronto 2026 / PetaScale, GitHub issue #6.
Folded 2026-10-15 23:01 to 2026-10-16 15:37 on the Compute Canada Killarney cluster (L40S and H100 GPUs, account aip-qmorris).

Code, benchmark and full record: `petadex/igem-toronto`, `resources/260806_issue6_purav_foldrun/` and
`notebooks/260806_issue6_purav_foldrun.ipynb` (Part II, Steps 7-16).

## Contents

| tier | model | sequence length | structures | hardware |
|---|---|---|---|---|
| `short_full/` | ESMFold2 (Full) | 40-149 aa | 123,940 | L40S |
| `mid_fast/` | ESMFold2-Fast | 150-1,740 aa | 1,667,482 | L40S |
| `tail_h100/` | ESMFold2-Fast | 1,741-2,500 aa | 2,386 | H100 |
| duplicates (in their twin's tier) | as twin | | 31 | |
| `afdb/` | AlphaFold DB v4 (copied, not folded) | any | 19,403 | |
| **total with a structure** | | | **1,813,242** | |
| no structure: longer than 2,500 aa | | 2,501-59,275 aa | 1,054 | |

`structure_index.tsv` lists every one of the 1,814,296 centroids: `orf_id`, `source` (`esmfold2_full`,
`esmfold2_fast`, `afdb` or `not_folded_over_2500aa`), `structure_key` and `metrics_key` (relative to
this folder) and `copied_from`.

- **Duplicates:** 31 centroids have exactly the same sequence as a lower-numbered centroid. They were
  folded once; the twin's structure and metrics were copied under the duplicate's ID (`copied_from` and
  `inputs/manifest.tsv` `folded_as` give the twin; the file contents still name the twin).
- **AFDB:** 19,430 centroids have an AlphaFold DB structure listed in the corpus (`centroid_afdb_s3_key`).
  For 19,403 the AFDB sequence is exactly the centroid's, and the AFDB model (`s3://petadex-protein-structures/af_db_v260701/structures/AF-<acc>-F1-model_v4.cif`)
  is copied to `afdb/structures/<NNN>/orf<ID>.cif.gz` with the accession in `copied_from` (every one
  checked against the centroid md5). For the other 27 the AFDB entry is a related protein, not the same
  sequence (1-16 substitutions, or a different length; `qc/afdb_inexact_refolded.tsv`), so those were
  folded with ESMFold2-Fast like the rest of `mid_fast/` (2026-10-17). One duplicate's twin is an AFDB
  centroid, so it has the twin's AFDB model.
  pLDDT is in the B-factor column as for ESMFold2. AlphaFold DB is CC-BY 4.0 (Varadi et al. 2024).
  There are no metrics JSONs for this tier.
- **Over 2,500 aa:** 1,054 centroids (`inputs/over_cap_2500.tsv`) are not folded; this is a hard cap.

## Inputs (provenance)

- Centroid IDs: Hugging Face `petadex/petadex-corpus-full`, `dims/clusters60.parquet`, column `centroid_orf_id`.
- Sequences: Hugging Face `petadex/petadex-corpus-centroids90`, `centroids90.parquet`, column
  `centroid_sequence` (every 60% centroid is also a 90% centroid), uppercased, trailing `*` removed.
  Verified: all 1,793,781 sequences of the main run (and the 27 folded on 2026-10-17) equal `centroid_sequence` and match `centroid_md5`;
  identical to `orfs.sequence` wherever checked.
- Corpora built 2026-10-15 by Dennis Zhu. Built by `make_fold_input_hf.py --level 60`.

## Model and settings

`esm` 3.4.1.post1 native `EsmFold2Model` (`biohub/ESMFold2` for `short_full`, `biohub/ESMFold2-Fast`
otherwise), ESMC-6B backbone in bf16, fused kernels, 32 diffusion sampling steps, 3 refinement loops,
1 diffusion sample, seed 0, single sequence (no MSA). Chunked pair operations above 1,000 aa
(chunk 64 on L40S, chunk 32 on H100). Chosen by a 17-configuration benchmark against PDB ground truth,
AlphaFold DB and a Full 100-step / 20-loop reference (notebook Steps 9-15): Fast loses confidence on
sequences under 150 aa, so those use Full.

## Layout and file formats

```
structure_index.tsv                      every centroid -> source and file (see Contents)
<tier>/structures/<NNN>/orf<ID>.cif.gz   all-atom mmCIF (gzip); per-residue pLDDT (0-100) in the B-factor column
<tier>/metrics/<NNN>/orf<ID>.json        {id, seq_len, num_sampling_steps, num_loops, backend, esmc_precision,
                                          kernel, chunk, fast, peak_vram_gb, wall_s, residues_per_s,
                                          mean_plddt, ptm, iptm, status}
<tier>/progress.shard<K>of<N>.csv        one row per fold: id, seq_len, wall_s, residues_per_s, peak_vram_gb,
                                          mean_plddt, ptm, iptm, status (same metrics, one file per shard)
inputs/                                  manifest.tsv (all 1,814,296 centroids: orf_id, c90/c60/c30 ids, origin,
                                          length, tier, skip_reason, md5, afdb_s3_key, folded_as), tier FASTAs,
                                          over_cap_2500.tsv
qc/                                      final completeness checks, QC samples, run record
```

`<NNN>` is the ORF number modulo 1000, zero-padded: `orf1041339001` is in `structures/001/`. To find
an ORF's file, look it up in `structure_index.tsv`.

## Quality notes

- Every input sequence has a structure and metrics file (checked per ORF); 0 sequence-level fold failures.
- QC on 3,000 structures: sequences match the input, heavy atoms complete, no chain gaps > 6 A. Some
  tail structures (1,741-2,500 aa) have stretched CA-CA bonds (4.4-5.3 A) in low-confidence loops.
- Median mean-pLDDT: `short_full` 67 (about half under 70; mostly short Logan stop-to-stop fragments),
  `mid_fast` 92, `tail_h100` 83. Filter or weight by pLDDT for downstream use.
- Reproducible: 80 structures overlapping the benchmark reproduce it at median TM-score 1.000.

## Cost

851.6 GPU-hours (839.2 L40S, 12.3 H100), 16 h 36 min wall-clock at up to 80 GPUs; 93.7 GB.
