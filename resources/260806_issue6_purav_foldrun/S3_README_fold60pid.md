# PETadex v1.1: predicted structures of the 60% identity cluster centroids

ESMFold2 all-atom structure predictions for the centroids of the 1,814,296 60%-identity clusters in
PETadex v1.1 (iGEM Toronto 2026 / PetaScale, GitHub issue #6). Folded 2026-10-15 23:01 to 2026-10-16
15:37 on the Compute Canada Killarney cluster (L40S and H100 GPUs, account aip-qmorris).

Code, benchmark and full record: `petadex/igem-toronto`, `resources/260806_issue6_purav_foldrun/` and
`notebooks/260806_issue6_purav_foldrun.ipynb` (Part II, Steps 7-16).

## Contents

| tier | model | sequence length | structures | hardware |
|---|---|---|---|---|
| `short_full/` | ESMFold2 (Full) | 40-149 aa | 123,940 | L40S |
| `mid_fast/` | ESMFold2-Fast | 150-1,740 aa | 1,667,455 | L40S |
| `tail_h100/` | ESMFold2-Fast | 1,741-2,500 aa | 2,386 | H100 |
| **total** | | | **1,793,781** | |

Not folded, by design (listed in `inputs/manifest.tsv` and `inputs/over_cap_2500.tsv`):
19,430 centroids that already have their own AlphaFold DB structure, 31 exact-duplicate sequences
(`folded_as` gives the ORF whose structure to use), and 1,054 sequences longer than the 2,500 aa hard cap.
1,793,781 + 19,430 + 31 + 1,054 = 1,814,296.

## Inputs (provenance)

- Centroid IDs: Hugging Face `petadex/petadex-corpus-full`, `dims/clusters60.parquet`, column `centroid_orf_id`.
- Sequences: Hugging Face `petadex/petadex-corpus-centroids90`, `centroids90.parquet`, column
  `centroid_sequence` (every 60% centroid is also a 90% centroid), uppercased, trailing `*` removed.
  Verified: all 1,793,781 folded sequences equal `centroid_sequence` and match `centroid_md5`;
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
an ORF's tier, look it up in `inputs/manifest.tsv`.

## Quality notes

- Every input sequence has a structure and metrics file (checked per ORF); 0 sequence-level fold failures.
- QC on 3,000 structures: sequences match the input, heavy atoms complete, no chain gaps > 6 A. Some
  tail structures (1,741-2,500 aa) have stretched CA-CA bonds (4.4-5.3 A) in low-confidence loops.
- Median mean-pLDDT: `short_full` 67 (about half under 70; mostly short Logan stop-to-stop fragments),
  `mid_fast` 92, `tail_h100` 83. Filter or weight by pLDDT for downstream use.
- Reproducible: 80 structures overlapping the benchmark reproduce it at median TM-score 1.000.

## Cost

851.6 GPU-hours (839.2 L40S, 12.3 H100), 16 h 36 min wall-clock at up to 80 GPUs; 93.7 GB.
