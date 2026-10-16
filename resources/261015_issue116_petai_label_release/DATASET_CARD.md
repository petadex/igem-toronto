---
pretty_name: PETai Activity Labels
language:
  - en
tags:
  - biology
  - proteins
  - enzymes
  - plastic-degradation
  - PETase
  - protein-engineering
task_categories:
  - tabular-regression
  - text-classification
size_categories:
  - 10K<n<100K
configs:
  - config_name: labels
    data_files:
      - split: train
        path: labels/train-*
  - config_name: sequences
    data_files:
      - split: train
        path: sequences/train-*
  - config_name: split
    data_files:
      - split: train
        path: split/train-*
  - config_name: clusters
    data_files:
      - split: train
        path: clusters/train-*
---

# PETai Activity Labels

13,092 measured activity and stability labels on 1,219 plastic-degrading enzyme
sequences, curated from 696 published references plus one in-house halo screen.
Built for iGEM Toronto 2026 (project PetaBite) as the supervised label set for
the [PETadex](https://huggingface.co/datasets/petadex/petadex-corpus-full)
corpus of 307,155,746 catalytic ORFs.

PETadex supplies candidates. This supplies the targets. No activity label exists
anywhere in the PETadex release, so the two are published separately and joined
on sequence MD5.

## Read this before training on the split

**The shipped `split` config is unusable for evaluation.** The test fold holds 7
negative labels on 5 of its 170 sequences. The in-house halo screen contributes
168 of the 200 negatives in the release and was allocated 178 to train, 12 to
val and 0 to test. No classification metric, calibration curve or enrichment
number computed on this test fold is interpretable.

A second defect compounds it: label provenance is inverted across the folds.
Train is 19.6% LLM-extracted and test is 62.9%, and the 4,898 LLM labels contain
zero negatives. A model trained here learns from human-curated measurements and
is scored on unvalidated machine extractions.

The split is shipped anyway because the homology-clean clustering underneath it
is sound and worth reusing. It is a correct 30%-identity partition with a
passing leakage check. Only the stratification is wrong. A corrected re-split
is planned as v2. Until then, treat `split` as a cluster assignment rather than
an evaluation protocol, and re-stratify it yourself.

## Configs

| Config | Rows | Description |
|---|---|---|
| `labels` | 13,092 | One measurement per row, keyed by `seq_md5` |
| `sequences` | 1,219 | One row per sequence, with cluster and split assignment |
| `split` | 1,137 | The 30%-identity fold assignment for the clean sequences |
| `clusters` | 1,219 | Cluster membership at 90%, 50% and 30% identity |

`raw/sequences.fasta` holds all 1,219 sequences, named by MD5.
`raw/summary.json` is the release summary emitted by the build.

### `labels` schema

| Column | Type | Notes |
|---|---|---|
| `seq_md5` | string | Join key, MD5 of the uppercase sequence |
| `source` | string | `activity` 6,655, `llm` 4,898, `kinetics` 1,342, `halo` 197 |
| `task` | string | `quantitative` 12,895, `binary` 197 |
| `variable` | string | What was measured, 99.95% present |
| `unit` | string | See the unit warning below |
| `substrate` | string | 92.16% present. PET 4,287, pNPB 1,078, BHET 822 |
| `value` | float | Not one target. See below |
| `is_negative` | bool | 200 of 13,092 (1.53%) |
| `ph` | float | 53.25% present |
| `temperature_c` | float | 29.74% present |
| `time_h` | float | 9.99% present |
| `ref` | string | 696 distinct references |
| `curator` | string | 89.75% present |

## `value` is not a single regression target

The column mixes incompatible units and spans -372.91 to 1.21e9 as a result. Do
not regress on it directly. The three collisions that matter:

- **Relative activity** is split across `%` (1,370), `percent` (917) and
  `% relative activity` (369). Same quantity, three spellings.
- **Melting data**, 904 labels in all, mixes `fraction of folded protein`
  (751), which are CD spectroscopy curve points rather than a melting
  temperature, with `°C` (148), `K` (3) and `mJ mg⁻¹` (2). Only the 151
  temperature-unit rows are actual Tm values.
- **Kinetics** are three separate targets, not one: `kcat` (559), `Km` (551)
  and `kcat/Km` (473).

Clean target choices, both well populated: binary halo activity for
classification, and PET monomer release for regression.

## Conditions are trainable but unavailable at inference

pH is on 53.25% of labels and temperature on 29.74%, with both on 3,398
(25.95%). This is far better coverage than the PETadex corpus, which has
`melting_point` and `thermostability` null on all 307M rows.

The asymmetry is the catch. A PETadex candidate mined from Logan has no measured
pH or temperature, so a model conditioned on them at training time has nothing
to condition on at inference. Condition dropout during training is the standard
fix.

## Overlap with the PETadex corpus is small

Joining on exact sequence MD5 against all 131 `orfs` shards of
`petadex/petadex-corpus-full`, 401 of 1,219 labelled sequences (32.90%) appear
in the corpus at all. Of those, 327 come from NCBI nr, 109 from PAZy and
**21 from Logan**.

Two consequences. Corpus-derived columns cannot be model features, since they
exist for only 216 of the 796 train sequences. And the 21 Logan hits are the
only direct evidence of how a model trained here behaves on the population it
is meant to screen.

This is exact-MD5 overlap and the label set is full of engineered point
variants, so it understates true homology coverage. An MMseqs2 search against
the 22,235 30%-identity centroids is the number that actually bears on
transfer, and is not yet run.

## Cluster-size asymmetry

Bootstrap confidence intervals must resample **clusters, not sequences**.

| Fold | Sequences | 30% clusters | Mean size | Max | Top 3 share |
|---|---|---|---|---|---|
| train | 796 | 70 | 11.37 | 124 | 43.84% |
| val | 171 | 63 | 2.71 | 16 | 23.39% |
| test | 170 | 63 | 2.70 | 14 | 22.94% |

Effective training diversity is nearer 70 families than 796 sequences.

## Fragments

82 of the 1,219 sequences are fragments or contained duplicates and carry no
split assignment. They hold 3,188 labels (24.35% of the release), including
1,301 PET measurements and 773 of the 904 melting labels (85.5%).

73 of those 82 name a `contained_in` parent that is already in a fold, so 3,150
of the 3,188 labels can be recovered by remapping onto the parent. Only 9
sequences with 38 labels are genuinely unrecoverable. The remap is not applied
in this version.

## Provenance and limitations

- 8,194 labels are human-curated, 4,898 are LLM-extracted from published PDFs.
  The LLM subset is unvalidated against a human-labelled gold set at release
  and contains zero negatives, which is itself evidence of extraction bias
  rather than of biology. Filter on `source != "llm"` for an evaluation set.
- The halo screen is in-house iGEM Toronto wet lab data under the condition
  `BHET25 @ 120h, active at z>3.0`. It is a BHET assay, not a PET polymer
  assay, so its negatives mean "no BHET hydrolysis detected" rather than
  "not a PETase".
- Values are transcribed from assays run under widely differing protocols.
  Cross-reference comparisons are not controlled.
- `value` has no uncertainty estimate.

## Licensing

**Not yet settled.** The underlying measurements are facts extracted from 696
published references, but the aggregation, the in-house halo screen and the
curation are iGEM Toronto 2026 work and no license has been chosen. Contact the
team before redistribution or commercial use.

## Citation

```bibtex
@misc{petai_activity_labels_2026,
  title  = {PETai Activity Labels},
  author = {{iGEM Toronto 2026}},
  year   = {2026},
  url    = {https://huggingface.co/datasets/petadex/petai-activity-labels}
}
```

Built 2026-10-14. Published 2026-10-15. See `igem-toronto` issue 116.
