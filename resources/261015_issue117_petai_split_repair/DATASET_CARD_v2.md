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
  - config_name: normalised
    data_files:
      - split: train
        path: normalised/train-*
---

# PETai Activity Labels

13,054 measured activity and stability labels on 1,137 plastic-degrading enzyme
sequences, curated from published references plus one in-house halo screen.
Built for iGEM Toronto 2026 (project PetaBite) as the supervised label set for
the [PETadex](https://huggingface.co/datasets/petadex/petadex-corpus-full)
corpus of 307,155,746 catalytic ORFs.

PETadex supplies candidates. This supplies the targets. No activity label exists
anywhere in the PETadex release, so the two are published separately and joined
on sequence MD5.

**This is v2.** v1 shipped a split that could not evaluate, and the fix turned
out to change the modelling target rather than just the folds. Read the next two
sections before using either version.

## What changed in v2

**The 73 contained fragments are remapped onto their parents.** v1 stranded 82
sequences with 3,188 labels, 24.4% of the release, because they carried labels
but no fold assignment. 73 of them name a `contained_in` parent already in a
fold, so 3,150 of those labels are now attached to the parent. Lineage is kept:
every moved label carries `orig_seq_md5` and `remapped`. The 9 sequences with no
parent, holding 38 labels, are dropped. All 200 negatives survive.

**The split unit is now a homology component, not a 30% cluster.** This is the
correction that matters, and it applies to v1 as well.

The release ships three clusterings, at 90%, 50% and 30% identity, and they are
not nested. 14 of the 318 50%-identity clusters straddle two different 30%
clusters, so they were produced by independent runs rather than a cascade.
Splitting on `cluster_30` alone therefore leaves sequences that share a
50%-identity cluster, and so are more similar than the 30% criterion implies,
sitting in different folds. **The v1 card claimed a correct 30%-identity
partition with a passing leakage check. That was true for `cluster_30` and
incomplete: four of v1's 50%-identity clusters span folds.**

v2 assigns connected components of the graph linking any two sequences that
share a cluster at any threshold. The result is clean at 30%, 50% and 90%
identity simultaneously: zero groups at any threshold span folds.

**The evaluated target is regression, not classification.** See below. This was
forced by the data, not chosen.

## Why binary activity cannot be evaluated here

v1's headline defect was that test held 7 negative labels on 5 of 170 sequences.
It cannot be fixed by re-stratifying, and measuring why is the main result of
v2.

All 200 negatives sit in 13 of the 196 30%-identity clusters, and 160 of them in
the four largest clusters in the set. Once the split unit is corrected to a
homology component, 546 sequences, 48% of the release, collapse into a single
component that holds:

- 181 of the 200 negatives
- 187 of the 197 halo binary labels
- 64% of all labels

Only 19 negatives and 10 halo labels sit anywhere else. At 48% of sequences that
component cannot be an evaluation fold, so it has to be train, and binary halo
activity therefore has no homology-clean test set at any fold size.

The component is not breakable cheaply. One 50%-identity cluster alone spans
five of the large 30% clusters and contributes 196 sequences, and the six
largest clusters are linked by a dense mesh of independent groups rather than a
chain with cut points. Forcing the clusterings to nest costs 211 sequences and
32.6% of the labels.

This is biology, not a curation defect. These sequences are all PET hydrolases
and cutinases, one tight homologous family, and the halo screen sampled exactly
that family. The label set simply does not contain homology-independent inactive
sequences in quantity.

**Consequence for users.** Train a binary activity head if useful, but any
number it produces is training-set performance. It cannot be reported as held
out. A stratified `cluster_30` split does reach 54 negative labels on 53 test
sequences, but it costs 13,494 cross-fold pairs at 50% identity or better across
430 sequences and 108 pairs at 90% or better across 37, so metrics on it are
inflated by near-duplicates.

## The split

Folds are homology components. The evaluated set is PET monomer release in
concentration units, restricted to human-curated labels, which also removes
v1's provenance inversion from every metric.

| Fold | Components | Sequences | Labels | Evaluated labels | Evaluated sequences | Halo |
|---|---|---|---|---|---|---|
| train | 168 | 956 | 11,700 | 1,564 | 123 | 197 |
| val | 4 | 130 | 629 | 90 | 27 | 0 |
| test | 12 | 51 | 725 | 89 | 16 | 0 |

Groups spanning folds, by threshold: 0 at 90%, 0 at 50%, 0 at 30%.

**The evaluated test fold is 89 labels on 16 sequences.** That is small, and it
is the honest ceiling for a homology-controlled evaluation of this label set.
The 18 components that carry evaluable labels pool 179 labels on 43 sequences,
so k-fold cross-validation over components gives the same homology guarantee
with roughly twice the evaluation data, and a bootstrap confidence interval over
components rather than sequences. That is the recommended protocol.

Resample components, never sequences. Train holds 168 components over 956
sequences and one of them is 546 sequences, so effective training diversity is
far below the sequence count.

## Named targets

`value` is not a single regression target. The release spells units 621
different ways. 15 of those are pure encoding duplicates, since the micro sign
U+00B5 and Greek mu U+03BC both occur, and NFKC folding plus dash unification
takes 621 to 590. **Do not lowercase unit strings**: `nm` and `nM`, `mm` and
`mM`, and `µm` and `µM` all occur in the data and mean different things.

The `normalised` config assigns each label to one of eight named targets with a
single canonical unit, or leaves it unassigned. 5,807 of 13,054 labels (44.5%)
are assigned.

| Target | Labels | Sequences | Canonical unit | train / val / test |
|---|---|---|---|---|
| `pet_release_conc` | 1,997 | 197 | µM | 1,812 / 93 / 92 |
| `relative_activity` | 1,208 | 82 | percent | 1,172 / 22 / 14 |
| `kcat` | 713 | 266 | s-1 | 542 / 100 / 71 |
| `kcat_over_km` | 535 | 230 | M-1 s-1 | 391 / 85 / 59 |
| `km` | 529 | 226 | µM | 373 / 96 / 60 |
| `weight_loss` | 490 | 91 | percent | 481 / 6 / 3 |
| `halo_binary` | 197 | 192 | z vs empty vector | 197 / 0 / 0 |
| `melting_temperature` | 138 | 68 | °C | 127 / 0 / 11 |

Deliberate exclusions:

- **`fraction of folded protein`** (751 labels) is CD spectroscopy curve points,
  not a melting temperature, and is not folded into `melting_temperature`.
- **`mmol` and `µmol`** are amounts rather than concentrations and need an assay
  volume, so they stay out of `pet_release_conc`.
- **`weight_loss` is separate from `relative_activity`** although both are
  percentages, because mass loss and activity are different quantities. v1's
  card treated the percent spellings as one quantity.
- **PET substrate variants** (`gfPET`, `pcPET`, `hcPET`, `3PET`, `PET-02mer` and
  others) are PET in some form but heterogeneous in morphology and oligomer
  length. Including them adds 129 human-curated labels, all of which fall in
  train, so they cannot improve evaluation. `pet_release_conc` requires the
  exact `PET` substrate token.

Rows carry a `flag` rather than being silently coerced:

| Flag | Count | Meaning |
|---|---|---|
| `scaled_unit` | 60 | The unit carries its own scale factor, such as `s-1 M-1 x 10-6` |
| `ambiguous_fraction` | 51 | A percent target whose value is in (0, 1], so it may be a fraction |
| `mw_required` | 27 | A `Km` in g/L or mg/mL, needing a molecular weight |
| `out_of_range` | 13 | A melting temperature outside 20 to 120 C |

**Correction to the v1 card.** It stated that the 151 temperature-unit rows are
actual Tm values. 13 of those 151 are physically impossible, including a 401.1
in the `°C` column and an 8.8 in the `K` column, so 138 are usable. The v1
card's fragment section also reported 1,301 recoverable PET measurements; the
correct figure is 1,281 for exact-token PET and 1,314 treating `substrate` as a
semicolon-separated list.

## Configs

| Config | Rows | Description |
|---|---|---|
| `labels` | 13,054 | One measurement per row, keyed by `seq_md5` |
| `sequences` | 1,137 | One row per sequence, with cluster, component and fold |
| `split` | 1,137 | Fold assignment with component and `cluster_30` |
| `clusters` | 1,137 | Cluster membership at 90%, 50% and 30% identity |
| `normalised` | 13,054 | `labels` plus `target`, `value_canonical`, `unit_canonical`, `flag` |

`raw/sequences.fasta` holds all 1,137 sequences, named by MD5.
`raw/summary.json` and `raw/targets_summary.json` are the build summaries.

`labels` keeps the v1 schema and adds two columns:

| Column | Type | Notes |
|---|---|---|
| `orig_seq_md5` | string | The sequence the label was recorded on, before the remap |
| `remapped` | bool | True for the 3,150 labels moved onto a parent |

## Conditions are trainable but unavailable at inference

pH is on 53.25% of labels and temperature on 29.74%, with both on 25.95%. This
is far better coverage than the PETadex corpus, which has `melting_point` and
`thermostability` null on all 307M rows.

The asymmetry is the catch. A PETadex candidate mined from Logan has no measured
pH or temperature, so a model conditioned on them at training time has nothing
to condition on at inference. Condition dropout during training is the standard
fix.

## Overlap with the PETadex corpus is small

Joining on exact sequence MD5 against all 131 `orfs` shards of
`petadex/petadex-corpus-full`, 401 of the 1,219 v1 sequences (32.90%) appear in
the corpus at all. Of those, 327 come from NCBI nr, 109 from PAZy and 21 from
Logan.

Two consequences. Corpus-derived columns cannot be model features, since they
exist for a minority of train sequences. And the 21 Logan hits are the only
direct evidence of how a model trained here behaves on the population it is
meant to screen.

This is exact-MD5 overlap and the label set is full of engineered point
variants, so it understates true homology coverage. An MMseqs2 search against
the 22,235 30%-identity centroids is the number that actually bears on transfer,
and is not yet run.

## Provenance and limitations

- Human-curated labels and LLM-extracted labels are mixed in `labels`. The LLM
  subset is unvalidated against a human-labelled gold set and contains zero
  negatives, which is evidence of extraction bias rather than of biology. v1's
  folds inverted this, with train at 19.6% LLM and test at 62.9%. v2 defines the
  evaluated set as human-curated only, so no LLM label enters a metric. Filter
  on `source != "llm"` for any evaluation of your own.
- The halo screen is in-house iGEM Toronto wet lab data under the condition
  `BHET25 @ 120h, active at z>3.0`. It is a BHET assay, not a PET polymer assay,
  so its negatives mean "no BHET hydrolysis detected" rather than "not a
  PETase". It is train-only in v2 for the homology reason above.
- Values are transcribed from assays run under widely differing protocols.
  Cross-reference comparisons are not controlled.
- `value` has no uncertainty estimate.
- The 30%, 50% and 90% clusterings disagree with each other as described above.
  Use the `component` column, not `cluster_30`, if you re-split.

## Licensing

**Not yet settled.** The underlying measurements are facts extracted from
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

Built 2026-10-14. v1 published 2026-10-15. v2 prepared 2026-10-15. See
`igem-toronto` issues 116 and 117.
