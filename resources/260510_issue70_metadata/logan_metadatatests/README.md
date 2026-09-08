# Logan bioplastics metadata benchmark

## Use these notebooks

- `03_build_bioplastics_metadata_variant.ipynb`
- `04_benchmark_bioplastics_metadata.ipynb`

These are the executed, reviewed notebooks for the study definition of metadata:
plastic or bioplastic associations, enzyme and evidence annotations, biochemical
protein properties, sequence-similarity candidates, study context, and source
provenance.

## Earlier generic experiment

Notebooks `01_build_logan_database_variants.ipynb` and
`02_benchmark_metadata_vs_no_metadata.ipynb` measure generic ENA/run metadata.
They are retained as a separate administrative-metadata experiment, but they do
not answer the bioplastics-metadata question and should not be cited as the
study's primary comparison.

## Interpretation

The corrected benchmark compares identical canonical Logan sequence rows:

- `no_bioplastics_metadata`: sequence identity, nucleotide sequence, length, GC
  fraction, and source-run accession.
- `bioplastics_metadata`: the same sequence rows plus PlasticDB and PAZy
  reference annotations, PETadex-derived classification rules, ProtParam
  biochemical properties, verified plastic-study context, DIAMOND candidate
  similarities, and provenance.

DIAMOND rows are a permissive similarity screen, not predictions or proof of
plastic-degradation function. The microplastic study runs establish experimental
context only.

## Principal evidence

- `results/bioplastics_validation.json`
- `results/bioplastics_reference_manifest.json`
- `results/plastic_classification_audit.csv`
- `results/homology_control_validation.json`
- `results/homology_run_summary.json`
- `results/bioplastics_shared_query_proofs.json`
- `results/bioplastics_query_plans.json`
- `results/bioplastics_capability_query_plans.json`
- `results/bioplastics_storage_summary.csv`
- `results/bioplastics_build_measurements.csv`
- `results/bioplastics_query_summary.csv`
- `results/bioplastics_metadata_capabilities.csv`