# PlasticDB + PAZy Metadata Analysis

Issue #70 — Metadata Analysis of different species under PAZy and PlasticDB.

Run `main.py` to reproduce all analysis outputs. Results are written to `/outputs`.

## Data

- `data/plasticdb_microorganisms.tsv` — PlasticDB full download (2,535 entries, 933 species, 70 plastic types, 1974–2025)
- `data/pazy_proteins.csv` — PAZy biochemically characterised enzymes (24 enzymes, 21 organisms, 9 plastic types)

## Scripts

- `scripts/data_loader.py` — load and clean both datasets
- `scripts/analysis.py` — taxonomic diversity, evidence scoring, gap analysis, cross-DB comparison
- `scripts/visualization.py` — Plotly figure builders
- `scripts/novel_discovery.py` — novelty potential scoring, phylogenetic gap detection
