
# Plastic Biodegradation Meta-Analysis Toolkit



---

## Data Sources

| Database | URL | Coverage | Entries |
|----------|-----|----------|---------|
| **PlasticDB** | [plasticdb.org](https://plasticdb.org) | All reported plastic-degrading microorganisms | ~875 species, ~2,535 entries |
| **PAZy** | [pazy.eu](https://www.pazy.eu) | Thoroughly biochemically characterised plastic-active enzymes | ~24 curated enzyme entries |


---

## Repository Structure

```
plastic-biodegradation-analysis/
├── README.md
├── requirements.txt
├── .gitignore
│
├── data/                          # Downloaded/cached database files
│   ├── plasticdb_microorganisms.tsv  ← PlasticDB TSV (auto-downloaded)
│   └── pazy_proteins.csv             ← PAZy curated enzymes (scraped/cached)
│
├── src/                           # Core analysis library
│   ├── __init__.py
│   ├── data_loader.py             # Load, clean, categorise both databases
│   ├── analysis.py                # Meta-analysis: diversity, gaps, novelty scoring
│   ├── visualization.py           # matplotlib + plotly charts (all figures)
│   └── novel_discovery.py         # Novel species discovery pipeline
│
├── notebooks/                     # Analysis notebooks (plain Python / Jupytext)
│   ├── 01_data_exploration.py     # Raw data structure, completeness, distributions
│   ├── 02_taxonomic_analysis.py   # Genus/species diversity, phylogenetic gaps
│   ├── 03_plastic_substrate_analysis.py  # Substrate coverage, co-occurrence, gaps
│   ├── 04_geographic_temporal_trends.py  # Geography & time trends
│   ├── 05_novel_species_discovery.py     # Full discovery pipeline + report
│   └── 06_cross_database_analysis.py     # PlasticDB vs PAZy comparison
│
├── outputs/
│   ├── figures/                   # Saved PNG / HTML charts
│   ├── reports/                   # CSV tables + JSON summaries
│   └── discovery_report.txt       # Text discovery report
│
├── generate_sample_results.py     # Run all analyses, save figures + reports
└── app.py                         # Streamlit interactive explorer
```

---

## Quick Start

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Generate all sample results & figures
```bash
python generate_sample_results.py
```
This downloads the latest data, runs all analyses, saves ~10 publication-quality figures to `outputs/figures/`, and prints the full discovery report.

### 3. Launch the interactive Streamlit app
```bash
streamlit run app.py --server.port 5000
```
Navigate to `http://localhost:5000` to explore all analyses interactively.

### 4. Run individual notebooks
Notebooks are plain Python files (Jupytext light format). To open as Jupyter notebooks:
```bash
pip install jupytext
jupytext --to notebook notebooks/01_data_exploration.py
jupyter notebook notebooks/01_data_exploration.ipynb
```
Or run directly:
```bash
python notebooks/01_data_exploration.py
```

---

## Analysis Modules

### `src/data_loader.py`
- `load_plasticdb()` — downloads and cleans the PlasticDB TSV with 22 fields
- `fetch_pazy_proteins()` — scrapes/loads PAZy characterised enzyme table
- `get_unique_organisms()` — one-row-per-organism aggregate with plastic breadth
- `get_plastic_summary()` — per-plastic statistics including evidence coverage
- `load_all()` — convenience loader returning all four DataFrames

### `src/analysis.py`
- `evidence_quality_score()` — 0–100 composite score (sequence, genbank, enzyme, analytical grade, not-extrapolated)
- `taxonomic_diversity()` — Shannon diversity, singleton analysis, top genera
- `temporal_trend_analysis()` — year-by-year growth, rolling averages, cumulative species
- `geographic_distribution()` — activity by country/region
- `isolation_environment_profile()` — soil, marine, compost, etc.
- `research_gap_analysis()` — gap priority score per plastic type
- `plastic_co_occurrence()` — organisms sharing multiple plastic substrates
- `cross_database_comparison()` — PlasticDB vs PAZy overlap analysis
- `compute_novelty_potential()` — multi-factor novelty scoring for every organism

### `src/visualization.py`
All chart functions accept `backend='plotly'` (interactive) or `backend='matplotlib'` (publication PNG).
- Plastic type distribution bar chart
- Temporal trends (dual-axis: entries/yr + cumulative species)
- Geographic activity horizontal bar chart
- Plastic co-occurrence heatmap
- Novelty potential bubble scatter
- Evidence quality donut chart
- Top 20 genera bar chart
- Plastic category sunburst
- Research gap priority chart
- PlasticDB vs PAZy coverage comparison

### `src/novel_discovery.py`
- `identify_phylogenetic_gaps()` — singleton genera with high discovery priority
- `underexplored_environments()` — isolation environments with high species count but low characterisation
- `plastic_specific_candidates()` — per-plastic candidate ranking for follow-up experiments
- `generate_discovery_report()` — comprehensive discovery summary dict
- `format_discovery_report_text()` — human-readable text report

---

## Key Findings (as of May 2026)

### Research Explosion
- **2022–2023** were peak years for plastic biodegradation publications
- Cumulative species count has grown ~6× since 2015
- PET research dominates the well-characterised literature (IsPETase, LCC, FAST-PETase)

### Taxonomic Landscape
- **875 unique species** across **200+ genera** reported in PlasticDB
- >50% of genera are **singletons** — one described species — indicating massive under-sampling
- *Pseudomonas*, *Bacillus*, *Aspergillus* dominate the count, but many rare clades show multi-plastic activity

### Substrate Coverage Gaps
- **PE, LDPE, HDPE, PP, PS, PVC** — the six most environmentally abundant plastics — have the **highest research gap scores**: few organisms, low sequence evidence, high extrapolation rates
- Biodegradable plastics (PHA, PHB, PCL, PLA) are the best characterised; commodity plastics remain the primary challenge

### Geographic Bias
- **India, Japan, South Korea, France** account for a disproportionate share of reports
- Large geographic regions (Africa, Middle East, Latin America) are barely represented — major biodiversity sampling opportunity

### Cross-Database Gap
- Only ~8 plastic types have entries in both PlasticDB *and* PAZy (the biochemically-verified set)
- Most PlasticDB organisms lack protein sequences in the database — priority targets for genomic/proteomic follow-up

---

## Novel Discovery Methodology

The novelty potential score combines four components:
| Component | Weight | Rationale |
|-----------|--------|-----------|
| Plastic breadth | 30% | Organisms degrading more plastic types are more versatile |
| Taxonomic rarity | 25% | Rare genera are phylogenetically undersampled |
| Recency | 20% | Recently reported organisms may represent new research fronts |
| Evidence gap | 25% | High breadth but low molecular evidence = highest experimental value |

The discovery pipeline additionally identifies:
- **Phylogenetic gaps** — genera with only one described species in the DB
- **Underexplored environments** — isolation sources with low characterisation rates
- **Hard plastic candidates** — organisms reported on PE/PP/PS/PVC worth follow-up sequencing

---

## Streamlit Application

The `app.py` Streamlit application provides 8 interactive pages:

| Page | Content |
|------|---------|
| Overview | Key metrics, top plastics, evidence quality, publication trends |
| Taxonomy | Genus/species rankings, multi-plastic generalists |
| Plastic Substrates | Sunburst, co-occurrence heatmap, per-plastic explorer |
| Geography & Time | Temporal growth, geographic bar chart, environment breakdown |
| Research Gaps | Gap priority scores, unstudied regions, evidence tables |
| Novel Discovery | Novelty scatter, phylogenetic gaps, environment gaps, hard-plastic candidates |
| Cross-Database | PlasticDB vs PAZy comparison, characterised enzyme table |
| Data Explorer | Filterable full dataset with CSV download |

---

## Extending the Toolkit

### Adding a new data source
Implement a loader function in `src/data_loader.py` following the same pattern as `load_plasticdb()`, then integrate the result into `load_all()`.

### Adding a new analysis
Add a function to `src/analysis.py` that accepts a DataFrame and returns a DataFrame or dict, then wire it into `run_full_analysis()`.

### Adding a new visualisation
Add a function to `src/visualization.py` following the `plot_*(df, backend='plotly')` signature, then add it to `save_all_figures()`.

---

