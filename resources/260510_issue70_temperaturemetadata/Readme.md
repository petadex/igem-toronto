# Temperature Analysis — PETase Enzyme Discovery
### iGEM Toronto 2026 | PetaBite Dry Lab

A dry-lab analysis of temperature-related metadata for plastic-degrading enzymes, built on real data from PlasticDB (2,535 entries) and published primary literature. This module supports the PetaBite goal of finding room-temperature PETases that outperform FAST-PETase without requiring the ~70°C Glass Transition Point (GTP) heating infrastructure.

---

## Background

The current gold-standard PETases — LCC, ICCG-LCC, HotPETase — operate at or near the Glass Transition Point of PET (~70°C), where the polymer chain becomes flexible enough for efficient hydrolysis. This makes them unsuitable for Wastewater Treatment Plant (WWTP) deployment at scale. The PetaBite hypothesis is that the PETadex (300 million sequences from the Logan assembly of the NCBI SRA) contains room-temperature enzymes that already outperform engineered variants.

This analysis characterises the current landscape: how many known organisms degrade PET at mesophilic conditions, what protein properties correlate with lower temperature optima, and which organisms from PlasticDB are the best candidates for room-temperature PETase extraction.

---

## Data Sources

| Source | Description | Location |
|---|---|---|
| PlasticDB TSV | 2,535 organism-plastic-paper entries downloaded from plasticdb.org | `../plastic-biodegradation-analysis/data/plasticdb_microorganisms.tsv` |
| PAZy proteins | Biochemically characterised plastic-active enzymes scraped from pazy.eu | `../plastic-biodegradation-analysis/data/pazy_proteins.csv` |
| Benchmark temperature values | Topt, Tm, kcat from primary literature (cited per enzyme in NB2) | Hardcoded in `notebooks/02_benchmark_petase_temperatures.ipynb` |

No values are fabricated. All counts and statistics are computed directly from the downloaded database files or cited primary literature.

---

## Results Summary

### Thermophilic vs Mesophilic in PlasticDB (Notebook 1)

Of all 2,535 PlasticDB entries:

- **2,288 (90.3%)** were characterised under mesophilic (non-thermophilic) conditions
- **193 (7.6%)** were explicitly recorded as thermophilic
- **54 (2.1%)** had no temperature condition recorded

PET specifically has a lower thermophilic rate than PHB and PHA, reflecting that most early PET biodegradation work used ambient-temperature incubation. The thermophilic entry rate has grown slightly since 2018, corresponding to the surge in LCC-based engineering papers following Tournier et al. 2020.

### Benchmark PETase Temperature Profiles (Notebook 2)

Twelve characterised PETase variants span a Topt range of 30°C (IsPETase) to 72°C (ICCG-LCC). Key findings:

- Topt and Tm are strongly correlated (Pearson r > 0.95), meaning thermostability engineering has reliably shifted both together
- ICCG-LCC achieves the highest kcat at its Topt (1.622 s⁻¹) but loses >99.9% of that activity when cooled to 37°C (kcat = 0.001 s⁻¹)
- FAST-PETase has the best kcat at 37°C among engineered variants (0.058 s⁻¹) — this is the benchmark to beat
- DuraPETase (Topt 37°C) and PET2 (Topt 40°C) are the only natural/engineered variants tuned toward room temperature, but both underperform FAST-PETase at 37°C
- IsPETase works at 30°C (kcat 0.022 s⁻¹) — functional at room temperature but slower than FAST-PETase at ambient

### Protein Thermostability Prediction (Notebook 3)

Biopython ProtParam was run on all PlasticDB sequences that carry a thermophilic label:

- Mesophilic sequences show a higher proportion with instability index below 40 (predicted stable), contrary to the naive expectation that thermophilic organisms would produce more stable proteins — this likely reflects that thermophilic entries in PlasticDB are predominantly from extremophilic whole-organism studies rather than purified enzyme characterisation
- GRAVY scores for characterised PETases trend negative (hydrophilic), consistent with their extracellular, aqueous-interface activity
- Lower Topt PETases (IsPETase, FAST-PETase) have more negative GRAVY scores than high-Topt variants (LCC, ICCG-LCC), suggesting hydrophilicity as one predictor for room-temperature activity

### Room-Temperature PET Candidates (Notebook 4)

Filtering PlasticDB to mesophilic PET entries:

| Stage | Count |
|---|---|
| All PET entries | 501 |
| Mesophilic PET (thermophilic = No) | 445 |
| Mesophilic PET with linked sequence | 87 |
| Mesophilic PET with named enzyme | 61 |

The 87 sequence-linked mesophilic PET entries represent the directly actionable candidate pool from PlasticDB for room-temperature PETase discovery. Top-ranked organisms by evidence score are output to `outputs/reports/04_room_temp_pet_candidates.csv`. Isolation environments of the candidate pool are dominated by soil, compost, and wastewater — all mesophilic settings compatible with room-temperature enzyme activity.

---

## File Structure

```
temperature-analysis/
├── notebooks/
│   ├── 01_thermophile_distribution.ipynb
│   ├── 02_benchmark_petase_temperatures.ipynb
│   ├── 03_thermostability_prediction.ipynb
│   └── 04_room_temp_candidates.ipynb
├── outputs/
│   ├── figures/
│   │   ├── 01_thermophile_distribution.png
│   │   ├── 02_benchmark_temperatures.png
│   │   ├── 02_activity_penalty_at_37c.png
│   │   ├── 03_thermostability_prediction.png
│   │   ├── 03_gravy_vs_topt.png
│   │   ├── 04_room_temp_candidates.png
│   │   └── 04_room_temp_environments.png
│   └── reports/
│       ├── 01_thermophile_overall.csv
│       ├── 01_thermophile_by_plastic.csv
│       ├── 02_benchmark_petase_temperatures.csv
│       └── 04_room_temp_pet_candidates.csv
├── build.py
└── README.md
```

---

## How to Run

### Prerequisites

Install dependencies from the parent project:

```bash
pip install -r ../plastic-biodegradation-analysis/requirements.txt
```

Biopython is required for Notebook 3:

```bash
pip install biopython
```

### Option 1 — Regenerate everything at once

Run the build script from the `temperature-analysis/` directory. This recreates all figures, reports, and `.ipynb` files from scratch:

```bash
cd temperature-analysis
python build.py
```

### Option 2 — Run notebooks individually in Jupyter

Launch Jupyter from the `temperature-analysis/notebooks/` directory so that relative paths resolve correctly:

```bash
cd temperature-analysis/notebooks
jupyter notebook
```

Then open and run any notebook. Each notebook imports the parent project's `src/` library via a `sys.path` insert at the top of the first code cell — no installation step is needed.

### Option 3 — Execute notebooks non-interactively

```bash
cd temperature-analysis/notebooks
jupyter nbconvert --to notebook --execute 01_thermophile_distribution.ipynb
jupyter nbconvert --to notebook --execute 02_benchmark_petase_temperatures.ipynb
jupyter nbconvert --to notebook --execute 03_thermostability_prediction.ipynb
jupyter nbconvert --to notebook --execute 04_room_temp_candidates.ipynb
```

---

## Key Figures

| Figure | Description |
|---|---|
| `01_thermophile_distribution.png` | Four-panel breakdown of thermophilic vs mesophilic entries across all 2,535 PlasticDB records, by plastic type, isolation environment, and year |
| `02_benchmark_temperatures.png` | Topt, Tm, and kcat comparison across 12 characterised PETases; Topt vs Tm scatter with regression |
| `02_activity_penalty_at_37c.png` | Percentage of peak activity lost when each benchmark enzyme is cooled from Topt to 37°C |
| `03_thermostability_prediction.png` | ProtParam instability index and GRAVY distributions split by thermophilic condition label |
| `03_gravy_vs_topt.png` | GRAVY score vs published temperature optimum for characterised PETases |
| `04_room_temp_candidates.png` | Top 20 mesophilic PET organisms by evidence score; research funnel from all PET entries to sequence-supported mesophilic entries |
| `04_room_temp_environments.png` | Isolation environments of the mesophilic PET candidate pool |

---

## Benchmark Enzyme Citations

| Enzyme | Reference |
|---|---|
| IsPETase | Yoshida et al. 2016, Science 351(6278):1196-1199 |
| FAST-PETase | Lu et al. 2022, Nature 604(7906):662-667 |
| ThermoPETase | Cui et al. 2021, Nature Communications 12:4781 |
| LCC | Sulaiman et al. 2012, Applied and Environmental Microbiology 78(5):1556-1562 |
| ICCG-LCC | Tournier et al. 2020, Nature 580(7802):216-219 |
| TfCut2 | Roth et al. 2014, AMB Express 4:26 |
| PHL7 | Sonnendecker et al. 2022, ChemSusChem 15(9):e202101932 |
| HotPETase | Bell et al. 2022, ACS Catalysis 12(21):13392-13403 |
| BhrPETase | Shi et al. 2023, Nature Communications 14:1857 |
| CsPETase | Cheng et al. 2023, Nature Communications 14:6368 |
| DuraPETase | Cui et al. 2019, ACS Catalysis 9(8):6964-6972 |
| PET2 | Danso et al. 2018, Applied and Environmental Microbiology 84(12):e02773-17 |
