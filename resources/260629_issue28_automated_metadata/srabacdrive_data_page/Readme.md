# Organism Profiles

Per-organism HTML profile pages for 941 unique organisms relevant to bioplastic and plastic biodegradation research.
Sourced from PlasticDB (932 organisms) and NCBI SRA bioplastic-specific searches (9 additional).
All data pulled from live public APIs. Nothing fabricated or estimated.

## What is in this folder

```
organism-profiles/
  build.py                    Full-featured build script (reference implementation)
  fast_fetch.py               Original chunked fetcher for SRA + BacDive + pages
  enhance.py                  Enhanced fetcher: genome, PubMed, SRA bioplastic expand, pages
  data/
    sra_stats.csv             NCBI SRA run counts and metadata   (932/932 organisms)
    bacdive_data.csv          BacDive physiological data          (932/932 organisms)
    genome_data.csv           NCBI Assembly genome stats          (932/932 queried, 673 with assembly)
    pubmed_counts.csv         PubMed plastic/total paper counts   (932/932 organisms)
    extra_organisms.csv       Organisms added via SRA bioplastic search
  notebooks/
    01_fetch_sra_data.ipynb           Interactive: fetch SRA data from NCBI
    02_fetch_bacdive_data.ipynb       Interactive: scrape BacDive public pages
    03_generate_organism_pages.ipynb  Interactive: generate HTML pages
  pages/
    index.html               Searchable master index of all 941 organisms
    <slug>.html              One HTML page per organism (940 pages)
```

## Data sources (all real, no fabrication)

| Variable | Source | Coverage |
|---|---|---|
| Plastics degraded, enzyme, year, evidence | PlasticDB TSV | 2,535 entries, 932 organisms |
| SRA run count, platforms, strategies, bases | NCBI SRA E-utilities (esearch + esummary) | 932/932; 690 with runs |
| Culture temp, pH, oxygen, morphology, isolation | BacDive public strain pages | 932/932; 474 found |
| Genome size, assembly level, contig N50, coverage | NCBI Assembly (esearch + esummary + meta XML) | 932/932; 673 with assembly |
| Plastic/bioplastic paper count, total papers | NCBI PubMed esearch | 932/932; 760 with plastic papers |

## Bioplastic relevance classification

Each plastic type in PlasticDB is classified as **bioplastic** (green) or **conventional** (blue):

**Bioplastics**: PHA, PHB, PHBV, PHO, PHBH, PHC, PHV, P3HP, P4HB, P34HB, PCL, PLA, PBS, PBSA,
PBAT, PBSeT, PVA, PEF, PEA, PPL, Ecovio-FT, and all blend variants.

**Conventional plastics**: PU, LDPE, HDPE, LLDPE, PE, PET, PETG, PP, PS, PVC, Nylon, PC, PES, NR, PEG.

Out of 932 organisms: **440 have at least one bioplastic entry** (tagged "Bioplastic research" or
"Bioplastic + Conventional" in the header and index).

## Charts on each page (Chart.js, CDN)

Every organism page includes four interactive charts derived from real data:

1. **Plastic degradation profile** (horizontal bar) -- entry count per plastic type; bioplastics green,
   conventional blue.
2. **Research timeline** (bar chart) -- PlasticDB entries by publication year.
3. **Evidence methods** (horizontal bar) -- breakdown of analytical evidence used (weight loss,
   spectrophotometry, HPLC, SEM/FTIR, CO2, clear zone, etc.).
4. **Enzyme families** (doughnut) -- cutinase, lipase, depolymerase, PETase, etc.

## Index features

`pages/index.html` provides:

- Live text filter by organism name
- Filter buttons: All / Bioplastic / Conventional only / SRA-expanded
- Columns: organism, relevance badge, plastics pills, PlasticDB entries, SRA runs, PubMed count,
  Genome assembly (yes/no), BacDive (yes/no), first research year

## Notable data points (real values)

| Organism | Genome | Assembly | SRA runs | PubMed plastic papers |
|---|---|---|---|---|
| Ideonella sakaiensis | 6.18 Mbp | Complete Genome | 14 | 170 |
| Cupriavidus necator | 8.37 Mbp | Complete Genome | - | 1,405 |
| Thermobifida fusca | 3.76 Mbp | Complete Genome | 52 | - |
| Aspergillus tubingensis | 35.48 Mbp | Complete Genome | 143 | 176 |
| Pseudomonas putida | - | - | - | 6,977 |
| Bacillus megaterium | 5.74 Mbp | Complete Genome | 10,580 | 2,771 |

## Regenerating / extending the data

From the `organism-profiles/` directory:

```bash
# Status overview
python enhance.py --mode status

# Fetch or refresh genome data (NCBI Assembly)
python enhance.py --mode genome --batch 80

# Fetch or refresh PubMed counts
python enhance.py --mode pubmed --batch 100

# Search SRA for new bioplastic-specific organisms
python enhance.py --mode sra-expand

# Regenerate all HTML pages from cached CSVs
python enhance.py --mode pages
```

All modes are idempotent and resume from the cache automatically.

## Rate limits

- NCBI SRA / Assembly / PubMed: 0.35-0.4 s sleep per request (stays under 3 req/s)
- BacDive: 0.5 s sleep per request

## Gotchas

- NCBI Assembly `esummary` uses `assemblystatus` (not `assemblylevel`) and buries genome
  size inside a `meta` XML blob under `<Stat category="total_length">`. The `gcpercent`
  field does not exist in the JSON summary; GC% requires downloading the assembly stats file.
- BacDive search endpoint is `/search?search=` (not `/?search=`).
- The `sra_stats.csv` platform/strategy columns are sampled from the top 5 SRA runs only;
  `sra_run_count` is exact (from esearch).
