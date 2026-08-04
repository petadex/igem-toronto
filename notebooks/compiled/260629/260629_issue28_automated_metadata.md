## Automated Metadata Propagation & BacDive Integration

Lead     : `Amar / @AmDeep` and `Dennis / @UofTDriv`

Issue    : [Github Issue #28](https://github.com/petadex/igem-toronto/issues/28) — _Automated Metadata_

Start    : `2026-06-29`

Complete : `2026-06-29`

Files    : `~/igem-toronto/files/260629_issue28_automated_metadata/` — local working directory (active analysis, intermediate outputs)

S3 files : `s3://petadex/automated_metadata/` — remote archive (final outputs, shared with team)

---

### Data Accessed: 2026-06-29
```bash
# commands used to pull input data
aws s3 cp s3://petadex/sra/petadex_biosamples.csv ./resources/260629_issue28_automated_metadata/data/
aws s3 cp s3://petadex/accessions/accessions_sequences.csv ./resources/260629_issue28_automated_metadata/data/
```

## Introduction

The PETadex database compiles millions of Logan sequence candidates (predominantly PETase-like enzymes) identified in environmental metagenomic datasets. However, these sequences lack critical ecological context—such as the temperature, pH, and growth media of the source environments. Downstream machine learning tasks, such as Dennis's protein thermostability predictor, depend heavily on these environmental labels.

To address this characterization gap, we propagate metadata from NCBI SRA libraries and BioSamples. Using BioSample accessions, we map sequencing runs to detected TaxIDs using NCBI SRAstat/Athena tables. We then query the BacDive API to fetch physiological and culture traits for these organisms. Finally, by joining the taxa relative abundances with BacDive traits, we perform statistical inference to calculate abundance-weighted environmental conditions (pH, temperature, salinity, etc.) for each sample, which are subsequently mapped back to the candidate sequences.

## Objectives

1. **SRA-to-Taxonomy Mapping**: Link the 6.8 million BioSamples from `petadex_biosamples.csv` to their corresponding SRA runs and the taxonomic IDs of detected organisms via AWS Athena queries on the SRA Metadata and SRA STAT tables.
2. **Organism Table & Web Integration**: Build a normalized `organism` database table in the SQL dump and create an Organism tab/page on the PETadex website linking to BacDive.
3. **BacDive Trait Harvesting**: Extract physiological and growth conditions (temperature, pH, growth media) from the BacDive database for all unique TaxIDs using deduplicated queries to the BacDive API.
4. **Weighted Metadata Inference**: Compute statistical abundance-weighted temperature and pH ranges/means for each BioSample to fill in blank metadata fields, enabling environmental annotation of PETase sequences.

---

## Materials and Methods

### System Initialization

- Compute environment: AWS Cloud EC2 instance, remote PostgreSQL read replica (`serratus-aurora` for Logan metadata).
- Software: R (v4.2.0) and Python (v3.11.0) with `pandas`, `duckdb` (v0.8.1), and `requests` libraries.
- Setup document: [AWS setup guide](https://github.com/orgs/petadex/projects/4) / Dry Lab AWS thread.

### Data Initialization

The input files are obtained from local resources extracted from the S3 bucket:

```bash
# Accessed: 2026-06-29
# Input datasets:
# - petadex_biosamples.csv: List of 6.81M unique BioSample SAMN IDs linked to PETadex runs.
# - accessions_sequences.csv: GenBank accession IDs and amino acid sequences for candidate enzymes.
```

### Processing Steps

The execution workflow is structured as follows:

```bash
# Step 1: Extract unique BioSample accessions and query SRA/SRAstat tables on Athena
# to join SAMN -> SRS -> SRR -> TaxID and counts (abundances).
Rscript scripts/extract_unique_samns.R petadex_biosamples.csv athena_input_runs.csv

# Step 2: Extract unique TaxIDs from the SRA mapping results to deduplicate BacDive API requests.
python scripts/extract_unique_taxids.py athena_output_taxa.csv unique_taxa.csv

# Step 3: Fetch physiological traits (temperature, pH, growth media) from BacDive API
# for all unique taxa, saving to a local cache.
python scripts/query_bacdive_api.py unique_taxa.csv local_bacdive_cache.csv

# Step 4: Left join taxonomic abundances per BioSample with BacDive trait data,
# and compute abundance-weighted mean, median, and IQR boundaries for temperature and pH.
Rscript scripts/calculate_weighted_traits.R athena_output_taxa.csv local_bacdive_cache.csv biosample_environmental_report.csv

# Step 5: Join the computed environmental metadata back to PETase sequences via run mappings
# and load into SQL database as 'organism' and 'biosample_metadata' tables.
psql -h petadex.ccz9y6yshbls.us-east-1.rds.amazonaws.com -U user -d petadex -f sql/load_metadata_tables.sql
```

All intermediate data and scripts are archived in S3.

```bash
# Uploaded: 2026-06-29
aws s3 cp biosample_environmental_report.csv s3://petadex/automated_metadata/
aws s3 cp local_bacdive_cache.csv s3://petadex/automated_metadata/
```

## Results & Discussion

**Key metric**: Target is to resolve metadata for over **95%** of the 6.8 million BioSamples linked to PETase sequences.

**Breakdown by category (Expected)**:
- BioSamples with direct NCBI metadata (lat/lon, country): 3.8M (approx 55%)
- BioSamples with BacDive-propagated environmental metadata: ~2.8M additional samples
- BioSamples unresolvable (lacking taxonomy or BacDive coverage): <5%

**Observations**:
- Preliminary testing indicates environmental propagation via metagenomic taxonomy (e.g. ocean water) matches known site conditions closely, whereas host-associated samples (e.g. gut metagenomes) yield narrow, constant ranges (37°C, pH ~7).
- Restricting BacDive queries to deduplicated unique taxa prevents rate-limiting blocks and minimizes API latency.
- Statistical inference using abundance-weighted IQR and trimmed means reduces the influence of noisy taxa detections and transient species.

**Follow-up questions**:
- How should we resolve edge cases where a metagenomic sample contains equal proportions of thermophilic and mesophilic species?
- Should we incorporate auxiliary physical databases (e.g. World Ocean Atlas, soil databases) for BioSamples that have coordinates but lack taxonomy/BacDive entries?
- How do we map predicted environmental temperatures directly to the downstream ESM-C thermostability model's input features?
