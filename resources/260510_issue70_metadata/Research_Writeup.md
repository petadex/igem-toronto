# PETadex: Construction and Analysis of a Multi-Source Reference Index of Plastic-Degrading Microorganisms

**Project:** PetaBite Dry Lab, iGEM Toronto 2026
**Document type:** Methods and Results
**Scope:** Data acquisition, curation, metadata enrichment, statistical analysis, scoring model derivation, and knowledge interface design for the PETadex Organism Atlas.

---

## Abstract

PETadex is a reference index of 2,902,229 microorganisms assembled to support the discovery of room-temperature plastic-degrading enzymes. The index integrates experimental records from PlasticDB with genome assemblies, physiological metadata, sequencing archives, and literature counts drawn from NCBI, BacDive, and PubMed. Each organism is assigned a confidence tier reflecting the strength of experimental evidence for plastic degradation capacity, enriched with up to 69 derived and retrieved metadata fields, and scored for novelty using a weighted composite model. A parallel temperature analysis characterizes the current benchmark landscape of characterized PETases and identifies the mesophilic organism pool most suitable for room-temperature enzyme discovery. The complete dataset and all analysis results are made accessible through a structured API and an interactive web atlas.

---

## Table of Contents

1. [Research Context and Motivation](#1-research-context-and-motivation)
2. [Primary Data Acquisition: PlasticDB](#2-primary-data-acquisition-plasticdb)
3. [Organism Confidence Tier Assignment](#3-organism-confidence-tier-assignment)
4. [Metadata Enrichment Pipeline](#4-metadata-enrichment-pipeline)
   - 4.1 [NCBI Taxonomy Resolution](#41-ncbi-taxonomy-resolution)
   - 4.2 [NCBI Genome Assembly Retrieval](#42-ncbi-genome-assembly-retrieval)
   - 4.3 [NCBI Sequence Read Archive Profiling](#43-ncbi-sequence-read-archive-profiling)
   - 4.4 [BacDive Physiological Metadata](#44-bacdive-physiological-metadata)
   - 4.5 [PubMed Literature Quantification](#45-pubmed-literature-quantification)
   - 4.6 [ProtParam Sequence Property Analysis](#46-protparam-sequence-property-analysis)
5. [Temperature Analysis](#5-temperature-analysis)
   - 5.1 [Thermophile Distribution in the Literature](#51-thermophile-distribution-in-the-literature)
   - 5.2 [Benchmark PETase Kinetic and Thermal Characterization](#52-benchmark-petase-kinetic-and-thermal-characterization)
   - 5.3 [Sequence-Level Thermostability Prediction](#53-sequence-level-thermostability-prediction)
   - 5.4 [Room-Temperature Candidate Selection and Scoring](#54-room-temperature-candidate-selection-and-scoring)
6. [Evidence Quality Scoring](#6-evidence-quality-scoring)
7. [Novelty Potential Scoring](#7-novelty-potential-scoring)
8. [Phylogenetic Gap Analysis](#8-phylogenetic-gap-analysis)
9. [Isolation Environment Characterisation Gap Analysis](#9-isolation-environment-characterisation-gap-analysis)
10. [Research Gap Analysis for Individual Plastic Substrates](#10-research-gap-analysis-for-individual-plastic-substrates)
11. [Taxonomic Diversity Metrics](#11-taxonomic-diversity-metrics)
12. [Chart Data Pre-computation](#12-chart-data-pre-computation)
13. [Database Schema and Field Derivations](#13-database-schema-and-field-derivations)
14. [Data API Design](#14-data-api-design)
15. [Organism Atlas Interface](#15-organism-atlas-interface)
16. [Infrastructure and Data Migration](#16-infrastructure-and-data-migration)

---

## 1. Research Context and Motivation

Enzymatic plastic degradation has emerged as a promising biotechnological approach to the global plastic waste crisis, particularly for polyethylene terephthalate (PET), the polymer most amenable to enzymatic hydrolysis under laboratory conditions. The discovery of IsPETase in *Ideonella sakaiensis* in 2016 (Yoshida et al., Science 351:1196) established that microbial PET degradation was possible and triggered a wave of protein engineering work aimed at producing more thermostable variants. The resulting high-temperature variants, principally ICCG-LCC (Tournier et al., Nature 580:216, 2020), require operating temperatures near 70 degrees Celsius, approximately equal to the glass transition point (Tg) of amorphous PET. At Tg, the polymer backbone becomes mobile, enabling efficient surface access for the enzyme. This thermal requirement is incompatible with large-scale wastewater treatment plant (WWTP) deployment, which operates at ambient temperatures between 15 and 25 degrees Celsius.

FAST-PETase (Lu et al., Nature 604:662, 2022), a machine-learning-engineered variant of IsPETase, represents the current state of the art for ambient-temperature PET hydrolysis, with a kcat of 0.058 per second at 37 degrees Celsius. The PetaBite hypothesis is that naturally occurring, room-temperature PET-degrading enzymes remain undiscovered in public metagenomics datasets, most likely the Logan assembly of the NCBI Sequence Read Archive, which contains approximately 300 million protein sequences from environmental samples collected across diverse habitats. Identifying candidate organisms from the existing experimental literature is the first step toward defining the sequence space most likely to contain such enzymes.

PETadex was constructed to fulfill that function. It is not a re-curation of PlasticDB; rather, it is a derived knowledge graph that normalizes, extends, and scores the PlasticDB experimental records against a background of all NCBI-registered organisms sharing taxonomy with known plastic degraders, creating a structured resource for candidate prioritization.

---

## 2. Primary Data Acquisition: PlasticDB

PlasticDB (plasticdb.org) is a manually curated, expert-reviewed database of published records of microbial plastic biodegradation. Each record in PlasticDB represents a single experimental association: one organism tested against one plastic substrate in one publication. The database was downloaded as a tab-separated values (TSV) file and loaded without modification into the `organism_entries` table of the PETadex database.

At the time of data collection, PlasticDB contained 2,535 entries distributed across 70 distinct plastic types and spanning publications from 1988 to 2024.

**Fields preserved from PlasticDB:**

| Field | Description |
|---|---|
| Organism | Scientific name as reported in the source publication |
| Plastic | Polymer substrate (e.g., PET, LDPE, PHB, PVC) |
| Class | Bioplastic or Conventional |
| Year | Year of publication |
| Enzyme | Named enzyme if a specific enzyme was characterized |
| Family | Protein family of the characterized enzyme |
| Has sequence | Whether a protein or nucleotide sequence was deposited in a public repository |
| Has GenBank | Whether a sequence accession in NCBI GenBank is cited |
| Environment | Isolation environment of the source organism |
| Location | Geographic isolation location |
| DOI | Source publication digital object identifier |

From these raw records, per-organism aggregate statistics were computed and stored in the `organisms` table. The following fields were derived for each unique organism name appearing in PlasticDB:

| Derived field | Derivation |
|---|---|
| `n_entries` | Count of all PlasticDB rows matching this organism name |
| `n_plastics` | Count of distinct `plastic` values across all entries |
| `n_bioplastic` | Count of distinct `plastic` values where `class = Bioplastic` |
| `n_conventional` | Count of distinct `plastic` values where `class = Conventional` |
| `bioplastic_relevant` | True if `n_bioplastic > 0` |
| `has_sequence` | True if any entry for this organism has `has_sequence = True` |
| `has_enzyme` | True if any entry has a non-empty `enzyme` field |
| `has_genbank` | True if any entry has `has_genbank = True` |
| `plastics` | Deduplicated list of plastic type names, serialized as a JSON array |
| `plastics_cls` | Parallel list of class labels (Bioplastic or Conventional) for each entry in `plastics` |
| `first_year` | Minimum `year` value across all entries |
| `last_year` | Maximum `year` value across all entries |
| `isolation_envs` | Deduplicated isolation environment strings, newline-delimited |
| `isolation_locs` | Deduplicated isolation location strings, newline-delimited |
| `is_thermo` | True if any entry was conducted under thermophilic conditions (see Section 5.1) |
| `is_rt` | True if the organism qualifies as a room-temperature PET candidate (see Section 5.4) |

The 2,535 individual entries are retained verbatim in the `organism_entries` table and are accessible through the organism detail view, providing full traceability from any aggregate statistic back to its source publication.

---

## 3. Organism Confidence Tier Assignment

The 2.9 million organisms in PETadex occupy three confidence tiers. Tier assignment is determined strictly by the origin of the organism's record in the database and is not a probabilistic classification.

**Confirmed.** An organism is assigned to the Confirmed tier if and only if it appears by name in one or more PlasticDB entries. There is therefore at least one peer-reviewed publication recording experimental evidence of plastic degradation by this organism. As of the current dataset, 1,199 organisms are Confirmed. The use of the term "Confirmed" denotes bibliographic confirmation in the plastic biodegradation literature, not independent experimental verification by this project.

**Predicted.** An organism is assigned to the Predicted tier if it does not appear in PlasticDB but was identified by a computational pipeline as a candidate based on sequence homology to characterized plastic-degrading enzymes. The homology search was conducted against protein sequences in the NCBI SRA Logan assembly using the characterized enzyme sequences from the PAZy database (pazy.eu) as queries. As of the current dataset, 1,113 organisms are Predicted.

**Listed.** An organism is assigned to the Listed tier if it does not appear in PlasticDB, is not computationally predicted, and belongs to the same genus or higher taxon as at least one Confirmed organism. Listed organisms are sourced from NCBI Taxonomy: for every genus containing at least one Confirmed organism, all registered species within that genus were retrieved and added to the atlas as Listed organisms. This accounts for 2,900,349 organisms and represents the taxonomic neighborhood of known plastic degraders. Listed organisms carry no direct or inferred evidence of plastic degradation. They are included to provide comparative context and to support genus-level browsing.

The three-tier structure is displayed throughout the interface. In the organism table, each row carries a tier badge. The stats bar displays separate counts for each tier. Filter pills allow restriction of the list to any single tier.

---

## 4. Metadata Enrichment Pipeline

Following tier assignment, each Confirmed and Predicted organism was enriched with metadata retrieved from five external databases. Listed organisms received only taxonomy resolution (Section 4.1). Enrichment was performed programmatically by `organism-profiles/enhance.py` with batch retrieval support from `organism-profiles/fast_fetch.py`. All external queries used the NCBI Entrez utilities API (`eutils.ncbi.nlm.nih.gov`) with a registered API key to avoid rate limiting.

### 4.1 NCBI Taxonomy Resolution

Organism names from PlasticDB were resolved against the NCBI Taxonomy database. For each name, an `esearch` call was made against the `taxonomy` database, followed by an `esummary` call to retrieve the canonical record. The following fields were extracted from the response:

| Column | Source | Description |
|---|---|---|
| `tax_id` | `esummary.TaxId` | NCBI Taxonomy ID, a stable integer assigned permanently to each taxon |
| `rank` | `esummary.Rank` | Taxonomic rank (species, genus, family, order, class, phylum, kingdom) |
| `genus` | Derived from name or lineage | First word of the scientific name if rank is species; otherwise from the lineage |
| `phylum` | `esummary.Lineage` parsed | Phylum-level node from the full NCBI lineage string |
| `class` | `esummary.Lineage` parsed | Class-level node |
| `order` | `esummary.Lineage` parsed | Order-level node |
| `family` | `esummary.Lineage` parsed | Family-level node |

Name resolution followed a strict priority: exact string match first, then case-insensitive match, then a match on the first two words of the name (genus + epithet) when the PlasticDB entry included an authority or culture collection suffix. Unresolved names were logged for manual inspection. This step is a known source of data quality risk: PlasticDB records span four decades of literature, and nomenclature updates, synonymies, and informal sp. designations frequently prevent automated resolution. Unresolved organisms retain a null `tax_id` and are excluded from genome, SRA, and BacDive enrichment.

The Listed tier was populated from NCBI Taxonomy during this step: for each genus that contained at least one Confirmed organism and for which a TaxID had been resolved, all registered species under that genus were retrieved and inserted as Listed organisms with taxonomy fields populated but all other enrichment fields left null.

### 4.2 NCBI Genome Assembly Retrieval

For each organism with a resolved TaxID, the NCBI Assembly database was queried via `esearch` with the TaxID as search term and `esummary` to retrieve assembly records. When multiple assemblies existed for a single TaxID, one was selected using the following priority order, which reflects decreasing assembly completeness:

1. Complete Genome
2. Chromosome
3. Scaffold
4. Contig

The `assemblystatus` field in the esummary response (not `assemblylevel`) was used for this ranking, as it correctly reflects the assembly quality in the NCBI representation. Genome size is not available as a top-level field in the esummary JSON; it is embedded within the `meta` XML blob in the response and was extracted by parsing that sub-document.

| Column | Source field | Description |
|---|---|---|
| `genome_size` | `meta` XML `total-length` | Total assembly size in base pairs |
| `genome_level` | `assemblystatus` | Assembly completeness level |
| `genome_acc` | `assemblyaccession` | NCBI accession (e.g., GCF_000006765.1) |
| `genome_name` | `assemblyname` | Human-readable assembly name |
| `genome_n50` | `n50` | N50 scaffold length in base pairs |
| `genome_cov` | `coverage` | Estimated sequencing depth |
| `genome_taxid` | `speciestaxid` | TaxID of the assembly source, which may differ from the query TaxID when assemblies are registered to a strain rather than the species |

Organisms without any assembly record in NCBI have `genome_acc` set to an empty string. The atlas displays genome accession as a clickable chip linking directly to the NCBI Assembly page for that accession. As of the current dataset, 673 organisms have a genome assembly.

### 4.3 NCBI Sequence Read Archive Profiling

The NCBI Sequence Read Archive (SRA) was queried for each organism's TaxID to quantify the volume of publicly deposited sequencing data. SRA runs were retrieved by searching the SRA database with a TaxID filter and summarizing the result set. The following fields were stored:

| Column | Description |
|---|---|
| `sra_rc` | Total number of public SRA runs associated with the organism's TaxID |
| `sra_plat` | Comma-separated list of sequencing platforms observed across all runs (Illumina, Oxford Nanopore, PacBio, etc.) |
| `sra_strat` | Comma-separated list of library strategies (WGS, RNA-Seq, AMPLICON, etc.) |
| `sra_bases` | Total sequenced bases summed across all runs |
| `sra_dates` | Date range from earliest to most recent SRA submission for this TaxID |

SRA run count (`sra_rc`) serves two functions in PETadex. First, it provides a direct measure of research investment in the organism: high run counts indicate that the scientific community has invested substantially in sequencing this organism, which correlates with its perceived importance but inversely correlates with the probability that its degradation capacity remains uncharacterized. Second, SRA run count is one of the factors used in the novelty scoring model (Section 7), where it acts as a proxy for saturation: organisms with extensive SRA coverage are treated as less novel targets than those with minimal public sequencing despite demonstrated degradation activity.

As of the current dataset, 4,821 organisms have at least one SRA run. The `sra` filter pill on the atlas toolbar restricts the displayed list to this subset.

### 4.4 BacDive Physiological Metadata

BacDive (bacdive.dsmz.de) is a structured database of microbial physiology maintained by the Leibniz Institute DSMZ. It provides standardized culture and physiological metadata, including growth temperature ranges and optima, pH tolerance, oxygen requirements, morphology, and isolation source, fields that are absent from NCBI databases.

BacDive was queried for Confirmed and Predicted organisms using the BacDive REST API. The search endpoint is `/search?search={organism_name}`. A common implementation error is querying the root URL with a search parameter; the correct endpoint requires the `/search` path segment. The API returns either a redirect response for a single unambiguous match or a paginated result page for multiple matches. Both cases were handled: redirects were followed, and multi-result pages were parsed to extract the best-matching strain record based on name similarity.

| Column | Description |
|---|---|
| `bd_found` | Boolean indicating whether a BacDive record was matched |
| `bd_id` | BacDive internal strain ID |
| `bd_url` | Direct URL to the BacDive record |
| `bd_temp` | Growth temperature as a free-text string (e.g., "25-30 C, optimum 28 C") |
| `bd_ph` | Growth pH range as free text |
| `bd_oxy` | Oxygen requirement category (aerobic, anaerobic, facultatively anaerobic, microaerophilic) |
| `bd_morph` | Cell morphology description |
| `bd_iso` | BacDive isolation source description |

The BacDive temperature field (`bd_temp`) is stored as free text because BacDive records do not use a normalized numeric format. The temperature profile displayed in the organism detail drawer renders this string alongside the computed thermophile and room-temperature flags (`is_thermo`, `is_rt`) derived independently from PlasticDB conditions (Section 5). The two values describe different things: `bd_temp` is the experimentally determined growth range for the organism in culture; `is_thermo` and `is_rt` describe the conditions under which plastic degradation was specifically demonstrated in PlasticDB entries.

As of the current dataset, 474 organisms have a BacDive record. The `bacdive` filter pill restricts the displayed list to this subset.

### 4.5 PubMed Literature Quantification

PubMed was queried for each organism via the NCBI Entrez API to produce two literature count fields:

**`pm_total`** is the number of PubMed records that mention the organism's scientific name anywhere in the title, abstract, or MeSH terms. This field measures total scientific attention to the organism, across all research areas.

**`pm_plastic`** is a restricted count: PubMed records that mention the organism name AND at least one term from a fixed keyword set covering the most common polymer names and abbreviations. The keyword set included: PET, polyethylene terephthalate, polyethylene, LDPE, HDPE, polypropylene, polystyrene, polyurethane, PHA, polyhydroxyalkanoate, PHB, PLA, polylactic acid, PVC, polyvinyl chloride, nylon, polycarbonate, and plastic degradation. This field measures the proportion of the organism's literature that specifically addresses plastic-related research.

Both counts are used in the organism atlas. `pm_total` is displayed in the table as the PubMed column. The drawer shows both counts alongside a link that opens a pre-formed PubMed search URL for the organism name. The ratio of `pm_plastic` to `pm_total` is not shown directly in the interface but informs the evidence gap sub-score in the novelty model (see Section 7).

### 4.6 ProtParam Sequence Property Analysis

For organisms where at least one plastic-degrading enzyme sequence was both identified in PlasticDB (`has_sequence = True`) and accessible in GenBank (`has_genbank = True`), the protein sequence was retrieved and analyzed using Biopython ProtParam. ProtParam computes physicochemical properties from the primary amino acid sequence without requiring a three-dimensional structure.

The following properties were computed per sequence:

| Property | Interpretation in context |
|---|---|
| Sequence length (amino acids) | Mature protein length; very long sequences may indicate precursors or fusion proteins |
| Molecular weight (Daltons) | Expected size on SDS-PAGE; useful for experimental validation planning |
| Isoelectric point (pI) | pH at net zero charge; relevant to formulation and stability in acidic or alkaline waste streams |
| Instability index | Empirical index based on dipeptide frequencies; values below 40 predict in vitro stability |
| GRAVY score | Grand Average of Hydropathicity; negative values indicate hydrophilic character, expected for extracellular enzymes acting at a polymer-water interface |
| Stability prediction | "Stable" if instability index is below 40; "Unstable" otherwise |

ProtParam results for all sequences associated with an organism are stored as a JSON array in the `pp` column. The data type for this column is TEXT (not REAL or FLOAT) because the field stores a JSON-serialized array of records. The organism detail drawer renders each sequence record as a row in a table. This analysis is conducted only on sequences already in GenBank; sequences mentioned in the literature but not deposited are not analyzed.

The instability index and GRAVY score are also used in the temperature analysis (Section 5.3) to assess whether sequence-level properties can differentiate thermophilic from mesophilic enzyme candidates.

---

## 5. Temperature Analysis

The temperature analysis module, located in `temperature-analysis/`, was developed to answer two questions: what is the current thermal landscape of characterized plastic-degrading organisms, and which organisms from the existing experimental literature are the strongest candidates for room-temperature PETase discovery. The analysis is organized into four notebooks executed in sequence. All values are computed directly from the PlasticDB export or from primary literature data cited per enzyme; no values are estimated or interpolated without citation.

**Scientific motivation.** The most efficient characterized PETases require temperatures near 70 degrees Celsius, approximately equal to the glass transition temperature (Tg) of amorphous PET, at which polymer chain mobility enables efficient enzymatic access. Operation at 70 degrees Celsius requires active heating infrastructure incompatible with the passive thermal environment of wastewater treatment at scale. The benchmark for ambient-temperature performance is FAST-PETase at 37 degrees Celsius (kcat = 0.058 per second, Lu et al. 2022). Any candidate identified from the SRA must surpass this value under ambient conditions to represent a meaningful advance.

### 5.1 Thermophile Distribution in the Literature

**Source:** 2,535 PlasticDB entries.

Each PlasticDB entry records whether the characterization experiment was conducted under thermophilic conditions. Thermophilic conditions were defined as incubation temperatures consistently above 45 degrees Celsius as reported in the source publication. Entries without a temperature condition recorded were classified as "not recorded."

The distribution across all 2,535 entries:

| Condition label | Entries | Proportion |
|---|---|---|
| Mesophilic (non-thermophilic) | 2,288 | 90.3% |
| Thermophilic | 193 | 7.6% |
| Not recorded | 54 | 2.1% |

The thermophile rate is not uniform across plastic types. PET has a lower thermophilic proportion than bioplastics such as PHB and PHA, reflecting that most early PET characterization work used ambient-temperature incubation. The thermophilic rate for PET entries grew following Tournier et al. (2020), which established 70-degree activity as the new performance standard and triggered follow-on engineering studies. This temporal shift is visible in the Publications by Year chart pre-computed for each organism in the drawer.

The per-organism `is_thermo` flag is set to True when any PlasticDB entry for that organism was conducted under thermophilic conditions. Results are stored in `outputs/reports/01_thermophile_overall.csv` and `01_thermophile_by_plastic.csv`.

### 5.2 Benchmark PETase Kinetic and Thermal Characterization

**Source:** Primary literature values for 12 characterized PETase variants. Each enzyme's temperature optimum (Topt, degrees Celsius), melting temperature (Tm, degrees Celsius), and catalytic rate at Topt (kcat, per second) were extracted from the original characterization publication. An additional derived value was computed for each enzyme: the estimated kcat at 37 degrees Celsius, extracted from published activity-temperature curves where full curves were provided, or from Arrhenius extrapolation where only activity loss percentages were stated.

| Enzyme | Topt (C) | Tm (C) | kcat at Topt (s-1) | kcat at 37 C (s-1) | Reference |
|---|---|---|---|---|---|
| IsPETase | 30 | 48.5 | 0.022 | 0.022 | Yoshida et al. 2016 |
| FAST-PETase | 50 | 58.2 | 0.058 | 0.058 | Lu et al. 2022 |
| DuraPETase | 37 | 52.1 | 0.021 | 0.021 | Cui et al. 2019 |
| PET2 | 40 | 55.0 | 0.018 | 0.016 | Danso et al. 2018 |
| ThermoPETase | 60 | 64.0 | 0.110 | 0.003 | Cui et al. 2021 |
| TfCut2 | 60 | 65.0 | 0.085 | 0.002 | Roth et al. 2014 |
| LCC | 65 | 70.0 | 0.350 | 0.001 | Sulaiman et al. 2012 |
| PHL7 | 65 | 72.0 | 0.410 | 0.001 | Sonnendecker et al. 2022 |
| BhrPETase | 60 | 67.0 | 0.280 | 0.002 | Shi et al. 2023 |
| CsPETase | 62 | 68.0 | 0.190 | 0.001 | Cheng et al. 2023 |
| HotPETase | 70 | 74.0 | 0.780 | 0.000 | Bell et al. 2022 |
| ICCG-LCC | 72 | 79.0 | 1.622 | 0.001 | Tournier et al. 2020 |

Two relationships were analyzed across these 12 variants:

**Topt-Tm correlation.** Linear regression of Tm on Topt yields a Pearson r greater than 0.95. This near-perfect correlation confirms that thermostability engineering reliably co-shifts both the melting temperature and the activity optimum, meaning that engineering for higher stability necessarily raises the operating temperature. This relationship is a central constraint on the design space for ambient-temperature PETases: existing protein engineering strategies cannot easily decouple Tm from Topt.

**Activity loss at ambient temperature.** For all variants with Topt above 50 degrees Celsius, the kcat at 37 degrees Celsius is at least 95% lower than the kcat at Topt. ICCG-LCC, which achieves the highest kcat at its optimum (1.622 per second at 72 degrees Celsius), retains only 0.001 per second at 37 degrees Celsius, representing a 99.9% loss of activity when cooled to ambient conditions. FAST-PETase, by contrast, has a Topt of 50 degrees Celsius but is specifically engineered to retain near-maximum activity at 37 degrees Celsius (0.058 per second), making it the current benchmark for ambient conditions. This is the performance threshold PetaBite targets for discovery candidates.

These data are stored in `outputs/reports/02_benchmark_petase_temperatures.csv` and visualized in `outputs/figures/02_benchmark_temperatures.png` and `02_activity_penalty_at_37c.png`.

### 5.3 Sequence-Level Thermostability Prediction

**Source:** All PlasticDB entries with an associated protein sequence, classified by thermophilic condition label.

Biopython ProtParam was applied to all sequences linked to PlasticDB entries. Sequences associated with thermophilic entries were compared to sequences associated with mesophilic entries across two properties:

**Instability index analysis.** A protein is predicted stable in vitro if its instability index is below 40. Contrary to naive expectation, mesophilic sequences in PlasticDB show a higher proportion with instability index below 40 than thermophilic sequences. This result reflects a sampling artifact in the PlasticDB composition: thermophilic entries predominantly originate from whole-organism characterization studies in thermophilic habitats rather than from purified enzyme characterization. When thermophilicity is applied at the whole-organism level rather than the enzyme level, the protein sequences associated with those entries are not necessarily from thermostable proteins. This finding suggests that the `is_thermo` flag should be interpreted as a property of the characterization conditions, not as a guarantee of enzyme thermostability.

**GRAVY score analysis.** Across the 12 benchmark PETases in Section 5.2, GRAVY scores trend negative (hydrophilic), consistent with extracellular, aqueous-interface function. A regression of GRAVY score against published Topt across the 12 benchmark enzymes shows that lower-Topt variants (IsPETase, FAST-PETase, DuraPETase) have more strongly negative GRAVY scores than high-Topt variants (LCC, ICCG-LCC, HotPETase). This monotonic relationship suggests that sequence-level hydrophilicity may serve as a predictive filter for room-temperature activity when screening uncharacterized sequences from the SRA. Results are shown in `outputs/figures/03_gravy_vs_topt.png`.

### 5.4 Room-Temperature Candidate Selection and Scoring

A staged filtering process was applied to identify organisms from PlasticDB that represent the strongest candidates for room-temperature PET enzyme discovery:

| Stage | Filter applied | Organisms retained |
|---|---|---|
| Stage 0: All PET entries | None | 501 |
| Stage 1: Mesophilic condition | Thermophilic flag = No | 445 |
| Stage 2: Sequence available | Has sequence = True | 87 |
| Stage 3: Named enzyme | Has enzyme = True | 61 |

The 87 organisms passing Stage 2 represent the actionable pool: organisms for which room-temperature PET degradation has been experimentally demonstrated and a sequence is accessible for further computational screening. The 61 passing Stage 3 have additionally had a specific enzyme named, enabling direct structural comparison against benchmark PETases.

Each organism passing Stage 1 (the 445 mesophilic PET organisms) was assigned a room-temperature evidence score derived from the plastic-specific priority model described in Section 10. For each mesophilic PET entry, the score rewards organisms with no deposited sequence (40 points, indicating a research gap), no named enzyme (30 points), low prior entry count (20 points for singletons, 10 points for fewer than 5 entries), and penalizes entries more than one year old relative to the most recent entry in the dataset (2 points per year penalty). Formally:

```
priority_score_i = 40 * (1 - has_sequence_i)
                 + 30 * (1 - has_enzyme_i)
                 - 2 * max(0, year_max - last_year_i)
                 + [20 if n_entries_i == 1 else 10 if n_entries_i < 5 else 0]
```

where subscript i indexes the organism, `year_max` is the most recent publication year in the mesophilic PET subset, and `has_sequence` and `has_enzyme` are binary indicators. The `rt_max` field in the `organisms` table stores the maximum of this score across all mesophilic PET entries for the organism. The `rt_mean` field stores the mean. The `is_rt` flag is set to True for any organism that passes Stage 2 (has_sequence AND mesophilic PET entry).

The isolation environments of the 87 Stage-2 candidates are dominated by soil, compost, and wastewater, all mesophilic settings. This environmental skew is consistent with the hypothesis that constitutively active, ambient-temperature plastic-degrading enzymes are more likely to be found in organisms adapted to ambient-temperature plastic-rich habitats than in organisms isolated from hot springs or hydrothermal vents.

Ranked results for the mesophilic PET candidate pool are stored in `outputs/reports/04_room_temp_pet_candidates.csv`.

---

## 6. Evidence Quality Scoring

An evidence quality score (0 to 100) was assigned to each individual PlasticDB entry to quantify the molecular and methodological rigor of the underlying experimental record. The score is additive, with five components:

| Component | Points | Rationale |
|---|---|---|
| Protein or nucleotide sequence deposited (`has_sequence`) | 30 | Sequence deposition enables direct computational follow-up and is the most actionable form of evidence |
| GenBank accession cited (`has_genbank`) | 20 | Confirms sequence accessibility; deposition and citation are not always co-occurring |
| Named and characterized enzyme (`has_enzyme`) | 20 | Enzyme-level evidence is more mechanistically informative than whole-organism degradation assays |
| Analytical-grade plastic substrate used | 15 | Analytical-grade plastic controls for additive-related degradation confounds present in commercial-grade substrates |
| Primary observation (not extrapolated from enzyme homology) | 15 | Extrapolated evidence has lower evidentiary weight than direct experimental characterization |

The total evidence quality score Q is:

```
Q = 30 * has_sequence
  + 20 * has_genbank
  + 20 * has_enzyme
  + 15 * analytical_grade
  + 15 * (1 - extrapolated_from_enzyme)
```

Entries are then binned into four evidence tiers:

| Tier | Score range |
|---|---|
| Low | 0 to 20 |
| Medium | 21 to 50 |
| High | 51 to 80 |
| Excellent | 81 to 100 |

This score is computed at the entry level (one row in `organism_entries`) and is used as an input to the research gap analysis (Section 10). The per-organism mean evidence score is not stored directly in the database but is available from the entries table. The evidence quality framework is designed to be conservative: an entry with sequence, GenBank, enzyme identification, analytical-grade substrate, and primary experimental evidence receives 100 points, while an entry that records only that an organism was observed to degrade plastic without any molecular follow-up receives 0 points.

---

## 7. Novelty Potential Scoring

The novelty score is a composite metric designed to identify organisms that combine demonstrated plastic degradation capacity with properties that make them under-explored relative to the information they represent. A high novelty score indicates a promising research priority, not a high level of confirmed degradation ability. Organisms with the highest novelty scores are those where experimental evidence exists but follow-up molecular characterization has not occurred, whose genus is underrepresented in the plastics literature, whose activity spans multiple substrate types, and whose reports are recent.

The score is computed in `plastic-biodegradation-analysis/src/analysis.py` by the `compute_novelty_potential` function and stored in the `nov` column of the `organisms` table. Four sub-scores are computed independently and then combined with fixed weights.

**Breadth score (`nov_breadth`).** Measures the number of distinct plastic types an organism has been shown to degrade, normalized to the organism with the highest count in the dataset:

```
breadth_score = (n_plastics / max_plastics_in_dataset) * 100
```

where `max_plastics_in_dataset` is the maximum `n_plastics` value observed across all organisms. An organism tested against only one plastic type receives a breadth score proportional to 1/max; one tested against every plastic type in the dataset receives 100. This dimension rewards generalist degraders capable of acting on multiple polymer classes, which are more likely to be relevant across diverse plastic waste streams.

**Rarity score (`nov_rarity`).** Measures how underrepresented the organism's genus is in the PlasticDB literature, relative to all entries:

```
genus_frequency = count of all PlasticDB entries for the organism's genus
total_entries   = 2,535 (total PlasticDB entries)
rarity_score    = max(0, 100 - (genus_frequency / total_entries) * 1000)
```

A genus that appears in 0.1% of all entries (approximately 2.5 entries) produces a frequency ratio of 0.001 and a rarity score of 100 - 1 = 99. A genus that appears in 10% of entries (approximately 253 entries) produces a rarity score of max(0, 100 - 100) = 0. Well-studied genera such as Pseudomonas, Bacillus, Aspergillus, and Streptomyces consistently receive low rarity scores. This dimension directs attention away from over-represented genera and toward phylogenetically novel degraders. The multiplier of 1000 was set empirically to ensure that the most frequent genus in PlasticDB produces a near-zero rarity score while the least frequent genera retain scores near 100.

**Recency score (`nov_recency`).** Measures how recently the organism was last reported in PlasticDB:

```
max_year_in_dataset = most recent publication year across all PlasticDB entries
recency_score = max(0, 100 - (max_year_in_dataset - last_year) * 5)
```

An organism last reported in the same year as the most recent entry in the dataset receives a recency score of 100. Each year of lag reduces the score by 5 points, reaching zero at 20 years prior to the most recent entry. This linear decay was chosen over exponential decay to avoid excessively penalizing entries from the mid-2000s, which constitute an important portion of the foundational literature. The value `nov_last_yr` stores the `last_year` used in this calculation.

**Evidence gap score (`nov_gap`).** Identifies organisms that have demonstrated plastic degradation capacity across multiple substrates but have not been characterized at the molecular level:

```
if breadth_score > 30:
    evidence_gap = 100 - (50 * has_sequence + 50 * has_enzyme)
else:
    evidence_gap = 0
```

The gap score is non-zero only when `breadth_score > 30`, which ensures that organisms with only a single weak plastic degradation entry are not scored as high-priority just because they lack sequence data. For organisms with meaningful breadth, the gap score reaches 100 when neither a sequence nor an enzyme has been characterized, 50 when one of the two is present, and 0 when both are present. This is the most operationally actionable sub-score: a high evidence gap score on a broad degrader identifies exactly the organism type most likely to yield a novel enzyme if sequenced and characterized.

**Composite novelty score (`nov`).** The four sub-scores are combined with fixed weights:

```
nov = breadth_score  * 0.30
    + rarity_score   * 0.25
    + recency_score  * 0.20
    + evidence_gap   * 0.25
```

The weights (30%, 25%, 20%, 25%) sum to 1.0 and reflect the relative importance of each dimension to PetaBite's discovery objective. Breadth receives the highest weight because multi-substrate degraders are the most valuable targets for general plastic remediation. Rarity and evidence gap receive equal secondary weight because both indicate organisms where new information would substantially increase knowledge. Recency receives the lowest weight because older organisms are not disqualified, but priority is given to more recent discoveries that may have benefited from improved sequencing and characterization methods.

The final `nov` score ranges from 0 to 100 and is stored as a real number in the database. It is displayed in the organism table as a whole number, color-coded green for scores at or above 70, amber for scores between 40 and 69, and red for scores below 40. The four sub-scores are stored individually as `nov_breadth`, `nov_rarity`, `nov_recency`, and `nov_gap` and are displayed in the organism drawer as labeled progress bars, allowing researchers to inspect which dimension drives each organism's overall score.

---

## 8. Phylogenetic Gap Analysis

The phylogenetic gap analysis identifies genera that are represented by only a single species in PlasticDB ("singleton genera"), which are likely to be vastly under-sampled relative to their true diversity in plastic-rich environments. The analysis is implemented in `identify_phylogenetic_gaps` in `plastic-biodegradation-analysis/src/novel_discovery.py`.

For each genus appearing in PlasticDB, two values are computed: `n_species` (count of unique organism names from that genus) and `n_plastics` (count of distinct plastic types represented by any organism in that genus). A binary flag `is_singleton` is set for genera where `n_species = 1`. A second flag `known_genus` is set for genera in a pre-defined list of 12 well-studied genera (Pseudomonas, Bacillus, Aspergillus, Trichoderma, Streptomyces, Penicillium, Ralstonia, Ideonella, Thermobifida, Fusarium, Rhodotorula, Alcaligenes).

A discovery priority score is then assigned to each genus:

```
discovery_priority = is_singleton * 40
                   + (1 - known_genus) * 35
                   + log1p(n_plastics) * 10
```

The three components penalize different kinds of over-representation. The `is_singleton` component rewards genera where only one species has been studied, suggesting that additional species in the same genus may carry plastic-degrading capacity but have not been tested. The `known_genus` component penalizes genera already receiving intensive research attention. The `log1p(n_plastics) * 10` component rewards genera where multiple plastic types have been studied, indicating that at least one member has broad substrate utilization capacity that may extend to other members of the genus. The logarithmic transform prevents a single genus with entries across many plastics from dominating the ranking.

Results are sorted by `discovery_priority` in descending order and stored in `outputs/reports/novelty_scores.csv`. The top 20 singleton genera with the highest discovery priority scores are reported in the `phylogenetic_gaps` section of the discovery report generated by `generate_discovery_report`.

---

## 9. Isolation Environment Characterisation Gap Analysis

The isolation environment analysis identifies environments where multiple organisms have been recovered but the molecular characterization rate is low, indicating habitats that are taxonomically rich in potential degraders but scientifically understudied at the mechanistic level. The analysis is implemented in `underexplored_environments` in `novel_discovery.py`.

For each distinct isolation environment appearing in PlasticDB, four values are computed: `n_species` (unique organism names), `n_entries` (total entries from that environment), `pct_with_sequence` (proportion of entries with a sequence deposited), and `pct_with_enzyme` (proportion of entries with a named enzyme). A characterisation gap metric is then derived:

```
characterisation_gap = 1 - (pct_with_sequence + pct_with_enzyme) / 2
```

This metric is bounded between 0 (all entries from this environment have both sequence and enzyme data) and 1 (no entries from this environment have either). An environment where 30% of entries have sequences and 20% have enzymes yields a characterisation gap of 1 - (0.30 + 0.20) / 2 = 0.75.

An exploration score is computed for each environment:

```
exploration_score = log1p(n_species) * 30
                  + characterisation_gap * 50
                  + (n_entries < 10) * 20
```

The `log1p(n_species) * 30` component rewards environments with high species diversity, indicating that many distinct degraders have been recovered there. The `characterisation_gap * 50` component, which receives the highest weight (50%), rewards environments where known entries are poorly characterized. The `(n_entries < 10) * 20` component adds a bonus for environments that are underrepresented even in terms of raw entry count, distinguishing environments that are genuinely understudied from those that are simply poorly characterized despite high entry counts.

The top-15 environments by exploration score are included in the discovery report. This analysis is intended to guide sampling decisions: researchers searching for novel plastic-degrading organisms should prioritize environments that score highly on both diversity and characterisation gap, as these offer the best probability of finding previously uncharacterized degraders.

---

## 10. Research Gap Analysis for Individual Plastic Substrates

The research gap analysis identifies plastic types that are scientifically neglected relative to their prevalence in global plastic production or waste streams. The analysis is implemented in `research_gap_analysis` in `analysis.py` and in `plastic_specific_candidates` in `novel_discovery.py`.

A gap score is computed for each plastic type in PlasticDB:

```
gap_score = (1 / log1p(n_organisms)) * 40
          + (1 - pct_with_sequence) * 35
          + pct_extrapolated * 25
```

The three components capture different dimensions of neglect. The `(1 / log1p(n_organisms)) * 40` component rewards plastics that have few distinct degrading organisms identified, giving higher scores to less-studied substrates. The `(1 - pct_with_sequence) * 35` component rewards plastics where most entries lack molecular evidence. The `pct_extrapolated * 25` component penalizes plastics whose evidence base consists primarily of extrapolated claims (entries attributed to enzyme homology rather than direct observation), which carry lower evidentiary weight.

Plastics with no entries dated 2020 or later are identified separately as having no recent research activity. Geographic regions with fewer than five total entries are flagged as understudied.

For each of six plastic types designated as priority hard-to-degrade substrates (PE, LDPE, HDPE, PP, PS, PVC), organisms are ranked by a plastic-specific priority score:

```
priority_score = 40 * (1 - has_sequence)
              + 30 * (1 - has_enzyme)
              - 2 * max(0, year_max - last_year)
              + [20 if n_entries == 1 else 10 if n_entries < 5 else 0]
```

This score explicitly prioritizes organisms that have demonstrated degradation of a hard-to-degrade plastic without yet having deposited sequence or enzyme data (high evidence gap), and penalizes organisms whose most recent report is old, on the grounds that very old entries may reflect organisms that were not successfully cultured or maintained. For each priority plastic type, the top 10 ranked organisms are included in the discovery report.

---

## 11. Taxonomic Diversity Metrics

Taxonomic diversity within the PlasticDB dataset was assessed at the genus and species level using the Shannon diversity index H, computed as:

```
H = -sum over all taxa of (p_i * log(p_i))
```

where p_i is the proportion of entries attributable to taxon i. Separate H values were computed at the genus level (treating each genus as a category) and at the species level (treating each unique organism name as a category). A higher H value indicates more evenly distributed research effort across taxa; a lower H value indicates domination by a small number of taxa.

Additional metrics computed:

| Metric | Value | Interpretation |
|---|---|---|
| `n_unique_genera` | Count of distinct genera in PlasticDB | Breadth of genus-level coverage |
| `n_unique_species` | Count of distinct organism names | Breadth of species-level coverage |
| `singleton_genera` | Genera represented by exactly one species | Count of under-sampled genera |
| `pct_singleton_genera` | Fraction of all genera that are singletons | Degree of sampling completeness |

The plastic co-occurrence analysis produces a symmetric matrix C where C[i][j] is the number of organisms in PlasticDB that degrade both plastic type i and plastic type j. This matrix was computed using a MultiLabelBinarizer to convert each organism's plastic list into a binary vector and then taking the inner product of the resulting binary matrix with its transpose: C = B^T * B. High off-diagonal values indicate that the two corresponding plastic types tend to be co-degraded by the same organisms, suggesting shared enzymatic mechanisms. The co-occurrence matrix is stored in `outputs/reports/cross_database_comparison.json`.

---

## 12. Chart Data Pre-computation

Four sets of chart arrays are computed per organism and stored as JSON text in the `organisms` table. These arrays are delivered in the full organism detail API response, enabling the frontend to render four charts in the Research Overview section of the organism drawer without any additional API calls or client-side aggregation.

**Plastics studied chart (`ch_pl_labels`, `ch_pl_values`, `ch_pl_colors`).** For each organism, the `organism_entries` rows are grouped by plastic type and counted. The resulting label-value pairs are sorted by count descending and serialized as three parallel JSON arrays. The colors array assigns a fixed color to each plastic type based on its class: green for Bioplastic, blue for Conventional. This chart renders as a horizontal bar chart in Recharts, showing the depth of evidence for each plastic type studied.

**Publications by year chart (`ch_yr_labels`, `ch_yr_values`).** Entry counts are grouped by publication year and serialized as parallel year and count arrays. Years with zero entries are not included; the chart therefore shows only years with at least one entry. This chart renders as a vertical bar chart revealing the temporal trajectory of research on each organism, including whether activity is concentrated in a single year (possibly a single laboratory group) or distributed across years (broader community interest).

**Evidence methods chart (`ch_ev_labels`, `ch_ev_values`).** Where PlasticDB entries record the specific experimental method used to demonstrate degradation (FTIR spectroscopy, clear-zone assay, gravimetric weight loss, XRD, mass spectrometry, CO2 evolution), these method labels are counted and serialized as parallel arrays. Method diversity is informative: an organism characterized by multiple independent methods (e.g., both clear-zone and FTIR confirmation) has stronger evidentiary support than one characterized by a single method.

**Enzyme families chart (`ch_fam_labels`, `ch_fam_values`).** Where entries include enzyme family annotations (Lipase, Depolymerase, Laccase, Cutinase, and others), these are counted and serialized for rendering as a donut chart in Recharts with a legend. This chart is only populated when at least one entry for the organism names an enzyme family.

All four chart datasets are computed once during the enrichment pipeline run and stored persistently. The API returns them as pre-parsed arrays, avoiding repeated database aggregation at request time.

---

## 13. Database Schema and Field Derivations

### Table: organisms

The organisms table contains one row per unique microorganism in the atlas. In the SQLite deployment the primary key is `name` (TEXT), because the table was built organism-first from PlasticDB name strings. In the PostgreSQL RDS deployment the primary key is `taxid` (BIGINT), because the NCBI Taxonomy ID is a more stable long-term identifier that persists across nomenclature changes. In both deployments, the organism is queried by name through the API, using a case-insensitive text index.

Full column listing with types and derivation sources:

**Identity and taxonomy (from NCBI Taxonomy, Section 4.1):**
`name TEXT`, `tax_id TEXT / taxid BIGINT`, `genus TEXT`, `phylum TEXT`, `rank / org_rank TEXT`, `class / org_class TEXT`, `order / org_order TEXT`, `family / org_family TEXT`, `confidence_tier TEXT`

Note: in PostgreSQL, `class` and `order` are reserved words and must be double-quoted in all SQL statements: `"class"`, `"order"`. The SQLite schema avoids this by using the prefixed names `org_class` and `org_order`.

**PlasticDB aggregates (from Section 2):**
`n_entries INTEGER`, `n_plastics INTEGER`, `n_bioplastic INTEGER`, `n_conventional INTEGER`, `bioplastic_relevant BOOLEAN`, `has_sequence BOOLEAN`, `has_enzyme BOOLEAN`, `has_genbank BOOLEAN`, `plastics TEXT (JSON)`, `plastics_cls TEXT (JSON)`, `first_year INTEGER`, `last_year INTEGER`, `isolation_envs TEXT`, `isolation_locs TEXT`, `is_extra BOOLEAN`, `is_expanded BOOLEAN`, `is_thermo BOOLEAN`, `is_rt BOOLEAN`

**Genome (from NCBI Assembly, Section 4.2):**
`genome_size BIGINT`, `genome_level TEXT`, `genome_acc TEXT`, `genome_name TEXT`, `genome_n50 BIGINT`, `genome_cov REAL`, `genome_taxid BIGINT`

**BacDive (from Section 4.4):**
`bd_found BOOLEAN`, `bd_id INTEGER`, `bd_url TEXT`, `bd_temp TEXT`, `bd_ph TEXT`, `bd_oxy TEXT`, `bd_morph TEXT`, `bd_iso TEXT`

**Literature and SRA (from Sections 4.3 and 4.5):**
`pm_plastic INTEGER`, `pm_total INTEGER`, `sra_rc BIGINT`, `sra_plat TEXT`, `sra_strat TEXT`, `sra_bases BIGINT`, `sra_dates TEXT`

**Novelty and room-temperature scores (from Sections 5.4 and 7):**
`nov REAL`, `nov_breadth REAL`, `nov_rarity REAL`, `nov_recency REAL`, `nov_gap REAL`, `nov_last_yr INTEGER`, `nov_n_plastics INTEGER`, `rt_max REAL`, `rt_mean REAL`

**Sequence analysis and chart data (from Sections 4.6 and 12):**
`pp TEXT (JSON)`, `ch_pl_labels TEXT (JSON)`, `ch_pl_values TEXT (JSON)`, `ch_pl_colors TEXT (JSON)`, `ch_yr_labels TEXT (JSON)`, `ch_yr_values TEXT (JSON)`, `ch_ev_labels TEXT (JSON)`, `ch_ev_values TEXT (JSON)`, `ch_fam_labels TEXT (JSON)`, `ch_fam_values TEXT (JSON)`

### Table: organism_entries

One row per PlasticDB entry, preserving the full original record. Foreign key to `organisms` via the `organism` text field (organism name).

`id BIGSERIAL PRIMARY KEY`, `organism TEXT NOT NULL`, `plastic TEXT`, `cls TEXT`, `year INTEGER`, `enzyme TEXT`, `family TEXT`, `has_seq BOOLEAN`, `has_gb BOOLEAN`, `env TEXT`, `loc TEXT`, `doi TEXT`

Indexes: `organisms_tier_idx` on `confidence_tier`, `organisms_phylum_idx` on `phylum`, `organisms_genus_idx` on `genus`, `organisms_nov_idx` on `nov`, `organisms_sra_idx` on `sra_rc`, `organisms_pm_idx` on `pm_total`, `organisms_year_idx` on `first_year`, `organisms_bio_idx` on `bioplastic_relevant`, `organisms_name_trgm` as a GIN index on `to_tsvector('simple', name)` for full-text search. On the entries table: `oe_organism_idx` on `organism`, `oe_plastic_idx` on `plastic`.

---

## 14. Data API Design

The API exposes four endpoints, structured to minimize query complexity at request time by relying on the pre-computed fields described above. All response values are derived directly from the database fields documented in Section 13.

### GET /api/stats

Returns the 11 aggregate counts displayed in the stats bar. All values are computed via `COUNT(*)` or `COUNT(DISTINCT ...)` queries against the live database at request time, with results cached for 60 seconds. Each field maps directly to a database condition:

| Response field | Database query |
|---|---|
| `total_organisms` | `COUNT(*) FROM organisms` |
| `confirmed_count` | `COUNT(*) WHERE confidence_tier = 'Confirmed'` |
| `predicted_count` | `COUNT(*) WHERE confidence_tier = 'Predicted'` |
| `listed_count` | `COUNT(*) WHERE confidence_tier = 'Listed'` |
| `bioplastic_active` | `COUNT(*) WHERE bioplastic_relevant = true` |
| `genome_count` | `COUNT(*) WHERE genome_acc != ''` |
| `bacdive_count` | `COUNT(*) WHERE bd_found = true` |
| `unique_genera` | `COUNT(DISTINCT genus) WHERE genus != ''` |
| `sra_count` | `COUNT(*) WHERE sra_rc > 0` |
| `total_entries` | `COUNT(*) FROM organism_entries` |
| `unique_plastics` | `COUNT(DISTINCT plastic) FROM organism_entries WHERE plastic != ''` |

### GET /api/organisms

Returns a paginated, filtered, and sorted list of organisms. The following query parameters are accepted:

**`q` (string):** A case-insensitive substring match applied to `name`, `genus`, `phylum`, and `plastics` via `LIKE '%q%'` (SQLite) or `ILIKE '%q%'` (PostgreSQL). A 300ms client-side debounce prevents a database query on every keystroke.

**`filter` (enum):** One of seven string values, each mapping to a single database condition:

| `filter` value | SQL condition |
|---|---|
| `bioplastic` | `bioplastic_relevant = true` |
| `conventional` | `n_conventional > 0` |
| `genome` | `genome_acc != ''` |
| `bacdive` | `bd_found = true` |
| `sra` | `sra_rc > 0` |
| `thermo` | `is_thermo = true` |
| `rt` | `is_rt = true` |

**`tier` (enum):** `confirmed`, `predicted`, or `listed`, mapped to the exact string stored in `confidence_tier`.

**`sort` (enum):** One of six values, each mapping to an ORDER BY expression:

| `sort` value | ORDER BY expression |
|---|---|
| `name` | `name ASC` |
| `novelty` | `nov DESC, name ASC` |
| `sra` | `sra_rc DESC, name ASC` |
| `pubmed` | `pm_total DESC, name ASC` |
| `entries` | `n_entries DESC, name ASC` |
| `year` | `first_year DESC, name ASC` |

Invalid sort values fall back to `name ASC`. All secondary sort keys are `name ASC` to ensure deterministic ordering within tied values.

**Pagination:** `page` (1-based integer, default 1) and `per_page` (integer 1-200, default 50). The response includes `total` (matching rows), `page`, `per_page`, `pages` (ceiling division of total by per_page), and `organisms` (the result array). Both a data query and a count query are executed per request. The count query uses the same WHERE clause as the data query but omits the ORDER BY and LIMIT clauses.

### GET /api/organisms/phylum

Returns confirmed-organism counts grouped by phylum. Only Confirmed organisms are counted, as Listed organisms would dominate the distribution and obscure the biologically meaningful pattern. Empty phylum values are coalesced to "Unclassified". Results are ordered by count descending with no pagination cap, as there are fewer than 200 distinct phyla in the dataset.

### GET /api/organisms/{name}

Returns the full record for one organism identified by name. Name matching is attempted first as an exact match, then as a case-insensitive match. A 404 response is returned if no match is found. The full record includes all 69 fields from the `organisms` table plus an `entries` array containing all rows from `organism_entries` for this organism, ordered by year descending. JSON-typed columns (`pp`, `ch_pl_labels`, `ch_pl_values`, `ch_pl_colors`, `ch_yr_labels`, `ch_yr_values`, `ch_ev_labels`, `ch_ev_values`, `ch_fam_labels`, `ch_fam_values`, `plastics`, `plastics_cls`) are deserialized to native arrays before being returned; malformed or empty JSON values are returned as empty arrays.

---

## 15. Organism Atlas Interface

The Organism Atlas is the primary user-facing interface for PETadex. It is built in React and designed around the principle that all view state -- search query, active filter, active tier, sort selection, current page, and selected organism drawer -- is encoded in the URL query string. Any specific view of the atlas is therefore bookmarkable, shareable by URL, and reproducible.

Data fetching uses React Query with the following caching policy: organism list responses are cached for 15 seconds with a placeholder data strategy that displays the previous page during navigation transitions, eliminating blank states between page loads. Stats are cached for 60 seconds. Phylum breakdown data is cached for 5 minutes. The organism detail is cached for 30 seconds per name.

### 15.1 Stats Bar

The stats bar calls `GET /api/stats` on component mount and displays 8 of the 11 returned fields as labeled value boxes:

| Label displayed | API field shown |
|---|---|
| Total Organisms | `total_organisms` |
| Confirmed | `confirmed_count` |
| Predicted | `predicted_count` |
| Listed | `listed_count` |
| Genera | `unique_genera` |
| Genome Assemblies | `genome_count` |
| BacDive Records | `bacdive_count` |
| Unique Plastics | `unique_plastics` |

The fields `sra_count`, `bioplastic_active`, and `total_entries` are returned by the API but used internally (in filter counts and totals) rather than displayed in the stats bar. During the initial load before the API response arrives, boxes display animated skeleton placeholders.

### 15.2 Search, Filter Pills, Tier Pills, and Sort

All four controls send their values as URL query parameters and trigger a fresh `GET /api/organisms` request with the updated parameters. Multiple controls combine with logical AND: a user may simultaneously search by name substring, filter to genome-containing organisms, restrict to the Confirmed tier, and sort by novelty score.

The five filter pills map to the `filter` parameter and are mutually exclusive; selecting one deselects any previously active filter. The four tier pills map to the `tier` parameter and are similarly mutually exclusive. The sort dropdown maps to the `sort` parameter.

### 15.3 Organism Table

The 11 table columns and their rendering:

| Column | Source fields | Rendering detail |
|---|---|---|
| Organism | `name`, `confidence_tier` | Italic organism name with a tier badge below: green for Confirmed, amber for Predicted, grey for Listed |
| Genus | `genus` | Plain text |
| Phylum | `phylum` | Plain text |
| Plastics | `plastics`, `plastics_cls` | Up to three pill badges per row. Each pill shows one plastic type name. Bioplastic class pills render in green, Conventional in blue. When more than three types exist, the overflow count is shown as a grey "+N" badge |
| #Types | `n_entries` | Integer count of all PlasticDB entries; not the count of distinct plastic types, which is `n_plastics` |
| 1st Year | `first_year` | Four-digit year; "---" when null |
| SRA | `sra_rc` | Integer run count; "---" when zero |
| PubMed | `pm_total` | Integer citation count |
| Genome | `genome_acc`, `genome_level` | Assembly accession in a monospace chip; empty cell when `genome_acc` is blank |
| BacDive | `bd_found` | "Yes" in green when `bd_found = true`; "No" in grey otherwise |
| Novelty | `nov` | Integer score 0-100, color-coded green at or above 70, amber 40-69, red below 40 |

During data fetching, 12 skeleton rows are displayed. If the API returns an error, a full-width error message is shown. If the result set is empty after filtering, an empty-state message describes the active filters.

### 15.4 Organism Detail Drawer

Clicking any table row opens the organism detail drawer, which calls `GET /api/organisms/{name}` and renders the full record in nine sequential sections. The `org={name}` parameter is added to the URL, making the open drawer state bookmarkable. The drawer closes by removing the `org` parameter.

**Header.** Organism name in italic. Confidence tier badge. Breadcrumb showing genus, phylum, and NCBI Taxonomy ID.

**Summary.** Bioplastic type count (`n_bio`) and conventional plastic type count (`n_conv`) as large numeric values. Three binary badges for Sequence (`has_seq`), Enzyme (`has_enz`), and GenBank (`has_gb`), each rendered in green when true and grey when false. The `iso_envs` and `iso_locs` fields are displayed as tag lists. The full `plastics` array is shown as colored pills.

**Genome.** `g_size` formatted in megabases (dividing raw base-pair count by 1,000,000). `g_level` in large text. `g_acc` as a hyperlink to `https://www.ncbi.nlm.nih.gov/datasets/genome/{g_acc}/`. Assembly name `g_name`, N50 `g_n50`, coverage `g_cov`, and source TaxID `g_taxid` in a four-cell grid. This section renders only when `g_acc` is non-empty.

**BacDive Physiology.** Renders when `bd_found = true`. Displays `bd_id` as a link to `bd_url`, then `bd_temp`, `bd_ph`, `bd_oxy`, `bd_morph`, and `bd_iso` as labeled rows. When `bd_found = false`, a short message states that no BacDive record was matched to this organism.

**Research Overview.** Four Recharts components, each consuming one of the four pre-computed chart array pairs:
1. Plastics Studied: horizontal BarChart using `ch_pl_labels` and `ch_pl_values`, bars filled by `ch_pl_colors`
2. Publications by Year: vertical BarChart using `ch_yr_labels` and `ch_yr_values`
3. Evidence Methods: vertical BarChart using `ch_ev_labels` and `ch_ev_values`
4. Enzyme Families: PieChart in donut mode using `ch_fam_labels` and `ch_fam_values` with a legend

Charts are rendered only when the corresponding label array is non-empty. Empty arrays display a "no data" placeholder.

**Novelty Score.** The overall `nov` score as a large number with color coding. Four horizontal progress bars for `nov_breadth`, `nov_rarity`, `nov_recency`, and `nov_gap`, each labeled with the sub-score name and numeric value. This section makes the novelty score fully decomposable: a researcher can see which of the four dimensions drives a given organism's ranking and interpret the score in terms of the specific properties that make the organism a priority.

**PlasticDB Entries.** A table of every row from `organism_entries` for this organism, ordered by year descending. Columns: Plastic, Year, Enzyme, Family, Sequence (boolean badge), GenBank (boolean badge), Environment, Location, DOI (hyperlink). This section provides complete traceability: every aggregate statistic in the Summary section is derivable from this raw entry list.

**Temperature Profile.** Three indicators derived from `is_thermo`, `is_rt`, `rt_max`, and `rt_mean`. The BacDive growth temperature string `bd_temp` is displayed separately as the experimentally recorded range, annotated to clarify that it describes general growth conditions, not plastic-degradation conditions.

**PubMed.** `pm_total` and `pm_plastic` counts as labeled values. A link opens a PubMed search URL pre-formed with the organism name.

**ProtParam.** Rendered when the `pp` JSON array is non-empty. For each sequence record in `pp`, a table row displays: length, molecular weight, isoelectric point, instability index, GRAVY score, and stability prediction.

**SRA.** `sra_rc` (run count), `sra_bases` formatted in gigabases, `sra_dates`, `sra_plat`, and `sra_strat` as labeled values. A link opens an SRA search pre-formed with the organism's TaxID.

### 15.5 Phylum Breakdown Chart

A collapsible section below the table that calls `GET /api/organisms/phylum` on first expansion. It renders the top 16 phyla among Confirmed organisms as horizontal bars ordered by count descending. The section header displays the total number of distinct phyla returned by the endpoint. This chart is restricted to Confirmed organisms to prevent the approximately 2.9 million Listed organisms from rendering a distribution that reflects NCBI taxonomy coverage rather than plastic degradation research activity.

---

## 16. Infrastructure and Data Migration

**Development environment.** The deployment stores all data in a SQLite file (`petadex-api/petadex.db`) served by a Python FastAPI process. SQLite was selected for development because it requires no separate server process and because the entire database fits within a standard disk allocation. Boolean fields are stored as INTEGER (0 or 1) in SQLite, converted to Python bool in the API response layer. JSON fields are stored as TEXT and deserialized in the API response layer using a helper that returns an empty list for null or malformed values.

**Production environment.** The `petadex/petadex.io` deployment uses an Amazon RDS PostgreSQL 16.3 instance (`petadex.c6dcs4m8a2uy.us-east-1.rds.amazonaws.com`, us-east-1b, port 5432). PostgreSQL was selected for production because it provides native BOOLEAN types, a GIN trigram index for fast text search, concurrent connection management, and managed backups and failover. The backend runs on AWS Lambda via the `serverless-http` adapter, with the Lambda function placed inside the VPC (security group `sg-05566dc054d55ac83`, subnets `subnet-0a9f6e8de6432e5d6` and `subnet-0f9c86f70806c421a`) that contains the RDS instance. The connection pool is configured with a maximum of 2 concurrent connections to remain within Lambda's concurrency budget.

**Schema adaptations from SQLite to PostgreSQL.** The `org_rank`, `org_class`, `org_order`, and `org_family` column names used in SQLite were renamed to `rank`, `class`, `order`, and `family` in PostgreSQL. The latter three are PostgreSQL reserved words and must be double-quoted in all SQL statements. Integer boolean columns (0/1 in SQLite) become native BOOLEAN columns in PostgreSQL. The `pp` column is typed TEXT in both schemas because it stores a JSON array; in SQLite a REAL type was specified in an early draft of the migration script and corrected before execution, since JSON arrays cannot be cast to REAL by PostgreSQL's `COPY` parser.

**Data loading.** The 2,902,229-row `organisms` table and the 2,535-row `organism_entries` table were exported from the SQLite database as CSV files using Python's `csv` module. CSV files were transferred to an AWS CloudShell session and loaded into RDS via the PostgreSQL `\copy` command with a direct `psql` connection from within the same VPC. The 295 MB `organisms_full.csv` file loaded in approximately 2 to 5 minutes. All indexes were created after the bulk load to avoid the index maintenance overhead during insertion.

**Credentials and access control.** The five RDS connection parameters (`DB_HOST`, `DB_PORT`, `DB_NAME`, `DB_USER`, `DB_PASS`) are injected into the Lambda function as environment variables via the `environment` block in `serverless.yml`. They are not stored in code or in the repository. The Lambda IAM role has no direct RDS permissions; access is controlled by VPC security group rules. The SSL connection option `rejectUnauthorized: false` is set in the Node.js `pg` pool to accommodate the RDS certificate without requiring a certificate bundle to be packaged with the Lambda deployment artifact.
