# Physiological Characterization and Genomic Sequencing Depth of Plastic-Degrading Microorganisms: An Analysis of BacDive and NCBI SRA Data in the Context of Bioplastic Enzyme Discovery

**Document type:** Research Analysis
**Subject:** BacDive physiological metadata and NCBI Sequence Read Archive data for the 932 confirmed and predicted plastic-degrading organisms in PETadex, and their collective implications for room-temperature bioplastic enzyme discovery.

---

## Abstract

Two secondary databases, the BacDive strain information system and the NCBI Sequence Read Archive, were queried systematically for all 932 confirmed and predicted plastic-degrading organisms in the PETadex atlas. BacDive returned physiological records for 474 organisms (50.9%), revealing a striking mesophilic dominance: 97.0% of organisms with temperature data grow optimally between 15 and 45 degrees Celsius, with a mean recorded growth temperature of 29.5 degrees Celsius across all matched records. Aerobic and obligately aerobic organisms account for 80.3% of all organisms with oxygen requirement data. Terrestrial soil environments are the single most represented isolation habitat, appearing in 124 of 207 environmentally sourced BacDive records. The NCBI Sequence Read Archive contains data for 690 of the 932 queried organisms (74.0%), encompassing 1,295,877 runs and 8.41 terabases of sequence data. The distribution of SRA runs is highly right-skewed, with a median of 19 runs per organism against a mean of 1,878, driven by a small number of clinically significant organisms that happen to also degrade bioplastics. Whole-genome sequencing constitutes the plurality of library strategy (542 organism-strategy pairs), and Illumina short-read sequencing accounts for 618 of all platform observations. Oxford Nanopore and PacBio long-read data are present for 150 and 115 organisms respectively, marking a measurable shift toward complete, closed genome assemblies within the plastic-degrading organism set. Cross-referencing BacDive and SRA coverage identifies 395 organisms with both physiological and genomic data, 79 with physiological data but no public sequencing, 295 with sequencing but no physiological record, and 163 with neither. This last group represents the highest-priority targets for experimental investment: organisms with confirmed plastic degradation capacity but essentially no secondary data to guide enzyme discovery. These findings, placed in the context of the novelty scoring model and the temperature analysis of characterized PETases, define the organism prioritization framework for PetaBite's search for room-temperature plastic-degrading enzymes.

---

## 1. Introduction

The identification of new plastic-degrading enzymes depends on understanding not only which organisms have been shown to degrade plastic substrates, but also the physiological conditions under which those organisms thrive and the depth of existing genomic data available to support enzyme discovery. Two publicly accessible databases address these needs from complementary angles.

BacDive, maintained by the Leibniz Institute DSMZ, is a structured repository of microbial physiology. It records experimentally determined growth parameters, including temperature optima and ranges, pH tolerance, oxygen requirements, and cell morphology, alongside isolation source and culture conditions. For any organism in this database, BacDive provides a rapid physiological profile without requiring primary literature access. In the context of plastic degradation research, BacDive data are particularly relevant because the conditions under which an enzyme functions in vivo are constrained by the physiology of its host organism. An aerobic mesophile isolated from soil is far more likely to produce an enzyme active at ambient temperature in an aerobic wastewater environment than an obligate anaerobe from a deep-sea thermal vent.

The NCBI Sequence Read Archive is the world's largest repository of high-throughput sequencing data, accumulating primary sequence output from published studies regardless of organism, tissue, or methodology. For microorganisms, SRA entries most commonly represent whole-genome sequencing, amplicon-based community surveys, or transcriptomic runs. The volume and type of SRA data associated with an organism reflects the cumulative investment of the research community in that organism's genomics. High SRA run counts signal a well-resourced area of study; low counts, for an organism with confirmed plastic degradation capacity, indicate that the genomic resources needed to identify and characterize plastic-degrading enzymes have not yet been generated.

This paper reports the results of a systematic query of both databases for the 932 confirmed and predicted plastic-degrading organisms in PETadex, analyzes the distributions and patterns in the retrieved data, and interprets those findings in the specific context of room-temperature bioplastic enzyme discovery.

---

## 2. Methods

### 2.1 Organism Cohort

The analysis cohort consisted of the 932 organisms classified as either Confirmed or Predicted in the PETadex confidence tier system. Confirmed organisms are those appearing by name in at least one peer-reviewed PlasticDB entry recording experimental evidence of plastic degradation. Predicted organisms were identified by sequence homology to characterized plastic-degrading enzymes in the PAZy database. Listed organisms (those belonging to genera of known degraders but lacking direct evidence) were excluded from the BacDive and SRA queries to keep the analysis focused on organisms with a defensible experimental or computational basis for plastic degradation capacity.

### 2.2 BacDive Query

Each organism name was submitted to the BacDive REST API using the search endpoint at `/search?search={organism_name}`. Exact name matches were preferred; where the API returned multiple candidates, the record with the highest name similarity score was selected. The following fields were extracted for each matched record: strain ID, record URL, growth temperature (as free text), growth pH range, oxygen requirement category, cell morphology, and isolation source using BacDive's hierarchical environment classification.

### 2.3 SRA Query

Each organism's NCBI Taxonomy ID was used to query the SRA database via the Entrez API. All public runs associated with the TaxID were retrieved and summarized to produce: total run count, sequencing platforms, library strategies, total sequenced bases, and submission date range. Organisms lacking an NCBI TaxID (due to failed name resolution) were excluded from SRA queries. The analysis was conducted on all 932 organisms for which TaxID resolution succeeded.

### 2.4 Temperature Parsing

BacDive temperature data are stored as free-text strings because BacDive does not enforce a numeric format. Temperature values were extracted from these strings using regular expression matching to identify all numeric substrings and filtering to values in the physiologically plausible range of 4 to 80 degrees Celsius. For organisms with multiple temperature values in their record (reflecting a range rather than a single optimum), the mean of all extracted values was used for binning. Three bins were applied: psychrophilic for mean temperatures below 15 degrees Celsius, mesophilic for 15 to 45 degrees Celsius, and thermophilic for values above 45 degrees Celsius. These thresholds follow the standard microbiological classification in use in the plastic degradation literature.

---

## 3. BacDive: Physiological Profiles of Plastic-Degrading Organisms

### 3.1 Database Coverage

Of the 932 organisms queried, 474 (50.9%) returned a BacDive record. The remaining 458 organisms (49.1%) had no match. The 50.9% coverage rate is meaningful: it implies that roughly half of all experimentally confirmed or predicted plastic-degrading organisms have been cultured and physiologically characterized to the degree necessary for formal database entry. The complementary implication is that the other half exist primarily as name records in the literature, without the culture-based physiology data that would allow confident prediction of the conditions under which their enzymes are active.

The organisms without a BacDive match are not necessarily less important as research targets. Their absence from BacDive may reflect either that they are poorly cultured (environmental isolates that have not been maintained in axenic culture), that they are taxonomically heterogeneous categories (genus-level or sp. designations that cannot be matched to a single type strain), or simply that they were described before BacDive's curation of their literature. For enzyme discovery purposes, the absence of a BacDive record shifts the burden of physiological inference to the PlasticDB experimental conditions and to genomic proximity to better-characterized relatives.

### 3.2 Growth Temperature Distribution

Temperature data were extracted for 435 of the 474 BacDive-matched organisms, representing 91.8% of matched records. The full temperature range across all extracted values spans 4.0 to 80.0 degrees Celsius, reflecting the full breadth of environments from which plastic-degrading organisms have been isolated and characterized.

The distribution is strongly mesophilic. Applying the three-bin classification:

| Temperature class | Organisms | Proportion |
|---|---|---|
| Psychrophilic (below 15 degrees Celsius) | 1 | 0.2% |
| Mesophilic (15 to 45 degrees Celsius) | 421 | 97.0% |
| Thermophilic (above 45 degrees Celsius) | 13 | 3.0% |
| Total with temperature data | 435 | 100% |

The mean growth temperature across all 435 organisms with data is 29.5 degrees Celsius. This single value carries significant weight for the PetaBite research program. The current state-of-the-art ambient-temperature PETase benchmark, FAST-PETase, was engineered to operate at 37 degrees Celsius, and the industrial case for room-temperature PET hydrolysis in wastewater treatment plants centers on organisms active between approximately 15 and 25 degrees Celsius. The finding that the mean growth temperature of the 435 characterized plastic-degrading organisms lies at 29.5 degrees Celsius, squarely within the mesophilic range and close to the ambient conditions of temperate-climate wastewater infrastructure, constitutes strong prior evidence that room-temperature plastic-degrading enzymes are not rare exceptions but rather the statistical norm within the existing characterized organism pool.

The 13 thermophilic organisms in the BacDive dataset represent 3.0% of those with temperature data, a figure that is lower than the 7.6% thermophilic fraction observed in the raw PlasticDB entry conditions (Section 5.1 of the main documentation). This discrepancy is expected: PlasticDB thermophile classifications reflect the temperature conditions of individual experiments, whereas BacDive thermophile classification reflects the organism's general growth physiology. Some organisms isolated from mesophilic environments may have been tested under thermophilic conditions in individual experiments without being thermophilic organisms in the BacDive sense. The lower rate observed in BacDive is therefore considered the more reliable estimate of true thermophile prevalence among plastic-degrading organisms.

The single psychrophilic organism in the dataset is notable as a potential source of room-temperature enzymes, as cold-adapted enzymes from psychrophiles are known to exhibit high catalytic rates at low temperatures as a structural compensation for reduced thermal energy, but the sample size is too small for generalization.

### 3.3 Oxygen Requirements

Of the 474 BacDive-matched organisms, 360 had an oxygen requirement recorded. The distribution is strongly aerobic:

| Oxygen category | Count | Proportion |
|---|---|---|
| Aerobe | 159 | 44.2% |
| Obligate aerobe | 130 | 36.1% |
| Facultative anaerobe | 52 | 14.4% |
| Anaerobe | 16 | 4.4% |
| Microaerophile | 2 | 0.6% |
| Obligate anaerobe | 1 | 0.3% |

Combining the aerobe and obligate aerobe categories, 289 organisms (80.3%) require or strongly prefer oxygen. Facultative anaerobes, which can grow with or without oxygen, constitute 14.4%. Strict anaerobes and microaerophiles together account for 5.3%.

The aerobic dominance has direct implications for bioplastic remediation system design. Wastewater treatment plants operate aerobic biological treatment stages as their primary treatment mechanism, and aerobic organisms are the natural occupants of composting facilities, soil bioreactors, and open-air treatment lagoons. The fact that 80.3% of characterized plastic-degrading organisms with oxygen data are aerobic or obligately aerobic means that conventional aerobic treatment infrastructure is broadly compatible with the metabolism of the organisms most likely to harbor useful plastic-degrading enzymes. It also means that enzyme assay conditions in any screening program based on these organisms should default to aerobic incubation.

The 14.4% facultative anaerobe fraction is significant for two reasons. First, facultative anaerobes are among the most environmentally resilient microorganisms, surviving in both aerobic and anoxic zones of treatment systems. Second, several of the highest-novelty organisms in PETadex are facultative anaerobes, including members of the genera Serratia and Enterobacter, which means their enzymes would need to function under aerobic conditions (as extracellular hydrolases acting at a polymer surface) even if the organism itself is capable of anaerobic growth.

### 3.4 Isolation Sources and Ecological Context

BacDive records isolation source information using a hierarchical classification system. Among the 474 matched organisms, the following top-level isolation categories were identified:

| Isolation category | Occurrences |
|---|---|
| Environmental | 207 |
| Terrestrial | 151 |
| Soil | 124 |
| Host (plant or animal) | 98 |
| Engineered environment | 88 |
| Aquatic | 61 |
| Marine | 25 |

Note that these categories are not mutually exclusive; BacDive uses hierarchical nested tags, so an organism from soil is tagged as Environmental, Terrestrial, and Soil simultaneously. The counts therefore reflect co-occurrence of tags rather than distinct organisms.

Soil is the single most specific environment appearing in 124 records, representing the plurality of the environmental isolation records. Terrestrial non-soil environments (forest litter, plant surfaces, decomposing organic matter) account for a further 27 records when the Terrestrial total is adjusted for the Soil subset. The 98 host-associated records cover organisms isolated from plant tissues, bovine rumen, and other animal hosts, environments that contain biopolymers structurally analogous to certain plastics and that represent a historically productive source of esterases and depolymerases.

The 88 engineered environment records are particularly relevant for wastewater treatment applications. Engineered environments in BacDive include industrial fermenters, wastewater treatment systems, activated sludge, and composting facilities. Organisms isolated from these environments have already demonstrated growth and activity under the conditions most relevant to large-scale plastic remediation. Their appearance in PlasticDB with confirmed plastic degradation records suggests that plastic degradation has already been selected for, at least partially, in engineered microbial communities.

Aquatic environments (61 records) include freshwater and sediment isolates. Marine organisms (25 records) are particularly relevant for ocean plastic remediation research, though they are underrepresented in PETadex relative to their likely ecological importance given the concentration of persistent synthetic plastics in marine systems.

### 3.5 Cell Morphology and Taxonomic Implications

Among the 474 BacDive-matched organisms, morphological data were present for 235. The dominant form observed was the Gram-negative rod (71 organisms, 30.2%), followed by Gram-positive organisms without further morphological specification (64, 27.2%), Gram-negative organisms without further specification (53, 22.6%), and Gram-positive rods (40, 17.0%). Coccoid forms were present in 7 records (3.0%).

The dominance of rod-shaped, Gram-negative bacteria among characterized plastic degraders reflects the broader taxonomic composition of the confirmed organism set, which is heavily weighted toward Proteobacteria genera such as Pseudomonas, Ralstonia, Ideonella, Serratia, and Acinetobacter. Gram-negative bacteria are generally characterized by an outer membrane that serves as a physical interface with the extracellular environment, including solid substrates such as plastic surfaces. The structural biology of this membrane interface may be relevant to the mechanism by which these organisms achieve contact with hydrophobic polymer surfaces, a recognized rate-limiting step in enzymatic plastic degradation.

The Gram-positive organisms in the set include Bacillus, Streptomyces, Arthrobacter, and Rhodococcus representatives, genera known for robust biosurfactant production and diverse secondary metabolite repertoires. Streptomyces sp. and Bacillus sp. collectively represent some of the broadest plastic substrate ranges in PETadex (16 and 20 plastic types respectively), suggesting that the metabolic versatility characteristic of these genera extends to the range of plastic polymers they can process.

### 3.6 The BacDive Coverage Gap and Its Research Implications

The 458 organisms without a BacDive match fall into two functionally distinct groups. The first group consists of named species that exist in the literature primarily from single studies, have not been deposited in culture collections, and have therefore never been formally strain-characterized for BacDive inclusion. Many of these represent precisely the high-novelty, low-characterization organisms that the PETadex novelty scoring model identifies as discovery priorities. Their absence from BacDive is not a sign of low importance; it is a sign that their characterization pipeline terminated at the publication stage without the follow-on culture and phenotyping work that BacDive requires.

The second group consists of genus-level or species-group designations (entries of the form "Pseudomonas sp." or "uncultured bacterium") that cannot be mapped to a single BacDive strain record, even when members of the genus are well-represented in BacDive. This is a structural limitation of name-based matching rather than a gap in physiological knowledge per se. For these organisms, the BacDive records of closest relatives within the genus provide a reasonable physiological approximation for the purposes of environmental and process design.

---

## 4. NCBI Sequence Read Archive: Genomic Data Depth for Plastic-Degrading Organisms

### 4.1 Coverage and Run Count Distribution

Of the 932 queried organisms, 690 (74.0%) had at least one public SRA run at the time of data collection. The remaining 242 organisms (26.0%) had no public sequencing data in SRA, implying that their genomic content is either available only in finished genome assemblies (NCBI Assembly rather than SRA), published in non-NCBI archives, or not yet sequenced and deposited.

The 690 organisms with SRA data collectively account for 1,295,877 public sequencing runs and 8.41 terabases of sequence data. The distribution of runs across organisms is highly right-skewed:

| Metric | Value |
|---|---|
| Total organisms with SRA runs | 690 |
| Total runs across all organisms | 1,295,877 |
| Total sequenced bases | 8.41 TB |
| Median runs per organism | 19 |
| Mean runs per organism | 1,878 |
| Maximum runs (single organism) | 613,710 |

The extreme divergence between the median (19) and the mean (1,878) reflects the presence of a small number of organisms with exceptional SRA coverage. The top three by run count are Escherichia coli (613,710 runs), Klebsiella pneumoniae (151,076 runs, represented twice in the dataset under two name variants), and Pseudomonas aeruginosa (93,879 runs). These organisms appear in PETadex because they have documented plastic degradation activity, predominantly involving bioplastics and polyurethane components, but their dominant SRA presence is driven by their status as major human pathogens subject to global genomic surveillance programs. Their inclusion in the plastic-degrading organism set is scientifically valid but creates a significant distributional artifact: after removing the top three organisms by run count, the mean across the remaining 687 organisms falls to approximately 659 runs, a value still elevated relative to the median but more representative of the sequencing depth available for organisms studied principally in the context of plastic degradation.

The more ecologically representative summary comes from the median: the typical plastic-degrading organism in PETadex has 19 public SRA runs. At median coverage depths for modern SRA submissions (100-fold genome coverage per WGS run), 19 runs represents robust but not exhaustive genomic characterization. For the purpose of enzyme discovery, 19 runs is generally sufficient to produce a high-quality draft genome assembly from which enzyme candidates can be identified by homology to characterized enzyme families.

### 4.2 Sequencing Platform Distribution

Illumina short-read sequencing is the dominant platform, appearing in 618 of all platform-organism associations. Oxford Nanopore long-read sequencing follows with 150 associations, and PacBio SMRT long-read sequencing with 115. DNBSEQ (BGI Genomics) accounts for 92 associations, reflecting the platform's adoption primarily in Asian research institutions. Older platforms (capillary sequencing, 454 pyrosequencing, SOLiD) are present in a small number of records reflecting legacy data.

The emergence of Oxford Nanopore and PacBio in the dataset is significant. These long-read platforms generate sequence reads of thousands to hundreds of thousands of base pairs, enabling de novo assembly of complete, circularized bacterial genomes without the fragmentation artifacts characteristic of short-read assembly. For enzyme discovery, complete genome assemblies are preferable to draft assemblies because they resolve repetitive regions that may flank or contain enzyme-encoding genes, and they enable accurate prediction of genomic context (neighboring genes, regulatory elements, operon structure) that is informative for understanding enzyme expression and secretion. The fact that 150 organisms have Oxford Nanopore data and 115 have PacBio data indicates that a meaningful subset of plastic-degrading organisms has already received the long-read sequencing treatment that produces the most tractable genomic resource for enzyme mining.

The co-occurrence of multiple platforms for a single organism is common: 89 organisms have both Illumina and Oxford Nanopore data, and 54 have both Illumina and PacBio, reflecting the common hybrid assembly approach in which long reads provide structural scaffolding and short reads provide base-level accuracy.

### 4.3 Library Strategy Distribution

Whole-genome sequencing (WGS) is the most prevalent library strategy, appearing in 542 organism-strategy associations. Amplicon sequencing, used for 16S rRNA gene surveys and environmental community profiling, accounts for 185 associations. RNA-Seq (transcriptomics) appears in 172 associations. Other strategies, including whole-genome amplification (WGA, 33), chromosome conformation capture (Hi-C, 8), and transposon sequencing (Tn-Seq, 4), account for the remainder.

The WGS dominance is favorable for enzyme discovery. A WGS library, submitted to SRA and assembled, provides the complete genomic sequence from which all protein-coding genes can be annotated and any with homology to plastic-degrading enzyme families can be identified. The 172 RNA-Seq associations are valuable for a different purpose: transcriptomic data from organisms grown in the presence of plastic substrates can identify which genes are upregulated in response to the polymer, providing evidence that specific genes are functionally relevant to plastic degradation rather than merely homologous to known enzymes.

The 185 amplicon-sequencing associations predominantly reflect 16S rRNA gene surveys in which the plastic-degrading organism was identified as part of an environmental community study rather than as a pure culture isolate. These records confirm the organism's presence in a given environment but provide no sequence data suitable for enzyme identification. They are therefore valuable for ecological mapping but not directly for the enzyme discovery pipeline.

### 4.4 Organisms Without SRA Data

The 242 organisms with no SRA runs constitute 26.0% of the queried cohort. These organisms fall into two categories. The first consists of organisms for which a genome sequence exists in NCBI Assembly but for which the raw sequencing reads were not deposited in SRA. Some older genome assemblies were submitted before mandatory raw data deposition became standard practice, and their sequences are accessible only as finished assemblies. These organisms do have exploitable genomic data, but it is not captured by an SRA run count query.

The second category consists of organisms for which no publicly accessible sequence data of any type exists. These organisms are known only from their PlasticDB experimental records: they have been shown to degrade plastic in a laboratory setting, a finding sufficient for inclusion in PETadex, but no genomic sequence has been deposited. For these organisms, the pathway from experimental observation to enzyme identification requires de novo genome sequencing, which is the most resource-intensive step in the discovery pipeline. They represent a priority target for sequencing investment, particularly the subset with high novelty scores and broad plastic substrate ranges.

### 4.5 Total Data Volume and Its Implications

The 8.41 terabases of sequencing data associated with plastic-degrading organisms in PETadex represents a substantial existing resource for enzyme discovery. However, the volume is very unevenly distributed. The top three organisms (Escherichia coli, Klebsiella pneumoniae, Pseudomonas aeruginosa) alone account for an estimated 90 to 95% of this total, given their dominant run counts and the high coverage depths typical of clinical pathogen surveillance sequencing. After excluding these three, the remaining 687 organisms collectively represent roughly 300 to 500 gigabases of sequencing data, which, at typical read lengths and depths for environmental bacterial sequencing, corresponds to approximately 300 to 600 complete or near-complete genome equivalents.

This volume is computationally tractable for systematic enzyme mining. A genome-scale homology search across 600 genomes using characterized plastic-degrading enzyme families as queries represents a manageable bioinformatics task and constitutes the most immediate practical application of the SRA data collected in PETadex.

---

## 5. Cross-Database Integration: Organisms with Both Physiological and Genomic Profiles

### 5.1 Coverage Overlap Analysis

The intersection of BacDive and SRA coverage identifies four distinct groups among the 932 queried organisms:

| Coverage category | Count | Proportion |
|---|---|---|
| BacDive record AND SRA runs | 395 | 42.4% |
| BacDive record only (no SRA runs) | 79 | 8.5% |
| SRA runs only (no BacDive record) | 295 | 31.6% |
| Neither BacDive nor SRA | 163 | 17.5% |

The 395 organisms with both BacDive physiology data and SRA sequencing data represent the most completely characterized subset of the plastic-degrading organism pool. For these organisms, there exists a direct connection between growth physiology (temperature, oxygen, pH tolerance) and the genomic sequence from which enzymes can be identified and recombinantly expressed. This connection is essential for rational enzyme selection: an enzyme from an organism that grows optimally at 28 degrees Celsius in aerobic soil conditions is more likely to be functional under wastewater treatment conditions than an enzyme from a thermophilic deep-sea organism, regardless of their sequence homology to characterized enzymes.

### 5.2 The Neither-BacDive-nor-SRA Group

The 163 organisms with neither a BacDive record nor any SRA sequencing data represent the most information-poor entries in the Confirmed tier. These organisms appear in PETadex because they are cited in PlasticDB experimental records with direct evidence of plastic degradation, but no physiological characterization has been deposited in BacDive and no sequencing data has been made publicly available. Their only retrievable metadata consists of what is recorded in the PlasticDB entry itself: the plastic substrate tested, the publication year, the isolation environment and location, and whether a sequence was reported (with or without deposition).

This group is not necessarily low-priority. Some of the highest-novelty organisms in PETadex fall into this category, including Actinomadura sp., which carries a novelty score of 52.0 with activity across 6 plastic types and zero SRA runs. Actinomadura is a genus of filamentous actinobacteria known for exceptional secondary metabolite production including polyketide enzymes and esterases, with no representatives in SRA and no BacDive record for this specific designation. The absence of data is itself informative: it signals that this organism has been studied in the specific context of plastic degradation but has not entered the broader genomics mainstream. Sequencing one representative member of Actinomadura sp. with confirmed plastic degradation capacity would constitute a genuine contribution to the genomic record.

### 5.3 Novelty Score Distribution Across Coverage Categories

The novelty scoring model (described in the main documentation) was designed independently of BacDive and SRA coverage status, using only PlasticDB-derived features. Nevertheless, the relationship between novelty score and coverage category reveals a consistent pattern: organisms with the highest novelty scores tend to have low SRA run counts despite broad plastic substrate ranges.

The top 10 organisms by novelty score and their SRA run counts:

| Organism | Novelty score | Plastic types | SRA runs |
|---|---|---|---|
| Phanerochaete chrysosporium | 68.3 | 11 | 238 |
| Bacillus sp. | 65.1 | 20 | 1,052 |
| Serratia marcescens | 63.3 | 7 | 6,719 |
| Stenotrophomonas sp. | 62.9 | 7 | 128 |
| Enterobacter sp. | 62.8 | 7 | 135 |
| Arthrobacter sp. | 58.1 | 7 | 302 |
| Penicillium sp. | 57.5 | 7 | 184 |
| Streptomyces sp. | 56.0 | 16 | 188 |
| Aspergillus niger | 54.4 | 12 | 3,101 |
| Fusarium sp. | 53.8 | 7 | 213 |

Seven of the top ten high-novelty organisms have fewer than 400 SRA runs each, despite demonstrating plastic degradation activity across 7 to 20 distinct polymer types. Streptomyces sp. is the most striking case: 16 plastic types spanning bioplastics and conventional polymers, a confirmed sequence deposit, and only 188 SRA runs. This combination of broad substrate range, molecular evidence, and limited public genomics represents the most precisely defined target profile for enzyme discovery investment.

Bacillus sp. with 20 plastic types and 1,052 runs is the most broadly active organism in the novelty top 10, but the genus-level designation (sp.) means that the sequencing data is distributed across multiple species that may or may not correspond to the specific strains tested in the PlasticDB entries. Resolving Bacillus sp. to species level would substantially increase the actionability of its SRA data.

Serratia marcescens with 6,719 runs is an outlier in this group: it has extensive SRA coverage, driven primarily by clinical surveillance sequencing of this opportunistic pathogen, but no deposited enzyme sequence despite activity across 7 plastic types. Its evidence gap sub-score is therefore high (50 out of 100 on the evidence gap dimension), contributing to its elevated overall novelty score despite abundant sequencing data. For Serratia marcescens, the bottleneck is not sequencing depth but enzyme attribution: the genomic data exists to identify candidate enzymes by homology, but no biochemical characterization linking a specific gene to the observed plastic degradation has been published.

---

## 6. BacDive and SRA Data in the Context of Bioplastic Enzyme Discovery

### 6.1 Temperature Data as a Selection Criterion

The finding that 97.0% of plastic-degrading organisms with BacDive temperature records are mesophilic fundamentally reframes the challenge of room-temperature enzyme discovery. The dominant narrative in the PETase engineering literature since Tournier et al. (2020) has treated ambient-temperature activity as a property that must be engineered into enzymes whose natural state is thermophilic, using protein engineering strategies such as directed evolution, disulfide bridge introduction, and stability-activity tradeoff optimization. The BacDive temperature data suggest a different framing: the vast majority of organisms known to degrade plastic are mesophilic, meaning that nature has already produced plastic-degrading enzymes active at ambient temperature, and these enzymes are simply undercharacterized relative to the engineered thermostable variants that dominate recent publications.

The mean growth temperature of 29.5 degrees Celsius across all 435 organisms with BacDive temperature data is particularly significant in the context of wastewater treatment. Activated sludge systems in temperate climates operate between 10 and 25 degrees Celsius seasonally. An enzyme from an organism with an optimal growth temperature of 29.5 degrees Celsius is likely to retain substantial activity across this range, whereas FAST-PETase, engineered for activity at 37 degrees Celsius, may perform less predictably under seasonal temperature variation in an open treatment system. Room-temperature enzyme discovery from the mesophilic majority is therefore both scientifically justified by the BacDive data and practically motivated by operational constraints.

### 6.2 Oxygen Requirements and Reactor Design

The aerobic dominance among plastic-degrading organisms (80.3% aerobic or obligately aerobic) has practical consequences for reactor design in bioremediation applications. Aerobic bioreactors, including rotating biological contactors, aerated lagoons, and activated sludge systems, are the standard configuration for municipal wastewater treatment and are technically mature across the relevant scale range. The compatibility of the dominant plastic-degrading organism physiology with aerobic reactor configurations means that no exotic reactor technology needs to be developed to support a biocatalytic plastic degradation process based on these organisms or their enzymes.

The 14.4% of organisms that are facultative anaerobes is also favorable: this group would function effectively in aerobic conditions while being capable of survival and activity during anaerobic periods caused by mixing failure, fouling, or process upsets. A plastic-degrading consortium including both obligate aerobes and facultative anaerobes would likely show greater operational robustness than one composed exclusively of strict aerobes.

The 4.7% of organisms classified as strict or obligate anaerobes are less immediately relevant to aerobic plastic degradation systems, but their enzymes may be valuable for specific applications such as anaerobic digestion enhancement, where adding bioplastic-degrading capacity to an anaerobic sludge digester would improve its performance on bioplastic-contaminated organic waste streams.

### 6.3 Soil Isolation Dominance and the Ambient-Temperature Hypothesis

The prominence of soil as the isolation environment in BacDive records (124 of 207 environmentally categorized records) reinforces the ambient-temperature hypothesis through an ecological argument. Soil microbiomes are exposed to ambient temperature plastic pollution; plastic waste deposited in soil undergoes colonization by soil microorganisms, and the organisms that persist on plastic surfaces are those capable of using the polymer as a carbon or energy source. The enrichment of soil environments in the BacDive records of plastic-degrading organisms reflects the selection pressure of ambient-temperature plastic exposure and implies that the enzymes produced by soil-derived degraders are adapted to ambient conditions by the ecological context of their isolation.

This argument applies with particular force to LDPE, HDPE, PE, PP, PS, and PVC, the conventional plastics that are most abundant in terrestrial plastic waste but also among the most poorly characterized in terms of degradation mechanism. The soil-derived organisms in PETadex with confirmed activity against these substrates, particularly those with BacDive temperature records confirming mesophilic growth, represent a concentrated source of candidate enzymes for these undercharacterized substrates.

### 6.4 SRA Platform Trends and the Long-Read Transition

The emergence of Oxford Nanopore and PacBio data for 150 and 115 plastic-degrading organisms respectively marks a transition in the type of genomic resource available for enzyme discovery. Short-read Illumina assemblies of bacterial genomes typically produce draft assemblies with 20 to 200 contiguous sequence segments, depending on genome complexity and repeat content. Long-read assemblies routinely produce complete, single-contig circular chromosomes for bacterial genomes, resolving the fragmentation that obscures gene clusters and operons in draft assemblies.

For plastic-degrading enzyme discovery, complete genomes offer a specific advantage: plastic-degrading enzyme genes are frequently clustered with accessory genes encoding surface-binding modules, secretion signal peptides, regulatory proteins, and auxiliary metabolic enzymes that together constitute a degradation pathway. Short-read assemblies often fragment these clusters across multiple contigs, making pathway-level analysis unreliable. Long-read assemblies preserve cluster integrity, enabling researchers to characterize not just the catalytic enzyme but the full biological context required to understand its natural function and predict its behavior under expression in heterologous hosts.

The 150 organisms with Oxford Nanopore data and 115 with PacBio data collectively represent organisms for which complete pathway-level genomic analysis is already feasible using existing public data.

### 6.5 The Transcriptomic Dimension

The 172 RNA-Seq associations in the SRA data deserve specific attention. Transcriptomic studies of plastic-degrading organisms typically measure gene expression changes when organisms are grown in the presence of the target polymer compared to a reference condition. Genes that are significantly upregulated in the plastic-containing condition are candidates for encoding plastic-degrading enzymes or the accessory proteins that support their function.

The availability of RNA-Seq data for 172 organism-study combinations means that, for a substantial fraction of the plastic-degrading organism set, expression data already exists that could be used to prioritize enzyme candidates identified by genomic homology. An enzyme that is homologous to a known PETase AND is strongly upregulated in the presence of plastic substrate has a much stronger case for direct experimental characterization than one identified by homology alone. Integrating SRA RNA-Seq data with genome-based enzyme candidate lists constitutes a high-value next step in the PETadex analysis pipeline that is enabled directly by the SRA run count data collected here.

---

## 7. Key Insights and Outcomes

The following findings emerge from the analysis of BacDive and SRA data for the 932 confirmed and predicted plastic-degrading organisms in PETadex. Each finding is grounded in the specific numbers recovered and has a direct implication for the PetaBite research direction.

**Finding 1: Room-temperature plastic degradation is the statistical norm, not the exception.**
Of the 435 organisms with BacDive temperature data, 421 (97.0%) are mesophilic, with a mean growth temperature of 29.5 degrees Celsius. This directly supports the hypothesis that room-temperature plastic-degrading enzymes are abundant in the uncharacterized genomic space and that the apparent rarity of ambient-temperature activity is an artifact of the research community's focus on thermophilic engineering following the 2020 ICCG-LCC publication.

**Finding 2: Aerobic metabolism dominates, aligning naturally with aerobic treatment infrastructure.**
80.3% of organisms with oxygen requirement data are aerobic or obligately aerobic. This alignment with the aerobic conditions of wastewater treatment plants means that enzyme expression systems and bioreactor designs based on these organisms can use standard aerobic configurations without modification.

**Finding 3: Soil is the principal natural source of ambient-temperature plastic degraders.**
Soil isolation dominates the BacDive environmental record, appearing in 124 of 207 environmentally categorized records. Organisms isolated from soil are adapted to the ambient temperature and aerobic conditions of plastic pollution in terrestrial environments, making them more ecologically appropriate sources of wastewater-relevant enzymes than organisms from thermal springs or deep-sea hydrothermal vents.

**Finding 4: Three-quarters of plastic-degrading organisms have public sequencing data, but coverage is deeply uneven.**
690 of 932 organisms (74.0%) have at least one SRA run. However, the median of 19 runs versus the mean of 1,878 runs, driven by clinical surveillance data for Escherichia coli and Klebsiella pneumoniae, reveals that most plastic-degrading organisms have modest but sufficient sequencing coverage, while a handful of medically significant organisms vastly dominate the run count statistics.

**Finding 5: Long-read sequencing is arriving for plastic-degrading organisms at a useful scale.**
150 organisms have Oxford Nanopore data and 115 have PacBio data. These complete or near-complete genome resources support pathway-level enzyme discovery, which is qualitatively more powerful than homology-based enzyme identification from draft assemblies.

**Finding 6: Transcriptomic data exists for 172 organism-study combinations and is underutilized.**
RNA-Seq data associated with plastic-degrading organisms in SRA represents a ready-made experimental filter for enzyme candidate prioritization. A systematic analysis integrating SRA RNA-Seq expression data with homology-based enzyme candidates from genome sequences would substantially improve the precision of enzyme discovery without requiring any new experimental work.

**Finding 7: The highest-novelty organisms are genomically undersequenced relative to their biological importance.**
Seven of the top ten organisms by novelty score have fewer than 400 SRA runs each, despite plastic substrate ranges spanning 7 to 20 polymer types. Stenotrophomonas sp. (7 plastics, 128 runs), Enterobacter sp. (7 plastics, 135 runs), and Streptomyces sp. (16 plastics, 188 runs) are the most acute examples: their combination of broad substrate range and limited public genomics makes them the clearest targets for de novo or expanded sequencing investment.

**Finding 8: 163 organisms have neither BacDive nor SRA data despite confirmed plastic degradation.**
This group, representing 17.5% of the queried cohort, constitutes the highest-priority targets for basic data generation. Their plastic degradation is experimentally documented, but no physiology and no sequence information is publicly available to guide enzyme identification. Each organism in this group that is sequenced and physiologically characterized adds a genuinely new data point to the global knowledge base of plastic-degrading enzymes.

---

## 8. Limitations

The BacDive match rate of 50.9% introduces a selection bias that must be acknowledged. The 474 matched organisms are systematically more likely to have been formally cultured, deposited in a culture collection, and studied repeatedly by multiple groups, because these are the prerequisites for BacDive inclusion. The physiological conclusions derived from the matched set (mesophilic dominance, aerobic dominance, soil isolation prevalence) may therefore over-represent well-characterized, easily cultured organisms and under-represent the environmental metagenome fraction of the plastic-degrading organism pool. The true temperature distribution of all plastic-degrading enzymatic capacity, including enzymes from uncultured organisms accessible only through metagenomics, may differ from the distribution observed here.

The SRA run count distribution is severely right-skewed due to the inclusion of clinical surveillance organisms (Escherichia coli, Klebsiella pneumoniae, Pseudomonas aeruginosa) that appear in PlasticDB for bioplastic degradation but whose SRA presence is driven by pathogen surveillance rather than plastic degradation research. Summary statistics (mean, total run count) computed over the full cohort are therefore not representative of the plastic-degradation-specific research investment. Median values and post-exclusion summaries are more informative for this purpose.

BacDive temperature data are extracted from free-text strings using regular expression parsing, which introduces a risk of systematic error if BacDive authors use non-standard formats. The thresholds used for temperature binning (15 and 45 degrees Celsius) follow standard microbiological convention but are not universally agreed upon; alternative thresholds would produce different classification counts.

Name-based matching between PlasticDB organism names and BacDive or SRA records is imperfect. Genus-level designations (sp.), synonym groups, and name changes between the time of a PlasticDB publication and the current NCBI taxonomy record can cause failures to match records that logically belong to the same organism. The 458 BacDive non-matches include an unknown fraction attributable to this name resolution failure rather than true absence from BacDive.

---

## 9. Future Directions

The analysis reported here identifies several specific follow-on studies that are technically feasible using the data already collected in PETadex and the public databases queried:

**Systematic RNA-Seq integration.** The 172 RNA-Seq SRA associations should be systematically processed to extract differential expression data for organisms grown in the presence of plastic substrates. Combining these expression profiles with genome-based enzyme candidate lists would produce a ranked list of candidate enzymes with both homology support and transcriptional induction evidence, which is a substantially higher-confidence candidate list than homology alone.

**Long-read assembly-based pathway mapping.** The 150 organisms with Oxford Nanopore data and 115 with PacBio data should be assembled to complete circular chromosomes (if not already available in NCBI Assembly) and analyzed for plastic-degrading enzyme gene clusters. The identification of conserved gene cluster architectures across phylogenetically diverse organisms would provide evidence for convergent evolution of plastic degradation pathways and help prioritize which genomic contexts are most likely to contain functional enzymes.

**Targeted sequencing of the 163 neither-BacDive-nor-SRA organisms.** A subset of these organisms, selected based on novelty score, plastic substrate range, and availability of living cultures, should be prioritized for whole-genome sequencing and submission to SRA, and for physiological characterization and BacDive submission. This step would convert the most information-poor high-novelty organisms into actionable enzyme discovery targets.

**Expanded BacDive-temperature-correlated enzyme screening.** The 421 mesophilic organisms with BacDive temperature data and at least one SRA run should be assembled into a focused screening cohort for room-temperature enzyme discovery. For each organism, the genome should be annotated for esterase, cutinase, lipase, laccase, and depolymerase family genes, and candidates identified by homology to characterized PETases should be prioritized for recombinant expression and kinetic characterization at 25 and 37 degrees Celsius.

**Integration of BacDive isolation environment data with plastic substrate geography.** The soil-dominated isolation profile of plastic-degrading organisms suggests that soil metagenomics from plastic pollution sites would be a productive sampling strategy for new enzyme discovery. Cross-referencing the BacDive isolation environment data with geographic coordinates and published plastic pollution maps would identify specific sampling locations most likely to yield new plastic-degrading organisms whose enzymes are adapted to the ambient temperatures and aerobic conditions of terrestrial plastic waste environments.
