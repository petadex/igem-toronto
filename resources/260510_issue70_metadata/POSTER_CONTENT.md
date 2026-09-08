# PETadex Organism Atlas

The PETadex Organism Atlas was developed as a searchable reference for plastic-degrading and plastic-associated microorganisms. Data from PlasticDB, NCBI taxonomy and genome assemblies, BacDive, SRA, PubMed, and protein analyses were integrated into a unified PostgreSQL dataset containing approximately 2.9 million organisms. The completed interface supports organism search, evidence and confidence filters, sorting, pagination, taxonomic summaries, and detailed organism profiles. Each profile combines plastics studied, publication history, genome metadata, environmental context, enzyme families, sequence availability, temperature traits, and a novelty score. The application was prepared for integration with the existing Gatsby frontend, Express API, AWS Lambda deployment, and AWS RDS database while preserving the established TaxID API and database schema.

## Key Outcomes

- Integrated approximately 2.9 million organism records and 2,535 PlasticDB entries.
- Classified records into confirmed, predicted, and listed evidence tiers.
- Connected research evidence with taxonomy, genome, BacDive, SRA, PubMed, and ProtParam metadata.
- Added interactive search, filtering, sorting, charts, pagination, and organism detail profiles.
- Added validation for JSON fields, duplicate names, orphaned entries, and transactional database reloads.

## Conclusion

The Atlas converts distributed biological and literature data into a single research interface for comparing plastic-associated organisms. The results show that evidence is concentrated in a relatively small group of confirmed organisms, while the broader taxonomy provides a large search space for identifying related strains, research gaps, and candidates for future biodegradation studies.

![PETadex Organism Atlas overview](images/organism-atlas-overview.jpg)

**Figure 1.** Searchable organism overview showing dataset totals, confidence tiers, filters, taxonomic fields, research indicators, and novelty scores.

![Enriched organism profile for Pseudomonas aeruginosa](images/organism-detail-profile.jpg)

**Figure 2.** Enriched profile for *Pseudomonas aeruginosa* combining plastic evidence, genome metadata, BacDive physiology, publication trends, enzyme families, novelty scoring, and PlasticDB records.