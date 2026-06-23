# Deliverables for iGEM Dry Lab Members
**Date:** June 21st, 2026\
**Author:** Thomas Quigley

## Expectations
**Deliverable for July 20th:** Each member will be an author of the PETadex paper if they can provide the following for their given projects:
1. Data on the SQL database.
2. Data is displayed on the website (talk to Sara).
3. A paragraph written on the methods.
4. A paragraph written on the results.
5. A figure (in SVG format).
6. A table of summary stats.
7. Thorough documentation (we are publishing this)

This is a hard deadline, please get things done by then.

## Projects
- **Adi**
  - LLM agent to extract activity data
  - AI side of the activity scraping experiment (see below)
  - Graph RAG (openvirome.com)
- **Alex**
  - Validated PDBs (experimental)
  - Predicted PDBs (ESM, AF)
  - Prediction and measurement quality metrics
  - Fold all 90% centroids using ESMFold2 (or 60%/30% centroids if too expensive)
    - With PETadex MSAs
    - With PTMs
- **Angela**
  - DeepLoc
  - SignalP
  - Phylogenetic search tool
    - 90% ancestral reconstructions
    - 60% ancestral reconstructions of 90% ARs
    - 30% ancestral reconstructions of 30% ARs
- **Oscar**
  - Final Atlas (base ESMC) + controls
    - 3 resolutions (30% for all domains, 60% for each 30% cluster, 90% for each 60% cluster, 100% for each 90% cluster)
  - Add metadata search for the Atlas (for all Logan + GenBank)
  - Add boolean for if there is experimental data
  - Integrate all the pages to point to the atlas
  - Compare base ESMC to PETadex-ESM for the atlas
- **Nino**
  - A list of every single cluster of accessory domains within the PETadex
  - Validate if the domain exists against TED, then TEDs annotations tools to find novel
  - A page that shows all the accessory domains (TED)
  - Link these domains to TED, InterPro, CATH, all external DBs
- **Lisa**
  - CATH domain page
  - Mechanism
  - Literature about the domain
  - Show all HMMs for domains
  - HMM logos for each component
  - Sequence logos for all sequences in the component
  - Wrangle people to finish CATH domains, or do yourself
- **Claire**
  - Annotate the sequences with as much information about the structure as possible
    - PTMs (PTM Tack)
    - Cofactor predictions
    - Metal binding site predictions
- **Amar/Denis**
  - SRA Stat data for all organisms
  - A page for each organism
  - Integrate BacDrive data for each organism page 
  - Link to BacDrive
- **Purav**
  - Fine-tune ESMC
  - Does the fine-tuned model improve predictability of the metadata/activity/structure
  - Metadata prediction / label propagation
  - ??
- **Owen**
  - Human side of the activity scraping experiment (see below)
  - Create schemas for the scraping and refine them
  - Management of iGEM team members for scraping
- **Special case: Sara**
  - Just need you to add everyones data to the website, you can skip deliverables other than making sure everyone has their data on the website.
- **Declan**
  - ??
- **David**
  - 


## Scraping Experiment
We are going to take as many papers as possible that have plastic-degrading enzyme activity data (~250 is my guess), and we are going to have iGEM members extract activity data from these tables. To supplement this, we are going to have an LLM agent also extract activity data from these tables, so we can see where they disagree.

### Protocol
1. Design schemas to extract the activity from these papers.
2. iGEM members will extract data from 10 papers (~40 hours). 
   1. **They cannot use LLMs. If they do they are banned from the paper.**
   2. Hopefully we will have overlap between members activity extractions. These will get us replicates for each datapoint (for much higher confidence).
3. LLM Agents will extract data from the papers. We will have the same number of models as replicates of humans.
4. Where the agents and the humans disagree, we need to figure out who is right.
5. At the end, we will have scraped ~250 papers and a confidence interval we have in the model to predict correctly.
6. If the confidence is high, we can apply this to all papers.
