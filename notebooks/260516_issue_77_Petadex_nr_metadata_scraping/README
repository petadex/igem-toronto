# Petadex nr metadata scraping pipeline (read me made using AI cause I could not be bothered lol sorry, also the comments in my code and some of the code itself like print statements and functions were also made with ai)

An automated pipeline designed to construct a comprehensive full-text literature corpus for domain-specific language model fine-tuning. The workflow maps biological database identifiers to academic literature, resolves cross-repository duplicates, and handles authenticated bulk retrieval across multiple publisher APIs and institutional proxies.

## Pipeline Architecture

The harvesting process is divided into five chronological phases. It begins with raw biological identifiers and concludes with a cleanly formatted, deduplicated file directory ready for machine learning ingestion.

### Phase 1: Database Identifier Mapping
Queries biological literature APIs to map protein accession numbers from local target databases to unique PubMed Identifiers (PMIDs).
* **`phase1_accession_mapper_basic.py`**: Executes chunked API requests to the NCBI Entrez API to map identifiers efficiently.
* **`phase1_accession_mapper_robust.py`**: An advanced implementation optimized for Google Colab environments featuring JSON payload parsing, offline pre-filtering of incompatible biological identifiers, and auto-healing capabilities for invalid sequences.

### Phase 2: Unstructured Metadata Resolution
* **`phase2_journal_mapper.py`**: Processes unstructured journal strings and cross-references publication metadata utilizing the Europe PMC API to compile a supplementary dataset of target literature. It utilizes a two-pass system (Exact Match and Rescue Pass) to maximize resolution success rates.

### Phase 3: Target Consolidation & DOI Resolution
* **`phase3_doi_resolution_engine.py`**: Acts as the central traffic controller. It merges the outputs of Phase 1 and Phase 2, performs set intersection to remove duplicates, and resolves the finalized PMIDs to Digital Object Identifiers (DOIs). It outputs the master `still_missing.csv` ledger used by all downstream extraction tools.

### Phase 4: Automated API & Direct Document Retrieval
A suite of targeted scripts designed to interface with publisher APIs and Open Access aggregators to retrieve full-text documents.
* **`phase4_elsevier_api_fetcher.py`**: Retrieves full-text XML from ScienceDirect utilizing the Elsevier Article Retrieval API. Includes a byte-size trapdoor to automatically reject phantom abstracts.
* **`phase4_wiley_api_fetcher.py`**: Interfaces with the Wiley Text and Data Mining (TDM) API to securely download full-text PDFs using specific institutional authentication headers.
* **`phase4_crossref_tdm_fetcher.py`**: Queries the CrossRef Metadata API to locate secondary TDM links prioritizing XML formats over PDF.
* **`phase4_epmc_dual_engine_scraper.py`**: Targets Europe PMC for Open Access PDFs using a sophisticated dual-engine approach. It defaults to a stealth HTTP/2 client and automatically degrades to a standard HTTP/1.1 client if the connection stream is unexpectedly closed.
* **`phase4_unpaywall_splash_piercer.py`**: Uses the Unpaywall API to identify Open Access links. If an endpoint routes to an HTML landing page rather than a direct file, it utilizes BeautifulSoup to parse the DOM and extract the hidden PDF source.
* **`phase4_direct_doi_cloudflare_scraper.py`**: Utilizes `curl_cffi` to spoof browser TLS fingerprints to bypass institutional Cloudflare challenges and extract embedded PDF assets directly from publisher domains.
* **`phase4_smart_cloudflare_spoofer.py`**: A specialized pre-processor that filters out unusable structural datasets (like PDB files) and isolates specific publishers into dedicated queues for authenticated Selenium processing.

### Phase 5: Institutional Proxy & UI Automation
Heavy-duty browser automation tools designed to handle documents hidden behind complex authentication walls or custom JavaScript repositories.
* **`phase5_custom_repository_scraper.py`**: Utilizes `undetected_chromedriver` and `PyAutoGUI` to navigate target repositories and simulate user interactions for PDF retrieval. Includes robust state recovery for uninterrupted batch processing.
* **`phase5_ezproxy_authenticated_fetcher.py`**: Automates retrieval through an institutional EZProxy portal. Features manual authentication handoff, adaptive human pacing (jitter), and scheduled system pauses to evade automated bot detection rate limits.

## Data Schema Summary
* **Total Aggregated Source Listings:** ~12,000
* **Unique Consolidated Targets:** 10,386 unique DOIs
* **Pipeline Extraction Success Rate:** 96.30% (10,002 valid documents secured) which can be found at https://petadexstorage.blob.core.windows.net/petadex-metadata-papers

## Dependencies
This pipeline requires standard data science libraries (`pandas`, `tqdm`, `beautifulsoup4`) alongside specialized networking and automation modules:
* `biopython`
* `curl_cffi`
* `undetected-chromedriver`
* `selenium`
* `PyAutoGUI`
