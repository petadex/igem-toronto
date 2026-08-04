# SRA / BioSample / BacDive hubs

## Goal

Give living sample provenance a place on the site: SRA runs, NCBI BioSamples, and organisms, using ingested `sra_metadata`, then attach Denis BacDive environmental means when available.

Related igem issue: [#28](https://github.com/petadex/igem-toronto/issues/28) (Automated Metadata).

## What shipped

Rough timeline Jul 2026. Branch: `feat/denis-bacdive-biosample-means` (stacked on the structure work for a clean upstream rebase).

1. **Hubs on top of `sra_metadata`**
   - `/biosamples`, `/sra/:acc`, `/biosample/:id`, `/organism/:name`
   - API under `/api/sra` (summary, run, biosample, organism search, pagination).
   - File-backed organism stats cache so search stays responsive.
   - Deep links from ORF provenance, curated metadata, identifier resolver, cluster organism stats, map popups.

2. **Denis BacDive CSV #3**
   - Per-BioSample means: optimum temperature and pH averages (plus counts).
   - Join key: `biosampleID` to BioSample id.
   - Loader / cache: `backend/src/lib/bacdiveMeansCache.js`
   - Example schema: `backend/data/biosample_bacdive_means.csv.example`
   - Env: `BACDIVE_MEANS_PATH` or `BACDIVE_MEANS_URL`
   - Status route: `GET /api/sra/bacdive/status`
   - UI: BacDive means panel on biosample pages (organism pages note availability). Replaces the old BacDrive stub when data is present.

3. **Docs**
   - `frontend/docs/sra-biosamples.md` (CSV #1 to #3 status table, cache invalidation).

## Key files (petadex.io)

- `backend/src/routes/sra.js`
- `backend/src/lib/bacdiveMeansCache.js`
- `backend/src/lib/sraStatsCache.js`
- `frontend/src/components/sra/SraShared.jsx`
- `frontend/src/pages/biosamples.js`
- `frontend/src/pages/biosample/[biosampleId].js`
- `frontend/src/pages/organism/[organismName].js`
- `frontend/src/pages/sra/[acc].js`

## Workflow

1. Confirm Denis’s S3 / CSV path and column names.
2. Prefer a local file or env URL (default remote URL was private / 403).
3. Parse once into `.cache/bacdive-biosample-means.json`.
4. Attach `bacdive` on biosample API responses.
5. Show a plain panel on the page when means exist; stay quiet when not.

## How to check

```bash
# Drop CSV then restart API
# BACDIVE_MEANS_PATH=./data/biosample_bacdive_means.csv
curl -s localhost:3001/api/sra/bacdive/status
```

- Open a BioSample that appears in CSV #3: means panel should render.
- BioSample missing from CSV: empty / unavailable, no crash.
- After CSV updates: delete `backend/.cache/bacdive-biosample-means.json` and restart.

## Outstanding

| CSV | Contents | Status |
|-----|----------|--------|
| #1 | Unique BacDive organisms in SRA stats | Not published |
| #2 | #1 filtered to organisms with environmental data | Not published |
| #2.5 | BioSample SRA rows trimmed to BacDive organisms from #2 | Not published |
| #3 | Per-BioSample T / pH means | Wired |

Also: keep organism hub cache in sync when Denis reloads `sra_metadata` (delete `backend/.cache/sra-organism-stats.json` and restart).
