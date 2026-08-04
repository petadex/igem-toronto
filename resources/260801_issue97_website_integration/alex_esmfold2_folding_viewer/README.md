# Alex ESMFold2 Folding Viewer

## Goal

Show predicted 3D folds (and confidence metrics) on ORF / family / cluster pages using Alex’s S3 layout under `s3://petadex-protein-structures/`, without browsing millions of files.

## What shipped

Rough timeline Jul to Aug 2026. Branches: earlier structure prototype commits, then `feat/alex-structure-msa-lanes`.

1. **Resolve API**
   - Prefer experimental PDB from `pdb_accessions` when present.
   - Else predicted CIF if the object is publicly readable (Range GET probe).
   - Routes under `/api/structure` (orf, accession, metrics, content proxy).

2. **Alex schema (confirmed)**
   - `{lane}/structures/orf{id}.cif`
   - `{lane}/metrics/orf{id}.json`
   - Demo lane (public today): `esmfold2-centroids/test2`
   - Production lane: `esmfold2-centroids/60pid` (still 403 as of Aug 2026)
   - MSA experimental lane: `esmfold2-centroids/60pid-msa` (path confirmed by Alex; still 403)

3. **Metrics**
   - Parse Alex JSON (mean_plddt, ptm, per-residue pLDDT, PAE).
   - CIF B-factors match pLDDT.
   - Legacy npy / NPZ still supported if pointed at.

4. **UI**
   - `FoldingViewer` via `StructurePanel` on ORF sequence, family centroid, 90% cluster, curated Structure tab.
   - Mol*, confidence table, PAE heatmap, SAE stub.
   - Base / MSA toggle (MSA URLs only appear when that lane is configured and the CIF exists).
   - Relabeled former “finetune” wording to **MSA** (Alex: MSA is the experimental group).

5. **CORS**
   - Bucket CORS is petadex.net-oriented. Localhost uses `GET /api/structure/content/orf/:orfId` as a proxy.

6. **Verified sample**
   - ORF `4981589` on `test2` (CIF + metrics). Listing the bucket is still 403, so other example IDs are unknown unless someone shares keys.

## Key files (petadex.io)

- `backend/src/routes/structure.js`
- `backend/src/lib/structureMetrics.js`
- `frontend/src/components/structure/FoldingViewer.jsx`
- `frontend/src/components/StructurePanel.js`
- `frontend/docs/protein-structures.md`
- `backend/.env.example` (`STRUCTURE_S3_LANE`, `STRUCTURE_S3_MSA_LANE`)

## Env (after Dennis opens ACL)

```bash
# STRUCTURE_S3_LANE=esmfold2-centroids/60pid
# STRUCTURE_S3_MSA_LANE=esmfold2-centroids/60pid-msa
```

Default today stays on `test2`. Leave `STRUCTURE_S3_MSA_LANE` empty until `60pid-msa` is readable (avoids probing a private prefix on every resolve).

## Workflow

1. Ask Alex for lane names and file layout.
2. Probe public GET on a known ORF key (see `../scripts/probe_structure_s3.sh`).
3. Wire resolve + metrics + Mol* proxy.
4. Relabel UI when experimental naming changes (finetune to MSA).
5. Ask Dennis for ACL when prefixes stay 403.

## How to check

- Backend: `GET /api/structure/orf/4981589` and `/api/structure/metrics/4981589`
- Site: `http://localhost:8000/sequence/orf/4981589` (with API on :3001)
- Run the probe script and expect `test2` = 206, `60pid` / `60pid-msa` = 403 until Dennis acts.

## Outstanding

| Item | Owner |
|------|--------|
| Public GET on `60pid` and `60pid-msa` | Dennis |
| Optional ListObjects on those prefixes | Dennis |
| Optional localhost CORS on the bucket | Dennis / Alex (proxy already works) |
| Exp. figures beside the viewer | Purav (full folds) |
| Flip env to production lanes | Sara, after ACL opens |
