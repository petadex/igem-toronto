# Sequence and structure viewers

## Goal

Show sequence features and 3D structure on ORF / curated sequence pages in a way that matches InterPro-style expectations: zoomable tracks, residue highlights, and a structure tab.

Related igem issues: [#54](https://github.com/petadex/igem-toronto/issues/54), [#65](https://github.com/petadex/igem-toronto/issues/65), [#56](https://github.com/petadex/igem-toronto/issues/56).

## What shipped

Rough timeline May to Jun 2026:

1. **Nightingale / sequence UI**
   - InterPro-style enrichment layout.
   - Hover, labels, light/dark contrast, zoom fixes.
   - Amino acid letter highlight on hover.

2. **Annotation API prototype**
   - `/api/sara-viewer` for domains, motifs, signal sequence by `orf_id`.
   - Prototype page with annotation source toggle (including neXtProt-oriented paths and `sara_*` SQL annotations).
   - Demo mode with always-on UniProt for testing.

3. **Mol* annotations**
   - Annotated 3D viewer with residue highlights.
   - Wired into the structure tab.
   - Callouts / pop-outs for some annotations.
   - Bulk select fixes for the viewing window.

4. **Catalytic domains feature viewer**
   - Shared Petadex catalytic-domains feature viewer across sequence pages (`feat(viewer): add Petadex catalytic-domains feature-viewer...`).

5. **Dev hygiene**
   - Static CSS imports for Mol* / feature-viewer to avoid HMR `removeChild` failures.
   - Dev 404 page improvements.

## Key files (petadex.io)

- `frontend/src/components/protein/ProteinViewer.js`
- Annotated viewer / Nightingale components under `frontend/src/components/protein/`
- `backend/src/routes/saraViewer.js` (and related annotation routes)
- Structure tab wiring on curated sequence / ORF pages

## Workflow

1. Confirm which annotation table or API the dry lab wants on a given page.
2. Resolve ORF or accession to features + optional structure URL.
3. Render sequence tracks and Mol* side by side (or tabbed).
4. Smoke-test light and dark mode, zoom, and selection.

## How to check

- Open a curated sequence Structure / Features tab and an ORF sequence page.
- Confirm tracks load, hover works, and structure (when available) rotates / zooms.
- Toggle annotation sources on the prototype page if it is still enabled in your checkout.

## Outstanding

- Keep annotation sources in sync when new dry lab tables land (PTMs, SignalP, etc.).
- Predicted ESMFold2 folds now go through the Folding Viewer path (see sibling folder); experimental PDB path stays on ProteinViewer / pdb_accessions.
