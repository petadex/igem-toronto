# CATH domain pages

## Goal

Give each curated CATH / Pfam-style domain a site page with literature, figures, HMMs, and a structure view, and link those pages to the atlas and enzyme pages.

Related igem issues: [#55](https://github.com/petadex/igem-toronto/issues/55) (Add CATH domain section), [#54](https://github.com/petadex/igem-toronto/issues/54) / [#65](https://github.com/petadex/igem-toronto/issues/65) (protein vis).

## What shipped

Rough timeline Apr to Jul 2026 on petadex.io:

1. **Skeleton + lit**
   - CATH section on the site, then filled Pfam narrative content and figures from lit review.
   - Atlas deep links that open filtered atlas views in a new tab.
   - Citation tooling so in-text references and a reference list stay consistent.

2. **Catalog tooling**
   - Profile HMM download and logos.
   - Domain architecture diagrams.
   - Site cross-links (atlas, enzymes, external DBs).

3. **Citation / validation fixes**
   - Citation-order validation for PF01425, PF01083, PF09995.
   - Cleaned duplicate InterPro refs and trailing commas in in-text URLs.
   - Full DOI canonicalization so papers like `10.1128/AEM...` do not collapse into one entry.
   - All 23 catalog profiles pass `build:validate-cath`.

4. **Inline Mol* (Jul 2026)**
   - Branch: `feat/cath-inline-molstar`
   - Replaced the molstar.org iframe with in-page `ProteinViewer`.
   - PDB dropdown when a domain lists multiple IDs.
   - Download PDB link to RCSB (`files.rcsb.org/download/{pdb}.pdb`).
   - Uses `accession={rcsbUrl}` so it works with upstream ProteinViewer (http URL path).

## Key files (petadex.io)

- `frontend/src/components/cath/CathDomainVisualizationPanel.js`
- CATH profile / catalog content under the frontend CATH content modules
- Citation helpers used by `build:validate-cath`

## Workflow

1. Curate or update a domain profile (text, figures, pdbIds, references).
2. Run catalog validation (`build:validate-cath` on the frontend).
3. Check the domain page: overview, HMM, architecture, structure block.
4. Confirm atlas / enzyme links open the right filters.

## How to check

- Open a CATH domain page that has `pdbIds`.
- Confirm the structure loads in-page (no molstar.org iframe).
- If multiple PDBs: switch the dropdown and confirm reload.
- Download PDB should open RCSB in a new tab.
- Domain with no `pdbIds` should still show the empty / placeholder state.

## Outstanding

- Keep profiles updated when Lisa / others finish more domain writeups.
- After Alex / Dennis open more experimental structures, consider linking PETadex predicted CIFs on domain pages where ORFs map cleanly (not done yet).
