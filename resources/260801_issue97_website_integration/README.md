# Sara: PETadex website integration

Index of website work that wires dry lab outputs onto petadex.net / [petadex.io](https://github.com/petadex/petadex.io).

**Workflow issue:** [petadex/igem-toronto#97](https://github.com/petadex/igem-toronto/issues/97)

**App code:** https://github.com/petadex/petadex.io (Sara fork: https://github.com/sarapr06/petadex.io)

Thomas’s dry lab note listed Sara as the person who puts other people’s data on the website. This folder is the record of how that was done, plus a few small helper scripts. It does not copy the website codebase.

## Themes

| Folder | What it covers | Older igem issues |
|--------|----------------|-------------------|
| [cath_domain_pages/](cath_domain_pages/) | CATH domain pages, lit, HMMs, inline Mol* | #55, #54, #65 |
| [sequence_and_structure_viewers/](sequence_and_structure_viewers/) | Nightingale, Mol* annotations, feature viewer | #54, #65, #56 |
| [phylo_tree_navigation/](phylo_tree_navigation/) | Family tree search and navigation tools | (website, 2026) |
| [alex_esmfold2_folding_viewer/](alex_esmfold2_folding_viewer/) | Alex predicted folds on the site | (website, 2026) |
| [sra_biosample_bacdive/](sra_biosample_bacdive/) | SRA / BioSample hubs and Denis BacDive means | #28 |
| [scripts/](scripts/) | Probe scripts (S3 ACL checks, etc.) | |

## How website work usually goes

1. Dry lab person ships data (S3 path, CSV, SQL table, or notebook result).
2. Confirm the public URL or API shape (and ACL if S3).
3. Wire resolve + UI on petadex.io.
4. Open a PR on the fork, rebase onto `petadex/petadex.io` when needed.
5. Log the result on the workflow issue thread.

## Current blockers (Aug 2026)

- Dennis: public GET on `esmfold2-centroids/60pid/` and `esmfold2-centroids/60pid-msa/` (same ACL as `test2`)
- Purav: full folds before Exp. figures beside the Folding Viewer
- Denis: BacDive CSV #1 / #2 / #2.5 not published yet (#3 is wired)

## Related open PRs / branches (fork)

- `feat/alex-structure-msa-lanes`
- `feat/denis-bacdive-biosample-means`
- `feat/cath-inline-molstar`
