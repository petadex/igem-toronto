# Phylo tree navigation

## Goal

Make family phylogenetic trees usable on the site: find a tip, see where it sits, inspect neighbors, narrow to a local region, and color leaves by metadata. Dry lab asked for navigable search on the existing tree panel (not a separate stack).

## What shipped

Rough timeline Jun to Jul 2026:

1. **Basics**
   - Read public S3 Newick / tree files without AWS credentials locally.
   - Shared viewer module, tree search, and links from search results.

2. **Viewer prototypes**
   - Compared interactive vs static viewers on a temporary prototype page, then reverted the comparison page once the production path was clear.

3. **Navigation tools on the live panel** (`showNavTools` on family / tree / prototype routes)
   - Path to root (traceback).
   - Nearby sequences (patristic distance / hop distance).
   - Local neighborhood: radius or k-NN, with clade toggle and “dim distant sequences”.
   - Metadata coloring (component, family_pid, organism, country, etc.).

4. **UX cleanup**
   - Removed “Fit to neighborhood”.
   - Dynamic k max and % presets.
   - Layman sidebar copy (“By tree distance”, “Closest N”, …).
   - Docs: `frontend/docs/phylo-tree-navigation.md`.

## Key files (petadex.io)

- `frontend/src/components/phyloTree/PhyloTreePanel.jsx`
- `frontend/src/components/phyloTree/PhyloTreeViewer.jsx`
- Topology helpers under `frontend/src/components/phyloTreePrototype/` (or merged equivalents)
- Family / tree pages that pass `showNavTools`

## Workflow

1. Confirm tree files exist for the family on the public S3 layout the API already uses.
2. Load Newick in the existing panel.
3. Add nav sidebar controls that only change highlight / dim / color (do not replace the renderer).
4. Document controls for non-experts.

## How to check

- Open `/family/:id` or `/tree/:id` (and `/phylo-tree-prototype` if still present) with `showNavTools`.
- Search or focus a tip: path-to-root should highlight.
- Neighbors list should update; radius / k-NN should dim outside sequences.
- Color by metadata and confirm the legend matches.

## Outstanding

- More metadata fields can be added as joins become available.
- Very large trees may need further performance work (not the main gap right now).
