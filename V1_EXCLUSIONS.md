# What comes out for V1

V2 (this repo) is the extended version. V1, published at correlate.cmm.se, is
the simpler build. This file lists what V1 leaves out, so the removal is a
known set rather than a hunt.

Every block below is marked in the source with a comment containing
`V2 ONLY`, so `grep -n "V2 ONLY" index.html app.js` finds them all.

`V1_EXCLUSIONS.md` lists what comes out. **`scripts/build_v1.py` performs the
removal** — run it instead of editing V1 by hand, then copy `web_data/` only
when the DepMap release changes.

## Dimensionality reduction (PCA / UMAP)

- `index.html`: the whole `#clbUmapSection` block, marked
  `===== V2 ONLY: Dimensionality Reduction (PCA / UMAP) =====`.
  It is already collapsed behind a button, so removing it leaves no gap.
- `app.js`: the `clbUmapToggle` click handler, and the `clbUmap*` methods
  (`clbUmapSelectMode`, `clbUmapCopySelected`, `clbUmapSelectInList`,
  `_resetUmap`, `_renderUmap*`).
- The UMAP gate panel is inside this block and goes with it.

## Gates

Two rectangles drawn on a plot, then the cell lines inside them compared.
Four separate places, all marked `V2 ONLY: gates`:

1. `index.html` scatter controls: the `control-box-compact` holding
   `#setGateABtn` / `#setGateBBtn` / `#compareGatesBtn` / `#clearGatesBtn` /
   `#gateStatus`.
2. `index.html` gene-effect toolbar: `#geGatesGroup`.
3. `index.html` gene-effect results: `#geGateComparePanel`.
4. `index.html` UMAP: `#clbUmapGatePanel` (inside the block above).

`app.js` side: `clearGEGates`, `compareGates`, `enrichrGateDiffGE`,
`enrichrGateMutations`, the `gateAnalyzeWithAI` handler, and the `drawrect`
dragmode branches that exist to draw a gate.

Note `clearGEGates()` is called from `_resetGEFilters()`, so either keep a
no-op stub or drop that call too.

## Expression basis (gene effect vs mRNA): NO LONGER EXCLUDED

Was V1-excluded until 2026-08-31, when Fredrik asked for mRNA correlation
in V1 too. The `#basisParams` block now ships in both builds; the guards
described below stay, they just never fire in either build any more.

`app.js` needs no edit. `_analysisBasis()` reads the radios with
`?.value || 'ge'`, and every other reference is null-guarded, so a build
without the markup simply always runs on gene effect.

## Not yet decided

Raised but not settled, so still in V2 and not marked:

- mutation analysis
- the cell line wiki
- Enrichr
- expression correlates

Ask before assuming any of these leaves.

## After removing

- The version badge, the `?v=` cache-bust on the `app.js` script tag, and the
  changelog are hand-edited in three places; see the memory note on versioning.
- Search for now-dead ids before shipping: a removed button whose handler
  remains is harmless, but a handler whose button is gone throws on load if it
  is not guarded with `?.`.
