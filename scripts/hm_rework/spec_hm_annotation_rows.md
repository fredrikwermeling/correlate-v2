# Heatmap rework: annotation-row order IS the hierarchy, plain sort controls, simpler clustering

App: /Users/fredrikwermeling/Documents/correlate_v2 (vanilla JS, single class in app.js, UI in index.html).
Serve for testing: `python3 scripts/export_audit/serve_capture.py` on port 8642 if not already running
(check with curl first). Open your own Chrome tab, and FIRST thing after every page load install the
de-throttle (occluded tabs are timer-throttled):
```js
const _st = window.setTimeout;
window.setTimeout = (fn, d, ...a) => _st(fn, Math.min(d || 0, 50), ...a);
window.requestAnimationFrame = (cb) => _st(() => cb(performance.now()), 16);
```
You MAY edit app.js and index.html. DO NOT commit, DO NOT push, DO NOT touch the version badge or the
changelog (the main session does that after review). Keep a working page at all times; `node --check app.js`
after every substantial edit.

## Why (user feedback, verbatim gist)

The user does not understand the current annotation-row design: "what is block?"; finds clustering as an
annotation row complicated, with too many settings; the cluster band colours add nothing at the default
auto-k of 2; wants clustering usable as the top layer with e.g. lineage as a downstream colour indicator,
AND lineage-first with clusters within each lineage; concludes "the easiest is to just define the order of
the annotation rows, and the hierarchy is defined according to order"; and wants each relevant row to offer
sorting by score high/low and low/high, by name, or by its own category levels (e.g. BRAF hotspot 0,1,2).

## Current model (read this code before designing)

- `_hmAnnRows`: up to 4 rows `{mode, gene, sortDir, block}`; modes: hotspot / fusion / cn / lineage /
  subtype / disease / cluster / ge / expr. Per-row buttons: ▦ (block toggle, one row at a time),
  ↕ (three-state sortDir null/'desc'/'asc'), ↑↓ (reorder), × (remove), plus "+ Add row" and a terminal
  "then by" select (`hmThenBy`: score/name).
- The OUTERMOST toggled row already drives the group strip / legend / drill / min-n
  (`_hmGroupSpecFromRows`, `_hmSortChainInfo`); other toggled rows sort within, top-down; rows below a
  cluster row or a continuous row are inert (`_hmSortInertReason`).
- Cell-line clustering exists TWICE: the "Cluster cell lines" checkbox + k select ("Auto (suggested)",
  2-12) + tree-detail select on the CLUSTERING toolbar row, mirrored both ways with a 'cluster'
  annotation row (max one). Cluster k-cut feeds groups "Cluster 1..k". Whole-cohort tree when outermost,
  per-innermost-block trees (g.clusterRoots) when below toggled rows.
- Key functions: `_hmResolveAnnotationRows`, `_hmSortChainInfo`, `_hmGroupSpecFromRows`,
  `_hmSortInertReason`, `_hmAnnRowNumericValueFor`, `_hmSortSummary`, `_hmBuildAndPaint` (ordering
  section), `_hmBuildAnnotation2`, `_hmViewState` / `_hmRestoreView` (save/reopen .json with migration
  from older formats), `_hmOpenModal` (defaults).
- Consumers of the ordering words/state that MUST stay correct: on-screen hint line (`hmHint`),
  exported figure caption (`_hmComposeCanvas` caption), CSV corner cell (`_hmExportCsv`), Methods text
  (`_methodsHeatmap`, `_hmSortSummary(d)`), AI export heatmap branch (~app.js:30800), save/reopen.

## Target design

### 1. Row order defines the hierarchy; the ▦ block toggle is GONE
- The list order of `_hmAnnRows` is the sort hierarchy. The TOPMOST row whose sort is ON forms the
  blocks (group strip + labels + legend + drill + min-n), the next sorted rows sort within those blocks,
  top-down, exactly as `_hmSortChainInfo` already does for the toggled chain. `hmThenBy` stays as the
  terminal tie-break. Rows with sort OFF are colour indicators only, wherever they sit.
- Remove the ▦ button and the `block` flag everywhere (state, UI, view-state; `_hmRestoreView` must
  migrate old saved views: a row saved with `block:true` moves to the FRONT of the row list with its
  sort on; other toggled rows keep their relative order behind it).
- The ↑↓ reorder arrows and × stay. The helper line under the rows is rewritten in plain words, e.g.:
  "Rows are applied top to bottom: the first row with Sort on splits the columns into blocks and rows
  below it sort within each block. Rows with Sort off just paint colours." (User is a general user; no
  jargon, no em-dashes, use / or commas.)

### 2. One per-row Sort select replaces ▦ + ↕
Replace the two icon buttons with ONE compact select per row, labelled by its options (keep row width
stable, no layout shift, reserve the space). Options depend on the row's mode:
- categorical (lineage / subtype / disease): Sort: off | Largest group first | Smallest group first |
  A to Z | Z to A | Score high to low | Score low to high   ("Score" = the block's median of the
  heatmap's own per-line score, the same clScore/medianScoreOf machinery the group labels use)
- alteration (hotspot / fusion / cn): Sort: off | Wild-type first (0, 1, 2) | Altered first (2, 1, 0) |
  Score high to low | Score low to high
- continuous (ge / expr): Sort: off | High to low | Low to high
- cluster: Sort: off | Tree order   (see 3; when off the row is inert and greyed like today's inert rows)
Internally keep `sortDir` for direction plus a new `sortKey` field; map old saved `sortDir` values to
the equivalent key+direction in `_hmRestoreView` (categorical desc = Largest first, asc = Smallest
first; alteration desc = Altered first, asc = Wild-type first; continuous desc = High to low).
The existing comparison code in the ordering section extends to the new keys (A to Z = localeCompare on
category label; Score = medianScoreOf per category). Blocks and within-block sorting use the same key
logic at both levels, as today.

### 3. Clustering: one place, fewer settings, no pointless colours
- REMOVE the "Cluster cell lines" checkbox, its k select and its tree-detail select from the CLUSTERING
  toolbar row (the "Cluster genes" checkbox + its tree-order/detail selects STAY, they are about gene
  rows). Cell-line clustering is now ONLY the "Cell-line clusters" annotation row.
- The cluster row carries exactly TWO controls: the Sort select (off | Tree order) and a small k select:
  "No split (tree only)" (default) | Auto | 2..8. Tree-detail for cell lines moves into the Settings
  panel next to the gene tree detail (both are power-user knobs).
- With "No split": the row orders columns by the dendrogram (whole-cohort when topmost, per-block when
  under a sorted row), the top dendrogram is drawn, but NO coloured band is painted for this row and NO
  Cluster 1..k groups are formed (so when topmost with no split, there is no group strip at all, plain
  tree order). With Auto/k: exactly today's k-cut behaviour (band, Cluster 1..k groups when topmost).
  This answers "k=2 colours add no information": by default there is no band, the tree shows the
  structure.
- Keep the auto-k silhouette machinery for the Auto option. Keep per-block trees for a cluster row
  below a sorted row (the user's "sort by lineage and within that the clusters").
- Migration: old saved views / restores with clusterCells true map to a cluster row with Sort on and
  k = the saved k (or Auto); the checkbox key disappears from new view-states but old keys must still
  restore.

### 4. Wording sweep
`_hmSortSummary`, the hint line, the caption, the CSV corner, Methods and the AI export must describe
the new model in the same words ("blocked by X (largest group first), then Y (wild-type first), then by
score"). Grep every consumer of `sortPlan` / `_hmSortSummary` / the word "blocks toggle" and update.
The word "blocks" may stay in DESCRIBING text (it reads fine: "split into blocks"), it is the ▦ control
that goes.

## Traps (from the project's history, do not rediscover these)
- Group-label plan is computed in layout BEFORE gridH is fixed; probe font and paint font must both come
  from `_hmSettings.labelFont`.
- Any change to per-column bands changes `groupExtra`/`ann2Extra` and must be audited in gridH, hit-test
  (`_hmHitTest*`), label painting, and `_hmComposeCanvas` (the export composer) together.
- Hidden-groups signature (`_hmHiddenGroupsSig`) must match the active grouping; clusters use 'cluster|k'.
- Synthetic groups need `orderedCellLines`, not just `cellLines` (CSV export reads it).
- State living outside form controls must survive `_hmOpenModal`'s reset AND ride `_hmViewState` /
  `_hmRestoreView`; set restored state AFTER the modal open.
- Never let controls appear and shove others aside: reserve space or keep them present but disabled.
- UI copy: plain words for a general user, no em-dashes (use / or commas).

## Browser smoke test (do all of these yourself, look at screenshots)
1. Fresh open: default = one untoggled TP53 row, plain score sort, no blocks, no cluster band. Caption
   and hint read correctly.
2. Add lineage row, Sort = Largest group first: lineage blocks with labels/legend; then set its Sort to
   A to Z and to Score high to low; blocks reorder accordingly.
3. Add Cell-line clusters row BELOW lineage (Sort = Tree order, No split): per-lineage dendrograms above
   the grid, columns tree-ordered within each lineage block, no cluster band.
4. Move the cluster row to the TOP (↑): whole-cohort tree order, NO group strip (no split), lineage row
   below now paints colours only. Then set k = Auto: Cluster 1..k blocks appear with the band and
   legend; lineage stays a colour band below.
5. TP53 hotspot row, Sort = Wild-type first: block order Wild-type, One copy, Both copies; Altered
   first reverses it.
6. min-n, "Hide lines without data", hidden groups (click a legend entry), drill, gates: all still work
   with the new model.
7. Export image (PNG) + .csv + Methods + Save view -> reopen the saved .json: everything reflects the
   new ordering words and state; the saved view round-trips; ALSO hand-write an OLD-format view state
   (block:true on a non-first row + clusterCells:true) and confirm migration.
8. Phone width (resize to 400px): the row controls do not overflow; nothing shoves other controls.

## Report
Summarize what changed (files, functions), what you tested with which screenshots (keep them in
scripts/hm_rework/out/), anything you deliberately did not do, and any open questions. DO NOT commit.
