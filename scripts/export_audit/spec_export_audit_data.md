# Export audit B: every DATA export (CSV, view files, Methods txt, AI) under varied filters

READ-ONLY AUDIT. You do not edit app.js or index.html. You test, capture, inspect, and report findings.
App: /Users/fredrikwermeling/Documents/correlate_v2 (vanilla JS; app object `window.app`).
OUT dir for everything you save: /Users/fredrikwermeling/Documents/correlate_v2/scripts/export_audit/out/

## Setup

Same as audit A (another agent may already run the server; share it, prefix your capture names `dataaudit_`):
check http://localhost:8642 with curl; if down, start
`python3 /Users/fredrikwermeling/Documents/correlate_v2/scripts/export_audit/serve_capture.py` in the
background. POST /save?name= APPENDS to the OUT dir/<name>; one unique name per capture.
Open your own tab. MANDATORY, FIRST thing after every page load, before ANY other interaction: de-throttle
the tab. Your tab is almost always occluded, and Chrome treats it as background: setTimeout clamps to
~once a minute and requestAnimationFrame never fires, so analyses appear to hang and you will misdiagnose
app bugs that are not there (this cost the 2026-08-20 run a lot of wall-clock):
```js
const _st = window.setTimeout;
window.setTimeout = (fn, d, ...a) => _st(fn, Math.min(d || 0, 50), ...a);
window.requestAnimationFrame = (cb) => _st(() => cb(performance.now()), 16);
```
Then install the capture shim (blob AND data:-URI downloads; the createObjectURL-only shim missed every
export that goes out as a data: URI):
```js
window._caps = [];
const _ocu = URL.createObjectURL.bind(URL);
URL.createObjectURL = (b) => { const u = _ocu(b); window._caps.push({blob: b, type: b.type, size: b.size}); return u; };
const _oclick = HTMLAnchorElement.prototype.click;
HTMLAnchorElement.prototype.click = function() {
  if (this.download && this.href?.startsWith('data:')) {
    try {
      const [meta, payload] = this.href.split(',');
      const mime = meta.slice(5).split(';')[0] || 'application/octet-stream';
      const bin = meta.includes('base64') ? atob(payload) : decodeURIComponent(payload);
      const arr = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) arr[i] = bin.charCodeAt(i);
      window._caps.push({blob: new Blob([arr], {type: mime}), type: mime, size: arr.length, download: this.download});
    } catch (e) { window._caps.push({download: this.download, err: String(e)}); }
    return; // suppress the real download
  }
  return _oclick.apply(this, arguments);
};
async function saveCap(name, idx = -1) {
  const c = window._caps.at(idx); if (!c?.blob) return c || null;
  await fetch('/save?name=' + encodeURIComponent(name), {method: 'POST', body: c.blob});
  return {type: c.type, size: c.size, download: c.download};
}
```
One export can push several caps: pick the one whose `type` matches what you asked for, not blindly the
last. AI exports are slow by nature (60-90 s, they build matrices over 1,208 lines): wait, don't diagnose.
For gzipped AI exports also patch `window.pako.gzip` to POST the pre-compression JSON (chunk any body over
4 MB into sequential POSTs to the same name).

## What to audit

Exercise each export twice: once in a plain state, once with real filters/settings applied (a lineage
filter, a hotspot filter, gates, a subset - whatever that view supports). Verify the exported file FOLLOWS
the state change.

CSV exports:
1. Correlation scatter: its .csv exports (points / data)
2. Selection inspect: ".csv, all genes" AND ".csv, shown genes per cell line" - in rest mode AND gate mode
3. Mutation analysis: results table CSV; the compare-modal Export CSV
4. Network view: correlations CSV, clusters CSV
5. Heatmap .csv (default AND the worked state: toggled block row, hidden group, silenced gene, min-n, gates
   - the corner cell and Group row must reflect ALL of it)
6. Cell Line Browser: "List on screen: CSV" and "List on screen: copy" (copy: no-error check only), plus
   "Genes: full profile" export if present
7. Cell line wiki: "Export gene info (.csv)"
8. UMAP/PCA export(s) if present
9. Oncoprint "Export data (CSV)" if reachable

Other data formats:
10. Heatmap "Save view" .json: save a worked state, validate the JSON shape, then RESTORE it through the
    real reopen path and confirm the drawing reproduces (column count, blocks, hidden groups)
11. Methods "Download .txt" from three different views (heatmap, mutation analysis, scatter): both sections
    present, filename pattern <view>_methods_DepMap26Q1_<date>.txt, content matches the on-screen state
12. Export for AI, standard: three sources (scatter, mutation, network) - capture the JSON, check
    whatIsInThisFile.present/absent honesty, context matches the view's settings, cellLineOrder consistent
    with matrices, extras._listScope present. (The selection/gate source was audited separately on
    2026-08-19 - skip it.)
13. Custom export for AI: one request block naming ~6 genes + a tissue cohort; confirm the parser applies it,
    the export honors gene list + cohort, and cohort replacement drops stale extras with a stated reason.

## What "as expected" means (check each explicitly)

- First-line/description convention: files open with a line describing their own scope, and it must match
  the CURRENT filters (change a filter, re-export, confirm the line follows).
- csvName convention: filenames carry DepMap release + date.
- Formula-injection guard: no cell begins with = + - or @ unescaped (spot-check text columns AND try to
  provoke it: a pasted-list label or gene that could start with a dash; plain negative NUMBERS are exempt
  by design and must stay numeric).
- Counts: row counts consistent with what the view shows (respecting stated caps; caps must be stated).
- Separators/quoting: fields containing commas or quotes survive a proper CSV parse (parse with python csv).
- Encoding: UTF-8, no mojibake in cell-line names.
- Gate-mode inspect CSVs must name gate A / gate B, not "selection"/"rest".
- The heatmap CSV column order must equal the drawn order.
- AI export: no section documented that is absent, no field claiming something false about the view (read
  context critically against the on-screen state, the way a blind reader would).

## Report

Findings table: export / state / VERDICT (ok, cosmetic, broken) / evidence file. Prose for every non-ok:
expected vs got, repro steps. Console errors anywhere = a finding. Captured files under the OUT dir with
self-explanatory names, listed. DO NOT fix anything. DO NOT commit. Close your tab when done; leave the
server running.
