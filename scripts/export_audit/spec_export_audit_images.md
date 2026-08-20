# Export audit A: every IMAGE export, every format, varied settings

READ-ONLY AUDIT. You do not edit app.js or index.html. You test, capture, inspect, and report findings.
App: /Users/fredrikwermeling/Documents/correlate_v2 (vanilla JS; app object `window.app`).
OUT dir for everything you save: /Users/fredrikwermeling/Documents/correlate_v2/scripts/export_audit/out/

## Setup

1. Server: check `curl -s -o /dev/null -w "%{http_code}" http://localhost:8642/index.html`. If not 200, start
   `python3 /Users/fredrikwermeling/Documents/correlate_v2/scripts/export_audit/serve_capture.py` in the
   background (it serves the repo AND accepts `POST /save?name=<file>` which APPENDS the body to
   the OUT dir/<file>). Another audit agent may be running against the same server: that is fine, use unique
   capture names prefixed `imgaudit_`.
2. Open your own Chrome tab on http://localhost:8642/index.html.
3. MANDATORY, FIRST thing after every page load, before ANY other interaction: de-throttle the tab.
   Your tab is almost always occluded (another agent's window or any other window covers it), and Chrome
   treats an occluded tab as background: setTimeout is clamped to ~once a minute and requestAnimationFrame
   NEVER fires. Without this patch, analyses and exports appear to hang forever and you will burn time
   misdiagnosing app bugs that are not there (this cost the 2026-08-20 run a lot of wall-clock).
   ```js
   const _st = window.setTimeout;
   window.setTimeout = (fn, d, ...a) => _st(fn, Math.min(d || 0, 50), ...a);
   window.requestAnimationFrame = (cb) => _st(() => cb(performance.now()), 16);
   ```
   If something still seems hung for minutes, suspect a slow-by-nature export (whole-page screenshots take
   30-60 s, AI exports 60-90 s) before suspecting a bug.
4. Automated Chrome cannot complete downloads. Install this shim once per page load, BEFORE exporting.
   It captures BOTH blob downloads and data:-URI downloads (PNG exports go out as data: URIs; the older
   record-only shim could not save them):
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
     const c = window._caps.at(idx);
     if (!c?.blob) return c || null;
     await fetch('/save?name=' + encodeURIComponent(name), {method: 'POST', body: c.blob});
     return {type: c.type, size: c.size, download: c.download};
   }
   ```
   Use ONE unique name per capture (the server appends; reusing a name corrupts the file). One export can
   push SEVERAL caps (intermediate SVG blob + final PNG): pick the cap whose `type` matches the format you
   asked for, not blindly the last one.
5. INSPECT captured images by Reading them (you can view PNG/JPG; SVG is text; PDF via Read pages). For TIFF
   check the magic bytes (II*/MM*) and a plausible size.

## What to audit (drive the real UI: open the dialog, set the controls, click the export button)

The shared "Export image..." dialog (per-view button) with variations, for EACH of these views:
1. Correlation scatter (open TP53 vs MDM2; also with a Skin lineage filter + a few highlighted lines)
2. Gene effect chart (a gene, tissue view; and a filtered variant)
3. Expression-correlate scatter (via mutation analysis > expression correlates), if it has Export image
4. Correlation analysis chart
5. Mutation analysis results (its image export)
6. Network (run Test Genes; export the chart alone AND "the whole page" option; with a highlight set)
7. Heatmap (default open; AND a worked state: blocked by a toggled row + cluster row + gates + hidden group,
   so the export must carry strips, bands, dendrograms, gates row, legends)
8. Cell Line Browser "Export image..." (whole-browser capture)
9. UMAP/PCA plot export if present

Per view, vary across the runs (not full cross-product; cover each option at least twice overall):
- Format: PNG, SVG, TIFF, PDF, PowerPoint (whichever the dialog offers per view)
- Size/preset and DPI/resolution options
- Background: white vs transparent where offered
- With and without an active filter banner (the banner must appear in the export when filters are active)

## What "looks good" means (check each explicitly)

- File is valid (PNG signature + plausible dimensions for the chosen size/DPI; SVG parses as XML and contains
  text elements for labels; PDF has pages; PPTX is a zip containing ppt/slides/slide1.xml).
- Read the PNGs and LOOK at them: complete legends (not clipped - a historic Plotly bug class here), axis
  labels present, filter banner on top when filters are active, gene names italic in the network, heatmap
  carries group strip + annotation bands + gates row + both legends + dendrograms, background actually
  transparent (check alpha) or white as chosen.
- Version/branding: exports read the version badge where they stamp it - if a stamped version appears, it
  must say v.88.85 (or the current badge).
- SVG: text must be real text (not curves), no clipPath left over hiding the legend (historic bug),
  fonts named.
- Consistency: the same view exported PNG vs SVG shows the same content.
- Copy buttons: skip actual clipboard verification (blocked in automation), but click them once each and
  confirm no console error.

## Report

A findings table: view / format / settings / VERDICT (ok, cosmetic, broken) / evidence (captured filename).
Then a prose section for every non-ok finding: what you expected, what you got, reproduction steps.
Console errors anywhere = a finding. Save every captured file under the OUT dir with self-explanatory names
and list them. DO NOT fix anything. DO NOT commit. Close your tab when done; leave the server running.
