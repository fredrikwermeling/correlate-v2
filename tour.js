// A guided tour of Correlate: one page per feature, each with a real chart
// drawn from the loaded data and a button that opens that view for real.
// Reads window.app only. No network calls; the current page is kept in
// localStorage so the tour reopens where it was left.
(function () {
    'use strict';

    const A = () => window.app || null;
    const STEP_KEY = 'correlateTourStep';
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const phone = () => window.innerWidth <= 640;
    const hasPlotly = () => typeof window.Plotly !== 'undefined';

    // ------------------------------------------------------------ data layer
    const D = {
        ready() {
            const a = A();
            return !!(a && a.metadata && a.metadata.cellLines && a.geneEffects && a.geneIndex);
        },
        ids() { return A().metadata.cellLines; },
        name(id) { try { return A().getCellLineName(id) || id; } catch (e) { return id; } },
        lineage(id) { return A().cellLineMetadata?.lineage?.[id] || ''; },
        row(gene) {
            const a = A();
            const gi = a.geneIndex.get(String(gene).toUpperCase());
            if (gi === undefined) return null;
            return a.geneEffects.subarray(gi * a.nCellLines, (gi + 1) * a.nCellLines);
        },
        ge(row, i) {
            if (!row) return NaN;
            const v = row[i];
            return (!isFinite(v) || v <= -990) ? NaN : v;
        },
        expr(gene, i) {
            try { const v = A().getExpressionValueByGEIndex(gene, i); return isFinite(v) ? v : NaN; }
            catch (e) { return NaN; }
        },
        hotspot(gene, id) {
            const g = A().mutations?.geneData?.[gene];
            return g ? (g.mutations?.[id] || 0) : 0;
        },
        lineId(name) {
            try {
                const map = A()._buildCellLineNameToIdMap();
                return map.get(String(name).toUpperCase()) || map.get(String(name).toUpperCase().replace(/[^A-Z0-9]/g, '')) || null;
            } catch (e) { return null; }
        }
    };
    const median = (arr) => {
        if (!arr.length) return NaN;
        const s = arr.slice().sort((a, b) => a - b);
        const m = s.length >> 1;
        return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
    };
    function pearson(x, y) {
        let n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (let i = 0; i < x.length; i++) {
            const a = D.ge(x, i), b = D.ge(y, i);
            if (!isFinite(a) || !isFinite(b)) continue;
            n++; sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
        }
        if (n < 3) return null;
        const den = Math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
        if (!(den > 0)) return null;
        return { r: (n * sxy - sx * sy) / den, n };
    }
    function scatterXY(g1, g2) {
        const r1 = D.row(g1), r2 = D.row(g2);
        const xs = [], ys = [];
        for (let i = 0; i < D.ids().length; i++) {
            const x = D.ge(r1, i), y = D.ge(r2, i);
            if (isFinite(x) && isFinite(y)) { xs.push(x); ys.push(y); }
        }
        return { xs, ys };
    }
    function tissueGroups(gene, useExpr, minN) {
        const ids = D.ids();
        const row = useExpr ? null : D.row(gene);
        if (!useExpr && !row) return [];
        const by = new Map();
        for (let i = 0; i < ids.length; i++) {
            const lin = D.lineage(ids[i]);
            if (!lin) continue;
            const v = useExpr ? D.expr(gene, i) : D.ge(row, i);
            if (!isFinite(v)) continue;
            if (!by.has(lin)) by.set(lin, []);
            by.get(lin).push(v);
        }
        return [...by.entries()].filter(e => e[1].length >= (minN || 12))
            .map(e => ({ name: e[0], vals: e[1], med: median(e[1]) }));
    }

    // ------------------------------------------------------- chart styling
    // The same colors and wording the app's own popouts use.
    const CFG = { displayModeBar: false, responsive: true, displaylogo: false };
    const plotH = () => (phone() ? 240 : 300);
    function layout(extra) {
        return Object.assign({
            height: plotH(),
            margin: { t: 44, r: 16, b: 48, l: phone() ? 48 : 60 },
            paper_bgcolor: '#ffffff', plot_bgcolor: '#fafafa', showlegend: false,
            font: { family: 'Arial, Helvetica, sans-serif', size: phone() ? 11 : 12, color: '#374151' },
            hovermode: false
        }, extra || {});
    }
    const title = (text, sub) => ({
        text: `<b>${esc(text)}</b>` + (sub ? `<br><span style="font-size:10px;color:#6b7280;">${esc(sub)}</span>` : ''),
        font: { size: phone() ? 12 : 14 }, x: 0.5, xanchor: 'center'
    });
    const axis = (t) => ({
        title: { text: t, font: { size: phone() ? 11 : 12 }, standoff: 6 },
        zeroline: true, zerolinecolor: '#000', zerolinewidth: 2, tickfont: { size: phone() ? 10 : 11 }
    });
    function regressionTrace(xs, ys) {
        const n = xs.length;
        let sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
        const den = n * sxx - sx * sx;
        if (!(Math.abs(den) > 1e-9)) return null;
        const slope = (n * sxy - sx * sy) / den, b = (sy - slope * sx) / n;
        const lo = Math.min.apply(null, xs), hi = Math.max.apply(null, xs);
        return { x: [lo, hi], y: [slope * lo + b, slope * hi + b], mode: 'lines', type: 'scatter',
            line: { color: '#6ba544', width: 3 }, hoverinfo: 'skip' };
    }
    const boxTraces = (groups) => groups.map(gp => ({
        type: 'box', name: `${gp.name} (n=${gp.vals.length})`, x: gp.vals,
        boxpoints: 'all', jitter: 0.35, pointpos: 0,
        marker: { color: 'rgba(80,80,80,0.5)', size: 4 }, line: { color: '#374151' },
        fillcolor: 'rgba(200,200,200,0.3)', hoverinfo: 'skip'
    }));
    const dots = { color: '#9ca3af', size: 6, opacity: 0.6 };

    // Open a pair's scatter with the hotspot overlay, and optionally the
    // hotspot filter, already set: the popout applies the preset itself
    // just before its first draw.
    function openPairWithHotspot(a, g1, g2, mutGene, filterLevel) {
        a._inspectPreset = { hotspotGene: mutGene, hotspotMode: 'color', filterGene: filterLevel ? mutGene : null, filterLevel: filterLevel || null };
        a.openInspectByGenes(g1, g2);
    }

    const P53_SET = ['TP53', 'MDM2', 'MDM4', 'CDKN1A', 'PPM1D', 'USP7'];
    const setText = () => P53_SET.filter(g => A().geneIndex.has(g)).join('\n');

    // ---------------------------------------------------------------- pages
    // Each page: a title, a few plain sentences, an optional chart drawn from
    // the data, and the action that opens the same thing for real.
    const PAGES = [
        {
            title: 'Welcome',
            body: () => `<p>Correlate is built on DepMap, the Broad Institute's map of what ${num(D.ids().length)} cancer cell lines
                depend on. Every gene was knocked out with CRISPR in every line, and the app lets you ask what those results mean.</p>
                <p>This tour shows one example per feature. Each page has a real chart from the data and a button that opens
                that view in the app, so you can try the same thing with your own genes.</p>
                <ul>
                    <li><b>Gene effect</b>, the score behind everything here.</li>
                    <li><b>Correlations and gene set analysis</b>, which genes rise and fall together.</li>
                    <li><b>Mutations</b>, whether a dependency follows a mutation.</li>
                    <li><b>Cell lines</b>, finding them and reading their pages.</li>
                    <li><b>Expression, heatmaps and drug response.</b></li>
                </ul>`
        },
        {
            title: 'The gene effect score',
            need: () => !!D.row('SOX10'),
            body: () => `<p>A gene effect score says how much a cell line needs a gene. <b>0</b> means knocking the gene out changed
                nothing; <b>about -1</b> is a typical essential gene; more negative is a stronger need; above 0 means the cells grew better without it.</p>
                <p>The chart shows SOX10 across tissues. Skin lines sit far to the left: melanoma cannot live without SOX10, and almost nothing else cares.
                A dependency that follows the tissue like this usually marks the gene that holds that lineage's identity.</p>`,
            plot: (div) => {
                const groups = tissueGroups('SOX10', false, 14).sort((a, b) => b.med - a.med);
                const keep = groups.slice(0, 4).concat(groups.slice(-3));
                return Plotly.newPlot(div, boxTraces(keep.sort((a, b) => b.med - a.med)), layout({
                    title: title('SOX10 Gene Effect by tissue'),
                    xaxis: axis('SOX10 Gene Effect'),
                    yaxis: { automargin: true, tickfont: { size: phone() ? 9 : 10 } },
                    margin: { t: 40, r: 16, b: 46, l: 4 }
                }), CFG);
            },
            action: { label: 'Open SOX10 in the Gene Effect view', run: (a) => a.openGeneEffectModal('SOX10', 'tissue') }
        },
        {
            title: 'A correlation between two genes',
            need: () => !!D.row('TP53') && !!D.row('MDM2'),
            body: () => {
                const st = pearson(D.row('TP53'), D.row('MDM2'));
                return `<p>Each dot is one cell line, placed by its TP53 score and its MDM2 score. The green line is the trend.
                    Here r = ${st ? st.r.toFixed(2) : '?'} across ${st ? num(st.n) : '?'} lines: the two run opposite ways.</p>
                    <p>The reason is biology. A line with working p53 needs MDM2 to keep p53 in check, so it dies without MDM2, while a line that has
                    already lost p53 does not care. That is the kind of relationship a correlation across many cell lines can reveal.</p>
                    <p>In the app you can open any gene pair this way, filter to one tissue, color the dots, and label the cell lines.</p>`;
            },
            plot: (div) => {
                const { xs, ys } = scatterXY('TP53', 'MDM2');
                const traces = [{ x: xs, y: ys, mode: 'markers', type: 'scatter', marker: dots, hoverinfo: 'skip' }];
                const rl = regressionTrace(xs, ys); if (rl) traces.push(rl);
                return Plotly.newPlot(div, traces, layout({
                    title: title('TP53 vs MDM2', `n = ${num(xs.length)} cell lines`),
                    xaxis: axis('TP53 Gene Effect'), yaxis: axis('MDM2 Gene Effect')
                }), CFG);
            },
            action: { label: 'Open TP53 vs MDM2 in the app', run: (a) => a.openInspectByGenes('TP53', 'MDM2') }
        },
        {
            title: 'The same pair, split by TP53 mutation',
            need: () => !!D.row('TP53') && !!D.row('MDM2') && !!A().mutations?.geneData?.TP53,
            body: () => `<p>Here the dots are colored by TP53 mutation status: grey lines are wild-type, blue carry a hotspot mutation on one copy,
                red on both. The pattern from the last page falls into two groups. The wild-type lines are the ones far down the MDM2 axis:
                they still have working p53 and need MDM2 to hold it back. The mutated lines sit near zero on both axes: with p53 already gone, MDM2 no longer matters.</p>
                <p>In any scatter you can color the dots this way with the Hotspot overlay, or go one step further and keep only the wild-type or only the mutated lines
                with the Hotspot filter, which recomputes the correlation on that group alone.</p>`,
            plot: (div) => {
                const ids = D.ids(), rx = D.row('TP53'), ry = D.row('MDM2');
                const g = { wt: { x: [], y: [] }, m1: { x: [], y: [] }, m2: { x: [], y: [] } };
                for (let i = 0; i < ids.length; i++) {
                    const x = D.ge(rx, i), y = D.ge(ry, i);
                    if (!isFinite(x) || !isFinite(y)) continue;
                    const lvl = D.hotspot('TP53', ids[i]);
                    const b = lvl >= 2 ? g.m2 : lvl === 1 ? g.m1 : g.wt;
                    b.x.push(x); b.y.push(y);
                }
                return Plotly.newPlot(div, [
                    { x: g.wt.x, y: g.wt.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: dots, name: `WT (n=${g.wt.x.length})` },
                    { x: g.m1.x, y: g.m1.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#3b82f6', size: 7, opacity: 0.8 }, name: `1 mut (n=${g.m1.x.length})` },
                    { x: g.m2.x, y: g.m2.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#dc2626', size: 7, opacity: 0.85 }, name: `2 mut (n=${g.m2.x.length})` }
                ], layout({
                    title: title('TP53 vs MDM2', 'colored by TP53 hotspot mutation'),
                    xaxis: axis('TP53 Gene Effect'), yaxis: axis('MDM2 Gene Effect'),
                    showlegend: true, legend: { x: 1, y: 0, xanchor: 'right', yanchor: 'bottom', bgcolor: 'rgba(255,255,255,0.85)', font: { size: phone() ? 9 : 10 } }
                }), CFG);
            },
            actions: [
                { label: 'Open with the TP53 overlay', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', null) },
                { label: 'Open only the TP53 wild-type lines', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', '0') },
                { label: 'Open only the TP53 mutated lines', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', '1+2') }
            ]
        },
        {
            title: 'Gene set analysis',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4,
            body: () => `<p>This is the app's main tool. Paste a set of genes, and the app correlates their gene effect profiles across all
                cell lines and draws the pairs above your cutoff as a network. Genes in one complex or pathway usually end up linked.</p>
                <p>The grid shows the p53 pathway genes against each other. Red is a positive correlation, blue is negative. MDM2 and MDM4, the two brakes on p53,
                move together, and both run opposite to TP53 itself.</p>
                <p>The button runs this exact set in the app. Lower the correlation cutoff and run again to see weaker links appear.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const z = genes.map(g1 => genes.map(g2 => g1 === g2 ? 1 : (pearson(D.row(g1), D.row(g2))?.r ?? 0)));
                return Plotly.newPlot(div, [{
                    type: 'heatmap', z, x: genes, y: genes, zmin: -1, zmax: 1,
                    colorscale: [[0, '#2166ac'], [0.5, '#f7f7f7'], [1, '#b2182b']],
                    colorbar: { thickness: 10, len: 0.9, tickfont: { size: 9 } }, hoverinfo: 'skip',
                    text: z.map(r => r.map(v => v.toFixed(2))), texttemplate: '%{text}', textfont: { size: phone() ? 9 : 11 }
                }], layout({
                    title: title('Correlation of gene effects, p53 pathway'),
                    xaxis: { tickfont: { size: phone() ? 9 : 11 } }, yaxis: { tickfont: { size: phone() ? 9 : 11 }, autorange: 'reversed' },
                    margin: { t: 44, r: 16, b: 40, l: phone() ? 54 : 64 }
                }), CFG);
            },
            action: {
                label: 'Run this gene set in the app',
                run: () => {
                    document.getElementById('modeGeneSetBtn')?.click();
                    const ta = document.getElementById('geneTextarea');
                    if (ta) { ta.value = setText(); ta.dispatchEvent(new Event('input', { bubbles: true })); }
                    setTimeout(() => document.getElementById('runAnalysis')?.click(), 150);
                }
            }
        },
        {
            title: 'Coloring a scatter by a mutation',
            need: () => !!D.row('BRAF') && !!D.row('MAPK1') && !!A().mutations?.geneData?.BRAF,
            body: () => `<p>The same scatter, now colored by BRAF mutation status: grey lines are wild-type, blue carry a hotspot mutation on one copy,
                red on both. The mutated lines sit low on both axes: they depend on BRAF and on MAPK1, the kinase at the end of the pathway BRAF drives.</p>
                <p>That is what a targeted drug is built on. A dependency that follows a mutation tells you the mutation is doing the driving.
                In any correlation view, pick a gene under Hotspot overlay to color the dots this way.</p>`,
            plot: (div) => {
                const ids = D.ids(), rx = D.row('BRAF'), ry = D.row('MAPK1');
                const g = { wt: { x: [], y: [] }, m1: { x: [], y: [] }, m2: { x: [], y: [] } };
                for (let i = 0; i < ids.length; i++) {
                    const x = D.ge(rx, i), y = D.ge(ry, i);
                    if (!isFinite(x) || !isFinite(y)) continue;
                    const lvl = D.hotspot('BRAF', ids[i]);
                    const b = lvl >= 2 ? g.m2 : lvl === 1 ? g.m1 : g.wt;
                    b.x.push(x); b.y.push(y);
                }
                return Plotly.newPlot(div, [
                    { x: g.wt.x, y: g.wt.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: dots },
                    { x: g.m1.x, y: g.m1.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#3b82f6', size: 7, opacity: 0.8 } },
                    { x: g.m2.x, y: g.m2.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#dc2626', size: 7, opacity: 0.85 } }
                ], layout({
                    title: title('BRAF vs MAPK1', 'colored by BRAF hotspot mutation'),
                    xaxis: axis('BRAF Gene Effect'), yaxis: axis('MAPK1 Gene Effect')
                }), CFG);
            },
            action: { label: 'Open this scatter with the overlay', run: (a) => openPairWithHotspot(a, 'BRAF', 'MAPK1', 'BRAF', null) }
        },
        {
            title: 'Mutation analysis',
            need: () => !!D.row('BRAF') && !!A().mutations?.geneData?.BRAF,
            body: () => `<p>Mutation analysis splits the panel by one mutation and asks which dependencies differ between the two groups.
                The chart is the simplest case: BRAF's own gene effect, by BRAF mutation status. Each dot is a cell line.</p>
                <p>Wild-type lines sit near 0. Mutated lines sit far to the left, because a mutant oncogene keeps its pathway switched on and the cells
                then cannot do without it. The analysis ranks every gene this way and lists the ones that differ most.</p>`,
            plot: (div) => {
                const ids = D.ids(), row = D.row('BRAF');
                const wt = [], m1 = [], m2 = [];
                for (let i = 0; i < ids.length; i++) {
                    const v = D.ge(row, i); if (!isFinite(v)) continue;
                    const lvl = D.hotspot('BRAF', ids[i]);
                    (lvl >= 2 ? m2 : lvl === 1 ? m1 : wt).push(v);
                }
                const jit = (base, n) => Array.from({ length: n }, () => base + (Math.random() - 0.5) * 0.5);
                return Plotly.newPlot(div, [
                    { x: wt, y: jit(0, wt.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#888888', size: 5, opacity: 0.7 } },
                    { x: m1, y: jit(1, m1.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#3b82f6', size: 6, opacity: 0.8 } },
                    { x: m2, y: jit(2, m2.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#dc2626', size: 6, opacity: 0.85 } }
                ], layout({
                    title: title('BRAF Gene Effect by BRAF mutation status'),
                    xaxis: axis('BRAF Gene Effect'),
                    yaxis: { tickmode: 'array', tickvals: [0, 1, 2], range: [-0.6, 2.6], automargin: true,
                        ticktext: [`WT (n=${wt.length})`, `1 mut (n=${m1.length})`, `2 mut (n=${m2.length})`], tickfont: { size: phone() ? 9 : 10 } },
                    margin: { t: 40, r: 16, b: 46, l: 4 }
                }), CFG);
            },
            action: {
                label: 'Run a BRAF mutation analysis',
                run: (a) => {
                    document.getElementById('modeMutationBtn')?.click();
                    setTimeout(() => {
                        const sel = document.getElementById('mutationHotspotSelect');
                        if (sel) {
                            if (![...sel.options].some(o => o.value === 'BRAF')) sel.add(new Option('BRAF', 'BRAF'));
                            sel.value = 'BRAF'; sel.dispatchEvent(new Event('change'));
                        }
                        setTimeout(() => a.runMutationAnalysis?.(), 300);
                    }, 300);
                }
            }
        },
        {
            title: 'The Cell Line Browser',
            body: () => `<p>The browser is where you find the cell lines that fit a project. Filter by tissue, subtype, disease, sex, a mutation, a fusion
                or a copy-number change, tick the lines you want, and send them on to a heatmap, an export, or a comparison of what they depend on.</p>
                <p>The chart shows how the panel is spread across tissues. Lung and blood cancers are the largest groups; some tissues have only a handful of lines,
                which matters when you compare groups.</p>`,
            plot: (div) => {
                const counts = new Map();
                D.ids().forEach(id => { const l = D.lineage(id); if (l) counts.set(l, (counts.get(l) || 0) + 1); });
                const rows = [...counts.entries()].sort((a, b) => a[1] - b[1]).slice(-12);
                return Plotly.newPlot(div, [{
                    type: 'bar', orientation: 'h', x: rows.map(r => r[1]), y: rows.map(r => r[0]),
                    marker: { color: 'rgba(122, 185, 80, 0.85)' }, hoverinfo: 'skip'
                }], layout({
                    title: title('Cell lines per tissue', 'the twelve largest groups'),
                    xaxis: { title: { text: 'Cell lines', font: { size: phone() ? 11 : 12 }, standoff: 6 }, tickfont: { size: phone() ? 10 : 11 } },
                    yaxis: { automargin: true, tickfont: { size: phone() ? 9 : 10 } },
                    margin: { t: 46, r: 16, b: 46, l: 4 }
                }), CFG);
            },
            action: { label: 'Open the Cell Line Browser', run: (a) => a.openCellLineBrowser() }
        },
        {
            title: 'A cell line\'s page',
            need: () => !!D.lineId('A375') && !!A().globalSignatures?.byCellLine,
            body: () => `<p>Every cell line has a page with everything the app knows about it: where it came from, its drivers, copy number, fusions,
                what it depends on, what it expresses, how it responds to drugs, and how to check that a stock really is that line.</p>
                <p>Charts like this one place the line among all the others. Here the bars are every cell line's ploidy, the average number of chromosome
                copies, and the red line is A375. Around 2 is a normal set; higher means the genome has been doubled at some point.</p>`,
            plot: (div) => {
                const id = D.lineId('A375'), sig = A().globalSignatures.byCellLine;
                const vals = []; D.ids().forEach(x => { const v = sig[x]?.Ploidy; if (typeof v === 'number') vals.push(v); });
                const mine = sig[id]?.Ploidy;
                const lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
                return Plotly.newPlot(div, [{
                    type: 'histogram', x: vals, marker: { color: '#9ca3af', line: { color: '#ffffff', width: 1 } },
                    xbins: { start: lo, end: hi, size: (hi - lo) / 30 }, hoverinfo: 'skip'
                }], layout({
                    title: title('Ploidy across the panel', 'A375 marked in red'),
                    xaxis: { title: { text: 'Ploidy', font: { size: phone() ? 11 : 12 }, standoff: 6 }, tickfont: { size: phone() ? 10 : 11 },
                        showgrid: false, zeroline: false, showline: true, linecolor: '#d1d5db' },
                    yaxis: { showgrid: false, zeroline: false, showticklabels: false, showline: false }, bargap: 0.15,
                    shapes: typeof mine === 'number' ? [{ type: 'line', x0: mine, x1: mine, y0: 0, y1: 1, yref: 'paper', line: { color: '#dc2626', width: 2 } }] : [],
                    margin: { t: 46, r: 16, b: 46, l: 20 }
                }), CFG);
            },
            action: { label: 'Open the A375 page', run: (a) => a.openCellLineWiki(D.lineId('A375')) }
        },
        {
            title: 'Expression',
            need: () => !!A().loadExpressionData,
            body: () => `<p>Alongside the CRISPR screen, the app carries each line's mRNA levels. Values are log2(TPM+1), so a step of 1 is a doubling and
                anything above 1 is clearly expressed. Every gene view and every scatter can be switched from gene effect to mRNA.</p>
                <p>The chart shows MITF, the melanocyte lineage factor: skin lines express it far above the rest. Comparing this with SOX10's dependency
                on the second page shows how expression and dependency tell the same story from two sides.</p>`,
            plot: async (div) => {
                if (!A().expressionLoaded) {
                    div.innerHTML = '<div class="tour-loading">Loading expression data...</div>';
                    try { await A().loadExpressionData(); } catch (e) { }
                    if (!A().expressionLoaded) { div.innerHTML = '<div class="tour-loading">Expression data is not available right now.</div>'; return; }
                    div.innerHTML = '';
                }
                const groups = tissueGroups('MITF', true, 14).sort((a, b) => a.med - b.med);
                const keep = groups.slice(0, 3).concat(groups.slice(-4));
                return Plotly.newPlot(div, boxTraces(keep.sort((a, b) => a.med - b.med)), layout({
                    title: title('MITF expression by tissue'),
                    xaxis: { title: { text: 'MITF mRNA, log2(TPM+1)', font: { size: phone() ? 11 : 12 }, standoff: 6 }, tickfont: { size: phone() ? 10 : 11 }, zeroline: false },
                    yaxis: { automargin: true, tickfont: { size: phone() ? 9 : 10 } },
                    margin: { t: 40, r: 16, b: 46, l: 4 }
                }), CFG);
            },
            action: { label: 'Open MITF expression by tissue', run: (a) => a.openGeneEffectModal('MITF', 'tissue', { dataType: 'expr' }) }
        },
        {
            title: 'The gene set heatmap',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4 && typeof A()._hmOpenModal === 'function',
            body: () => `<p>The heatmap is the wide view: many genes across many cell lines at once. Each row is a gene, each column a cell line, and the
                color is that line's score for that gene, z-scored against all lines so genes with different ranges share one scale.</p>
                <p>Here the p53 pathway genes are shown across the skin lines. Orange is a dependency, purple means the knockout helped growth.
                In the app you can group and sort the columns by tissue, mutation or any other annotation, and cluster them.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const ids = D.ids();
                const stats = genes.map(g => {
                    const row = D.row(g); let n = 0, s = 0, ss = 0;
                    for (let i = 0; i < ids.length; i++) { const v = D.ge(row, i); if (isFinite(v)) { n++; s += v; ss += v * v; } }
                    const m = s / n, sd = Math.sqrt(Math.max(ss / n - m * m, 1e-9));
                    return { row, m, sd };
                });
                const cols = [];
                for (let i = 0; i < ids.length; i++) if (D.lineage(ids[i]) === 'Skin') cols.push(i);
                const use = cols.slice(0, 45);
                const z = stats.map(st => use.map(i => { const v = D.ge(st.row, i); return isFinite(v) ? (v - st.m) / st.sd : null; }));
                return Plotly.newPlot(div, [{
                    type: 'heatmap', z, y: genes, x: use.map(i => D.name(ids[i])), zmin: -2.5, zmax: 2.5,
                    colorscale: [[0, '#e66101'], [0.5, '#f7f7f7'], [1, '#5e3c99']], hoverinfo: 'skip',
                    colorbar: { thickness: 10, len: 0.9, tickfont: { size: 9 }, title: { text: 'z', font: { size: 10 } } }
                }], layout({
                    title: title('p53 pathway genes across skin lines', 'gene effect, z-scored per gene'),
                    xaxis: { showticklabels: false, title: { text: `${use.length} skin cell lines`, font: { size: phone() ? 11 : 12 }, standoff: 6 } },
                    yaxis: { tickfont: { size: phone() ? 9 : 11 }, autorange: 'reversed' },
                    margin: { t: 44, r: 16, b: 40, l: phone() ? 54 : 64 }
                }), CFG);
            },
            action: {
                label: 'Open this set in the heatmap',
                run: (a) => {
                    a._hmOpenModal();
                    setTimeout(() => {
                        const p = document.getElementById('hmPreset'), g = document.getElementById('hmGenes');
                        if (p && [...p.options].some(o => o.value === 'custom')) p.value = 'custom';
                        if (g) g.value = setText();
                        try { a._hmSyncPresetUI?.(); a._hmRedraw?.(); } catch (e) { }
                    }, 400);
                }
            }
        },
        {
            title: 'Drug response',
            need: () => !!D.lineId('A375') && Array.isArray(A().drugResponse?.compounds)
                && A().drugResponse.compounds.some(c => /vemurafenib/i.test(c.name || '')),
            body: () => `<p>Each cell line's page also carries how it responded to a panel of drugs in the PRISM screen. The score is an AUC from 0 to 1:
                0 means the drug killed the cells at every dose, 1 means they grew as if untreated.</p>
                <p>The chart shows vemurafenib, the BRAF V600E inhibitor, across the panel, with A375 marked. A375 carries BRAF V600E and sits among the
                most sensitive lines, which is exactly what the mutation predicts. Killing in a dish is not clinical response, but it is where the story starts.</p>`,
            plot: (div) => {
                const id = D.lineId('A375');
                const cp = A().drugResponse.compounds.find(c => /vemurafenib/i.test(c.name || ''));
                const vals = []; D.ids().forEach(x => { const v = cp.auc?.[x]; if (typeof v === 'number') vals.push(v); });
                const mine = cp.auc?.[id];
                const lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
                return Plotly.newPlot(div, [{
                    type: 'histogram', x: vals, marker: { color: '#9ca3af', line: { color: '#ffffff', width: 1 } },
                    xbins: { start: lo, end: hi, size: (hi - lo) / 30 }, hoverinfo: 'skip'
                }], layout({
                    title: title('Vemurafenib response across the panel', 'A375 marked in red'),
                    xaxis: { title: { text: 'AUC (0 = all cells killed, 1 = no effect)', font: { size: phone() ? 10 : 12 }, standoff: 6 }, tickfont: { size: phone() ? 10 : 11 },
                        showgrid: false, zeroline: false, showline: true, linecolor: '#d1d5db' },
                    yaxis: { showgrid: false, zeroline: false, showticklabels: false, showline: false }, bargap: 0.15,
                    shapes: typeof mine === 'number' ? [{ type: 'line', x0: mine, x1: mine, y0: 0, y1: 1, yref: 'paper', line: { color: '#dc2626', width: 2 } }] : [],
                    margin: { t: 46, r: 16, b: 50, l: 20 }
                }), CFG);
            },
            action: { label: 'Open the A375 page', run: (a) => a.openCellLineWiki(D.lineId('A375')) }
        },
        {
            title: 'Where to go next',
            body: () => `<p>That is the whole app in outline. Three good ways to start:</p>
                <ul>
                    <li><b>Run the p53 example</b> as a gene set analysis, then swap in a set of your own.</li>
                    <li><b>Open the Cell Line Browser</b>, filter to your tissue, and open a few pages.</li>
                    <li><b>Open a correlation</b> for two genes you know, and color it by a mutation.</li>
                </ul>
                <p>Every popout has an Export button for figures and data, and Export for AI packages the view for a language model to read.
                The How it works page has the longer explanation of the statistics.</p>`,
            actions: [
                { label: 'Run the p53 example', run: () => PAGES.find(p => p.title === 'Gene set analysis').action.run(A()) },
                { label: 'Open the Cell Line Browser', run: (a) => a.openCellLineBrowser() },
                { label: 'Open TP53 vs MDM2', run: (a) => a.openInspectByGenes('TP53', 'MDM2') }
            ]
        }
    ];

    // ------------------------------------------------------------------ UI
    const CSS = `
#tourModal .modal { max-width: 720px; }
#tourModal .modal-body { overflow-y: auto; padding: 14px 20px; }
#tourModal .tour-step { font-size: 11px; color: #6b7280; font-weight: 400; margin-left: 10px; }
#tourModal .tour-text { font-size: 13px; line-height: 1.6; color: #374151; }
#tourModal .tour-text p { margin: 0 0 10px; }
#tourModal .tour-text ul { margin: 0 0 10px 18px; padding: 0; }
#tourModal .tour-text li { margin-bottom: 4px; }
#tourModal .tour-chart { min-height: 120px; margin: 4px 0 12px; border: 1px solid var(--gray-200); border-radius: 6px; overflow: hidden; }
#tourModal .tour-loading { padding: 40px 12px; text-align: center; color: #6b7280; font-size: 12px; }
#tourModal .tour-actions { display: flex; flex-wrap: wrap; gap: 8px; margin: 4px 0 6px; }
#tourModal .tour-actions .btn { min-height: 40px; }
#tourModal .modal-footer { display: flex; align-items: center; justify-content: space-between; gap: 10px; padding: 10px 20px; border-top: 1px solid var(--gray-200); flex-shrink: 0; }
#tourModal .tour-dots { display: flex; gap: 5px; flex-wrap: wrap; justify-content: center; }
#tourModal .tour-dot { width: 8px; height: 8px; border-radius: 50%; background: #d1d5db; border: none; padding: 0; cursor: pointer; }
#tourModal .tour-dot.on { background: var(--green-600); }
#tourModal .tour-dot.done { background: #9ecf82; }
#tourModal .tour-nav { min-width: 84px; min-height: 40px; }
@media (max-width: 640px) {
  #tourModal .tour-text { font-size: 14px; }
  #tourModal .tour-nav { min-width: 72px; min-height: 44px; font-size: 14px; }
  #tourModal .tour-actions .btn { flex: 1 1 100%; min-height: 44px; font-size: 14px; }
}`;

    let pages = [], step = 0, built = false;

    function build() {
        if (built) return;
        built = true;
        const style = document.createElement('style');
        style.textContent = CSS;
        document.head.appendChild(style);
        const wrap = document.createElement('div');
        wrap.className = 'modal-overlay';
        wrap.id = 'tourModal';
        wrap.style.zIndex = '1500';
        wrap.innerHTML = `
            <div class="modal" role="dialog" aria-labelledby="tourTitle">
                <div class="modal-header">
                    <h3 id="tourTitle" style="margin:0;"></h3>
                    <button type="button" class="modal-close" id="tourClose" aria-label="Close the tour">&times;</button>
                </div>
                <div class="modal-body">
                    <div class="tour-text" id="tourText"></div>
                    <div class="tour-chart" id="tourChart" style="display:none;"></div>
                    <div class="tour-actions" id="tourActions"></div>
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-outline btn-sm tour-nav" id="tourBack">Back</button>
                    <div class="tour-dots" id="tourDots"></div>
                    <button type="button" class="btn btn-success btn-sm tour-nav" id="tourNext">Next</button>
                </div>
            </div>`;
        document.body.appendChild(wrap);
        wrap.addEventListener('click', (e) => { if (e.target === wrap) close(); });
        document.getElementById('tourClose').addEventListener('click', close);
        document.getElementById('tourBack').addEventListener('click', () => show(step - 1));
        document.getElementById('tourNext').addEventListener('click', () => {
            if (step >= pages.length - 1) { close(); saveStep(0); return; }
            show(step + 1);
        });
        document.addEventListener('keydown', (e) => {
            if (!isOpen()) return;
            if (e.key === 'Escape') close();
            else if (e.key === 'ArrowRight') { e.preventDefault(); if (step < pages.length - 1) show(step + 1); }
            else if (e.key === 'ArrowLeft') { e.preventDefault(); show(step - 1); }
        });
    }

    const isOpen = () => document.getElementById('tourModal')?.classList.contains('active');
    function saveStep(i) { try { localStorage.setItem(STEP_KEY, String(i)); } catch (e) { } }
    function loadStep() { try { return parseInt(localStorage.getItem(STEP_KEY), 10) || 0; } catch (e) { return 0; } }

    function purgeChart() {
        const div = document.getElementById('tourChart');
        if (div && hasPlotly()) { try { Plotly.purge(div); } catch (e) { } }
        if (div) { div.innerHTML = ''; div.style.display = 'none'; }
    }

    function actionButton(act) {
        const b = document.createElement('button');
        b.type = 'button';
        b.className = 'btn btn-outline btn-sm';
        b.textContent = act.label;
        b.addEventListener('click', () => {
            close();
            try { act.run(A()); } catch (e) { console.warn('Tour action failed:', e); }
        });
        return b;
    }

    function show(i) {
        step = Math.max(0, Math.min(pages.length - 1, i));
        saveStep(step);
        const p = pages[step];
        document.getElementById('tourTitle').innerHTML = `Tour: ${esc(p.title)} <span class="tour-step">${step + 1} of ${pages.length}</span>`;
        document.getElementById('tourText').innerHTML = p.body();
        purgeChart();
        const chart = document.getElementById('tourChart');
        if (p.plot && hasPlotly()) {
            chart.style.display = 'block';
            chart.innerHTML = '<div class="tour-loading">Drawing...</div>';
            const mine = step;
            Promise.resolve().then(() => {
                if (mine !== step) return;
                chart.innerHTML = '';
                return p.plot(chart);
            }).catch(e => { console.warn('Tour chart failed:', e); chart.innerHTML = '<div class="tour-loading">This chart could not be drawn.</div>'; });
        }
        const acts = document.getElementById('tourActions');
        acts.innerHTML = '';
        (p.actions || (p.action ? [p.action] : [])).forEach(a => acts.appendChild(actionButton(a)));
        document.getElementById('tourBack').style.visibility = step === 0 ? 'hidden' : 'visible';
        document.getElementById('tourNext').textContent = step === pages.length - 1 ? 'Finish' : 'Next';
        const dots = document.getElementById('tourDots');
        dots.innerHTML = '';
        pages.forEach((pg, k) => {
            const d = document.createElement('button');
            d.type = 'button';
            d.className = 'tour-dot' + (k === step ? ' on' : k < step ? ' done' : '');
            d.title = pg.title;
            d.setAttribute('aria-label', `Go to page ${k + 1}, ${pg.title}`);
            d.addEventListener('click', () => show(k));
            dots.appendChild(d);
        });
        document.querySelector('#tourModal .modal-body').scrollTop = 0;
    }

    function open() {
        if (!D.ready()) {
            A()?.showCopyNotification?.('The tour opens once the data has loaded.');
            return;
        }
        build();
        pages = PAGES.filter(p => { try { return !p.need || p.need(); } catch (e) { return false; } });
        document.getElementById('tourModal').classList.add('active');
        const saved = loadStep();
        show(saved > 0 && saved < pages.length ? saved : 0);
    }
    function close() {
        purgeChart();
        document.getElementById('tourModal')?.classList.remove('active');
    }

    window.CorrelateTour = { open, close, isOpen };
})();
