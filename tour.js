// Tour of the app: one scrolling page that presents Correlate the way a
// product page would, with real screenshots of the app (tour_img/*.webp,
// captured from the live views) beside short feature text, and a button
// under each that opens the same view with the same example loaded.
// Reads window.app only. Counts come from the loaded data where available.
(function () {
    'use strict';

    const A = () => window.app || null;
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const gi = (g) => `<i>${esc(g)}</i>`;
    const IMG = 'tour_img/';
    // The pictures are cached by URL; the script's own ?v= tag renews them.
    const IMG_V = (() => { const m = (document.currentScript?.src || [...document.scripts].map(x => x.src).find(x => /tour\.js/.test(x)) || '').match(/[?&]v=([^&]+)/); return m ? '?v=' + m[1] : ''; })();

    // ------------------------------------------------------------ data
    const D = {
        nLines() { try { return A().metadata.cellLines.length; } catch (e) { return 1208; } },
        nGenes() { try { return A().nGenes || A().metadata.genes.length; } catch (e) { return 18531; } },
        nExprLines() { try { return A().expressionMetadata?.nCellLines || 1718; } catch (e) { return 1718; } },
        lineId(name) {
            try {
                const map = A()._buildCellLineNameToIdMap();
                return map.get(String(name).toUpperCase()) || null;
            } catch (e) { return null; }
        },
    };

    // ------------------------------------------------------------ actions
    // Every action reproduces the state shown in the screenshot beside it.
    const HERO_SET = ['TP53', 'MDM2', 'MDM4', 'PPM1D', 'USP7', 'CDKN1A', 'TP53BP1', 'CHEK2', 'ATM', 'TSC1', 'TSC2'];
    const setVal = (id, v) => {
        const el = document.getElementById(id); if (!el) return;
        el.value = v;
        el.dispatchEvent(new Event('input', { bubbles: true }));
        el.dispatchEvent(new Event('change', { bubbles: true }));
    };
    function runExampleSet(o) {
        o = o || {};
        const a = A();
        document.getElementById('modeGeneSetBtn')?.click();
        const basis = document.querySelector(`input[name="analysisBasis"][value="${o.basis || 'ge'}"]`);
        if (basis && !basis.checked) { basis.checked = true; basis.dispatchEvent(new Event('change', { bubbles: true })); }
        if (o.cutoff) setVal('correlationCutoff', o.cutoff);
        const ta = document.getElementById('geneTextarea');
        if (ta) { ta.value = HERO_SET.filter(g => a.geneIndex.has(g)).join('\n'); ta.dispatchEvent(new Event('input', { bubbles: true })); }
        const before = a.results;
        setTimeout(() => document.getElementById('runAnalysis')?.click(), 150);
        if (!o.colorGE) return;
        let tries = 0;
        const arm = () => {
            if (a.results && a.results !== before && a.results.success) {
                const cb = document.getElementById('colorByGeneEffect');
                if (cb && !cb.checked) { cb.checked = true; cb.dispatchEvent(new Event('change', { bubbles: true })); }
                return;
            }
            if (tries++ < 80) setTimeout(arm, 250);
        };
        setTimeout(arm, 600);
    }
    const wait = (ms) => new Promise(r => setTimeout(r, ms));
    const tick = (id) => { const cb = document.getElementById(id); if (cb && !cb.checked) cb.click(); };
    const until = async (test, max) => { for (let i = 0; i < (max || 80); i++) { if (test()) return true; await wait(250); } return false; };
    // The "Test Genes" set with pasted screen statistics, nodes coloured by
    // LFC, genes without a correlation shown too.
    async function runStatsExample(a, withStats) {
        document.getElementById('modeGeneSetBtn')?.click();
        [...document.querySelectorAll('.input-tab')].find(b => /with stats/i.test(b.textContent))?.click();
        await wait(300);
        a.loadTestGenesWithStats();
        setVal('correlationCutoff', 0.4);
        const before = a.results;
        await wait(200);
        document.getElementById('runAnalysis')?.click();
        await until(() => a.results && a.results !== before && a.results.success);
        if (withStats) { tick('showUncorrelatedGenes'); await wait(800); tick('colorByStats'); }
    }
    async function openMelanoma(a, hotspotGene) {
        a.openCellLineBrowser(); await wait(500);
        setVal('clbTissueFilter', 'Skin'); await wait(500);
        setVal('clbSubtypeFilter', 'Melanoma'); await wait(500);
        setVal('clbHotspotFilter', hotspotGene || ''); await wait(700);
    }
    async function compareMelanoma(a) {
        await openMelanoma(a, 'BRAF');
        a._resetGEInspectScope();
        a._geInspectScope = 'all';
        a._runBusy({ button: 'clbInspectGEBtn', label: 'Comparing' }, () => a.inspectSelectionGE());
    }
    async function pcaMelanoma(a) {
        await openMelanoma(a, '');
        [...document.querySelectorAll('#clbModalCard button')].find(b => /PCA \/ UMAP/.test(b.textContent))?.click();
        await wait(800);
        setVal('clbUmapMethod', 'pca'); setVal('clbUmapDataType', 'expr');
        tick('clbUmapShowLoadings');
        const before = a._clbUmapData;
        document.getElementById('clbUmapRunBtn')?.click();
        await until(() => a._clbUmapData && a._clbUmapData !== before, 120);
        await wait(500);
        const cg = document.getElementById('clbUmapColorGene');
        if (cg) { cg.value = 'MITF'; cg.dispatchEvent(new Event('input', { bubbles: true })); }
        setVal('clbUmapColorGeneType', 'expr');
        document.getElementById('clbUmapColorGeneBtn')?.click();
    }
    async function heatmapMelanoma(a) {
        await openMelanoma(a, '');
        document.getElementById('clbHeatmapBtn')?.click();
        await wait(800);
        setVal('hmLineage', 'Skin'); await wait(400);
        setVal('hmSubtype', 'Melanoma'); await wait(400);
        const rows = a._hmAnnRows || [];
        const mut = rows.find(r => r.mode === 'hotspot');
        if (mut) { mut.gene = 'TP53'; mut.sortKey = 'cat-asc'; a._hmRenderAnnRowsBlock?.(); }
        setVal('hmPreset', 'tp53_members');
    }
    async function runMutationExample(a) {
        document.getElementById('modeMutationBtn')?.click();
        await wait(700);
        setVal('mutationHotspotSelect', 'BRAF');
        setVal('lineageFilter', 'Skin');
        await wait(500);
        const before = a.mutationResults;
        document.getElementById('runMutationAnalysisBtn')?.click();
        await until(() => a.mutationResults && a.mutationResults !== before, 120);
        await wait(400);
        a.showGeneEffectDistribution('BRAF');
    }
    const ACTIONS = {
        hero: () => runExampleSet({ cutoff: 0.5 }),
        heroStats: (a) => runStatsExample(a, true),
        heroMrna: () => runExampleSet({ cutoff: 0.3, basis: 'expr' }),
        pair: (a) => a.openInspectByGenes('TP53', 'MDM2'),
        eda2r: (a) => { a._inspectPreset = { colorByGene: 'EDA2R', colorByGeneType: 'expr' }; a.openInspectByGenes('TP53', 'MDM2'); },
        mutation: (a) => runMutationExample(a),
        compare: (a) => compareMelanoma(a),
        pca: (a) => pcaMelanoma(a),
        heatmap: (a) => heatmapMelanoma(a),
        myc: (a) => a.openGeneBrowser('MYC'),
        kras: (a) => a.openGeneBrowser('KRAS'),
        tsc1: (a) => a.openGeneBrowser('TSC1'),
        clb: (a) => openMelanoma(a, 'BRAF'),
        wiki: (a) => a.openCellLineWiki(D.lineId('A375') || 'ACH-000219'),
        stats: () => { const m = document.getElementById('infographicModal'); if (m) m.style.display = 'flex'; },
    };

    // ------------------------------------------------------------ page
    const shot = (file, alt, act, o) => {
        o = o || {};
        // A figure is a picture; the one way to open its view is the button
        // under it, so there is never a hidden click target and a button for
        // the same thing.
        return `<figure class="tour-shot${o.cls ? ' ' + o.cls : ''}">
            <div class="tour-shot-bar"><span></span><span></span><span></span></div>
            <img src="${IMG}${file}.webp${IMG_V}" alt="${esc(alt)}" loading="${o.eager ? 'eager' : 'lazy'}" decoding="async">
            ${o.cap ? `<figcaption>${o.cap}</figcaption>` : ''}
            ${act && o.btn ? `<div class="tour-shot-act">${tryBtn(act, o.btn)}</div>` : ''}
        </figure>`;
    };
    const tryBtn = (act, label) => `<button type="button" class="tour-try" data-act="${act}">${label}<span aria-hidden="true"> &rarr;</span></button>`;
    const checks = (items) => `<ul class="tour-checks">${items.map(t => `<li>${t}</li>`).join('')}</ul>`;

    function pageHTML() {
        const nL = num(D.nLines()), nE = num(D.nExprLines()), nG = num(D.nGenes()), nSet = HERO_SET.length;
        return `
<header class="tour-top">
    <div class="tour-wrap tour-top-in">
        <div class="tour-brand"><b>Correlate</b><span>Tour of the app</span></div>
        <nav class="tour-nav" aria-label="Sections">
            <a href="#tour-network">Network</a>
            <a href="#tour-corr">Correlation plot</a>
            <a href="#tour-cells">Cell Line Browser</a>
            <a href="#tour-gene">Gene Browser</a>
            <a href="#tour-more">More</a>
        </nav>
        <button type="button" class="tour-close" id="tourClose" aria-label="Close the tour">&times;</button>
    </div>
</header>

<section class="tour-hero">
    <div class="tour-wrap">
        <h1>Correlate brings together CRISPR knockout screens, mRNA expression, mutations, copy number and drug response from DepMap, adds gene and cell line information from many other databases, and connects them.</h1>
        <p class="tour-sub">Use it to:</p>
        <ul class="tour-hero-list">
            <li>Find connections within your own gene list, for example hits from a CRISPR screen or a transcriptomics experiment, from unbiased correlations of CRISPR gene effect or mRNA across the cell lines.</li>
            <li>Find cell lines with specific mutations and mutation patterns, and compare omics data between groups of cell lines.</li>
            <li>Get extensive information about any gene: dependency, expression, mutations, protein complexes, drugs and literature, gathered from DepMap, UniProt, Cellosaurus and other sources.</li>
            <li>Export any view for an AI and discuss the data with it.</li>
        </ul>
    </div>
    <div class="tour-wrap tour-hero-shot">
        ${shot('net_p53', 'A network built from eleven p53 pathway and TSC genes', 'hero', { eager: true, cap: `A network from a list of ${nSet} genes: the p53 pathway and the TSC complex. Blue edges join genes that cells need together, red edges join opposites: cells that need ${gi('MDM2')} tolerate losing ${gi('TP53')}.`, btn: 'Try it: build this network' })}
    </div>
</section>

<section class="tour-stats">
    <div class="tour-wrap tour-stats-in">
        <div><b>${nL}</b><span>cell lines with DepMap CRISPR gene effect</span></div>
        <div><b>${nE}</b><span>cell lines with DepMap mRNA expression</span></div>
        <div><b>${nG}</b><span>genes</span></div>
        <div><b>10+</b><span>wiki pages per cell line, from DepMap, Cellosaurus and other sources</span></div>
    </div>
</section>

<section class="tour-sec" id="tour-network">
    <div class="tour-wrap tour-feature">
        <div class="tour-text">
            <p class="tour-eyebrow">1 &middot; Network</p>
            <h2>Paste a gene set, get a network</h2>
            <p>Correlate compares every pair of genes in your set across ${nL} cell lines. Two genes that matter in the same cell lines are joined by a line, called an edge. Blue means they go together, red means they go opposite ways, thick means a strong connection. Genes that end up joined usually work together, in the same complex or pathway.</p>
            ${checks([
                'Paste your own screen results next to the gene list (LFC and FDR) and the network shows them: the nodes take the screen colours, and genes without a connection are shown too.',
                'Build the network from the CRISPR data or from the mRNA data. A badge on the network says which.',
                'Expand the set to find new partners, send a group to Enrichr, open the network full page.',
            ])}
        </div>
        <div class="tour-stack">
            ${shot('net_stats', 'A list of 22 cancer genes with screen statistics pasted alongside, nodes coloured by LFC', 'heroStats', { cap: 'A list of 22 cancer genes with results from a screen pasted alongside: the nodes take the screen colours (LFC), and every gene in the list is shown.', btn: 'Build it with the screen results' })}
            ${shot('net_mrna', 'The same eleven genes built on mRNA expression', 'heroMrna', { cap: `The ${nSet} genes from the top of the page built on mRNA expression instead: a different wiring, because genes that are expressed together are not the same as genes that are needed together.`, btn: 'Build this set on mRNA' })}
        </div>
    </div>
</section>

<section class="tour-sec tour-alt" id="tour-corr">
    <div class="tour-wrap">
        <div class="tour-text tour-text-wide">
            <p class="tour-eyebrow">2 &middot; Correlation plot</p>
            <h2>Double-click an edge to see the data behind it</h2>
            <p>One dot per cell line, one gene against the other, with r and p. You can also open the plot on its own: type any two genes.</p>
        </div>
        <div class="tour-wide-shot">
            ${shot('edge_to_plot', 'Double-clicking the edge between TP53 and MDM2 opens their correlation plot', 'pair', { cap: `The edge between ${gi('TP53')} and ${gi('MDM2')}, and the plot behind it: ${nL} cell lines, r −0.72. Cell lines that tolerate losing ${gi('TP53')} depend on ${gi('MDM2')}.`, btn: 'Open TP53 vs MDM2' })}
        </div>
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <h3>Then ask why</h3>
                ${checks([
                    'Filter to a cancer type or a mutation. r and p follow the filter.',
                    `Colour the dots by tissue, by a mutation, or by a third gene. <b>Suggest</b> finds the gene that explains the pattern best.`,
                    'Draw gates around two groups of cell lines and <b>Compare</b> them: genes, tissues and mutations that differ, in the same view the Cell Line Browser uses.',
                    '<b>Export for AI</b> sends the plot with the data behind it, so you can ask an AI what could explain the pattern.',
                ])}
            </div>
            ${shot('scatter_eda2r', 'TP53 against MDM2 gene effect, coloured by EDA2R mRNA', 'eda2r', { cap: `The same plot coloured by ${gi('EDA2R')} mRNA, a p53 target gene: the cell lines that depend on ${gi('MDM2')} are the ones with active p53.`, btn: 'Open it coloured by EDA2R' })}
        </div>
    </div>
</section>

<section class="tour-sec" id="tour-cells">
    <div class="tour-wrap">
        <div class="tour-feature">
            <div class="tour-text">
                <p class="tour-eyebrow">3 &middot; Cell Line Browser</p>
                <h2>Find the right cell lines</h2>
                <p>Narrow the ${nL} cell lines down by cancer type, mutation, fusion, copy number, sex, drug response, or how much they express or need any gene. Every cell line has its own page.</p>
                ${checks([
                    `Example: melanoma cell lines with mutated ${gi('BRAF')}. The count updates as you go.`,
                    'From the list: compare the selection with the rest, a PCA or UMAP, a heatmap of any gene set.',
                    'Download a cell line’s full tables, or the list on screen, as CSV.',
                ])}
            </div>
            ${shot('clb_mel', 'The Cell Line Browser filtered to melanoma cell lines with mutated BRAF, A375 selected', 'clb', { cap: `Melanoma cell lines with mutated ${gi('BRAF')}, A375 opened in the side panel.`, btn: 'Open the browser with these filters' })}
        </div>
        <div class="tour-two tour-two-tight">
            ${shot('volcano_mel', 'Volcano plot comparing the 50 BRAF-mutant melanoma cell lines with every other cell line', 'compare', { cap: `Compare: the 50 selected cell lines against the other 1,158. ${gi('BRAF')}, ${gi('SOX10')} and ${gi('MAPK1')} are their dependencies.`, btn: 'Compare these cell lines with the rest' })}
            ${shot('pca_mitf', 'PCA of the melanoma cell lines on mRNA, coloured by MITF expression, with loadings', 'pca', { cap: `PCA of the melanoma cell lines on mRNA, coloured by ${gi('MITF')} expression. The loadings name the genes behind the axes.`, btn: 'Open this PCA' })}
        </div>
        <div class="tour-wide-shot">
            ${shot('heatmap_tp53', 'Heatmap of the TP53 pathway members across the melanoma cell lines, with a TP53 mutation row', 'heatmap', { cap: `Heatmap: the ${gi('TP53')} pathway members across the melanoma cell lines, sorted by ${gi('TP53')} mutation. Wild-type cell lines depend on ${gi('MDM2')} and ${gi('USP7')}; mutant ones do not.`, btn: 'Open this heatmap' })}
        </div>
    </div>
</section>

<section class="tour-sec tour-alt" id="tour-gene">
    <div class="tour-wrap">
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <p class="tour-eyebrow">4 &middot; Gene Browser</p>
                <h2>One gene, one page</h2>
                <p>Type a gene. Its page says how much cells need it and which cancers need it most, how it is expressed, which genes follow it, and which complexes and pathways it belongs to. Save it as a PDF.</p>
                ${checks([
                    'A verdict for dependency and for expression, ranked against every gene.',
                    'Co-dependencies and co-expression, with complex partners ringed.',
                    'Mutations, fusions, copy number, drug response, structure and sequence.',
                ])}
            </div>
            ${shot('gb_myc', 'The Gene Browser page for MYC: description, chromosome position and records', 'myc', { cap: `${gi('MYC')} in the Gene Browser: description, chromosome position and the records it links to.`, btn: 'Open MYC in the Gene Browser' })}
        </div>
        <div class="tour-wide-shot">
            ${shot('gb_kras_panels', 'Dependency and expression sections for KRAS', 'kras', { cap: `Dependency and expression, here for ${gi('KRAS')}: where the gene sits among all genes, the most dependent tissues, and the highest and lowest expressing ones.`, btn: 'Open KRAS in the Gene Browser' })}
        </div>
        <div class="tour-feature">
            <div class="tour-text">
                <h3>Does the complex act as one unit?</h3>
                ${checks([
                    `Protein complexes from the EBI Complex Portal, each member tagged with its co-dependency.`,
                    `${gi('TSC2')} follows ${gi('TSC1')} at r 0.90. ${gi('TBC1D7')}, the third member, only at 0.25.`,
                ])}
            </div>
            <div class="tour-stack">
                ${shot('gb_codep', 'Co-dependencies of TSC1 as a table and a fan figure', 'tsc1', { cap: `Co-dependencies of ${gi('TSC1')}, complex partners ringed.`, btn: 'Open TSC1 in the Gene Browser' })}
                ${shot('gb_complex', 'The TSC1-TSC2 complex card with co-dependency per member', 'tsc1', { cap: `The ${gi('TSC1')}-${gi('TSC2')} complex card.` })}
            </div>
        </div>
    </div>
</section>

<section class="tour-sec" id="tour-more">
    <div class="tour-wrap">
        <p class="tour-eyebrow">More</p>
        <div class="tour-feature">
            <div class="tour-text">
                <h3>Mutation analysis</h3>
                <p>Wild-type against mutated for any gene, every gene ranked by the difference in gene effect or mRNA. Any gene in the table opens as a strip plot.</p>
            </div>
            ${shot('strip_braf', 'BRAF gene effect in skin cell lines by BRAF mutation status', 'mutation', { cap: `${gi('BRAF')} gene effect in skin cell lines, split by ${gi('BRAF')} mutation: the mutated cell lines depend on it, the wild-type ones do not.`, btn: 'Open BRAF mutation analysis in skin' })}
        </div>
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <h3>Take it with you</h3>
                <p>Under every view: Export image (PNG, SVG, PDF, PowerPoint at print size), Methods and CSV. An exported PNG or SVG reopens in Correlate with its analysis. Every figure on this page is such an export.</p>
                <p><b>Export for AI</b> writes the data behind the view into one file made for a language model. Upload it to ChatGPT, Claude or Gemini and discuss the results: which of these dependencies has a drug, which cell line fits my experiment, what links these genes?</p>
            </div>
            ${shot('export_dialog', 'The Export image dialog', null, { cls: 'tour-shot-dialog', cap: 'Export image: what to export, format, size in cm, resolution, background.' })}
        </div>
        <div class="tour-chips">
            <div><b>Cell line pages</b>ten sections per cell line: origin, drivers, genome, dependencies, drug response</div>
            <div><b>Compilation</b>a grid of correlation plots, one per pair</div>
            <div><b>Threshold scan</b>how many edges the set keeps as the cutoff moves</div>
            <div><b>Correlation matrix</b>of the set, on gene effect or mRNA</div>
            <div><b>Gates</b>draw a cohort on any plot and carry it on</div>
            <div><b>Green Listed</b>library genes and curated gene classes</div>
        </div>
        <div class="tour-cta tour-cta-center">${tryBtn('hero', 'Build the example network')}${tryBtn('wiki', 'Open the A375 page')}${tryBtn('stats', 'The statistics behind the app')}</div>
        <p class="tour-foot">Data: DepMap 26Q1 (Broad Institute), with EBI Complex Portal, Reactome, PRISM and the other sources listed under Reference Data. Correlate is developed in the Wermeling lab.</p>
    </div>
</section>`;
    }

    // ------------------------------------------------------------ style
    const CSS = `
#tourPage { position: fixed; inset: 0; z-index: 1500; overflow-y: auto; overflow-x: hidden; background: #fff; color: #111827; display: none; -webkit-overflow-scrolling: touch; }
#tourPage.active { display: block; }
#tourPage * { box-sizing: border-box; }
#tourPage .tour-wrap { max-width: 1120px; margin: 0 auto; padding: 0 24px; }
#tourPage .tour-top { position: sticky; top: 0; z-index: 5; background: rgba(255,255,255,.94); backdrop-filter: blur(8px); border-bottom: 1px solid #e5e7eb; }
#tourPage .tour-top-in { display: flex; align-items: center; gap: 18px; height: 54px; }
#tourPage .tour-brand { display: flex; align-items: baseline; gap: 10px; white-space: nowrap; }
#tourPage .tour-brand b { font-size: 20px; color: #3b5e23; font-weight: 800; letter-spacing: -.01em; }
#tourPage .tour-brand span { font-size: 13px; color: #6b7280; }
#tourPage .tour-nav { display: flex; gap: 2px; margin-left: auto; overflow-x: auto; scrollbar-width: none; }
#tourPage .tour-nav::-webkit-scrollbar { display: none; }
#tourPage .tour-nav a { font-size: 13.5px; color: #374151; text-decoration: none; padding: 7px 10px; border-radius: 6px; white-space: nowrap; }
#tourPage .tour-nav a:hover { background: #f3f4f6; }
#tourPage .tour-nav a.on { color: #3b5e23; background: #eef6e6; font-weight: 600; }
#tourPage .tour-close { flex: none; width: 38px; height: 38px; border: 1px solid #e5e7eb; background: #fff; border-radius: 50%; font-size: 24px; line-height: 1; color: #374151; cursor: pointer; }
#tourPage .tour-close:hover { background: #f3f4f6; }
#tourPage .tour-eyebrow { font-size: 12px; font-weight: 700; letter-spacing: .12em; text-transform: uppercase; color: #5d9239; margin: 0 0 10px; }
#tourPage .tour-hero { background: linear-gradient(160deg, #1f3a12 0%, #3b5e23 55%, #4f7f30 100%); color: #fff; padding: 64px 0 0; position: relative; display: flow-root; }
#tourPage .tour-hero .tour-eyebrow { color: #b9dd9a; }
#tourPage .tour-hero h1 { font-size: 34px; line-height: 1.25; font-weight: 700; letter-spacing: -.01em; margin: 0 0 18px; max-width: 980px; }
#tourPage .tour-sub { font-size: 21px; line-height: 1.6; color: #e3efd8; max-width: 800px; margin: 0 0 6px; }
#tourPage .tour-hero-list { list-style: none; margin: 0 0 26px; padding: 0; max-width: 800px; }
#tourPage .tour-hero-list li { font-size: 20px; line-height: 1.5; color: #fff; padding: 4px 0 4px 26px; position: relative; }
#tourPage .tour-hero-list li::before { content: ''; position: absolute; left: 4px; top: 13px; width: 9px; height: 9px; border-radius: 50%; background: #b9dd9a; }
#tourPage .tour-cta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }
#tourPage .tour-cta-center { justify-content: center; margin: 28px 0 18px; }
#tourPage .tour-try { font: inherit; font-size: 15px; font-weight: 700; color: #fff; background: #5d9239; border: 1px solid #5d9239; border-radius: 8px; padding: 11px 18px; cursor: pointer; box-shadow: 0 1px 2px rgba(0,0,0,.08); }
#tourPage .tour-try:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-hero .tour-try { background: #fff; color: #2f4d1a; border-color: #fff; }
#tourPage .tour-hero .tour-try:hover { background: #eef6e6; }
#tourPage .tour-hero-shot { margin-top: 44px; margin-bottom: -90px; position: relative; z-index: 1; }
#tourPage .tour-stats { background: #fff; padding: 120px 0 10px; }
#tourPage .tour-stats-in { display: grid; grid-template-columns: repeat(4, 1fr); gap: 18px; border-top: 1px solid #e5e7eb; border-bottom: 1px solid #e5e7eb; padding: 22px 0; }
#tourPage .tour-stats-in div { text-align: center; }
#tourPage .tour-stats-in b { display: block; font-size: 32px; font-weight: 800; color: #3b5e23; letter-spacing: -.02em; line-height: 1.1; }
#tourPage .tour-stats-in span { font-size: 13px; color: #6b7280; }
#tourPage .tour-sec { padding: 72px 0; }
#tourPage .tour-alt { background: #f6faf2; }
#tourPage .tour-feature { display: grid; grid-template-columns: minmax(0, 5fr) minmax(0, 7fr); gap: 48px; align-items: center; }
#tourPage .tour-feature + .tour-feature, #tourPage .tour-two + .tour-feature, #tourPage .tour-feature + .tour-two, #tourPage .tour-two + .tour-wide-shot, #tourPage .tour-wide-shot + .tour-feature { margin-top: 56px; }
#tourPage .tour-feature + .tour-chips { margin-top: 56px; }
#tourPage .tour-flip .tour-text { order: 2; }
#tourPage .tour-text h2 { font-size: 32px; line-height: 1.15; font-weight: 800; letter-spacing: -.02em; margin: 0 0 14px; color: #111827; }
#tourPage .tour-text h3 { font-size: 21px; font-weight: 800; margin: 0 0 12px; color: #111827; letter-spacing: -.01em; }
#tourPage .tour-text p { font-size: 16.5px; line-height: 1.6; color: #374151; margin: 0 0 14px; }
#tourPage .tour-text-wide { max-width: 820px; margin-bottom: 30px; }
#tourPage .tour-checks { list-style: none; margin: 0; padding: 0; }
#tourPage .tour-checks li { position: relative; padding-left: 28px; font-size: 15px; line-height: 1.55; color: #374151; margin: 0 0 10px; }
#tourPage .tour-checks li::before { content: ''; position: absolute; left: 0; top: 4px; width: 18px; height: 18px; border-radius: 50%; background: #eef6e6; border: 1px solid #b9dd9a; }
#tourPage .tour-checks li::after { content: ''; position: absolute; left: 6px; top: 8px; width: 5px; height: 8px; border: solid #3b5e23; border-width: 0 2px 2px 0; transform: rotate(45deg); }
#tourPage .tour-shot { margin: 0; background: #fff; border: 1px solid #d9dee5; border-radius: 12px; overflow: hidden; box-shadow: 0 14px 40px rgba(17,24,39,.16), 0 2px 6px rgba(17,24,39,.06); }
#tourPage .tour-shot-bar { height: 22px; background: #f3f4f6; border-bottom: 1px solid #e5e7eb; display: flex; align-items: center; gap: 5px; padding: 0 10px; }
#tourPage .tour-shot-bar span { width: 8px; height: 8px; border-radius: 50%; background: #d1d5db; }
#tourPage .tour-shot img { display: block; width: 100%; height: auto; }
#tourPage .tour-shot figcaption { font-size: 13px; line-height: 1.45; color: #6b7280; padding: 9px 12px 10px; border-top: 1px solid #eef0f3; background: #fff; }
#tourPage .tour-shot-act { padding: 2px 12px 12px; background: #fff; }
#tourPage .tour-shot-act .tour-try { width: 100%; }
#tourPage .tour-hero-shot { max-width: 600px; }
#tourPage .tour-hero-shot .tour-shot { box-shadow: 0 30px 70px rgba(0,0,0,.35); }
#tourPage .tour-hero .tour-shot-act .tour-try { background: #5d9239; color: #fff; border-color: #5d9239; }
#tourPage .tour-hero .tour-shot-act .tour-try:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-stack { display: flex; flex-direction: column; gap: 18px; }
#tourPage .tour-shot-dialog { max-width: 420px; margin: 0 auto; }
#tourPage .tour-two { display: grid; grid-template-columns: 1fr 1fr; gap: 36px; align-items: start; }
#tourPage .tour-two-tight { gap: 24px; }
/* Paired figures share one frame height and keep their buttons on one
   line, whatever the pictures' own proportions. */
#tourPage .tour-two .tour-shot { display: flex; flex-direction: column; }
#tourPage .tour-two .tour-shot img { aspect-ratio: 5 / 6; object-fit: contain; padding: 6px; }
#tourPage .tour-two .tour-shot figcaption { flex: 1; }
#tourPage .tour-stack .tour-shot img { aspect-ratio: 4 / 3; object-fit: contain; padding: 6px; }
#tourPage .tour-card-fig h3 { font-size: 19px; font-weight: 800; margin: 18px 0 8px; color: #111827; }
#tourPage .tour-card-fig p { font-size: 15px; line-height: 1.55; color: #374151; margin: 0 0 14px; }
#tourPage .tour-wide-shot { margin-top: 40px; }
#tourPage .tour-cards { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-top: 6px; }
#tourPage .tour-card { background: #fff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 14px 16px; }
#tourPage .tour-card b { display: block; font-size: 15px; color: #3b5e23; margin-bottom: 4px; }
#tourPage .tour-card span { font-size: 13.5px; line-height: 1.5; color: #4b5563; }
#tourPage .tour-more { padding: 56px 0 40px; border-top: 1px solid #e5e7eb; }
#tourPage .tour-chips { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; }
#tourPage .tour-chips div { background: #f9fafb; border: 1px solid #e5e7eb; border-radius: 10px; padding: 12px 14px; font-size: 13px; line-height: 1.45; color: #4b5563; }
#tourPage .tour-chips b { display: block; color: #111827; font-size: 14px; margin-bottom: 2px; }
#tourPage .tour-foot { text-align: center; font-size: 12.5px; color: #9ca3af; margin: 10px 0 0; }
@media (max-width: 900px) {
  #tourPage .tour-feature, #tourPage .tour-two, #tourPage .tour-cards { grid-template-columns: 1fr; }
  #tourPage .tour-feature { gap: 26px; }
  #tourPage .tour-flip .tour-text { order: 0; }
  #tourPage .tour-stats-in { grid-template-columns: 1fr 1fr; }
  #tourPage .tour-chips { grid-template-columns: 1fr 1fr; }
  #tourPage .tour-hero h1 { font-size: 26px; }
  #tourPage .tour-text h2 { font-size: 26px; }
  #tourPage .tour-sec { padding: 48px 0; }
}
@media (max-width: 640px) {
  #tourPage .tour-wrap { padding: 0 16px; }
  #tourPage .tour-brand span, #tourPage .tour-nav { display: none; }
  #tourPage .tour-close { margin-left: auto; }
  #tourPage .tour-hero { padding-top: 40px; }
  #tourPage .tour-hero h1 { font-size: 22px; }
  #tourPage .tour-sub { font-size: 16px; }
  #tourPage .tour-hero-shot { margin-bottom: -40px; }
  #tourPage .tour-stats { padding-top: 64px; }
  #tourPage .tour-stats-in b { font-size: 24px; }
  #tourPage .tour-chips { grid-template-columns: 1fr; }
  #tourPage .tour-try { width: 100%; text-align: center; min-height: 44px; }
}`;

    // ------------------------------------------------------------ shell
    let built = false, io = null;
    function build() {
        if (built) return;
        built = true;
        const style = document.createElement('style');
        style.textContent = CSS;
        document.head.appendChild(style);
        const page = document.createElement('div');
        page.id = 'tourPage';
        page.setAttribute('role', 'dialog');
        page.setAttribute('aria-label', 'Tour of the app');
        document.body.appendChild(page);
        page.addEventListener('click', (e) => {
            const act = e.target.closest('[data-act]');
            if (act) { e.preventDefault(); run(act.getAttribute('data-act')); return; }
            const link = e.target.closest('.tour-nav a');
            if (link) {
                e.preventDefault();
                const t = page.querySelector(link.getAttribute('href'));
                if (t) page.scrollTo({ top: t.offsetTop - 54, behavior: 'smooth' });
            }
        });
        document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && isOpen()) close(); });
    }
    const isOpen = () => !!document.getElementById('tourPage')?.classList.contains('active');

    function render() {
        const page = document.getElementById('tourPage');
        page.innerHTML = pageHTML();
        page.querySelector('#tourClose').addEventListener('click', close);
        // Section links follow the scroll.
        if (io) io.disconnect();
        const links = [...page.querySelectorAll('.tour-nav a')];
        io = new IntersectionObserver((entries) => {
            entries.forEach(en => {
                if (!en.isIntersecting) return;
                links.forEach(l => l.classList.toggle('on', l.getAttribute('href') === '#' + en.target.id));
            });
        }, { root: page, rootMargin: '-40% 0px -55% 0px', threshold: 0 });
        links.forEach(l => { const t = page.querySelector(l.getAttribute('href')); if (t) io.observe(t); });
    }

    function run(name) {
        const fn = ACTIONS[name];
        if (!fn || !A()) return;
        close();
        try { const r = fn(A()); if (r && r.catch) r.catch(e => console.warn('Tour action failed:', e)); }
        catch (e) { console.warn('Tour action failed:', e); }
    }

    function open() {
        build();
        render();
        const page = document.getElementById('tourPage');
        page.classList.add('active');
        page.scrollTop = 0;
    }
    function close() {
        document.getElementById('tourPage')?.classList.remove('active');
    }

    window.CorrelateTour = { open, close, isOpen };
})();
