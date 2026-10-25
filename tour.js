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
    const HERO_SET = ['TP53', 'MDM2', 'MDM4', 'CDKN1A', 'PPM1D', 'USP7', 'MTOR', 'RPTOR', 'MLST8', 'RHEB', 'TSC1', 'TSC2',
        'KRAS', 'RAF1', 'MAPK1', 'MAP2K1', 'SHOC2', 'NRAS', 'SOX10', 'MITF', 'TFAP2A', 'PAX3'];
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
    function openBrowserFiltered(a, tissue, hotspotGene) {
        a.openCellLineBrowser();
        setTimeout(() => {
            if (tissue) setVal('clbTissueFilter', tissue);
            if (hotspotGene) setVal('clbHotspotFilter', hotspotGene);
        }, 400);
    }
    function compareSkinBraf(a) {
        openBrowserFiltered(a, 'Skin', 'BRAF');
        setTimeout(() => {
            a._resetGEInspectScope();
            a._geInspectScope = 'lineage';
            a._runBusy({ button: 'clbInspectGEBtn', label: 'Comparing' }, () => a.inspectSelectionGE());
        }, 1600);
    }
    function runMutationExample() {
        document.getElementById('modeMutationBtn')?.click();
        setTimeout(() => {
            setVal('mutationHotspotSelect', 'BRAF');
            setVal('lineageFilter', 'Skin');
            setTimeout(() => document.getElementById('runMutationAnalysisBtn')?.click(), 500);
        }, 700);
    }
    const ACTIONS = {
        hero: () => runExampleSet({ cutoff: 0.3, colorGE: true }),
        heroMrna: () => runExampleSet({ cutoff: 0.3, basis: 'expr' }),
        tyr: (a) => { a._inspectPreset = { colorByGene: 'TYR', colorByGeneType: 'expr' }; a.openInspectByGenes('SOX10', 'MITF'); },
        braf: (a) => { a._inspectPreset = { hotspotGene: 'BRAF', hotspotMode: 'color' }; a.openInspectByGenes('BRAF', 'MAPK1'); },
        mutation: () => runMutationExample(),
        compare: (a) => compareSkinBraf(a),
        kras: (a) => a.openGeneBrowser('KRAS'),
        tsc1: (a) => a.openGeneBrowser('TSC1'),
        clb: (a) => openBrowserFiltered(a, 'Skin', 'BRAF'),
        wiki: (a) => a.openCellLineWiki(D.lineId('A375') || 'ACH-000219'),
        stats: () => { const m = document.getElementById('infographicModal'); if (m) m.style.display = 'flex'; },
    };

    // ------------------------------------------------------------ page
    const shot = (file, alt, act, o) => {
        o = o || {};
        return `<figure class="tour-shot${o.cls ? ' ' + o.cls : ''}"${act ? ` data-act="${act}" title="Open this view in the app"` : ''}>
            <div class="tour-shot-bar"><span></span><span></span><span></span></div>
            <img src="${IMG}${file}.webp" alt="${esc(alt)}" loading="${o.eager ? 'eager' : 'lazy'}" decoding="async">
            ${o.cap ? `<figcaption>${o.cap}</figcaption>` : ''}
        </figure>`;
    };
    const tryBtn = (act, label, primary) => `<button type="button" class="tour-try${primary === false ? ' tour-try-2' : ''}" data-act="${act}">${label}<span aria-hidden="true"> &rarr;</span></button>`;
    const checks = (items) => `<ul class="tour-checks">${items.map(t => `<li>${t}</li>`).join('')}</ul>`;

    function pageHTML() {
        const nL = num(D.nLines()), nE = num(D.nExprLines()), nG = num(D.nGenes()), nSet = HERO_SET.length;
        return `
<header class="tour-top">
    <div class="tour-wrap tour-top-in">
        <div class="tour-brand"><b>Correlate</b><span>Tour of the app</span></div>
        <nav class="tour-nav" aria-label="Sections">
            <a href="#tour-network">Network</a>
            <a href="#tour-links">Correlations</a>
            <a href="#tour-compare">Compare groups</a>
            <a href="#tour-gene">Gene Browser</a>
            <a href="#tour-cells">Cell lines</a>
            <a href="#tour-export">Export</a>
        </nav>
        <button type="button" class="tour-close" id="tourClose" aria-label="Close the tour">&times;</button>
    </div>
</header>

<section class="tour-hero">
    <div class="tour-wrap">
        <p class="tour-eyebrow">Human cancer cell lines, DepMap 26Q1</p>
        <h1>From a gene list to a network of dependencies.</h1>
        <p class="tour-sub">Correlate finds the genes whose CRISPR knockout effects, or mRNA levels, rise and fall together across the DepMap cancer cell lines, and turns that into networks, plots and pages you can take straight into a presentation or a manuscript.</p>
        <div class="tour-cta">
            ${tryBtn('hero', 'Build the example network')}
            ${tryBtn('kras', 'Open KRAS in the Gene Browser', false)}
        </div>
    </div>
    <div class="tour-wrap tour-hero-shot">
        ${shot('hero_network', 'The example gene set as a network, coloured by gene effect', 'hero', { eager: true })}
    </div>
</section>

<section class="tour-stats">
    <div class="tour-wrap tour-stats-in">
        <div><b>${nL}</b><span>cell lines with CRISPR gene effect</span></div>
        <div><b>${nE}</b><span>cell lines with mRNA expression</span></div>
        <div><b>${nG}</b><span>genes</span></div>
        <div><b>10+</b><span>pages of wiki per cell line</span></div>
    </div>
</section>

<section class="tour-sec" id="tour-network">
    <div class="tour-wrap tour-feature">
        <div class="tour-text">
            <p class="tour-eyebrow">01 &middot; The network</p>
            <h2>Paste a gene set, get a network</h2>
            <p>Every pair of genes in your set is compared across the cell lines. Two genes whose profiles correlate above your cutoff are connected by an edge, a line between them: blue for a positive correlation, red for a negative one, thicker for a stronger one. Genes that end up in a tightly connected group usually belong to the same complex or pathway. A single edge to a gene outside such a group points to a partner you did not expect.</p>
            ${checks([
                `Build the network on CRISPR gene effect or on mRNA expression; a badge on the network says which one. The same gene set often gives two different networks, and the difference is informative: on mRNA (picture) the broadly expressed genes ${gi('MTOR')}, ${gi('RPTOR')}, ${gi('TSC1')}, ${gi('TSC2')} and ${gi('USP7')} merge into one block, while the melanocyte genes stay together.`,
                'Colour the nodes by gene effect, by expression, or by your own screen statistics: paste LFC and FDR next to the gene list and your screen sits on top of the DepMap network.',
                'Expand the set to find new partners, send any cluster to Enrichr, open a full-page view for large sets, and drag genes where you want them before you export.',
            ])}
            <div class="tour-cta">${tryBtn('heroMrna', 'Try it: build the example set on mRNA')}${tryBtn('hero', 'Build it on gene effect', false)}</div>
        </div>
        ${shot('net_mrna', 'The example set built on mRNA expression', 'heroMrna', { cap: `The same ${nSet} genes built on mRNA expression. Built on gene effect (top of the page) they split into the p53, mTORC1, MAPK and melanocyte groups.` })}
    </div>
</section>

<section class="tour-sec tour-alt" id="tour-links">
    <div class="tour-wrap tour-feature tour-flip">
        <div class="tour-text">
            <p class="tour-eyebrow">02 &middot; Behind every edge</p>
            <h2>Every edge in the network is a plot you can open</h2>
            <p>Double-click an edge, the line between two genes, and you see the cell lines behind it: one dot per cell line, one gene against the other, with r, p and the regression line. This is where you judge whether a correlation is real, driven by one cancer type, or explained by a third gene.</p>
            ${checks([
                'Filter to a cancer type, subtype or disease, to a mutation, a fusion or a copy-number change, or any combination. The statistics update with the filter.',
                `Colour by tissue, by the mutation status of any gene, or by the expression or gene effect of a third gene. <b>Suggest</b> ranks the genes that explain the pattern best; here it points at ${gi('TYR')}.`,
                'Highlight cell lines by name or group, draw gates to define cohorts, compare tissues side by side, and open any gene’s effect across all cell lines by tissue.',
            ])}
            <div class="tour-cta">${tryBtn('tyr', 'Try it: SOX10 vs MITF coloured by TYR')}${tryBtn('braf', 'BRAF vs MAPK1 by mutation', false)}</div>
        </div>
        <div class="tour-stack">
            ${shot('scatter_tyr', 'SOX10 against MITF gene effect, coloured by TYR mRNA', 'tyr', { cap: `${gi('SOX10')} vs ${gi('MITF')} gene effect, coloured by ${gi('TYR')} mRNA: the cell lines that depend on both are the ones that express the pigmentation gene.` })}
            ${shot('scatter_braf', 'BRAF against MAPK1 gene effect with BRAF-mutant cell lines coloured', 'braf', { cls: 'tour-shot-sm', cap: `${gi('BRAF')} vs ${gi('MAPK1')}, with ${gi('BRAF')}-mutant cell lines in colour.` })}
        </div>
    </div>
</section>

<section class="tour-sec" id="tour-compare">
    <div class="tour-wrap">
        <div class="tour-text tour-text-wide">
            <p class="tour-eyebrow">03 &middot; Compare two groups</p>
            <h2>Which genes separate two groups of cell lines?</h2>
            <p>Two entrances to the same question. <b>Mutation analysis</b> splits the panel into wild-type and mutated for a gene you choose and ranks every gene by the difference in gene effect or expression. The <b>Cell Line Browser</b> compares any selection you can filter to against the rest, as volcano plots. Both use Welch’s t-test per gene, and both let you restrict the comparison to one lineage so a tissue difference is not mistaken for a mutation effect.</p>
        </div>
        <div class="tour-two">
            <div class="tour-card-fig">
                ${shot('mutation_table', 'Mutation analysis table for BRAF in skin cell lines', 'mutation')}
                <h3>Mutation analysis: ${gi('BRAF')} in skin</h3>
                <p>Wild-type, all mutated and double-mutant cell lines in separate column groups, with the differences and p-values beside them. Search the table, send the list to Enrichr or a heatmap, open any gene’s plot.</p>
                ${tryBtn('mutation', 'Try it: BRAF mutation analysis in skin')}
            </div>
            <div class="tour-card-fig">
                ${shot('volcano', 'Volcano plots comparing BRAF-mutant skin cell lines with the other skin lines', 'compare')}
                <h3>Selection compare: 50 ${gi('BRAF')}-mutant skin cell lines vs the other 25</h3>
                <p>Gene effect and mRNA side by side, q-values from Benjamini-Hochberg, the top genes labelled and listed. Adjust the thresholds and the plots redraw; the lists go straight to a network or Enrichr.</p>
                ${tryBtn('compare', 'Try it: compare BRAF-mutant skin lines')}
            </div>
        </div>
    </div>
</section>

<section class="tour-sec tour-alt" id="tour-gene">
    <div class="tour-wrap">
        <div class="tour-text tour-text-wide">
            <p class="tour-eyebrow">04 &middot; The Gene Browser</p>
            <h2>One gene, everything the data says about it</h2>
            <p>Type a gene and get its page: what it is, how essential it is and where, how it is expressed, which genes co-depend and co-express with it, which protein complexes and pathways it belongs to, its mutations, fusions, copy number, structure and sequence. The page exports as a PDF of five pages or more, and as a file for a language model.</p>
        </div>
        <div class="tour-feature">
            <div class="tour-text">
                <h3>A summary written from the data</h3>
                ${checks([
                    `A verdict for dependency and for expression, ranked against every gene: common essential or selective, expressed everywhere or restricted to a few cancers. ${gi('KRAS')} is essential in half of the screened lines, most strongly in pancreas and bowel, and hotspot-mutated cell lines depend on it most.`,
                    'The most dependent tissues and cell lines, the highest and lowest expressing tissues, and the distributions behind them, each with a Copy figure button.',
                    'Hotspot mutations by position and cancer type, fusions, copy number, and drug response from PRISM where the gene is a known target.',
                ])}
                <div class="tour-cta">${tryBtn('kras', 'Try it: open KRAS')}</div>
            </div>
            ${shot('gb_top', 'The Gene Browser page for KRAS: description, ideogram and summary', 'kras', { cap: `${gi('KRAS')} in the Gene Browser: description, chromosome position and the summary line.` })}
        </div>
        <div class="tour-two tour-two-tight">
            ${shot('gb_dependency', 'Dependency section for KRAS with the gene effect distribution', 'kras', { cap: 'Dependency: the gene effect distribution with the gene marked, the most dependent tissues and lines.' })}
            ${shot('gb_expression', 'Expression section for KRAS with the mRNA distribution', 'kras', { cap: 'Expression: the same for mRNA, with the highest and lowest tissues.' })}
        </div>
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <h3>Does the complex act as one unit?</h3>
                ${checks([
                    'Co-dependencies: the genes whose knockout effect follows this gene’s across the lines, as a table and a figure, with known partners from the same complex or pathway ringed.',
                    'Co-expression as a heatmap of every pair, with genes on the same chromosome arm flagged so an amplicon is not mistaken for shared function.',
                    `Protein complexes from the EBI Complex Portal and Reactome pathways, each member tagged with its co-dependency. ${gi('TSC2')} follows ${gi('TSC1')} at r 0.90; ${gi('TBC1D7')}, the third member, at 0.25.`,
                ])}
                <div class="tour-cta">${tryBtn('tsc1', 'Try it: open TSC1 and its complex')}</div>
            </div>
            <div class="tour-stack">
                ${shot('gb_codep', 'Co-dependencies of TSC1 as a table and a fan figure', 'tsc1')}
                ${shot('gb_complex', 'The TSC1-TSC2 complex card with co-dependency per member', 'tsc1', { cls: 'tour-shot-sm' })}
            </div>
        </div>
        <div class="tour-wide-shot">
            ${shot('gb_coexpr', 'Co-expression table and heatmap for TSC1', 'tsc1', { cap: `Co-expression with ${gi('TSC1')}: the heatmap shows how the correlated genes relate to each other; names in amber sit within 10 Mb on the same chromosome.` })}
        </div>
    </div>
</section>

<section class="tour-sec" id="tour-cells">
    <div class="tour-wrap tour-feature">
        <div class="tour-text">
            <p class="tour-eyebrow">05 &middot; Cell lines</p>
            <h2>Find the right cell lines, then read their pages</h2>
            <p>The Cell Line Browser filters the panel by tissue, subtype and disease, by sex and loss of Y, by mutations, fusions and copy number, by expression or dependency of any gene, and by drug response, in any combination. Every cell line has a page of ten sections or more assembled from DepMap and other sources.</p>
            ${checks([
                `Combine filters: skin, ${gi('BRAF')} mutated, female, sensitive to vemurafenib. The count updates as you go, and the list is one click from a heatmap, PCA or UMAP.`,
                'Each cell line opens with a plain-language summary: origin, candidate drivers, derived cell lines, genome state. The wiki adds provenance, STR, signatures, dependencies, expression and drug response.',
                'Download a cell line’s full gene effect, expression and mutation tables, or the list on screen, as CSV.',
            ])}
            <div class="tour-cta">${tryBtn('clb', 'Try it: skin cell lines with mutated BRAF')}${tryBtn('wiki', 'Open the A375 page', false)}</div>
        </div>
        <div class="tour-stack">
            ${shot('clb', 'The Cell Line Browser filtered to skin cell lines with mutated BRAF, A375 selected', 'clb', { cap: `Skin cell lines with mutated ${gi('BRAF')}, with A375 opened in the side panel.` })}
            ${shot('wiki', 'The wiki page for A375', 'wiki', { cls: 'tour-shot-sm', cap: 'The A375 wiki page: executive summary, provenance, classification and onwards.' })}
        </div>
    </div>
</section>

<section class="tour-sec tour-alt" id="tour-export">
    <div class="tour-wrap tour-feature tour-flip">
        <div class="tour-text">
            <p class="tour-eyebrow">06 &middot; Take it with you</p>
            <h2>Every view leaves with you</h2>
            <p>The same buttons sit under every network, plot, table and page.</p>
            <div class="tour-cards">
                <div class="tour-card"><b>Export image</b><span>PNG, SVG, TIFF, PDF or PowerPoint at the print size you choose. PNG and SVG carry their settings, so opening the image in Correlate later brings the analysis back.</span></div>
                <div class="tour-card"><b>Export for AI</b><span>One file with the data behind the view, written for a language model. Upload it and ask: which of these dependencies has a drug? Which cell line should I use for this experiment?</span></div>
                <div class="tour-card"><b>Methods</b><span>How the numbers on screen were calculated, as a paragraph for a methods section and as a plain walk-through.</span></div>
                <div class="tour-card"><b>Tables</b><span>CSV of every result; every gene for a cell line; every cell line for a gene.</span></div>
            </div>
        </div>
        ${shot('export_dialog', 'The Export image dialog', null, { cls: 'tour-shot-dialog', cap: 'Export image: what to export, format, size in cm, resolution, background.' })}
    </div>
</section>

<section class="tour-sec tour-more">
    <div class="tour-wrap">
        <p class="tour-eyebrow">Also in Correlate</p>
        <div class="tour-chips">
            <div><b>Heatmap</b>a gene set across the cell lines, clustered, with annotation rows</div>
            <div><b>Compilation</b>a grid of correlation plots or gene effect panels, one per pair or gene</div>
            <div><b>Threshold scan</b>how many edges your set keeps as the cutoff moves</div>
            <div><b>PCA and UMAP</b>of cell lines or genes, with gene classes and loadings</div>
            <div><b>Correlation matrix</b>of the set, on gene effect or mRNA</div>
            <div><b>Gates</b>draw a cohort on any plot and carry it to the next view</div>
            <div><b>Open saved figure</b>drop an exported PNG or SVG back in and continue</div>
            <div><b>Green Listed</b>library genes and curated gene classes, offline</div>
        </div>
        <div class="tour-cta tour-cta-center">${tryBtn('hero', 'Build the example network')}${tryBtn('stats', 'The statistics behind the app', false)}</div>
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
#tourPage .tour-hero h1 { font-size: 46px; line-height: 1.08; font-weight: 800; letter-spacing: -.02em; margin: 0 0 18px; max-width: 760px; }
#tourPage .tour-sub { font-size: 18px; line-height: 1.6; color: #e3efd8; max-width: 720px; margin: 0 0 26px; }
#tourPage .tour-cta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }
#tourPage .tour-cta-center { justify-content: center; margin: 28px 0 18px; }
#tourPage .tour-try { font: inherit; font-size: 15px; font-weight: 700; color: #fff; background: #5d9239; border: 1px solid #5d9239; border-radius: 8px; padding: 11px 18px; cursor: pointer; box-shadow: 0 1px 2px rgba(0,0,0,.08); }
#tourPage .tour-try:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-try-2 { background: transparent; color: #3b5e23; border-color: #a3c98a; font-weight: 600; box-shadow: none; }
#tourPage .tour-try-2:hover { background: #eef6e6; border-color: #5d9239; }
#tourPage .tour-hero .tour-try { background: #fff; color: #2f4d1a; border-color: #fff; }
#tourPage .tour-hero .tour-try:hover { background: #eef6e6; }
#tourPage .tour-hero .tour-try-2 { background: transparent; color: #fff; border-color: rgba(255,255,255,.55); }
#tourPage .tour-hero .tour-try-2:hover { background: rgba(255,255,255,.12); }
#tourPage .tour-hero-shot { margin-top: 44px; margin-bottom: -90px; position: relative; z-index: 1; }
#tourPage .tour-stats { background: #fff; padding: 120px 0 10px; }
#tourPage .tour-stats-in { display: grid; grid-template-columns: repeat(4, 1fr); gap: 18px; border-top: 1px solid #e5e7eb; border-bottom: 1px solid #e5e7eb; padding: 22px 0; }
#tourPage .tour-stats-in div { text-align: center; }
#tourPage .tour-stats-in b { display: block; font-size: 32px; font-weight: 800; color: #3b5e23; letter-spacing: -.02em; line-height: 1.1; }
#tourPage .tour-stats-in span { font-size: 13px; color: #6b7280; }
#tourPage .tour-sec { padding: 72px 0; }
#tourPage .tour-alt { background: #f6faf2; }
#tourPage .tour-feature { display: grid; grid-template-columns: minmax(0, 5fr) minmax(0, 7fr); gap: 48px; align-items: center; }
#tourPage .tour-feature + .tour-feature, #tourPage .tour-two + .tour-feature, #tourPage .tour-feature + .tour-two, #tourPage .tour-two + .tour-wide-shot { margin-top: 56px; }
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
#tourPage .tour-shot[data-act] { cursor: pointer; transition: transform .18s, box-shadow .18s; }
#tourPage .tour-shot[data-act]:hover { transform: translateY(-2px); box-shadow: 0 20px 48px rgba(17,24,39,.22), 0 2px 6px rgba(17,24,39,.08); }
#tourPage .tour-shot-bar { height: 22px; background: #f3f4f6; border-bottom: 1px solid #e5e7eb; display: flex; align-items: center; gap: 5px; padding: 0 10px; }
#tourPage .tour-shot-bar span { width: 8px; height: 8px; border-radius: 50%; background: #d1d5db; }
#tourPage .tour-shot img { display: block; width: 100%; height: auto; }
#tourPage .tour-shot figcaption { font-size: 13px; line-height: 1.45; color: #6b7280; padding: 9px 12px 10px; border-top: 1px solid #eef0f3; background: #fff; }
#tourPage .tour-hero-shot .tour-shot { box-shadow: 0 30px 70px rgba(0,0,0,.35); }
#tourPage .tour-stack { display: flex; flex-direction: column; gap: 18px; }
#tourPage .tour-shot-sm { width: 78%; margin-left: auto; }
#tourPage .tour-flip .tour-shot-sm { margin-left: 0; margin-right: auto; }
#tourPage .tour-shot-dialog { max-width: 420px; margin: 0 auto; }
#tourPage .tour-two { display: grid; grid-template-columns: 1fr 1fr; gap: 36px; align-items: start; }
#tourPage .tour-two-tight { gap: 24px; }
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
  #tourPage .tour-shot-sm { width: 100%; }
  #tourPage .tour-stats-in { grid-template-columns: 1fr 1fr; }
  #tourPage .tour-chips { grid-template-columns: 1fr 1fr; }
  #tourPage .tour-hero h1 { font-size: 34px; }
  #tourPage .tour-text h2 { font-size: 26px; }
  #tourPage .tour-sec { padding: 48px 0; }
}
@media (max-width: 640px) {
  #tourPage .tour-wrap { padding: 0 16px; }
  #tourPage .tour-brand span, #tourPage .tour-nav { display: none; }
  #tourPage .tour-close { margin-left: auto; }
  #tourPage .tour-hero { padding-top: 40px; }
  #tourPage .tour-hero h1 { font-size: 28px; }
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
