/* Correlate Quest: an arcade quiz whose questions are generated from the data
   the app already holds. Everything here reads window.app; the only writes are
   the two "open in app" hand-offs, which close the quiz first. */
(function () {
    'use strict';

    // ---------------------------------------------------------------- palette
    const PAL = {
        bg: '#1a1c2c',
        panel: '#262b44',
        ink: '#f9f9f9',
        dim: '#8b93b8',
        yellow: '#ffd23f',
        green: '#3fbf7f',
        red: '#ff5e5b',
        blue: '#5bc0eb'
    };

    const TIMER_MS = 20000;
    const START_LIVES = 3;
    const SCORE_KEY = 'correlateQuizScores';
    const SOUND_KEY = 'correlateQuizSound';
    const MAX_SCORES = 10;

    // ------------------------------------------------------------------- css
    const CSS = `
#cq-root{position:fixed; inset:0; z-index:20000; background:${PAL.bg}; color:${PAL.ink};
  font-family:'Press Start 2P','Courier New',monospace; overflow-y:auto; -webkit-overflow-scrolling:touch;
  display:block; letter-spacing:0;}
#cq-root *{box-sizing:border-box; font-family:inherit;}
#cq-root canvas{image-rendering:pixelated; image-rendering:crisp-edges; display:block;}
.cq-shell{max-width:640px; margin:0 auto; padding:56px 12px 40px; min-height:100%;}
.cq-shell.cq-center{display:flex; flex-direction:column; justify-content:center;}
.cq-corner{position:fixed; top:6px; width:44px; height:44px; z-index:2; background:${PAL.panel};
  color:${PAL.ink}; border:4px solid ${PAL.ink}; font-size:12px; line-height:1; cursor:pointer; padding:0;}
.cq-corner:active{transform:translate(2px,2px);}
#cq-close{right:8px;}
#cq-sound{right:60px; font-size:11px;}
.cq-title{font-size:19px; line-height:1.9; color:${PAL.yellow}; text-align:center;
  text-shadow:4px 4px 0 ${PAL.red}; margin:18px 0 10px;}
.cq-tag{font-size:8px; line-height:1.8; color:${PAL.dim}; text-align:center; margin:0 0 22px;}
.cq-blink{font-size:12px; color:${PAL.green}; text-align:center; margin:18px 0; animation:cq-blink 1.1s steps(1) infinite;}
@keyframes cq-blink{50%{opacity:0;}}
.cq-btn{display:block; width:100%; min-height:48px; margin:0 0 12px; padding:12px 10px;
  background:${PAL.panel}; color:${PAL.ink}; border:4px solid ${PAL.ink}; box-shadow:5px 5px 0 rgba(0,0,0,0.55);
  font-size:12px; line-height:1.6; text-align:center; cursor:pointer; overflow-wrap:anywhere;}
.cq-btn:active{transform:translate(3px,3px); box-shadow:2px 2px 0 rgba(0,0,0,0.55);}
.cq-btn.cq-go{background:${PAL.green}; color:${PAL.bg}; border-color:${PAL.ink};}
.cq-btn.cq-sel{background:${PAL.blue}; color:${PAL.bg};}
.cq-btn.cq-quiet{background:transparent; font-size:10px; min-height:40px;}
.cq-ans{text-align:left; font-size:16px; line-height:1.5; display:flex; gap:10px; align-items:center;}
.cq-ans .cq-key{flex:0 0 auto; color:${PAL.yellow}; font-size:12px;}
.cq-ans.cq-right{background:${PAL.green}; color:${PAL.bg}; border-color:${PAL.ink};}
.cq-ans.cq-right .cq-key{color:${PAL.bg};}
.cq-ans.cq-wrong{background:${PAL.red}; color:${PAL.bg};}
.cq-ans.cq-wrong .cq-key{color:${PAL.bg};}
.cq-ans[disabled]{cursor:default;}
.cq-hud{display:flex; align-items:center; justify-content:space-between; gap:8px; font-size:10px;
  border:4px solid ${PAL.panel}; padding:8px; margin-bottom:10px; flex-wrap:wrap;}
.cq-hud b{color:${PAL.yellow}; font-weight:normal;}
.cq-hearts{color:${PAL.red}; letter-spacing:2px;}
.cq-bar{height:14px; border:4px solid ${PAL.ink}; padding:2px; margin-bottom:14px; background:${PAL.bg};}
.cq-bar > i{display:block; height:100%; background:${PAL.green}; width:100%;}
.cq-bar > i.cq-low{background:${PAL.red};}
.cq-cat{display:inline-block; font-size:9px; background:${PAL.blue}; color:${PAL.bg};
  padding:6px 8px; margin-bottom:12px;}
.cq-q{font-size:11px; line-height:1.9; margin:0 0 16px;}
.cq-quote{font-size:10px; line-height:2; color:${PAL.ink}; background:${PAL.panel};
  border-left:6px solid ${PAL.yellow}; padding:10px; margin:0 0 16px;}
.cq-fig{margin:0 0 16px;}
.cq-plot{background:#fff; border:4px solid ${PAL.ink}; padding:6px 6px 2px; margin:0 0 16px;}
.cq-plot, .cq-plot *{font-family:Arial,Helvetica,sans-serif;}
.cq-plot .js-plotly-plot, .cq-plot .plot-container{width:100% !important;}
.cq-exp{border:4px solid ${PAL.yellow}; padding:12px; margin:4px 0 12px; cursor:pointer;}
.cq-exp p{margin:0 0 10px; font-family:'Open Sans',system-ui,sans-serif; font-size:14px; line-height:1.6;}
.cq-exp .cq-verdict{font-size:11px; line-height:1.6; margin-bottom:10px;}
.cq-exp .cq-verdict.ok{color:${PAL.green};} .cq-exp .cq-verdict.no{color:${PAL.red};}
.cq-link{background:none; border:none; color:${PAL.blue}; font-size:9px; line-height:1.8;
  text-decoration:underline; cursor:pointer; padding:6px 0; display:inline-block;}
.cq-big{font-size:26px; color:${PAL.yellow}; text-align:center; line-height:1.6; margin:16px 0;}
.cq-stat{font-size:11px; line-height:2.2; text-align:center; color:${PAL.ink};}
.cq-stat span{color:${PAL.yellow};}
.cq-h2{font-size:13px; color:${PAL.green}; text-align:center; margin:22px 0 14px; line-height:1.6;}
.cq-note{font-size:8px; line-height:2; color:${PAL.dim}; text-align:center; margin:12px 0;}
.cq-how{font-size:10px; line-height:2.1; color:${PAL.ink}; border:4px solid ${PAL.panel}; padding:12px; margin-bottom:14px;}
.cq-table{width:100%; border-collapse:collapse; font-size:10px;}
.cq-table td{padding:8px 4px; border-bottom:3px solid ${PAL.panel}; line-height:1.5;}
.cq-table td.cq-r{text-align:right; color:${PAL.yellow};}
.cq-table tr.cq-me td{color:${PAL.green};}
.cq-table td.cq-co{font-size:8px; color:${PAL.dim}; overflow-wrap:anywhere;}
.cq-list{max-height:44vh; overflow-y:auto; border:4px solid ${PAL.panel}; padding:8px; margin-bottom:16px;}
.cq-list .cq-btn{margin-bottom:8px; font-size:11px; min-height:44px; text-align:left;}
.cq-row{display:flex; gap:10px;}
.cq-row .cq-btn{flex:1 1 0; min-width:0;}
.cq-slots{display:flex; gap:12px; justify-content:center; margin:16px 0;}
.cq-slot{text-align:center;}
.cq-slot .cq-ltr{font-size:30px; color:${PAL.yellow}; border:4px solid ${PAL.ink}; width:56px; height:64px;
  line-height:56px; background:${PAL.panel};}
.cq-slot button{width:56px; height:44px; margin-top:6px; background:${PAL.panel}; color:${PAL.ink};
  border:4px solid ${PAL.ink}; font-size:11px; cursor:pointer; padding:0;}
.cq-slot button:first-of-type{margin-top:8px;}
.cq-spin{width:48px; height:48px; margin:40px auto; background:${PAL.yellow}; animation:cq-spin 0.8s steps(4) infinite;}
@keyframes cq-spin{0%{transform:rotate(0);}100%{transform:rotate(360deg);}}
.cq-shake{animation:cq-shake 0.3s steps(2) 2;}
@keyframes cq-shake{0%{transform:translateX(0);}25%{transform:translateX(-8px);}75%{transform:translateX(8px);}100%{transform:translateX(0);}}
.cq-flash{animation:cq-flash 0.3s steps(1) 2;}
@keyframes cq-flash{50%{background:#4a1f28;}}
@media (min-width:641px){
  .cq-shell{padding:64px 16px 48px;}
  .cq-title{font-size:30px;}
  .cq-tag{font-size:11px;}
  .cq-q{font-size:12px;}
  .cq-quote{font-size:11px;}
  .cq-hud{font-size:11px;}
}
@media (prefers-reduced-motion:reduce){
  .cq-blink,.cq-shake,.cq-flash,.cq-spin{animation:none !important;}
}
/* The tour bar lives outside #cq-root: it has to stay readable while the app's
   own modals are open, so it sits above them but below the quiz itself. */
#cq-tour{position:fixed; left:0; right:0; bottom:0; z-index:19000; display:none;
  background:${PAL.bg}; color:${PAL.ink}; border-top:4px solid ${PAL.ink};
  font-family:'Press Start 2P','Courier New',monospace; padding:8px 10px 10px;
  max-height:34vh; overflow-y:auto; -webkit-overflow-scrolling:touch;
  box-shadow:0 -6px 0 rgba(0,0,0,0.45);}
#cq-tour *{box-sizing:border-box; font-family:inherit;}
#cq-tour.cq-collapsed{max-height:none;}
#cq-tour.cq-collapsed .cq-tour-body{display:none;}
.cq-tour-head{display:flex; align-items:center; gap:8px;}
.cq-tour-n{font-size:9px; line-height:1.6; color:${PAL.yellow}; flex:1 1 auto; overflow-wrap:anywhere;}
.cq-tour-pts{font-size:9px; line-height:1.6; color:${PAL.green};}
.cq-tour-chev{width:36px; height:30px; flex:0 0 auto; background:${PAL.panel}; color:${PAL.ink};
  border:3px solid ${PAL.ink}; font-size:10px; cursor:pointer; padding:0;}
.cq-tour-body{padding-top:8px;}
.cq-tour-title{font-size:11px; line-height:1.7; color:${PAL.blue}; margin-bottom:7px;}
#cq-tour .cq-tour-do{font-family:'Open Sans',system-ui,sans-serif; font-size:13px; line-height:1.5; margin-bottom:5px;}
#cq-tour .cq-tour-why{font-family:'Open Sans',system-ui,sans-serif; font-size:12px; line-height:1.5;
  color:${PAL.dim}; margin-bottom:7px;}
#cq-tour .cq-tour-hint{font-family:'Open Sans',system-ui,sans-serif; font-size:12px; line-height:1.5;
  color:${PAL.yellow}; margin-bottom:7px;}
#cq-tour .cq-tour-hint:empty{margin:0;}
.cq-tour-btns{display:flex; gap:6px; margin-bottom:8px;}
.cq-tour-b{flex:1 1 0; min-width:0; min-height:38px; background:${PAL.panel}; color:${PAL.ink};
  border:3px solid ${PAL.ink}; font-size:8px; line-height:1.5; cursor:pointer; padding:4px 2px;}
.cq-tour-b:active{transform:translate(2px,2px);}
.cq-tour-prog{display:flex; gap:4px; flex-wrap:wrap;}
.cq-sq{display:block; width:10px; height:10px; background:${PAL.panel}; border:2px solid ${PAL.dim};}
.cq-sq-on{background:${PAL.green}; border-color:${PAL.green};}
.cq-sq-now{border-color:${PAL.yellow};}
.cq-tour-done{font-size:14px; line-height:1.6; color:${PAL.green}; text-align:center; padding:16px 0;}
@media (min-width:641px){
  #cq-tour{left:auto; right:16px; bottom:16px; width:360px; max-height:62vh; border:4px solid ${PAL.ink};}
  .cq-tour-title{font-size:12px;}
  .cq-tour-b{font-size:9px;}
}`;

    // --------------------------------------------------------------- helpers
    const A = () => window.app || null;
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const ri = (n) => Math.floor(Math.random() * n);
    const pick = (a) => a[ri(a.length)];
    function shuffle(a) {
        const r = a.slice();
        for (let i = r.length - 1; i > 0; i--) { const j = ri(i + 1); const t = r[i]; r[i] = r[j]; r[j] = t; }
        return r;
    }
    const reduced = () => {
        try { return window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return false; }
    };

    // ------------------------------------------------------------ data layer
    const D = {
        ready() {
            const a = A();
            return !!(a && a.metadata && a.metadata.cellLines && a.geneEffects && a.geneIndex);
        },
        ids() { return A().metadata.cellLines; },
        n() { return A().nCellLines || A().metadata.cellLines.length; },
        name(id) { try { return A().getCellLineName(id) || id; } catch (e) { return id; } },
        lineage(id) { return A().cellLineMetadata?.lineage?.[id] || ''; },
        subtype(id) {
            const m = A().cellLineMetadata || {};
            return m.oncotreeSubtype?.[id] || m.subtype?.[id] || '';
        },
        row(gene) {
            const a = A();
            const gi = a.geneIndex.get(String(gene).toUpperCase());
            if (gi === undefined) return null;
            return a.geneEffects.subarray(gi * a.nCellLines, (gi + 1) * a.nCellLines);
        },
        // Missing gene effect arrives as NaN from the loader, but old files used
        // a -999 sentinel; treat both as no data.
        ge(row, i) {
            if (!row) return NaN;
            const v = row[i];
            return (!isFinite(v) || v <= -990) ? NaN : v;
        },
        geAt(gene, i) { return D.ge(D.row(gene), i); },
        expr(gene, i) {
            const a = A();
            try {
                const v = a.getExpressionValueByGEIndex(gene, i);
                return isFinite(v) ? v : NaN;
            } catch (e) { return NaN; }
        },
        hotspot(gene, id) {
            const g = A().mutations?.geneData?.[gene];
            return g ? (g.mutations?.[id] || 0) : 0;
        },
        damaging(id) {
            const m = A()._damagingCountByCL;
            const v = m && m.get ? m.get(id) : undefined;
            return (typeof v === 'number') ? v : null;
        },
        fusions(id) {
            const f = A().clinicalFusions?.byCellLine?.[id];
            return Array.isArray(f) ? f : [];
        },
        compounds() {
            const c = A().drugResponse?.compounds;
            return Array.isArray(c) ? c : [];
        },
        pathways() { try { return A()._WIKI_PATHWAYS ? A()._WIKI_PATHWAYS() : {}; } catch (e) { return {}; } },
        hallmarks() { try { return A()._WIKI_SUBTYPE_HALLMARKS ? A()._WIKI_SUBTYPE_HALLMARKS() : {}; } catch (e) { return {}; } },
        summary(id) {
            try {
                // Block tags become spaces first: textContent alone runs the
                // last word of one line into the first word of the next.
                const html = String(A()._cellLineSummaryText(id))
                    .replace(/<\/(div|p|li|h\d)>|<br\s*\/?>/gi, ' ');
                const d = document.createElement('div');
                d.innerHTML = html;
                return (d.textContent || '').replace(/\s+/g, ' ').trim();
            } catch (e) { return ''; }
        },
        lineageCounts() {
            const out = new Map();
            const ids = D.ids();
            for (let i = 0; i < ids.length; i++) {
                const l = D.lineage(ids[i]);
                if (!l) continue;
                out.set(l, (out.get(l) || 0) + 1);
            }
            return out;
        }
    };

    // Lineage transcription factors, paralog dependencies and drug targets
    // that a cell line can genuinely live or die by. Without them the strong
    // dependency question lands on genes nobody has heard of.
    const SELECTIVE = ['SOX10', 'MITF', 'PAX8', 'IRF4', 'POU2AF1', 'SPI1', 'GATA1', 'GATA2', 'GATA3',
        'GATA6', 'TAL1', 'LMO2', 'LYL1', 'MYB', 'MEF2C', 'CEBPA', 'RUNX1', 'HNF1A', 'HNF1B', 'HNF4A',
        'FOXA1', 'FOXA2', 'TP63', 'SOX2', 'NKX2-1', 'ASCL1', 'NEUROD1', 'POU2F3', 'ESR1', 'AR',
        'TFAP2C', 'SPDEF', 'ELF3', 'KLF5', 'GRHL2', 'CDX2', 'TCF7L2', 'CTNNB1', 'ZEB1', 'TEAD1',
        'YAP1', 'WWTR1', 'VGLL1', 'MYCN', 'MYC', 'MCL1', 'BCL2L1', 'BCL2', 'CCND1', 'CDK4', 'CDK6',
        'MDM2', 'WRN', 'SMARCA2', 'SMARCA4', 'ARID1B', 'EP300', 'CREBBP', 'PTPN11', 'SHOC2', 'RAF1',
        'SOS1', 'EGFR', 'ERBB2', 'ERBB3', 'MET', 'ALK', 'FGFR1', 'FGFR2', 'FGFR3', 'FLT3', 'JAK1',
        'JAK2', 'KIT', 'PDGFRA', 'RET', 'ABL1', 'BRAF', 'KRAS', 'NRAS', 'HRAS', 'PIK3CA', 'RPTOR',
        'RICTOR', 'XPO1', 'SALL4', 'TFAP2A', 'PRDM1', 'ETV1', 'ETV4', 'ETV5', 'ERG', 'FLI1', 'EWSR1',
        'NUTM1', 'RARA', 'PAX3', 'PAX7', 'FOXO1', 'SMARCB1', 'VHL', 'HIF1A', 'EPAS1', 'NFE2L2',
        'KEAP1', 'SREBF1', 'SLC7A11', 'GPX4', 'ADSL', 'UMPS', 'DHODH'];

    // Genes worth naming in a question: the curated well-known list plus the
    // lineage factors and drug targets above. Round 1 drew on every gene in
    // the pathway panels, which put unrecognisable symbols in the options.
    let NOTABLE = null;
    function notableGenes() {
        if (NOTABLE) return NOTABLE;
        const set = new Set();
        WELL_KNOWN_GENES.forEach(g => set.add(g));
        SELECTIVE.forEach(g => set.add(g));
        const a = A();
        NOTABLE = [...set].filter(g => !a._isPolymorphicLocus?.(g) && a.geneIndex.has(g.toUpperCase()));
        return NOTABLE;
    }
    // The tighter list, for questions where the gene itself is the point.
    let KNOWN_GENES = null;
    function knownGenes() {
        if (KNOWN_GENES) return KNOWN_GENES;
        const a = A();
        KNOWN_GENES = WELL_KNOWN_GENES.filter(g => a.geneIndex.has(g.toUpperCase()));
        return KNOWN_GENES;
    }

    // ------------------------------------------------------ curated subjects
    // Round 1 asked about genes and cell lines nobody recognizes. Every data
    // question now draws its subject from these lists, so the answer teaches
    // something a student would want to carry away.
    const WELL_KNOWN_GENES = ['TP53', 'KRAS', 'NRAS', 'HRAS', 'BRAF', 'EGFR', 'ERBB2', 'MYC', 'MYCN',
        'PIK3CA', 'PTEN', 'AKT1', 'MTOR', 'RB1', 'CDKN2A', 'CDK4', 'CDK6', 'CCND1', 'MDM2', 'MDM4',
        'APC', 'CTNNB1', 'VHL', 'BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'MLH1', 'MSH2', 'MSH6',
        'IDH1', 'IDH2', 'FLT3', 'KIT', 'ABL1', 'BCR', 'JAK2', 'STAT3', 'BCL2', 'MCL1', 'NOTCH1',
        'SMAD4', 'NF1', 'NF2', 'STK11', 'KEAP1', 'NFE2L2', 'ARID1A', 'SMARCA4', 'SMARCB1', 'EZH2',
        'KMT2D', 'CREBBP', 'TERT', 'AR', 'ESR1', 'GATA3', 'SOX2', 'SOX10', 'MITF', 'PAX5', 'IRF4',
        'MYB', 'WT1', 'RUNX1', 'CEBPA', 'NPM1', 'DNMT3A', 'TET2', 'ASXL1', 'SF3B1', 'ALK', 'RET',
        'ROS1', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'PDGFRA', 'NTRK1', 'CDX2', 'HNF1A', 'FOXA1',
        'TEAD1', 'YAP1', 'MAPK1', 'MAP2K1', 'RAF1', 'PARP1', 'TOP2A', 'TUBB', 'PCNA', 'MCM2'];

    // Display name first, then the spellings DepMap might have stored it under.
    // Anything that does not resolve is simply never used.
    const WELL_KNOWN_LINES = [
        ['K562'], ['A375'], ['MCF7'], ['MDA-MB-231'], ['HCT116'], ['A549'], ['HepG2'], ['PANC-1'],
        ['U2OS'], ['Jurkat'], ['Raji'], ['THP-1'], ['MOLM-13'], ['786-O'], ['LNCaP', 'LNCAPCLONEFGC'],
        ['SK-N-SH'], ['SH-SY5Y'], ['NCI-H1975'], ['HCC827'], ['SK-BR-3'],
        ['DU145'], ['PC-3'], ['T47D'], ['HT-29'], ['SW620'], ['LoVo'], ['MIA PaCa-2'], ['Capan-1'], ['Saos-2'], ['IMR-32'], ['Kelly'], ['NCI-H460'], ['NCI-H358'],
        ['Calu-6'], ['Calu-1'], ['OVCAR-3', 'NIHOVCAR3'], ['SK-OV-3'], ['Karpas-422'], ['OCI-LY19'],
        ['MV4-11'], ['Kasumi-1'], ['NB4'], ['Caki-1'], ['SK-MEL-2'], ['SK-MEL-30'], ['WM-266-4'],
        ['Malme-3M'], ['RPMI-8226'], ['MM.1S'], ['KMS-11'], ['HuH-7'], ['SK-HEP-1'], ['SNU-449'],
        ['AGS'], ['MKN45'], ['KATO III'], ['HT-1080'], ['RD'], ['A-673'], ['TC-71'], ['SK-ES-1'],
        ['LN-229'], ['U-251 MG'], ['T98G'], ['D283 Med'], ['Daoy'], ['SK-N-BE(2)', 'SKNBE2'],
        ['Y79'], ['WERI-Rb-1'], ['SW48'], ['MELJUSO'], ['COLO800']
    ];

    const FAMOUS_FUSIONS = ['BCR-ABL1', 'EWSR1-FLI1', 'EML4-ALK', 'PML-RARA', 'TMPRSS2-ERG',
        'KMT2A-MLLT3', 'RUNX1-RUNX1T1', 'CBFB-MYH11', 'ETV6-RUNX1', 'FGFR3-TACC3', 'NPM1-ALK',
        'SS18-SSX1', 'PAX3-FOXO1'];

    const WELL_KNOWN_DRUGS = ['imatinib', 'dasatinib', 'nilotinib', 'ponatinib', 'vemurafenib',
        'dabrafenib', 'encorafenib', 'trametinib', 'cobimetinib', 'binimetinib', 'erlotinib',
        'gefitinib', 'osimertinib', 'afatinib', 'lapatinib', 'neratinib', 'palbociclib',
        'ribociclib', 'abemaciclib', 'olaparib', 'niraparib', 'talazoparib', 'venetoclax',
        'navitoclax', 'bortezomib', 'carfilzomib', 'doxorubicin', 'paclitaxel', 'docetaxel',
        'carboplatin', 'gemcitabine', 'cytarabine', 'fluorouracil', '5-fluorouracil',
        'methotrexate', 'etoposide', 'irinotecan', 'topotecan', 'vincristine', 'everolimus',
        'temsirolimus', 'crizotinib', 'alectinib', 'lorlatinib', 'regorafenib', 'cabozantinib',
        'alpelisib', 'ruxolitinib', 'gilteritinib', 'enzalutamide', 'fulvestrant', 'idasanutlin',
        'azacitidine', 'decitabine', 'lenalidomide', 'tazemetostat', 'belinostat', 'romidepsin'];

    // A hotspot frequency question is only worth asking where the answer is
    // also the textbook answer. Anything else is trivia about this panel.
    const TEXTBOOK_HOTSPOT_TISSUE = {
        BRAF: 'Skin', KRAS: 'Pancreas', VHL: 'Kidney', APC: 'Bowel', EGFR: 'Lung',
        IDH1: 'CNS/Brain', FLT3: 'Myeloid', NOTCH1: 'Lymphoid', CTNNB1: 'Liver',
        AR: 'Prostate', PIK3CA: 'Breast', TP53: 'Esophagus/Stomach'
    };

    // Why a pair of gene effects moves together, in words. Used by the
    // correlation questions instead of a generic "they correlate".
    const PAIR_WHY = {
        'TP53|MDM2': 'TP53 and MDM2 pull opposite ways: a cell line with working p53 needs MDM2 to keep it in check, so it dies without MDM2, while a p53-mutant line does not care.',
        'TP53|MDM4': 'MDM4 is the second brake on p53, so lines with working p53 lean on it the same way they lean on MDM2.',
        'TP53|CDKN1A': 'CDKN1A (p21) is switched on by p53, so the two track the same working pathway.',
        'KRAS|RAF1': 'RAF1 sits directly downstream of KRAS, so a line driven by KRAS also needs the kinase that carries its signal on.',
        'BRAF|MAPK1': 'MAPK1 (ERK2) is the end of the pathway BRAF drives, so BRAF-driven lines need both.',
        'BRAF|SOX10': 'Melanoma lines are driven by BRAF and kept in their lineage by SOX10, so the two dependencies show up in the same cells.',
        'CCND1|CDK4': 'Cyclin D1 works only with CDK4, so a line that needs one needs the other.',
        'CCND1|CDK6': 'Cyclin D1 partners CDK6 as well as CDK4, so the pair rises and falls together.',
        'RB1|CDK4': 'RB1 is the brake CDK4 releases. Lines that lost RB1 no longer need CDK4, so the two run opposite.',
        'CDKN2A|CDK4': 'CDKN2A (p16) blocks CDK4. Losing p16 leaves the line leaning on CDK4 instead.',
        'TSC1|TSC2': 'TSC1 and TSC2 only work as a pair, so knocking out either has the same consequence.',
        'MYC|MAX': 'MYC cannot bind DNA without MAX, so the two are needed by the same cell lines.',
        'PAX5|EBF1': 'PAX5 and EBF1 together hold a B cell in its identity, so B-cell lines depend on both.',
        'SDHA|SDHB': 'SDHA and SDHB are two subunits of the same respiratory complex.',
        'BRCA1|BARD1': 'BRCA1 is only stable as a pair with BARD1, so the two behave as one gene.',
        'MLH1|PMS2': 'MLH1 and PMS2 form one mismatch repair complex, so they share a fate.',
        'MSH2|MSH6': 'MSH2 and MSH6 form one mismatch repair complex, so they share a fate.'
    };
    const pairWhy = (g1, g2) => PAIR_WHY[g1 + '|' + g2] || PAIR_WHY[g2 + '|' + g1] || '';

    // Genes with a dependency that follows the lineage, used for the box plots.
    const LINEAGE_GENES = ['SOX10', 'MITF', 'PAX5', 'IRF4', 'MYB', 'CTNNB1', 'HNF1A', 'FOXA1',
        'GATA3', 'AR', 'ESR1', 'CDX2', 'SOX2', 'NKX2-1', 'TEAD1', 'KLF5', 'SPI1', 'GATA1',
        'HNF4A', 'TP63', 'POU2AF1', 'MYCN'];

    let LINE_IDS = null;
    function resolvedLines() {
        if (LINE_IDS) return LINE_IDS;
        LINE_IDS = new Map();
        let map;
        try { map = A()._buildCellLineNameToIdMap(); } catch (e) { map = new Map(); }
        for (const entry of WELL_KNOWN_LINES) {
            const shown = entry[0];
            for (const cand of entry) {
                // DepMap stores some names stripped of punctuation, so the
                // plain name, the de-hyphenated name and the letters-and-digits
                // form are all tried before the line is given up on.
                const tries = [cand, cand.replace(/[-\s]/g, ''), cand.replace(/[^A-Za-z0-9]/g, '')];
                let id = null;
                for (const t of tries) {
                    const hit = map.get(t.toUpperCase());
                    if (hit) { id = hit; break; }
                }
                if (id) { LINE_IDS.set(shown, id); break; }
            }
        }
        return LINE_IDS;
    }
    const lineId = (name) => resolvedLines().get(name) || null;
    let FAMOUS_ID_SET = null;
    function famousIds() {
        if (!FAMOUS_ID_SET) FAMOUS_ID_SET = new Set([...resolvedLines().values()]);
        return FAMOUS_ID_SET;
    }

    // Compound names arrive in capitals; sentence case reads better in a
    // question and matches how the app's own drug tables show them.
    const drugName = (n) => String(n || '').replace(/[A-Z][A-Z0-9-]*/g,
        w => w.charAt(0) + w.slice(1).toLowerCase());

    let KNOWN_DRUGS = null;
    function knownCompounds() {
        if (KNOWN_DRUGS) return KNOWN_DRUGS;
        const want = new Set(WELL_KNOWN_DRUGS.map(d => d.toUpperCase().replace(/[^A-Z0-9]/g, '')));
        KNOWN_DRUGS = D.compounds().filter(c =>
            want.has(String(c.name || '').toUpperCase().replace(/[^A-Z0-9]/g, '')));
        return KNOWN_DRUGS;
    }

    // ==================================================== hand-written bank
    // Facts a biologist or student would actually want: how to read the data
    // types in this app, the big cancer genes, the famous lines, and what the
    // numbers mean. One fact per question, plain words, no trick questions.
    // cat = category, o = options, a = index of the true one, e = explanation,
    // level 1 basic / 2 intermediate / 3 advanced, line = famous cell line to
    // offer a wiki link for, tissue = cohort this question suits.
    const CONCEPTS = [
        // ------------------------------------------------ READING GENE EFFECT
        {
            id: 'ge1', cat: 'READING GENE EFFECT', level: 1,
            q: 'What does a CRISPR knockout screen measure?',
            o: ['How much a cell line needs each gene in order to grow',
                'How much RNA each gene makes in the cell line',
                'How many copies of each gene the cell line carries',
                'How quickly the cell line takes up a drug'],
            a: 0,
            e: 'Every gene is knocked out in turn and the cells are grown on. Genes the cells cannot do without fall away, and that drop is the gene effect score.'
        },
        {
            id: 'ge2', cat: 'READING GENE EFFECT', level: 1,
            q: 'On the gene effect scale, what does a score of 0 mean?',
            o: ['Losing the gene made no difference to growth',
                'The gene was not measured in that cell line',
                'The cell line cannot survive without the gene',
                'The gene is switched off in that cell line'],
            a: 0,
            e: 'The scale is anchored so that 0 is the behavior of a gene the cells do not need. No data is shown as missing, not as zero.'
        },
        {
            id: 'ge3', cat: 'READING GENE EFFECT', level: 1,
            q: 'A gene effect score of about -1 means what?',
            o: ['The gene is needed about as much as a typical essential gene',
                'The gene is needed about half as much as an average gene',
                'The cells grow better once the gene is gone',
                'The gene was measured in only half of the cell lines'],
            a: 0,
            e: 'The scale is set by common essential genes, which sit near -1. That makes -1 a useful landmark: below it is a real dependency.'
        },
        {
            id: 'ge4', cat: 'READING GENE EFFECT', level: 2,
            q: 'A gene effect score above 0, say +0.4, means what?',
            o: ['The cells grew a little better without the gene',
                'The gene was knocked out twice over',
                'The gene is essential in that cell line',
                'The measurement failed and should be ignored'],
            a: 0,
            e: 'A positive score means losing the gene helped growth. Tumor suppressors often look like this in the lines that still carry them.'
        },
        {
            id: 'ge5', cat: 'READING GENE EFFECT', level: 1,
            q: 'Which is the stronger dependency, a gene effect of -0.3 or one of -1.8?',
            o: ['-1.8, because more negative means a stronger need',
                '-0.3, because it is closer to zero',
                'They describe equally strong needs',
                'Neither, both are above the essential range'],
            a: 0,
            e: 'The more negative the score, the harder the cells were hit by losing that gene. -1.8 is a strong dependency, -0.3 is mild.'
        },
        {
            id: 'ge6', cat: 'READING GENE EFFECT', level: 1,
            q: 'What does it mean to call a gene "pan-essential"?',
            o: ['Almost every cell line needs it, whatever the cancer',
                'It is needed by exactly one kind of cancer',
                'It is expressed in every tissue of the body',
                'It is mutated in most human tumors'],
            a: 0,
            e: 'Pan-essential genes run the basic machinery of a cell, so knocking one out kills nearly any line. They are poor drug targets for that reason.'
        },
        {
            id: 'ge7', cat: 'READING GENE EFFECT', level: 1,
            q: 'Which group of genes is pan-essential in almost every cell line?',
            o: ['Ribosomal proteins', 'Olfactory receptors', 'Keratins', 'Antibody genes'],
            a: 0,
            e: 'Ribosomes, the proteasome and RNA polymerase are needed by any growing cell, so their genes come out strongly negative everywhere.'
        },
        {
            id: 'ge8', cat: 'READING GENE EFFECT', level: 2,
            q: 'Why is a selective dependency more interesting for drug discovery than a pan-essential gene?',
            o: ['It kills some tumors while leaving normal cells alone',
                'It is easier to make an antibody against',
                'It gives a larger drop in the screen',
                'It is always found on the same chromosome'],
            a: 0,
            e: 'A drug against a pan-essential gene poisons healthy tissue too. A selective dependency marks a weakness that only the tumor has.'
        },
        {
            id: 'ge9', cat: 'READING GENE EFFECT', level: 1,
            q: 'Which institute runs DepMap, the project behind this data?',
            o: ['The Broad Institute', 'The NCBI', 'The EMBL-EBI', 'Cold Spring Harbor Laboratory'],
            a: 0,
            e: 'DepMap is run at the Broad Institute in Cambridge, Massachusetts. It screens hundreds of cancer cell lines and releases the results openly.'
        },
        {
            id: 'ge10', cat: 'READING GENE EFFECT', level: 2,
            q: 'What does a z-score of a gene effect, taken against all cell lines, tell you?',
            o: ['How far this line sits from the average line, in standard deviations',
                'How many cell lines carry a mutation in that gene',
                'How well the gene was measured in that line',
                'How strongly the gene is expressed in that line'],
            a: 0,
            e: 'A z-score puts one line in the context of the whole panel. It answers "is this unusual?", which the raw score on its own cannot.'
        },
        {
            id: 'ge11', cat: 'READING GENE EFFECT', level: 2,
            q: 'A gene has a gene effect z-score of -3 in one cell line. What does that say?',
            o: ['That line needs the gene far more than the average line does',
                'The gene is deleted in that line',
                'The gene is essential in every line',
                'The gene is expressed three times higher there'],
            a: 0,
            e: 'Three standard deviations below the panel average is a strong outlier. This is exactly the shape of a selective dependency worth chasing.'
        },
        {
            id: 'ge12', cat: 'READING GENE EFFECT', level: 3,
            q: 'Why is comparing gene effects within one tissue different from comparing across the whole panel?',
            o: ['Tissue of origin drives many dependencies, so a panel-wide difference can just be lineage',
                'Gene effects are measured on a different scale inside a tissue',
                'Only one tissue at a time was screened',
                'Correlations cannot be calculated on fewer than 500 lines'],
            a: 0,
            e: 'Blood lines and solid tumor lines differ in thousands of genes. Staying inside one lineage asks a sharper question: what differs between lines of the same kind.'
        },
        {
            id: 'ge13', cat: 'READING GENE EFFECT', level: 2,
            q: 'What is actually counted at the end of a pooled CRISPR screen?',
            o: ['How much each guide RNA has grown or shrunk in the pool',
                'How much protein each gene made',
                'How many cells changed shape',
                'How many mutations appeared in each gene'],
            a: 0,
            e: 'Each cell carries one guide. Sequencing the pool before and after growth shows which guides were lost, and those name the genes the cells needed.'
        },
        {
            id: 'ge14', cat: 'READING GENE EFFECT', level: 3,
            q: 'Chronos, the method behind these scores, corrects for one big artifact of CRISPR screens. Which?',
            o: ['Cutting an amplified region many times is toxic on its own',
                'Guide RNAs bind the wrong gene',
                'Cells grow at different speeds',
                'Some genes have no guides designed against them'],
            a: 0,
            e: 'A gene present in many copies is cut many times, and the DNA damage alone slows the cells. Without that correction, amplified regions look falsely essential.'
        },
        {
            id: 'ge15', cat: 'READING GENE EFFECT', level: 2,
            q: 'Why does a dependency map need hundreds of cell lines rather than a few?',
            o: ['A dependency only looks selective when many different backgrounds are compared',
                'The screen fails in fewer than 100 lines',
                'Each line can only be screened for a few genes',
                'Statistics need at least 500 samples to work'],
            a: 0,
            e: 'One line tells you what that line needs. Hundreds of lines tell you which needs follow a mutation, a lineage or a pathway, which is the useful part.'
        },
        {
            id: 'ge16', cat: 'READING GENE EFFECT', level: 1,
            q: 'A cell line has a gene effect near 0 for a gene. What is the plainest reading?',
            o: ['The cell line grows just as well without that gene',
                'The gene is not present in that cell line',
                'The gene is mutated in that cell line',
                'The gene is needed only in later passages'],
            a: 0,
            e: 'Near zero means no measurable change in growth. Many genes look like this in most lines, which is why the interesting ones stand out.'
        },
        {
            id: 'ge17', cat: 'READING GENE EFFECT', level: 2,
            q: 'A gene is strongly expressed in a cell line but its gene effect is 0. What does that show?',
            o: ['Expression and dependency are separate measurements',
                'The expression measurement must be wrong',
                'The gene must be mutated',
                'The screen missed that gene'],
            a: 0,
            e: 'Plenty of genes are made in quantity and still not needed, often because another gene covers for them. Being present is not the same as being required.'
        },

        // -------------------------------------------------------- CANCER GENES
        {
            id: 'cg1', cat: 'CANCER GENES', level: 1,
            q: 'What is an oncogene?',
            o: ['A gene that drives cancer when it gains activity',
                'A gene that drives cancer when it is lost',
                'A gene found only in tumor cells',
                'A gene that repairs damaged DNA'],
            a: 0,
            e: 'Oncogenes are accelerators. One activating change in one copy is usually enough, which is why hotspot mutations cluster at the same few positions.'
        },
        {
            id: 'cg2', cat: 'CANCER GENES', level: 1,
            q: 'What is a tumor suppressor gene?',
            o: ['A gene that restrains growth and lets cancer through when it is lost',
                'A gene that drives growth when it is switched on too hard',
                'A gene that only works in normal tissue',
                'A gene that makes cells resistant to drugs'],
            a: 0,
            e: 'Tumor suppressors are brakes. Both copies usually have to be broken, so they are hit by stop codons, frameshifts and deletions rather than hotspots.'
        },
        {
            id: 'cg3', cat: 'CANCER GENES', level: 1,
            q: 'Which gene is mutated in more human cancers than any other?',
            o: ['TP53', 'BRCA1', 'ALK', 'NOTCH1'],
            a: 0,
            e: 'TP53 is mutated in roughly half of all human tumors. Losing p53 removes the response that would otherwise stop or kill a damaged cell.'
        },
        {
            id: 'cg4', cat: 'CANCER GENES', level: 1,
            q: 'The TP53 gene makes which protein?',
            o: ['p53', 'p16', 'p21', 'p110'],
            a: 0,
            e: 'p53 senses DNA damage and other stress, then halts the cell cycle or triggers cell death. It is often called the guardian of the genome.'
        },
        {
            id: 'cg5', cat: 'CANCER GENES', level: 2,
            q: 'About what share of pancreatic ductal cancers carry a KRAS mutation?',
            o: ['About 95 percent', 'About 50 percent', 'About 20 percent', 'About 5 percent'],
            a: 0,
            e: 'KRAS mutation is almost the definition of pancreatic ductal cancer. It is also common in bowel and lung cancer, though at lower rates.'
        },
        {
            id: 'cg6', cat: 'CANCER GENES', level: 1, tissue: 'Skin',
            q: 'BRAF V600E is found in about half of which cancer?',
            o: ['Melanoma', 'Pancreatic cancer', 'Prostate cancer', 'Acute myeloid leukemia'],
            a: 0,
            e: 'The V600E change locks BRAF on, driving the MAPK pathway. It is the reason BRAF inhibitors were developed for melanoma.'
        },
        {
            id: 'cg7', cat: 'CANCER GENES', level: 1, tissue: 'Lung',
            q: 'Activating EGFR mutations are most typical of which cancer?',
            o: ['Lung adenocarcinoma', 'Bowel cancer', 'Ovarian cancer', 'Chronic myeloid leukemia'],
            a: 0,
            e: 'EGFR mutations in lung adenocarcinoma made the first targeted lung drugs possible, and they are more common in patients who never smoked.'
        },
        {
            id: 'cg8', cat: 'CANCER GENES', level: 1,
            q: 'HER2, the target of trastuzumab, is made by which gene?',
            o: ['ERBB2', 'EGFR', 'ESR1', 'ERBB4'],
            a: 0,
            e: 'HER2 and ERBB2 are two names for the same thing. The gene is amplified in about one breast cancer in five.'
        },
        {
            id: 'cg9', cat: 'CANCER GENES', level: 2,
            q: 'What kind of protein does MYC make?',
            o: ['A transcription factor that drives growth programs',
                'A receptor on the cell surface',
                'A DNA repair enzyme',
                'A pump that removes drugs from the cell'],
            a: 0,
            e: 'MYC switches on the genes a cell needs to grow and divide. It is amplified or overdriven in a large share of cancers, and is hard to drug directly.'
        },
        {
            id: 'cg10', cat: 'CANCER GENES', level: 2,
            q: 'What does the RB1 protein normally do?',
            o: ['Holds the cell back from entering S phase',
                'Cuts damaged DNA out of the genome',
                'Carries growth signals from the surface to the nucleus',
                'Pumps calcium out of the cell'],
            a: 0,
            e: 'RB1 binds E2F and keeps the cell cycle shut until CDK4 and CDK6 release it. Losing RB1 removes that brake for good.'
        },
        {
            id: 'cg11', cat: 'CANCER GENES', level: 2,
            q: 'PTEN loss switches on which pathway?',
            o: ['PI3K and AKT', 'MAPK and ERK', 'Wnt and beta-catenin', 'Hedgehog'],
            a: 0,
            e: 'PTEN removes the phosphate that PI3K adds. Without PTEN the AKT signal stays on, pushing growth and survival.'
        },
        {
            id: 'cg12', cat: 'CANCER GENES', level: 3,
            q: 'The CDKN2A locus encodes two different proteins. Which pair?',
            o: ['p16INK4a and p14ARF', 'p53 and p21', 'p27 and p57', 'p110 and p85'],
            a: 0,
            e: 'Two reading frames share the locus. p16 blocks CDK4 and CDK6, p14ARF protects p53, so one deletion disables both major brakes at once.'
        },
        {
            id: 'cg13', cat: 'CANCER GENES', level: 2,
            q: 'What does MDM2 do to p53?',
            o: ['Tags it for destruction, keeping its level low',
                'Cuts it into an active fragment',
                'Carries it into the nucleus',
                'Locks it onto DNA'],
            a: 0,
            e: 'MDM2 is an E3 ligase that marks p53 for the proteasome. Cell lines with working p53 therefore lean on MDM2 to survive.'
        },
        {
            id: 'cg14', cat: 'CANCER GENES', level: 1, tissue: 'Bowel',
            q: 'Loss of APC is the classic first step in which cancer?',
            o: ['Bowel cancer', 'Kidney cancer', 'Melanoma', 'Myeloma'],
            a: 0,
            e: 'APC normally holds beta-catenin down. Losing it turns on the Wnt program, and that is where most colorectal tumors begin.'
        },
        {
            id: 'cg15', cat: 'CANCER GENES', level: 2, tissue: 'Kidney',
            q: 'Loss of VHL is the hallmark of which cancer?',
            o: ['Clear cell kidney cancer', 'Small cell lung cancer', 'Gastric cancer', 'Burkitt lymphoma'],
            a: 0,
            e: 'VHL normally destroys HIF when oxygen is plentiful. Without VHL the cell behaves as if starved of oxygen and builds new blood vessels.'
        },
        {
            id: 'cg16', cat: 'CANCER GENES', level: 2,
            q: 'Tumors that have lost BRCA1 or BRCA2 are unusually sensitive to which drug class?',
            o: ['PARP inhibitors', 'MEK inhibitors', 'Proteasome inhibitors', 'Aromatase inhibitors'],
            a: 0,
            e: 'BRCA loss removes one DNA repair route. Blocking PARP takes away the backup as well, and the cell cannot cope with both gone.'
        },
        {
            id: 'cg17', cat: 'CANCER GENES', level: 2,
            q: 'What does synthetic lethality mean?',
            o: ['Two changes are each survivable, but together they kill the cell',
                'A drug is lethal only at synthetic doses',
                'A gene is lethal when knocked out in every cell',
                'Two drugs given together are toxic to the patient'],
            a: 0,
            e: 'It is the logic behind PARP inhibitors in BRCA-mutant tumors: the mutation is already there, and the drug supplies the second hit.'
        },
        {
            id: 'cg18', cat: 'CANCER GENES', level: 3,
            q: 'Mutant IDH1 and IDH2 produce an unusual product. What is it?',
            o: ['2-hydroxyglutarate, an oncometabolite that blocks demethylases',
                'Extra lactate that acidifies the tumor',
                'A truncated protein that cannot fold',
                'Reactive oxygen that damages DNA'],
            a: 0,
            e: 'The mutant enzyme makes 2-hydroxyglutarate instead of its normal product. That metabolite jams the enzymes that remove methyl marks, freezing the cell in an immature state.'
        },
        {
            id: 'cg19', cat: 'CANCER GENES', level: 3,
            q: 'What does NF1 normally do?',
            o: ['Switches RAS off, acting as a brake on the pathway',
                'Switches RAS on in response to growth factors',
                'Repairs double strand breaks',
                'Controls the spindle checkpoint'],
            a: 0,
            e: 'NF1 is a GTPase activating protein for RAS. Losing it leaves RAS in its active form, which is why NF1 loss behaves much like a RAS mutation.'
        },

        // --------------------------------------------- FUSIONS AND CHROMOSOMES
        {
            id: 'fu1', cat: 'FUSIONS AND CHROMOSOMES', level: 1,
            q: 'What is a gene fusion?',
            o: ['Two genes joined into one after a chromosome break',
                'A gene present in twice the normal number of copies',
                'Two identical copies of a gene side by side',
                'A gene copied from RNA back into DNA'],
            a: 0,
            e: 'A break and rejoin puts the front of one gene onto the back of another. The hybrid protein often has activity the cell never asked for.'
        },
        {
            id: 'fu2', cat: 'FUSIONS AND CHROMOSOMES', level: 1,
            q: 'The BCR-ABL1 fusion is made by a swap between which two chromosomes?',
            o: ['9 and 22', '8 and 14', '15 and 17', '11 and 22'],
            a: 0,
            e: 'The t(9;22) swap makes the small marker chromosome first seen in Philadelphia, and puts the ABL1 kinase under the control of BCR.'
        },
        {
            id: 'fu3', cat: 'FUSIONS AND CHROMOSOMES', level: 1, tissue: 'Myeloid',
            q: 'The Philadelphia chromosome is the hallmark of which disease?',
            o: ['Chronic myeloid leukemia', 'Ewing sarcoma', 'Neuroblastoma', 'Multiple myeloma'],
            a: 0,
            e: 'It was the first chromosome change tied to a specific cancer, and it made chronic myeloid leukemia the first disease treated with a targeted kinase inhibitor.'
        },
        {
            id: 'fu4', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'The EWSR1-FLI1 fusion drives which tumor?',
            o: ['Ewing sarcoma', 'Synovial sarcoma', 'Rhabdomyosarcoma', 'Osteosarcoma'],
            a: 0,
            e: 'EWSR1-FLI1 acts as a rogue transcription factor. Ewing sarcoma lines depend on it heavily, which makes it a textbook example of fusion addiction.'
        },
        {
            id: 'fu5', cat: 'FUSIONS AND CHROMOSOMES', level: 2, tissue: 'Lung',
            q: 'The EML4-ALK fusion is found in which cancer?',
            o: ['Lung adenocarcinoma', 'Prostate cancer', 'Acute myeloid leukemia', 'Gastric cancer'],
            a: 0,
            e: 'EML4-ALK appears in a small share of lung adenocarcinomas, and those tumors respond well to ALK inhibitors such as crizotinib and alectinib.'
        },
        {
            id: 'fu6', cat: 'FUSIONS AND CHROMOSOMES', level: 2, tissue: 'Myeloid',
            q: 'The PML-RARA fusion defines which leukemia?',
            o: ['Acute promyelocytic leukemia', 'Chronic lymphocytic leukemia',
                'Hairy cell leukemia', 'Chronic myeloid leukemia'],
            a: 0,
            e: 'The fusion blocks myeloid cells from maturing. Recognising it turned a rapidly fatal leukemia into one of the most curable.'
        },
        {
            id: 'fu7', cat: 'FUSIONS AND CHROMOSOMES', level: 2, tissue: 'Myeloid',
            q: 'Which treatment works especially well in PML-RARA leukemia?',
            o: ['All-trans retinoic acid', 'Imatinib', 'Trastuzumab', 'Olaparib'],
            a: 0,
            e: 'Retinoic acid pushes the blocked cells to mature instead of killing them, and combined with arsenic it cures most patients without standard chemotherapy.'
        },
        {
            id: 'fu8', cat: 'FUSIONS AND CHROMOSOMES', level: 2, tissue: 'Prostate',
            q: 'The TMPRSS2-ERG fusion is common in which cancer?',
            o: ['Prostate cancer', 'Breast cancer', 'Bowel cancer', 'Melanoma'],
            a: 0,
            e: 'About half of prostate tumors carry it. The androgen-driven TMPRSS2 promoter is placed in front of ERG, so hormone signals now drive a transcription factor.'
        },
        {
            id: 'fu9', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'What does whole genome doubling mean?',
            o: ['The entire chromosome set was duplicated at some point',
                'One chromosome arm was duplicated',
                'The DNA was sequenced twice',
                'Every gene is expressed twice as strongly'],
            a: 0,
            e: 'A cell that fails to divide properly ends up with two copies of everything. It is common in tumors and tends to be followed by more chromosome chaos.'
        },
        {
            id: 'fu10', cat: 'FUSIONS AND CHROMOSOMES', level: 1,
            q: 'What is aneuploidy?',
            o: ['Having an abnormal number of chromosomes or chromosome arms',
                'Having no chromosomes at all in some cells',
                'Carrying two mutations in the same gene',
                'Having unusually long telomeres'],
            a: 0,
            e: 'Most cancer cells are aneuploid. The imbalance changes the dose of hundreds of genes at once, quite apart from any single mutation.'
        },
        {
            id: 'fu11', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'Microsatellite instability is caused by the loss of which system?',
            o: ['DNA mismatch repair', 'Nucleotide excision repair',
                'Homologous recombination', 'The spindle checkpoint'],
            a: 0,
            e: 'Mismatch repair fixes the slippage that happens in short repeats. Without it those repeats change length all over the genome, which is what the MSI test detects.'
        },
        {
            id: 'fu12', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'Loss of which genes most often causes microsatellite instability?',
            o: ['MLH1 and MSH2', 'BRCA1 and BRCA2', 'ATM and ATR', 'TP53 and RB1'],
            a: 0,
            e: 'MLH1, MSH2, MSH6 and PMS2 make up the mismatch repair machinery. MLH1 is often silenced by methylation rather than mutated.'
        },
        {
            id: 'fu13', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'Why does a driver fusion make an attractive drug target?',
            o: ['It exists only in the tumor cells, not in healthy tissue',
                'It is always found on the cell surface',
                'It is easier to sequence than a point mutation',
                'It never becomes resistant to treatment'],
            a: 0,
            e: 'A fusion protein has no counterpart in normal cells, so a drug against it can hit the tumor and leave the rest of the body alone.'
        },
        {
            id: 'fu14', cat: 'FUSIONS AND CHROMOSOMES', level: 3, tissue: 'Myeloid',
            q: 'RUNX1-RUNX1T1, from the t(8;21) swap, is found in which disease?',
            o: ['Acute myeloid leukemia', 'Ewing sarcoma', 'Follicular lymphoma', 'Neuroblastoma'],
            a: 0,
            e: 'It is one of the core binding factor leukemias. The fusion blocks the normal RUNX1 program that myeloid cells need in order to mature.'
        },
        {
            id: 'fu15', cat: 'FUSIONS AND CHROMOSOMES', level: 3,
            q: 'The PAX3-FOXO1 fusion marks which tumor?',
            o: ['Alveolar rhabdomyosarcoma', 'Ewing sarcoma', 'Synovial sarcoma', 'Osteosarcoma'],
            a: 0,
            e: 'It is used to separate alveolar rhabdomyosarcoma from the embryonal form, which has no such fusion and behaves better.'
        },
        {
            id: 'fu16', cat: 'FUSIONS AND CHROMOSOMES', level: 2,
            q: 'What separates a driver alteration from a passenger one?',
            o: ['A driver contributes to the cancer growing, a passenger came along by chance',
                'A driver is inherited, a passenger is acquired',
                'A driver is found in DNA, a passenger only in RNA',
                'A driver is always a fusion, a passenger always a point mutation'],
            a: 0,
            e: 'Tumors carry many changes and most do nothing. Telling the few that matter from the rest is much of what cancer genomics is for.'
        },

        // --------------------------------------------------- DRUGS AND TARGETS
        {
            id: 'dr1', cat: 'DRUGS AND TARGETS', level: 1, tissue: 'Myeloid',
            q: 'What does imatinib target?',
            o: ['The BCR-ABL1 kinase', 'The EGFR receptor', 'The proteasome', 'Microtubules'],
            a: 0,
            e: 'Imatinib was the first drug designed against a specific fusion protein, and it turned chronic myeloid leukemia into a manageable disease.'
        },
        {
            id: 'dr2', cat: 'DRUGS AND TARGETS', level: 1, tissue: 'Skin',
            q: 'What does vemurafenib target?',
            o: ['BRAF carrying the V600E change', 'KRAS carrying G12C', 'MEK', 'ALK fusions'],
            a: 0,
            e: 'Vemurafenib binds the mutant BRAF kinase. It shrinks melanomas quickly, although resistance through the same pathway usually follows.'
        },
        {
            id: 'dr3', cat: 'DRUGS AND TARGETS', level: 1, tissue: 'Lung',
            q: 'Erlotinib and osimertinib both target which protein?',
            o: ['EGFR', 'HER2', 'ALK', 'MET'],
            a: 0,
            e: 'Both block the EGFR kinase. Osimertinib is the later drug and reaches mutations that made the earlier ones stop working.'
        },
        {
            id: 'dr4', cat: 'DRUGS AND TARGETS', level: 1, tissue: 'Breast',
            q: 'What does trastuzumab target?',
            o: ['HER2 on the cell surface', 'The estrogen receptor', 'PARP', 'CDK4 and CDK6'],
            a: 0,
            e: 'Trastuzumab is an antibody against HER2, so it only helps tumors that carry extra copies of the ERBB2 gene.'
        },
        {
            id: 'dr5', cat: 'DRUGS AND TARGETS', level: 2, tissue: 'Breast',
            q: 'What does palbociclib target?',
            o: ['CDK4 and CDK6', 'CDK1 and CDK2', 'PARP1', 'MEK1 and MEK2'],
            a: 0,
            e: 'Blocking CDK4 and CDK6 keeps RB1 in its restraining form, so the cell cannot start a new round of division. It needs RB1 to be intact to work.'
        },
        {
            id: 'dr6', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'What does olaparib target?',
            o: ['PARP', 'BRCA1', 'ATM', 'The proteasome'],
            a: 0,
            e: 'Olaparib blocks and traps PARP on DNA. Cells with working BRCA repair the damage, cells without it cannot.'
        },
        {
            id: 'dr7', cat: 'DRUGS AND TARGETS', level: 2, tissue: 'Lymphoid',
            q: 'What does venetoclax target?',
            o: ['BCL2', 'MCL1', 'BTK', 'The proteasome'],
            a: 0,
            e: 'Venetoclax frees the death signal that BCL2 was holding back. Lymphoid tumors that lean on BCL2 are strikingly sensitive to it.'
        },
        {
            id: 'dr8', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'What does trametinib target?',
            o: ['MEK, one step below BRAF', 'BRAF itself', 'ERK, one step below MEK', 'RAS'],
            a: 0,
            e: 'Because MEK sits under BRAF in the same pathway, MEK and BRAF inhibitors are often given together to delay resistance.'
        },
        {
            id: 'dr9', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'In a drug screen, an AUC close to 0 for a cell line means what?',
            o: ['The drug killed the cells at almost every dose tested',
                'The drug had no effect at any dose',
                'The drug was never tested on that line',
                'The cells grew faster with the drug'],
            a: 0,
            e: 'AUC is the area under the dose response curve. Low means the curve is pressed to the floor, so the cells died across the dose range.'
        },
        {
            id: 'dr10', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'In a drug screen, an AUC close to 1 for a cell line means what?',
            o: ['The cells carried on growing as if untreated',
                'The drug killed the cells at every dose',
                'The measurement failed',
                'The drug worked only at the highest dose'],
            a: 0,
            e: 'An AUC near 1 means the survival curve stayed flat at the top. That line is resistant to the compound over the doses tested.'
        },
        {
            id: 'dr11', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'Why does killing a cell line in a dish not predict clinical response on its own?',
            o: ['A patient adds drug exposure, immune response and tumor environment to the picture',
                'Cell lines are always more resistant than tumors',
                'Dishes cannot measure cell death',
                'Only mouse experiments count as evidence'],
            a: 0,
            e: 'A dish measures one thing well: whether the cells need the target. Getting the drug there safely, and holding the response, are separate problems.'
        },
        {
            id: 'dr12', cat: 'DRUGS AND TARGETS', level: 2, tissue: 'Lung',
            q: 'Crizotinib is used against tumors carrying which alteration?',
            o: ['ALK fusions', 'BRAF V600E', 'BRCA1 loss', 'KRAS G12D'],
            a: 0,
            e: 'Crizotinib blocks the ALK kinase, and it also hits ROS1 and MET, so it is used across several fusion-driven lung cancers.'
        },
        {
            id: 'dr13', cat: 'DRUGS AND TARGETS', level: 2, tissue: 'Plasma Cell',
            q: 'What does bortezomib block?',
            o: ['The proteasome', 'The ribosome', 'Topoisomerase II', 'The spindle'],
            a: 0,
            e: 'Myeloma cells make huge amounts of antibody and depend on clearing the misfolded surplus, so blocking the proteasome hits them hardest.'
        },
        {
            id: 'dr14', cat: 'DRUGS AND TARGETS', level: 1,
            q: 'How does cisplatin damage cancer cells?',
            o: ['It cross-links DNA so the strands cannot be copied',
                'It blocks the ribosome',
                'It stops microtubules from forming',
                'It blocks a growth factor receptor'],
            a: 0,
            e: 'Platinum drugs bind DNA directly and stall replication. They are not targeted, which is why they work broadly and cause broad side effects.'
        },
        {
            id: 'dr15', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'How does paclitaxel work?',
            o: ['It locks microtubules in place so the cell cannot divide',
                'It breaks microtubules apart',
                'It blocks DNA synthesis',
                'It blocks the proteasome'],
            a: 0,
            e: 'Microtubules have to build and break down for a cell to pull its chromosomes apart. Paclitaxel freezes them, and the cell arrests in mitosis.'
        },
        {
            id: 'dr16', cat: 'DRUGS AND TARGETS', level: 3, tissue: 'Lung',
            q: 'Why does osimertinib still work when erlotinib has stopped?',
            o: ['It binds EGFR even when the T790M resistance change is present',
                'It is given at a much higher dose',
                'It blocks MEK instead of EGFR',
                'It is an antibody rather than a small molecule'],
            a: 0,
            e: 'T790M is the usual escape from first generation EGFR inhibitors. Osimertinib was designed to bind the changed pocket, and it also spares normal EGFR better.'
        },
        {
            id: 'dr17', cat: 'DRUGS AND TARGETS', level: 2,
            q: 'Gemcitabine and 5-fluorouracil belong to which drug class?',
            o: ['Antimetabolites that starve DNA synthesis',
                'Topoisomerase inhibitors', 'Alkylating agents', 'Kinase inhibitors'],
            a: 0,
            e: 'They imitate the building blocks of DNA and RNA. The cell takes them up, and replication stalls on the faulty parts.'
        },

        // --------------------------------------------- EXPRESSION AND OTHER DATA
        {
            id: 'ex1', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'What does TPM measure in RNA sequencing?',
            o: ['Transcripts per million, corrected for gene length and read depth',
                'The total number of reads for a gene',
                'The share of cells expressing a gene',
                'How many copies of a gene are in the genome'],
            a: 0,
            e: 'Raw read counts favour long genes and deeply sequenced samples. TPM removes both, so two genes or two samples can be compared.'
        },
        {
            id: 'ex2', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'Why is expression usually shown as log2(TPM + 1)?',
            o: ['It compresses a huge range, and the plus one keeps zeros usable',
                'It converts RNA into protein units',
                'It removes the need for normalization',
                'It makes every gene sum to one'],
            a: 0,
            e: 'Expression spans several orders of magnitude. Taking logs makes that readable, and adding one avoids taking the log of zero.'
        },
        {
            id: 'ex3', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'On a log2 expression scale, a difference of 1 means what?',
            o: ['Twice as much RNA', 'Ten times as much RNA',
                'One extra transcript per cell', 'One percent more RNA'],
            a: 0,
            e: 'Each step of 1 on a log2 scale is a doubling. A gap of 3 is therefore an eightfold difference.'
        },
        {
            id: 'ex4', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'What does an expression z-score against all cell lines tell you?',
            o: ['How unusual this line is compared with the whole panel',
                'How much protein the gene makes',
                'How many reads the gene received',
                'Whether the gene is mutated'],
            a: 0,
            e: 'It says where a line sits in the distribution rather than giving an absolute level, which is what you want when comparing across genes.'
        },
        {
            id: 'ex5', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'In DepMap’s relative copy number scale, what does a value of 1.0 mean?',
            o: ['The normal two copies', 'One copy, so half of normal',
                'One extra copy above normal', 'The gene is deleted'],
            a: 0,
            e: 'The scale is relative to the sample’s own ploidy, so 1.0 is normal, below 1 means loss, and above 1 means gain.'
        },
        {
            id: 'ex6', cat: 'EXPRESSION AND OTHER DATA', level: 3,
            q: 'On that relative copy number scale, what does a value near 2.0 suggest?',
            o: ['A real amplification, roughly double the normal dose',
                'Exactly two copies, which is normal',
                'A homozygous deletion',
                'A sequencing error'],
            a: 0,
            e: 'Because 1.0 is already the normal two copies, 2.0 means about four, and high values like this are where amplified oncogenes such as MYC or ERBB2 show up.'
        },
        {
            id: 'ex7', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'What is a hotspot mutation?',
            o: ['A change seen again and again at the same position, usually activating',
                'Any mutation in a cancer gene',
                'A mutation that breaks the protein',
                'A mutation inherited from a parent'],
            a: 0,
            e: 'BRAF V600E and KRAS G12D are hotspots. The same position turning up in tumor after tumor is strong evidence that it does something.'
        },
        {
            id: 'ex8', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'What is a damaging mutation?',
            o: ['One that breaks the protein, such as a stop codon or frameshift',
                'One that appears at the same spot in many tumors',
                'One that raises the protein level',
                'One found only in blood cancers'],
            a: 0,
            e: 'Damaging changes are how tumor suppressors are lost. They can happen anywhere in the gene, unlike the tight hotspots of oncogenes.'
        },
        {
            id: 'ex9', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'What is STR profiling used for?',
            o: ['Confirming a cell line really is the line it claims to be',
                'Measuring how fast a cell line grows',
                'Counting mutations in a cell line',
                'Finding which drugs a line resists'],
            a: 0,
            e: 'Short tandem repeats give each line a fingerprint. Comparing it against the reference catches cross-contamination and mixed-up stocks.'
        },
        {
            id: 'ex10', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'What is Cellosaurus?',
            o: ['A catalog of cell lines and what is known about each one',
                'A database of protein structures',
                'A gene expression atlas',
                'A drug screening platform'],
            a: 0,
            e: 'It records origin, identity, publications and, importantly, warnings about lines known to be misidentified or contaminated.'
        },
        {
            id: 'ex11', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'What is Oncotree?',
            o: ['A standard tree of cancer type names used to classify samples',
                'A tool for drawing phylogenies of tumor clones',
                'A list of cancer driver genes',
                'A method for clustering gene expression'],
            a: 0,
            e: 'It lets everyone use the same words for the same disease, so a lung adenocarcinoma in one dataset can be matched to one in another.'
        },
        {
            id: 'ex12', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'Why do cell line panels always record the tissue of origin?',
            o: ['Lineage shapes both expression and dependencies more than almost anything else',
                'It is needed to order the cells from a supplier',
                'It determines the growth medium price',
                'Statistics packages require a label'],
            a: 0,
            e: 'A blood line and a lung line differ in thousands of genes. Ignoring lineage is the fastest way to mistake a tissue difference for a biological finding.'
        },
        {
            id: 'ex13', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'Which measurement tells you a gene has been amplified?',
            o: ['Copy number', 'Expression', 'Gene effect', 'Mutation burden'],
            a: 0,
            e: 'Copy number counts the DNA. An amplified gene is usually expressed more as well, but the two are measured separately and can disagree.'
        },
        {
            id: 'ex14', cat: 'EXPRESSION AND OTHER DATA', level: 1,
            q: 'What does a cell line’s mutation burden count?',
            o: ['How many mutations that line carries in total',
                'How many cancer genes are mutated',
                'How many chromosomes are abnormal',
                'How many drugs the line resists'],
            a: 0,
            e: 'A high burden usually points to a broken repair system, such as mismatch repair loss or a faulty POLE proofreading domain.'
        },
        {
            id: 'ex15', cat: 'EXPRESSION AND OTHER DATA', level: 2,
            q: 'A cell line is described as MSI-high. What does that imply about its mutation burden?',
            o: ['It will be high, because mismatch repair has failed',
                'It will be low, because repair is working',
                'It says nothing about mutation burden',
                'It will be exactly average'],
            a: 0,
            e: 'MSI comes from mismatch repair loss, and that same loss lets mutations accumulate across the genome, so the two travel together.'
        },

        // ----------------------------------------------------------- STATISTICS
        {
            id: 'st1', cat: 'STATISTICS', level: 1,
            q: 'What range can a Pearson correlation r take?',
            o: ['-1 to 1', '0 to 1', '0 to 100', 'Any number at all'],
            a: 0,
            e: 'The sign says which way the relation runs and the size says how tight it is. Plus or minus 1 means the points sit exactly on a line.'
        },
        {
            id: 'st2', cat: 'STATISTICS', level: 1,
            q: 'What does r = 0 mean?',
            o: ['There is no straight-line relation between the two',
                'The two are exactly equal',
                'One is always double the other',
                'The measurement failed'],
            a: 0,
            e: 'It rules out a linear trend, not any relation at all. A neat U shape can give r near zero while being highly structured.'
        },
        {
            id: 'st3', cat: 'STATISTICS', level: 1,
            q: 'What does r = -0.9 describe?',
            o: ['A strong relation where one goes up as the other goes down',
                'A weak relation with a slight downward tilt',
                'No relation at all',
                'A strong relation where both rise together'],
            a: 0,
            e: 'The minus sign is the direction and the 0.9 is the strength, so this is about as tight an inverse relation as biology usually shows.'
        },
        {
            id: 'st4', cat: 'STATISTICS', level: 2,
            q: 'What does a p-value tell you?',
            o: ['How surprising this result would be if there were no real relation',
                'The chance that the finding is true',
                'How large the effect is',
                'How many samples were measured'],
            a: 0,
            e: 'It is a statement about the data under an assumption of nothing going on. It is not the probability that your hypothesis is right.'
        },
        {
            id: 'st5', cat: 'STATISTICS', level: 3,
            q: 'A correlation has a tiny r but a very small p-value. What is going on?',
            o: ['The relation is real but weak, and n is large',
                'The relation is strong and certain',
                'The calculation must be wrong',
                'The two variables are identical'],
            a: 0,
            e: 'With a thousand cell lines even a slight tilt clears significance. Always read the effect size next to the p-value, never instead of it.'
        },
        {
            id: 'st6', cat: 'STATISTICS', level: 1,
            q: 'Why does the number of cell lines behind a correlation matter?',
            o: ['A correlation from few lines can appear by chance alone',
                'More lines always give a higher r',
                'Fewer lines make the result more precise',
                'The number of lines does not affect the result'],
            a: 0,
            e: 'With only ten points a striking pattern is easy to get by luck. The same r across hundreds of lines is a very different claim.'
        },
        {
            id: 'st7', cat: 'STATISTICS', level: 2,
            q: 'Two genes have strongly correlated gene effects. What is the safest reading?',
            o: ['They may work in the same process, or share a common cause',
                'One gene controls the other directly',
                'They sit next to each other on a chromosome',
                'Knocking out one will replace the other'],
            a: 0,
            e: 'Co-dependency is a strong hint of a shared complex or pathway, but the correlation itself does not say which way, or whether anything connects them at all.'
        },
        {
            id: 'st8', cat: 'STATISTICS', level: 2,
            q: 'On a heatmap, what does a row z-score show?',
            o: ['Where each cell line sits relative to that gene’s own average',
                'The raw value for that gene',
                'How the gene compares with other genes',
                'How many lines have data for that gene'],
            a: 0,
            e: 'Scaling each row separately lets genes with different ranges share one color scale. It also means colors cannot be compared between rows as absolute values.'
        },
        {
            id: 'st9', cat: 'STATISTICS', level: 2,
            q: 'What does the slope of a regression line tell you?',
            o: ['How much y changes for a one unit change in x',
                'How tightly the points hug the line',
                'Whether the relation is statistically significant',
                'How many points were used'],
            a: 0,
            e: 'Slope is about size, r is about tightness. Two plots can share a slope and look completely different in how scattered they are.'
        },
        {
            id: 'st10', cat: 'STATISTICS', level: 3,
            q: 'Which of r and the regression slope is free of units?',
            o: ['r, which is why it can be compared across any two measurements',
                'The slope, because it is a ratio',
                'Both are free of units',
                'Neither is free of units'],
            a: 0,
            e: 'The slope carries the units of y over x, so it changes if you rescale an axis. r does not, which makes it the fair comparison across gene pairs.'
        },
        {
            id: 'st11', cat: 'STATISTICS', level: 2,
            q: 'What can a single extreme outlier do to a Pearson correlation?',
            o: ['Create or destroy the correlation on its own',
                'Nothing, r is resistant to outliers',
                'Only lower r, never raise it',
                'Change the p-value but not r'],
            a: 0,
            e: 'One far-out point can drag the line to it. Looking at the scatter plot before trusting an r is the cheapest check there is.'
        },
        {
            id: 'st12', cat: 'STATISTICS', level: 1,
            q: 'What is the median of a set of values?',
            o: ['The middle value when they are put in order',
                'The average of all the values',
                'The most common value',
                'The difference between the largest and smallest'],
            a: 0,
            e: 'The median ignores how extreme the extremes are, so it describes a skewed group better than the mean does.'
        },
        {
            id: 'st13', cat: 'STATISTICS', level: 3,
            q: 'You test 18,000 genes at p below 0.05 with nothing real going on. Roughly how many will pass?',
            o: ['About 900', 'About 90', 'About 9', 'None, because nothing is real'],
            a: 0,
            e: 'Five percent of 18,000 is 900 false hits. This is why genome-wide work corrects for multiple testing, for example with a false discovery rate.'
        },
        {
            id: 'st14', cat: 'STATISTICS', level: 2,
            q: 'Two genes correlate at r = 0.85 across a thousand cell lines. Which statement is safest?',
            o: ['Their dependencies rise and fall together across the panel',
                'One gene must regulate the other',
                'Both genes are essential',
                'They will correlate just as strongly inside any single tissue'],
            a: 0,
            e: 'A panel-wide correlation can be driven by lineage. Checking whether it survives inside one tissue is the usual next step.'
        },
        {
            id: 'st15', cat: 'STATISTICS', level: 2,
            q: 'What does an effect size add that a p-value does not?',
            o: ['How big the difference is, not just whether it exists',
                'How many samples were used',
                'Whether the test was one-sided',
                'Whether the data are normally distributed'],
            a: 0,
            e: 'A p-value answers "could this be noise?". An effect size answers "would anyone care?", and only the second one guides an experiment.'
        },

        // -------------------------------------------------------- FAMOUS CELL LINES
        {
            id: 'cl1', cat: 'FAMOUS CELL LINES', level: 1, line: 'K562', tissue: 'Myeloid',
            q: 'Which cell line carries the BCR-ABL1 fusion?',
            o: ['K562', 'A375', 'MCF7', 'HCT116'],
            a: 0,
            e: 'K562 came from a patient with chronic myeloid leukemia in blast crisis, and it is the standard BCR-ABL1 model.'
        },
        {
            id: 'cl2', cat: 'FAMOUS CELL LINES', level: 2, line: 'K562', tissue: 'Myeloid',
            q: 'K562 was derived from a patient with which disease?',
            o: ['Chronic myeloid leukemia', 'Burkitt lymphoma',
                'Acute promyelocytic leukemia', 'Multiple myeloma'],
            a: 0,
            e: 'It is one of the oldest human leukemia lines and is used for everything from BCR-ABL1 work to natural killer cell assays.'
        },
        {
            id: 'cl3', cat: 'FAMOUS CELL LINES', level: 1, line: 'A375', tissue: 'Skin',
            q: 'A375 is a melanoma line. Which mutation does it carry?',
            o: ['BRAF V600E', 'KRAS G12D', 'EGFR T790M', 'TP53 R175H'],
            a: 0,
            e: 'A375 is the usual test bed for BRAF and MEK inhibitors, because its growth depends on that pathway being switched on.'
        },
        {
            id: 'cl4', cat: 'FAMOUS CELL LINES', level: 1, line: 'MCF7', tissue: 'Breast',
            q: 'MCF7 is a model of which cancer?',
            o: ['Estrogen receptor positive breast cancer', 'Ovarian cancer',
                'Lung adenocarcinoma', 'Bowel cancer'],
            a: 0,
            e: 'MCF7 keeps its estrogen receptor and grows in response to estrogen, which made it the workhorse for hormone therapy research.'
        },
        {
            id: 'cl5', cat: 'FAMOUS CELL LINES', level: 2, line: 'MDA-MB-231', tissue: 'Breast',
            q: 'MDA-MB-231 is a model of which kind of breast cancer?',
            o: ['Triple negative breast cancer', 'HER2 amplified breast cancer',
                'Estrogen receptor positive breast cancer', 'Inherited BRCA1 breast cancer'],
            a: 0,
            e: 'It lacks the estrogen and progesterone receptors and does not overexpress HER2, so none of the targeted breast therapies apply to it.'
        },
        {
            id: 'cl6', cat: 'FAMOUS CELL LINES', level: 2, line: 'HCT116', tissue: 'Bowel',
            q: 'HCT116 is a colorectal line best known for which feature?',
            o: ['Microsatellite instability from mismatch repair loss',
                'A BRCA1 deletion', 'HER2 amplification', 'Loss of RB1'],
            a: 0,
            e: 'HCT116 has defective mismatch repair, so it is the standard MSI model. It also keeps working p53, unlike most colorectal lines.'
        },
        {
            id: 'cl7', cat: 'FAMOUS CELL LINES', level: 3, line: 'HCT116', tissue: 'Bowel',
            q: 'Which KRAS change does HCT116 carry?',
            o: ['G13D', 'G12C', 'Q61L', 'A146T'],
            a: 0,
            e: 'HCT116 is heterozygous for KRAS G13D, and matched clones with the mutant allele removed are widely used to study RAS dependency.'
        },
        {
            id: 'cl8', cat: 'FAMOUS CELL LINES', level: 1, line: 'A549', tissue: 'Lung',
            q: 'A549 came from which tissue?',
            o: ['Lung', 'Liver', 'Pancreas', 'Bone'],
            a: 0,
            e: 'A549 is a lung adenocarcinoma line carrying KRAS G12S, and it also has lost the KEAP1 brake on the NRF2 stress program.'
        },
        {
            id: 'cl9', cat: 'FAMOUS CELL LINES', level: 1, line: 'HepG2', tissue: 'Liver',
            q: 'HepG2 came from which organ?',
            o: ['Liver', 'Kidney', 'Stomach', 'Prostate'],
            a: 0,
            e: 'HepG2 keeps many liver functions, including making albumin and clotting factors, so it is used for metabolism and toxicity work.'
        },
        {
            id: 'cl10', cat: 'FAMOUS CELL LINES', level: 1, line: 'PANC-1', tissue: 'Pancreas',
            q: 'PANC-1 came from which organ?',
            o: ['Pancreas', 'Lung', 'Breast', 'Bowel'],
            a: 0,
            e: 'Like almost every pancreatic ductal line, PANC-1 carries a mutant KRAS, in its case G12D.'
        },
        {
            id: 'cl11', cat: 'FAMOUS CELL LINES', level: 1, line: 'U2OS', tissue: 'Bone',
            q: 'U2OS is a model of which cancer?',
            o: ['Osteosarcoma', 'Neuroblastoma', 'Glioblastoma', 'Melanoma'],
            a: 0,
            e: 'U2OS is a bone tumor line that keeps working p53, which makes it a favorite for DNA damage and cell imaging experiments.'
        },
        {
            id: 'cl12', cat: 'FAMOUS CELL LINES', level: 2, line: 'Jurkat', tissue: 'Lymphoid',
            q: 'Jurkat is a model of which cancer?',
            o: ['T cell acute lymphoblastic leukemia', 'B cell lymphoma',
                'Acute myeloid leukemia', 'Multiple myeloma'],
            a: 0,
            e: 'Jurkat cells signal like a T cell, so much of what is known about the T cell receptor pathway was worked out in them.'
        },
        {
            id: 'cl13', cat: 'FAMOUS CELL LINES', level: 2, line: 'Raji', tissue: 'Lymphoid',
            q: 'Raji is a model of which cancer?',
            o: ['Burkitt lymphoma', 'Hodgkin lymphoma', 'Myeloma', 'T cell leukemia'],
            a: 0,
            e: 'Raji carries the MYC translocation that defines Burkitt lymphoma, and it also carries Epstein-Barr virus.'
        },
        {
            id: 'cl14', cat: 'FAMOUS CELL LINES', level: 2, line: 'THP-1', tissue: 'Myeloid',
            q: 'THP-1 is a model of which cancer?',
            o: ['Acute myeloid leukemia of monocytic type', 'Chronic myeloid leukemia',
                'Burkitt lymphoma', 'Small cell lung cancer'],
            a: 0,
            e: 'THP-1 can be pushed to behave like a macrophage, which is why it turns up so often in immunology and inflammation papers.'
        },
        {
            id: 'cl15', cat: 'FAMOUS CELL LINES', level: 3, line: 'MOLM-13', tissue: 'Myeloid',
            q: 'MOLM-13 is an acute myeloid leukemia line driven by which alteration?',
            o: ['An FLT3 internal tandem duplication', 'BCR-ABL1',
                'A BRAF V600E mutation', 'PML-RARA'],
            a: 0,
            e: 'The FLT3-ITD keeps the receptor switched on without a ligand. MOLM-13 is therefore very sensitive to FLT3 inhibitors such as gilteritinib.'
        },
        {
            id: 'cl16', cat: 'FAMOUS CELL LINES', level: 2, line: '786-O', tissue: 'Kidney',
            q: 'The kidney cancer line 786-O has lost which gene?',
            o: ['VHL', 'PTEN', 'RB1', 'SMAD4'],
            a: 0,
            e: 'With VHL gone, HIF is never destroyed, so the cells behave as if short of oxygen. 786-O is the standard model for that biology.'
        },
        {
            id: 'cl17', cat: 'FAMOUS CELL LINES', level: 2, line: 'LNCaP', tissue: 'Prostate',
            q: 'LNCaP is a prostate cancer line that depends on which protein?',
            o: ['The androgen receptor', 'The estrogen receptor', 'HER2', 'ALK'],
            a: 0,
            e: 'LNCaP grows in response to androgen and is used to study hormone therapy and how resistance to it develops.'
        },
        {
            id: 'cl18', cat: 'FAMOUS CELL LINES', level: 2, line: 'SH-SY5Y', tissue: 'Peripheral Nervous System',
            q: 'SH-SY5Y and SK-N-SH are models of which cancer?',
            o: ['Neuroblastoma', 'Glioblastoma', 'Medulloblastoma', 'Retinoblastoma'],
            a: 0,
            e: 'SH-SY5Y is a subclone of SK-N-SH. Both can be pushed towards a neuron-like state, so they are used well beyond cancer research.'
        },
        {
            id: 'cl19', cat: 'FAMOUS CELL LINES', level: 3, line: 'NCI-H1975', tissue: 'Lung',
            q: 'The lung line NCI-H1975 carries which EGFR change?',
            o: ['L858R together with the T790M resistance change', 'An exon 19 deletion',
                'An EGFR amplification only', 'No EGFR change at all'],
            a: 0,
            e: 'Because it carries T790M, NCI-H1975 resists erlotinib and gefitinib but still responds to osimertinib. It is the standard resistance model.'
        },
        {
            id: 'cl20', cat: 'FAMOUS CELL LINES', level: 3, line: 'HCC827', tissue: 'Lung',
            q: 'The lung line HCC827 carries which EGFR change?',
            o: ['A deletion in exon 19', 'The T790M resistance change',
                'A KRAS G12C mutation', 'An ALK fusion'],
            a: 0,
            e: 'The exon 19 deletion leaves HCC827 addicted to EGFR, so it is exquisitely sensitive to EGFR inhibitors.'
        },
        {
            id: 'cl21', cat: 'FAMOUS CELL LINES', level: 2, line: 'SK-BR-3', tissue: 'Breast',
            q: 'SK-BR-3 is a breast cancer line known for what?',
            o: ['A strong HER2 amplification', 'Loss of BRCA1',
                'Being triple negative', 'Carrying an ALK fusion'],
            a: 0,
            e: 'SK-BR-3 overexpresses HER2 heavily, which makes it the usual test line for trastuzumab and other HER2 drugs.'
        },
        {
            id: 'cl22', cat: 'FAMOUS CELL LINES', level: 1, line: 'A375', tissue: 'Skin',
            q: 'Which of these is a melanoma cell line?',
            o: ['A375', 'K562', 'A549', 'HCT116'],
            a: 0,
            e: 'A375 is skin, K562 is leukemia, A549 is lung and HCT116 is bowel. Knowing a line’s tissue is the first step in reading any result from it.'
        },
        {
            id: 'cl23', cat: 'FAMOUS CELL LINES', level: 3, line: 'NB4', tissue: 'Myeloid',
            q: 'The leukemia line NB4 carries which fusion?',
            o: ['PML-RARA', 'BCR-ABL1', 'RUNX1-RUNX1T1', 'ETV6-RUNX1'],
            a: 0,
            e: 'NB4 is the standard acute promyelocytic leukemia model, and it matures rather than dies when given retinoic acid.'
        },
        {
            id: 'cl24', cat: 'FAMOUS CELL LINES', level: 3, line: 'HT-29', tissue: 'Bowel',
            q: 'The colorectal line HT-29 carries which driver mutation?',
            o: ['BRAF V600E', 'KRAS G12D', 'EGFR L858R', 'ALK fusion'],
            a: 0,
            e: 'BRAF V600E bowel cancers behave differently from BRAF melanomas, and HT-29 is the line used to show why BRAF inhibitors alone do less there.'
        },
        {
            id: 'cl25', cat: 'FAMOUS CELL LINES', level: 2, line: 'U-251 MG', tissue: 'CNS/Brain',
            q: 'U-251 MG is a model of which cancer?',
            o: ['Glioblastoma', 'Medulloblastoma', 'Neuroblastoma', 'Meningioma'],
            a: 0,
            e: 'U-251 MG and LN-229 are the usual glioblastoma lines. The older U-87 MG stock is known to be misidentified, so many labs avoid it.'
        },
        {
            id: 'cl26', cat: 'FAMOUS CELL LINES', level: 3, line: 'A-673', tissue: 'Bone',
            q: 'A-673 is a model of which tumor?',
            o: ['Ewing sarcoma', 'Osteosarcoma', 'Rhabdomyosarcoma', 'Chondrosarcoma'],
            a: 0,
            e: 'A-673 carries EWSR1-FLI1 and cannot grow without it, which is why it appears in almost every study of fusion dependency.'
        },

        // ------------------------------------------------------ HOW THE APP WORKS
        {
            id: 'ap1', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What does a gene set analysis do in this app?',
            o: ['Correlates the gene effect profiles of your genes across all cell lines',
                'Looks up which of your genes are mutated',
                'Ranks your genes by how often they are published',
                'Aligns the protein sequences of your genes'],
            a: 0,
            e: 'Every gene has one score per cell line. Comparing those profiles shows which genes rise and fall together, which usually means a shared pathway or complex.'
        },
        {
            id: 'ap2', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What does the correlation cutoff control?',
            o: ['How strong a correlation has to be before a link is drawn',
                'How many cell lines are included',
                'How many genes may be entered',
                'How thick the lines in the network are'],
            a: 0,
            e: 'Raise it and only the tightest relations survive. Lower it and more appears, including more that is coincidence.'
        },
        {
            id: 'ap3', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What is the Cell Line Browser for?',
            o: ['Finding and filtering the cell lines that fit a project',
                'Ordering cell lines from a supplier',
                'Uploading your own screening data',
                'Drawing pathway diagrams'],
            a: 0,
            e: 'You can filter by tissue, subtype, mutation, fusion and more, then open any line’s page to see what is known about it.'
        },
        {
            id: 'ap4', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What does a cell line’s wiki page show?',
            o: ['Everything the app knows about that line, on one page',
                'The published papers that used the line',
                'A protocol for growing the line',
                'The price and supplier of the line'],
            a: 0,
            e: 'Origin, drivers, copy number, fusions, drug responses and dependencies are gathered in one place, so you do not have to hunt through separate views.'
        },
        {
            id: 'ap5', cat: 'HOW THE APP WORKS', level: 2,
            q: 'In a scatter plot here, what does coloring by a hotspot mutation show?',
            o: ['Which points carry a mutation in the chosen gene',
                'Which points have the highest expression',
                'Which points come from the same tissue',
                'Which points were measured most reliably'],
            a: 0,
            e: 'It answers whether a relation is driven by the mutated lines. If the mutants sit apart from the rest, the mutation is part of the story.'
        },
        {
            id: 'ap6', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What does the gene set heatmap show?',
            o: ['Many genes across many cell lines at once, as a grid of colors',
                'One gene across the tissues',
                'The correlation between two genes',
                'Drug responses for one cell line'],
            a: 0,
            e: 'It is the wide view: patterns that a single scatter cannot show, such as a block of lines sharing a whole set of dependencies.'
        },
        {
            id: 'ap7', cat: 'HOW THE APP WORKS', level: 2,
            q: 'What is Export for AI?',
            o: ['A file holding the current view’s data and settings for a language model to read',
                'A way to send your results to the authors',
                'An automatic figure caption generator',
                'A machine learning model that predicts dependencies'],
            a: 0,
            e: 'It packages what is on screen, with the numbers behind it, so an assistant can answer questions about that exact view instead of guessing.'
        },
        {
            id: 'ap8', cat: 'HOW THE APP WORKS', level: 2,
            q: 'What does the mutation analysis rank?',
            o: ['Genes whose dependency differs most between mutated and wild-type lines',
                'Genes that are mutated most often',
                'Cell lines with the most mutations',
                'Drugs that work best on mutated lines'],
            a: 0,
            e: 'Split the panel by a mutation, then ask which dependencies change. That is how a genotype is tied to a weakness worth drugging.'
        },
        {
            id: 'ap9', cat: 'HOW THE APP WORKS', level: 2,
            q: 'In the network view, what does an edge between two genes mean?',
            o: ['Their gene effect profiles correlate above the cutoff',
                'They bind each other physically',
                'They sit close together on a chromosome',
                'They are mutated in the same cell lines'],
            a: 0,
            e: 'The edge is a statement about the screen, not about protein contacts. Genes in one complex usually do get an edge, which is why the map is useful.'
        },
        {
            id: 'ap10', cat: 'HOW THE APP WORKS', level: 1,
            q: 'What does the gene effect view by tissue show?',
            o: ['One gene’s dependency across every tissue in the panel',
                'Every gene’s dependency in one tissue',
                'How often a gene is mutated per tissue',
                'How much RNA a tissue makes in total'],
            a: 0,
            e: 'It is the fastest way to see whether a dependency is general or belongs to one lineage, which is usually the first question worth asking.'
        },
        {
            id: 'ap11', cat: 'HOW THE APP WORKS', level: 1,
            q: 'Where does the cell line data in this app come from?',
            o: ['DepMap, the Broad Institute’s cancer dependency project',
                'The Human Protein Atlas', 'The 1000 Genomes Project', 'UniProt'],
            a: 0,
            e: 'Gene effect, expression, copy number, mutations and drug response all come from DepMap releases, with cell line identity notes from Cellosaurus.'
        },
        {
            id: 'ap12', cat: 'HOW THE APP WORKS', level: 2,
            q: 'Why does the app let you restrict an analysis to one tissue before running it?',
            o: ['So a correlation is not just a difference between lineages',
                'To make the calculation finish faster',
                'Because gene effects are only valid within a tissue',
                'Because other tissues have no data'],
            a: 0,
            e: 'Panel-wide correlations often reflect blood versus solid rather than biology. Rerunning inside one tissue is the standard sanity check.'
        }
    ];

    // -------------------------------------------------------- correlation pool
    function pearson(x, y) {
        const n = x.length;
        let c = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (let i = 0; i < n; i++) {
            const a = x[i], b = y[i];
            if (!isFinite(a) || !isFinite(b) || a <= -990 || b <= -990) continue;
            c++; sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
        }
        if (c < 200) return null;
        const num0 = sxy - sx * sy / c;
        const den = Math.sqrt((sxx - sx * sx / c) * (syy - sy * sy / c));
        if (!(den > 0)) return null;
        return { r: num0 / den, n: c };
    }

    // -------------------------------------------------------------- figures
    function fitCanvas(canvas, w, h) {
        const dpr = Math.min(window.devicePixelRatio || 1, 2);
        canvas.width = Math.round(w * dpr);
        canvas.height = Math.round(h * dpr);
        canvas.style.width = w + 'px';
        canvas.style.height = h + 'px';
        const ctx = canvas.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.imageSmoothingEnabled = false;
        return ctx;
    }
    const pixFont = (px) => px + "px 'Press Start 2P', monospace";
    const clip = (s, n) => (s.length > n ? s.slice(0, n - 1) + '.' : s);

    // Horizontal bars: cell line and gene names are too wide to sit under
    // vertical bars at phone width.
    function barFigure(items, opts) {
        const o = opts || {};
        return function (canvas, w) {
            const rowH = 26, gap = 8, top = o.title ? 26 : 8;
            const h = top + items.length * (rowH + gap) + 8;
            const ctx = fitCanvas(canvas, w, h);
            ctx.fillStyle = PAL.panel; ctx.fillRect(0, 0, w, h);
            if (o.title) {
                ctx.fillStyle = PAL.dim; ctx.font = pixFont(8); ctx.textBaseline = 'middle';
                ctx.fillText(clip(o.title, Math.floor(w / 8)), 8, 14);
            }
            const labW = Math.min(Math.max(96, w * 0.34), 150);
            const fmt = o.fmt || (v => String(v));
            let valW = 0;
            items.forEach(it => { valW = Math.max(valW, fmt(it.value).length * 8); });
            // Values get their own right-hand column: drawing them at the end of
            // a bar puts them on top of the next label when the bar is long.
            const x0 = labW + 6, x1 = Math.max(x0 + 20, w - 12 - valW);
            let lo = 0, hi = 0;
            items.forEach(it => { lo = Math.min(lo, it.value); hi = Math.max(hi, it.value); });
            if (hi === lo) hi = lo + 1;
            const span = hi - lo;
            const zx = x0 + (0 - lo) / span * (x1 - x0);
            ctx.fillStyle = PAL.dim; ctx.fillRect(Math.round(zx), top, 3, h - top - 6);
            items.forEach((it, i) => {
                const y = top + i * (rowH + gap);
                ctx.font = pixFont(8); ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
                ctx.fillStyle = it.hi ? PAL.yellow : PAL.ink;
                ctx.fillText(clip(String(it.label), Math.floor(labW / 8)), 6, y + rowH / 2);
                const vx = x0 + (it.value - lo) / span * (x1 - x0);
                const bx = Math.min(zx, vx), bw = Math.max(Math.abs(vx - zx), 3);
                ctx.fillStyle = it.hi ? PAL.green : PAL.blue;
                ctx.fillRect(Math.round(bx), y + 4, Math.round(bw), rowH - 8);
                ctx.fillStyle = PAL.ink; ctx.font = pixFont(8); ctx.textAlign = 'right';
                ctx.fillText(fmt(it.value), w - 6, y + rowH / 2);
                ctx.textAlign = 'left';
            });
        };
    }

    // --------------------------------------------------------------- state
    const S = {
        opened: false, root: null, screen: 'title',
        cohort: null, nQ: 10,
        bank: [], qi: 0, score: 0, shown: 0, streak: 0, best: 0,
        lives: START_LIVES, correct: 0, answered: false, chosen: -1, gained: 0,
        deadline: 0, left: TIMER_MS, raf: 0, running: false,
        sound: false, actx: null, figure: null, figureEl: null, plotDiv: null, lastEntry: null, wired: false
    };

    // --------------------------------------------------------------- sound
    function beep(kind) {
        if (!S.sound) return;
        try {
            if (!S.actx) S.actx = new (window.AudioContext || window.webkitAudioContext)();
            const ctx = S.actx;
            if (ctx.state === 'suspended') ctx.resume();
            const spec = {
                correct: [[660, 0], [880, 0.08], [1320, 0.16]],
                wrong: [[220, 0], [160, 0.12]],
                tick: [[1200, 0]],
                over: [[440, 0], [330, 0.14], [220, 0.28], [110, 0.42]],
                blip: [[520, 0]],
                fanfare: [[523, 0], [659, 0.08], [784, 0.16], [1047, 0.26]]
            }[kind] || [[440, 0]];
            spec.forEach(([f, t]) => {
                const osc = ctx.createOscillator(), g = ctx.createGain();
                osc.type = 'square'; osc.frequency.value = f;
                g.gain.value = 0.05;
                osc.connect(g); g.connect(ctx.destination);
                const t0 = ctx.currentTime + t;
                g.gain.setValueAtTime(0.05, t0);
                g.gain.exponentialRampToValueAtTime(0.0001, t0 + 0.1);
                osc.start(t0); osc.stop(t0 + 0.12);
            });
        } catch (e) { /* audio is a nicety, never a failure */ }
    }

    // ---------------------------------------------------------- high scores
    function loadScores() {
        try {
            const raw = localStorage.getItem(SCORE_KEY);
            const arr = raw ? JSON.parse(raw) : [];
            return Array.isArray(arr) ? arr.filter(e => e && typeof e.score === 'number') : [];
        } catch (e) { return []; }
    }
    function saveScores(list) {
        try { localStorage.setItem(SCORE_KEY, JSON.stringify(list.slice(0, MAX_SCORES))); } catch (e) { }
    }
    function qualifies(score) {
        if (score <= 0) return false;
        const list = loadScores();
        if (list.length < MAX_SCORES) return true;
        return score > list[list.length - 1].score;
    }
    function addScore(entry) {
        const list = loadScores();
        list.push(entry);
        list.sort((a, b) => b.score - a.score);
        saveScores(list);
    }

    // ================================================== question generation
    // Each builder returns null when it cannot find a fair question; the bank
    // builder simply moves on to another template.
    // Questions name a cell line only from the recognisable list. A narrow
    // cohort can hold too few of those, and then the whole cohort is used.
    function ctxCellLines(cohort) {
        return (cohort.knownIdx && cohort.knownIdx.length >= 6) ? cohort.knownIdx : cohort.idx;
    }

    function tryTimes(n, fn) {
        for (let t = 0; t < n; t++) {
            const r = fn();
            if (r) return r;
        }
        return null;
    }

    function options4(correctLabel, wrongLabels) {
        const opts = [{ label: correctLabel, correct: true }]
            .concat(wrongLabels.slice(0, 3).map(l => ({ label: l, correct: false })));
        if (opts.length < 4) return null;
        return shuffle(opts);
    }

    function distinctNames(labels) {
        const seen = new Set();
        for (const l of labels) {
            const k = String(l).toUpperCase();
            if (seen.has(k)) return false;
            seen.add(k);
        }
        return true;
    }

    const TEMPLATES = [];
    const T = (id, tag, needs, build) => TEMPLATES.push({ id, tag, needs, build });

    // 1. WHO AM I ----------------------------------------------------------
    T('whoami', 'WHO AM I', () => true, (c) => tryTimes(200, () => {
        const pool = ctxCellLines(c);
        const i = pick(pool);
        const id = D.ids()[i];
        if (c.usedCL.has(id)) return null;
        const name = D.name(id);
        let text = D.summary(id);
        if (!text || text.length < 160) return null;
        const variants = [A().cellLineMetadata?.cellLineName?.[id], name, id]
            .filter(Boolean).map(String)
            .concat([String(name).replace(/[^A-Za-z0-9]/g, '')])
            .filter(v => v.length >= 2)
            .sort((a, b) => b.length - a.length);
        variants.forEach(v => {
            const rx = new RegExp(v.replace(/[.*+?^${}()|[\]\\-]/g, '\\$&'), 'gi');
            text = text.replace(rx, '???');
        });
        // The summary can open with a data warning banner. The identity
        // sentence is the one that starts with the masked name, so the quote
        // starts there and anything before it is dropped.
        const start = text.indexOf('??? is ');
        if (start < 0) return null;
        text = text.slice(start);
        if (text.length < 160) return null;
        if (text.length > 360) {
            const cut = text.lastIndexOf('. ', 360);
            text = text.slice(0, cut > 180 ? cut + 1 : 360);
        }
        const lin = D.lineage(id);
        const same = pool.map(k => D.ids()[k])
            .filter(x => x !== id && D.lineage(x) === lin && D.name(x) !== name);
        const other = pool.map(k => D.ids()[k]).filter(x => x !== id);
        const src = same.length >= 3 ? same : other;
        const wrong = shuffle(src).slice(0, 3).map(D.name);
        if (wrong.length < 3 || !distinctNames([name].concat(wrong))) return null;
        const opts = options4(name, wrong);
        if (!opts) return null;
        return {
            tag: 'WHO AM I', text: 'Which cell line is being described?', quote: text, options: opts,
            explain: `It is ${name}, ${lin ? 'a ' + lin.toLowerCase() + ' line' : 'one of the lines in the panel'}. ${D.subtype(id) ? 'Its type is ' + D.subtype(id) + '.' : ''}`.trim(),
            usedCL: [id],
            openInApp: { label: 'Open the ' + name + ' wiki', run: (a) => a.openCellLineWiki(id) }
        };
    }));

    // 2. WHERE FROM --------------------------------------------------------
    T('where', 'WHERE FROM', (c) => c.key === 'all', (c) => tryTimes(200, () => {
        const i = pick(ctxCellLines(c));
        const id = D.ids()[i];
        if (c.usedCL.has(id)) return null;
        const lin = D.lineage(id);
        if (!lin) return null;
        const others = c.lineages.map(l => l.name).filter(n => n !== lin);
        if (others.length < 3) return null;
        const opts = options4(lin, shuffle(others));
        if (!opts) return null;
        const nSame = c.lineageCount.get(lin) || 0;
        return {
            tag: 'WHERE FROM', text: `Which tissue does ${D.name(id)} come from?`, options: opts,
            explain: `${D.name(id)} is ${lin.toLowerCase()}. The panel holds ${num(nSame)} cell lines from that tissue.`,
            usedCL: [id],
            openInApp: { label: 'Open the ' + D.name(id) + ' wiki', run: (a) => a.openCellLineWiki(id) }
        };
    }));

    // 3. TOP DEPENDENCY ----------------------------------------------------
    T('dep', 'GENE EFFECT', () => true, (c) => tryTimes(60, () => {
        const i = pick(ctxCellLines(c));
        const id = D.ids()[i];
        if (c.usedCL.has(id)) return null;
        const a = A(), nCL = a.nCellLines;
        // The answer has to be specific to this line, not a gene every line needs.
        let hit = null, looked = 0;
        for (const gene of shuffle(notableGenes())) {
            if (looked >= 30 || hit) break;
            if (c.usedGene.has(gene)) continue;
            const row = D.row(gene);
            const v0 = D.ge(row, i);
            if (!isFinite(v0) || v0 > -1.1) continue;
            looked++;
            let seen = 0, mild = 0;
            for (let k = 0; k < nCL; k++) {
                const v = D.ge(row, k);
                if (!isFinite(v)) continue;
                seen++; if (v > -0.5) mild++;
            }
            if (seen >= 300 && mild / seen >= 0.7) hit = { gene, value: v0 };
        }
        if (!hit) return null;
        const wrong = [];
        for (const g of shuffle(notableGenes())) {
            if (wrong.length >= 3) break;
            if (g === hit.gene || c.usedGene.has(g)) continue;
            const v = D.geAt(g, i);
            if (isFinite(v) && v >= -0.2) wrong.push({ gene: g, value: v });
        }
        if (wrong.length < 3) return null;
        const opts = options4(hit.gene, wrong.map(w => w.gene));
        if (!opts) return null;
        const bars = shuffle([hit].concat(wrong)).map(x => ({
            label: x.gene, value: Math.round(x.value * 100) / 100, hi: x.gene === hit.gene
        }));
        return {
            tag: 'GENE EFFECT', text: `Which of these genes does ${D.name(id)} need most to grow?`,
            options: opts,
            explain: `${D.name(id)} depends on ${hit.gene} (gene effect ${hit.value.toFixed(2)}). A score near 0 means the cell line does not need the gene, and anything below -1 is a strong need.`,
            figure: barFigure(bars, { title: 'GENE EFFECT IN ' + D.name(id), fmt: v => v.toFixed(2) }),
            figureBefore: false,
            usedCL: [id], usedGene: [hit.gene].concat(wrong.map(w => w.gene)),
            openInApp: { label: 'See ' + hit.gene + ' across tissues', run: (a2) => a2.openGeneEffectModal(hit.gene, 'tissue') }
        };
    }));

    // 4. HIGHER EXPRESSION -------------------------------------------------
    T('expr', 'EXPRESSION', () => A().expressionLoaded, (c) => tryTimes(120, () => {
        const gene = pick(notableGenes());
        if (c.usedGene.has(gene)) return null;
        const pool = ctxCellLines(c);
        let hiI = -1, hiV = -Infinity, loI = -1, loV = Infinity;
        const sample = pool.length > 400 ? shuffle(pool).slice(0, 400) : pool;
        for (const i of sample) {
            const v = D.expr(gene, i);
            if (!isFinite(v)) continue;
            if (v > hiV) { hiV = v; hiI = i; }
            if (v < loV) { loV = v; loI = i; }
        }
        if (hiI < 0 || loI < 0 || hiI === loI) return null;
        if (hiV - loV < 3) return null;
        const hiId = D.ids()[hiI], loId = D.ids()[loI];
        if (c.usedCL.has(hiId) || c.usedCL.has(loId)) return null;
        const hiName = D.name(hiId), loName = D.name(loId);
        if (!distinctNames([hiName, loName])) return null;
        const opts = options4(hiName, [loName, 'About the same', 'Neither expresses it']);
        if (!opts) return null;
        return {
            tag: 'EXPRESSION', text: `Which cell line makes more ${gene} RNA, ${hiName} or ${loName}?`,
            options: opts,
            explain: `${hiName} sits at ${hiV.toFixed(1)} and ${loName} at ${loV.toFixed(1)} on a log2 scale, so ${hiName} carries roughly ${Math.round(Math.pow(2, hiV - loV))} times more ${gene} RNA.`,
            figure: barFigure([
                { label: hiName, value: Math.round(hiV * 10) / 10, hi: true },
                { label: loName, value: Math.round(loV * 10) / 10 }
            ], { title: gene + ' RNA, log2', fmt: v => v.toFixed(1) }),
            figureBefore: false,
            usedCL: [hiId, loId], usedGene: [gene]
        };
    }));

    // 7. HOW MANY ----------------------------------------------------------
    T('howmany', 'THE PANEL', (c) => c.lineages.length >= 1, (c) => tryTimes(60, () => {
        const l = pick(c.lineages);
        if (c.usedLineage.has(l.name)) return null;
        const truth = l.ids.length;
        const wrong = new Set();
        let guard = 0;
        while (wrong.size < 3 && guard++ < 60) {
            const f = 0.25 + Math.random() * 2.2;
            const v = Math.max(4, Math.round(truth * f));
            if (Math.abs(v - truth) / truth < 0.25) continue;
            if ([...wrong].some(w => Math.abs(w - v) / Math.max(w, v) < 0.25)) continue;
            wrong.add(v);
        }
        if (wrong.size < 3) return null;
        const opts = options4(num(truth), [...wrong].map(num));
        if (!opts) return null;
        return {
            tag: 'THE PANEL', text: `How many ${l.name.toLowerCase()} cell lines are in this panel?`,
            options: opts,
            explain: `The panel holds ${num(truth)} ${l.name.toLowerCase()} cell lines out of ${num(D.n())} in total.`,
            usedLineage: [l.name]
        };
    }));

    // 8. HALLMARK GENE -----------------------------------------------------
    T('hallmark', 'HALLMARK', () => Object.keys(D.hallmarks()).length >= 4, (c) => tryTimes(120, () => {
        const kb = D.hallmarks();
        const names = Object.keys(kb).filter(k => (kb[k].lookFor || []).length);
        const inCohort = names.filter(k => c.subtypes.has(k));
        const src = inCohort.length >= 1 ? inCohort : names;
        const sub = pick(src);
        const entry = kb[sub];
        const known = new Set(notableGenes());
        const good = (entry.lookFor || []).filter(g => !c.usedGene.has(g) && known.has(g));
        if (!good.length) return null;
        const gene = pick(good);
        const mine = new Set(entry.lookFor || []);
        const expected = String(entry.expected || '');
        const wrong = [];
        for (const other of shuffle(names)) {
            if (wrong.length >= 3 || other === sub) continue;
            for (const g of shuffle(kb[other].lookFor || [])) {
                if (wrong.length >= 3) break;
                if (mine.has(g) || wrong.indexOf(g) >= 0 || g === gene) continue;
                if (expected.indexOf(g) >= 0) continue;
                if (c.usedGene.has(g) || !known.has(g)) continue;
                wrong.push(g);
            }
        }
        if (wrong.length < 3) return null;
        const opts = options4(gene, wrong);
        if (!opts) return null;
        const first = expected.split('. ')[0];
        return {
            tag: 'HALLMARK', text: `Which gene is a hallmark of ${sub}?`,
            options: opts,
            explain: `${gene} is one of the genes to look for in ${sub}. ${first ? first + '.' : ''}`.trim().replace(/\.\.$/, '.'),
            usedGene: [gene].concat(wrong),
            openInApp: { label: 'See ' + gene + ' across tissues', run: (a2) => a2.openGeneEffectModal(gene, 'tissue') }
        };
    }));

    // 9. PATHWAY -----------------------------------------------------------
    T('pathway', 'PATHWAY', () => Object.keys(D.pathways()).length >= 4, (c) => tryTimes(120, () => {
        const pw = D.pathways();
        const names = Object.keys(pw).filter(k => (pw[k].genes || []).length);
        if (names.length < 4) return null;
        const right = pick(names);
        const mine = pw[right].genes || [];
        const known = new Set(notableGenes());
        const gene = pick(mine.filter(g => !c.usedGene.has(g) && known.has(g)) || []);
        if (!gene) return null;
        const wrong = shuffle(names).filter(n => n !== right && (pw[n].genes || []).indexOf(gene) < 0).slice(0, 3);
        if (wrong.length < 3) return null;
        const opts = options4(right, wrong);
        if (!opts) return null;
        const note = String(pw[right].note || '').split('. ')[0];
        return {
            tag: 'PATHWAY', text: `${gene} belongs to which pathway?`,
            options: opts,
            explain: `${gene} sits in ${right}. ${note ? note + '.' : ''}`.trim().replace(/\.\.$/, '.'),
            usedGene: [gene]
        };
    }));

    // 10. DRUG -------------------------------------------------------------
    T('drug', 'DRUG SCREEN', () => knownCompounds().length >= 6, (c) => tryTimes(80, () => {
        const i = pick(ctxCellLines(c));
        const id = D.ids()[i];
        if (c.usedCL.has(id)) return null;
        const comps = knownCompounds();
        const killers = [], sparers = [];
        for (const cp of comps) {
            const v = cp.auc?.[id];
            if (typeof v !== 'number') continue;
            if (v <= 0.4) killers.push({ cp, v });
            else if (v >= 0.8) sparers.push({ cp, v });
        }
        if (!killers.length || sparers.length < 3) return null;
        const hit = pick(killers);
        const wrong = shuffle(sparers).slice(0, 3);
        const names = [drugName(hit.cp.name)].concat(wrong.map(w => drugName(w.cp.name)));
        if (!distinctNames(names)) return null;
        const opts = options4(drugName(hit.cp.name), wrong.map(w => drugName(w.cp.name)));
        if (!opts) return null;
        const bars = shuffle([hit].concat(wrong)).map(x => ({
            label: drugName(x.cp.name), value: Math.round(x.v * 100) / 100, hi: x.cp.name === hit.cp.name
        }));
        const moa = hit.cp.moa || hit.cp.target || '';
        return {
            tag: 'DRUG SCREEN', text: `Which compound kills ${D.name(id)} best in the drug screen?`,
            options: opts,
            explain: `${drugName(hit.cp.name)} leaves ${D.name(id)} at ${hit.v.toFixed(2)} on a scale where 1 means the cells are untouched.${moa ? ' It is ' + (/^[aeiou]/i.test(moa) ? 'an ' : 'a ') + moa + (hit.cp.target ? ', aimed at ' + hit.cp.target : '') + '.' : ''}`,
            figure: barFigure(bars, { title: 'AUC IN ' + D.name(id) + ', LOWER KILLS MORE', fmt: v => v.toFixed(2) }),
            figureBefore: false,
            usedCL: [id],
            openInApp: { label: 'Open the ' + D.name(id) + ' wiki', run: (a2) => a2.openCellLineWiki(id) }
        };
    }));

    // 11. FUSION -----------------------------------------------------------
    T('fusion', 'FUSION', (c) => c.fusionLines.length >= 1 && c.noFusionLines.length >= 3, (c) => tryTimes(80, () => {
        const hitId = pick(c.fusionLines);
        if (c.usedCL.has(hitId)) return null;
        const famous = new Set(FAMOUS_FUSIONS);
        const calls = D.fusions(hitId).filter(f => famous.has(f.fusion));
        if (!calls.length) return null;
        const call = pick(calls);
        const others = shuffle(c.noFusionLines.filter(x => !c.usedCL.has(x))).slice(0, 3);
        if (others.length < 3) return null;
        const hitName = D.name(hitId);
        const wrongNames = others.map(D.name);
        if (!distinctNames([hitName].concat(wrongNames))) return null;
        const opts = options4(hitName, wrongNames);
        if (!opts) return null;
        const ctx = A().clinicalFusions?.fusionData?.[call.fusion]?.diseaseContext || '';
        return {
            tag: 'FUSION', text: `Which cell line carries the ${call.fusion} fusion?`,
            options: opts,
            explain: `${hitName} carries ${call.fusion}. ${ctx ? 'That fusion is the mark of ' + ctx + '.' : 'It is a ' + (D.lineage(hitId) || 'cancer').toLowerCase() + ' line.'}`,
            usedCL: [hitId].concat(others),
            openInApp: { label: 'Open the ' + hitName + ' wiki', run: (a2) => a2.openCellLineWiki(hitId) }
        };
    }));

    // 12. MUTATION BURDEN --------------------------------------------------
    T('burden', 'MUTATION LOAD', (c) => c.burden.length >= 12, (c) => tryTimes(80, () => {
        const top = c.burden[ri(Math.min(30, c.burden.length))];
        if (!top || c.usedCL.has(top.id) || top.n < 60) return null;
        const low = c.burden.filter(b => b.n > 0 && b.n <= top.n / 3 && !c.usedCL.has(b.id));
        if (low.length < 3) return null;
        const wrong = shuffle(low).slice(0, 3);
        const names = [D.name(top.id)].concat(wrong.map(w => D.name(w.id)));
        if (!distinctNames(names)) return null;
        const opts = options4(D.name(top.id), wrong.map(w => D.name(w.id)));
        if (!opts) return null;
        const bars = shuffle([top].concat(wrong)).map(b => ({
            label: D.name(b.id), value: b.n, hi: b.id === top.id
        }));
        return {
            tag: 'MUTATION LOAD', text: 'Which of these cell lines carries the most damaging mutations?',
            options: opts,
            explain: `${D.name(top.id)} carries ${num(top.n)} damaging mutations, well above the others here. A high load usually means broken DNA repair.`,
            figure: barFigure(bars, { title: 'DAMAGING MUTATIONS', fmt: v => num(v) }),
            figureBefore: false,
            usedCL: [top.id].concat(wrong.map(w => w.id))
        };
    }));


    // ============================================== app-style Plotly figures
    // Most questions are answered by reading a chart, and the chart is the one
    // the app itself would draw: same colors, same axis wording, same shapes.
    // Rendered with Plotly into the question card, before the answer.
    const hasPlotly = () => typeof window.Plotly !== 'undefined';
    const PLOT_CFG = { displayModeBar: false, responsive: true, displaylogo: false };
    const plotH = () => (window.innerWidth <= 640 ? 230 : 265);
    function plotLayout(extra) {
        const small = window.innerWidth <= 640;
        return Object.assign({
            height: plotH(),
            margin: { t: 42, r: 14, b: 46, l: small ? 46 : 56 },
            paper_bgcolor: '#ffffff',
            plot_bgcolor: '#fafafa',
            showlegend: false,
            font: { family: 'Arial, Helvetica, sans-serif', size: small ? 11 : 12, color: '#374151' },
            hovermode: false
        }, extra || {});
    }
    // The app writes its chart titles as bold text above the plot area.
    const plotTitle = (text, sub) => ({
        text: `<b>${text}</b>` + (sub ? `<br><span style="font-size:10px;color:#6b7280;">${sub}</span>` : ''),
        font: { size: window.innerWidth <= 640 ? 12 : 14 }, x: 0.5, xanchor: 'center'
    });
    const median = (arr) => {
        if (!arr.length) return NaN;
        const s = arr.slice().sort((a, b) => a - b);
        const m = s.length >> 1;
        return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
    };
    const mean = (arr) => arr.length ? arr.reduce((a, b) => a + b, 0) / arr.length : NaN;

    // ---------------------------------------------------------- pair pool
    // Every well-known gene against every other, bucketed by how strong the
    // relation is. Computed once per session, in about a tenth of a second.
    let PAIR_POOL = null;
    function pairPool() {
        if (PAIR_POOL) return PAIR_POOL;
        PAIR_POOL = { strong: [], weak: [], none: [], neg: [] };
        const genes = knownGenes().filter(g => D.row(g));
        const rows = genes.map(g => D.row(g));
        for (let i = 0; i < genes.length; i++) {
            for (let j = i + 1; j < genes.length; j++) {
                const st = pearson(rows[i], rows[j]);
                if (!st) continue;
                const rec = { g1: genes[i], g2: genes[j], r: st.r, n: st.n };
                if (st.r >= 0.5) PAIR_POOL.strong.push(rec);
                else if (st.r >= 0.25 && st.r < 0.5) PAIR_POOL.weak.push(rec);
                else if (Math.abs(st.r) < 0.1) PAIR_POOL.none.push(rec);
                else if (st.r <= -0.3) PAIR_POOL.neg.push(rec);
            }
        }
        // A pair we can explain properly comes first in every bucket.
        const score = (p) => (pairWhy(p.g1, p.g2) ? 1 : 0);
        Object.keys(PAIR_POOL).forEach(k => {
            PAIR_POOL[k] = shuffle(PAIR_POOL[k]).sort((a, b) => score(b) - score(a)).slice(0, 60);
        });
        return PAIR_POOL;
    }

    function scatterXY(g1, g2, cohortIdx) {
        const r1 = D.row(g1), r2 = D.row(g2);
        const xs = [], ys = [], ids = [];
        const idsAll = D.ids();
        const pool = cohortIdx && cohortIdx.length >= 150 ? cohortIdx : idsAll.map((_, i) => i);
        for (const i of pool) {
            const x = D.ge(r1, i), y = D.ge(r2, i);
            if (isFinite(x) && isFinite(y)) { xs.push(x); ys.push(y); ids.push(idsAll[i]); }
        }
        return { xs, ys, ids };
    }

    // Least squares line, drawn in the app's green.
    function regressionTrace(xs, ys) {
        const n = xs.length;
        let sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
        const den = n * sxx - sx * sx;
        if (!(Math.abs(den) > 1e-9)) return null;
        const slope = (n * sxy - sx * sy) / den;
        const intercept = (sy - slope * sx) / n;
        const lo = Math.min.apply(null, xs), hi = Math.max.apply(null, xs);
        return {
            x: [lo, hi], y: [slope * lo + intercept, slope * hi + intercept],
            mode: 'lines', type: 'scatter', line: { color: '#6ba544', width: 3 },
            hoverinfo: 'skip', showlegend: false
        };
    }

    const zeroLineAxis = (title) => ({
        title: { text: title, font: { size: window.innerWidth <= 640 ? 11 : 12 }, standoff: 6 },
        zeroline: true, zerolinecolor: '#000', zerolinewidth: 2,
        tickfont: { size: window.innerWidth <= 640 ? 10 : 11 }
    });

    const FIG = [];
    const FT = (id, tag, needs, build) => FIG.push({ id, tag, needs, build, isFigure: true });

    // F1. Two gene effects on a scatter, as the app's Correlation view draws it.
    FT('figCorr', 'CORRELATION', () => hasPlotly() && pairPool().strong.length >= 1,
        (c) => tryTimes(30, () => {
            const pool = pairPool();
            const buckets = [];
            if (pool.strong.length) buckets.push(['Strong positive correlation', pool.strong]);
            if (pool.weak.length) buckets.push(['Weak positive correlation', pool.weak]);
            if (pool.none.length) buckets.push(['No correlation', pool.none]);
            if (pool.neg.length) buckets.push(['Negative correlation', pool.neg]);
            if (!buckets.length) return null;
            const [label, list] = pick(buckets);
            const idx = ri(list.length);
            const p = list[idx];
            if (c.usedGene.has(p.g1) || c.usedGene.has(p.g2)) return null;
            list.splice(idx, 1);
            const { xs, ys } = scatterXY(p.g1, p.g2, null);
            if (xs.length < 200) return null;
            const why = pairWhy(p.g1, p.g2);
            const askR = Math.random() < 0.4;
            const BANDS = ['-1.0 to -0.5', '-0.5 to -0.2', '-0.2 to 0.2', '0.2 to 0.5', '0.5 to 1.0'];
            const bandOf = (r) => r < -0.5 ? BANDS[0] : r < -0.2 ? BANDS[1] : r < 0.2 ? BANDS[2] : r < 0.5 ? BANDS[3] : BANDS[4];
            const answer = askR ? bandOf(p.r) : label;
            const wrongs = askR
                ? shuffle(BANDS.filter(b => b !== answer)).slice(0, 3)
                : shuffle(['Strong positive correlation', 'Weak positive correlation',
                    'No correlation', 'Negative correlation'].filter(o => o !== answer)).slice(0, 3);
            const opts = options4(answer, wrongs);
            if (!opts) return null;
            return {
                tag: 'CORRELATION',
                text: askR
                    ? `Each dot is one cell line. Roughly what is the correlation r between these two gene effects?`
                    : `Each dot is one cell line. What does this scatter show?`,
                options: opts,
                explain: `r = ${p.r.toFixed(2)} across ${num(p.n)} cell lines. ${why || 'Genes whose dependencies track each other this way usually sit in the same complex or pathway.'}`,
                plot: (div) => {
                    const traces = [{
                        x: xs, y: ys, mode: 'markers', type: 'scatter',
                        marker: { color: '#9ca3af', size: 6, opacity: 0.6 }, hoverinfo: 'skip'
                    }];
                    const rl = regressionTrace(xs, ys);
                    if (rl) traces.push(rl);
                    return Plotly.newPlot(div, traces, plotLayout({
                        title: plotTitle(`${p.g1} vs ${p.g2}`, `n = ${num(p.n)} cell lines`),
                        xaxis: zeroLineAxis(`${p.g1} Gene Effect`),
                        yaxis: zeroLineAxis(`${p.g2} Gene Effect`)
                    }), PLOT_CFG);
                },
                usedGene: [p.g1, p.g2],
                openInApp: {
                    label: `Open ${p.g1} vs ${p.g2} in the app`,
                    run: (a) => a.openInspectByGenes(p.g1, p.g2)
                }
            };
        }));

    // F2. The same scatter, colored by a hotspot mutation the way the app's
    // Hotspot Overlay does: grey wild-type, blue one copy, red two.
    const HOTSPOT_TRIPLES = [
        { x: 'BRAF', y: 'MAPK1', m: 'BRAF' },
        { x: 'KRAS', y: 'RAF1', m: 'KRAS' },
        { x: 'TP53', y: 'MDM2', m: 'TP53' },
        { x: 'BRAF', y: 'SOX10', m: 'BRAF' },
        { x: 'NRAS', y: 'MAPK1', m: 'NRAS' }
    ];
    FT('figHotspot', 'HOTSPOT OVERLAY', () => hasPlotly() && !!A().mutations?.geneData,
        (c) => tryTimes(20, () => {
            const t = pick(HOTSPOT_TRIPLES);
            if (c.usedGene.has(t.y) || !A().mutations.geneData[t.m]) return null;
            if (!D.row(t.x) || !D.row(t.y)) return null;
            const ids = D.ids();
            const g = { wt: { x: [], y: [] }, m1: { x: [], y: [] }, m2: { x: [], y: [] } };
            const yWT = [], yMut = [];
            const rx = D.row(t.x), ry = D.row(t.y);
            for (let i = 0; i < ids.length; i++) {
                const xv = D.ge(rx, i), yv = D.ge(ry, i);
                if (!isFinite(xv) || !isFinite(yv)) continue;
                const lvl = D.hotspot(t.m, ids[i]);
                const b = lvl >= 2 ? g.m2 : lvl === 1 ? g.m1 : g.wt;
                b.x.push(xv); b.y.push(yv);
                (lvl >= 1 ? yMut : yWT).push(yv);
            }
            if (yMut.length < 12 || yWT.length < 40) return null;
            const mMut = mean(yMut), mWT = mean(yWT);
            const diff = mWT - mMut;
            const answer = Math.abs(diff) < 0.25 ? 'About the same'
                : (diff > 0 ? `Lines with a ${t.m} mutation` : `Lines with wild-type ${t.m}`);
            const opts = options4(answer, shuffle([`Lines with a ${t.m} mutation`,
                `Lines with wild-type ${t.m}`, 'About the same',
                'The plot cannot answer that'].filter(o => o !== answer)).slice(0, 3));
            if (!opts) return null;
            const why = pairWhy(t.x, t.y);
            return {
                tag: 'HOTSPOT OVERLAY',
                text: `Grey dots are wild-type for ${t.m}, colored dots carry a hotspot mutation. Which group depends more on ${t.y}?`,
                options: opts,
                explain: `Mean ${t.y} gene effect is ${mMut.toFixed(2)} in the ${num(yMut.length)} mutated lines and ${mWT.toFixed(2)} in the ${num(yWT.length)} wild-type lines. ${why || 'A dependency that follows the mutation is the pattern a targeted drug is built on.'}`,
                plot: (div) => {
                    const traces = [
                        { x: g.wt.x, y: g.wt.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#9ca3af', size: 6, opacity: 0.6 }, name: 'WT' },
                        { x: g.m1.x, y: g.m1.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#3b82f6', size: 7, opacity: 0.8 }, name: '1 mut' },
                        { x: g.m2.x, y: g.m2.y, mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#dc2626', size: 7, opacity: 0.85 }, name: '2 mut' }
                    ];
                    return Plotly.newPlot(div, traces, plotLayout({
                        title: plotTitle(`${t.x} vs ${t.y}`, `colored by ${t.m} hotspot mutation`),
                        xaxis: zeroLineAxis(`${t.x} Gene Effect`),
                        yaxis: zeroLineAxis(`${t.y} Gene Effect`)
                    }), PLOT_CFG);
                },
                usedGene: [t.x, t.y],
                openInApp: {
                    label: `Open ${t.x} vs ${t.y} in the app`,
                    run: (a) => a.openInspectByGenes(t.x, t.y)
                }
            };
        }));

    // Per-tissue values for one gene, either dependency or expression.
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
        return [...by.entries()]
            .filter(e => e[1].length >= (minN || 12))
            .map(e => ({ name: e[0], vals: e[1], med: median(e[1]) }));
    }

    // Rows the app draws as horizontal boxes with every point shown.
    function boxTraces(groups) {
        return groups.map(gp => ({
            type: 'box', name: `${gp.name} (n=${gp.vals.length})`,
            x: gp.vals, boxpoints: 'all', jitter: 0.35, pointpos: 0,
            marker: { color: 'rgba(80,80,80,0.5)', size: 4 },
            line: { color: '#374151' }, fillcolor: 'rgba(200,200,200,0.3)',
            hoverinfo: 'skip'
        }));
    }

    // F3. Gene effect by tissue, the app's Gene Effect view.
    FT('figTissueDep', 'GENE EFFECT BY TISSUE', () => hasPlotly(), (c) => tryTimes(40, () => {
        const gene = pick(LINEAGE_GENES.filter(g => A().geneIndex.has(g) && !c.usedGene.has(g)));
        if (!gene) return null;
        const groups = tissueGroups(gene, false, 14);
        if (groups.length < 5) return null;
        const sorted = groups.slice().sort((a, b) => a.med - b.med);
        const best = sorted[0];
        if (!(best.med < -0.35)) return null;
        // The answer has to be clearly the lowest, or the plot is unfair.
        if (!(sorted[1].med - best.med > 0.25)) return null;
        const others = shuffle(sorted.slice(1)).slice(0, 4);
        if (others.length < 3) return null;
        const show = shuffle([best].concat(others));
        const opts = options4(best.name, others.slice(0, 3).map(o => o.name));
        if (!opts) return null;
        return {
            tag: 'GENE EFFECT BY TISSUE',
            text: `Which tissue depends most on ${gene}?`,
            options: opts,
            explain: `${best.name} lines sit at a median gene effect of ${best.med.toFixed(2)} for ${gene}, well below the other tissues here. A dependency that follows the tissue usually means the gene holds that lineage's identity.`,
            plot: (div) => Plotly.newPlot(div,
                boxTraces(show.slice().sort((a, b) => b.med - a.med)),
                plotLayout({
                    title: plotTitle(`${gene} Gene Effect by tissue`),
                    xaxis: zeroLineAxis(`${gene} Gene Effect`),
                    yaxis: { automargin: true, tickfont: { size: window.innerWidth <= 640 ? 9 : 10 } },
                    margin: { t: 40, r: 14, b: 44, l: 4 }
                }), PLOT_CFG),
            usedGene: [gene],
            openInApp: {
                label: `See ${gene} across tissues in the app`,
                run: (a) => a.openGeneEffectModal(gene, 'tissue')
            }
        };
    }));

    // F4. The same plot for expression.
    FT('figTissueExpr', 'EXPRESSION BY TISSUE', () => hasPlotly() && A().expressionLoaded,
        (c) => tryTimes(40, () => {
            const gene = pick(LINEAGE_GENES.filter(g => A().geneIndex.has(g) && !c.usedGene.has(g)));
            if (!gene) return null;
            const groups = tissueGroups(gene, true, 14);
            if (groups.length < 5) return null;
            const sorted = groups.slice().sort((a, b) => b.med - a.med);
            const best = sorted[0];
            if (!(sorted[0].med - sorted[1].med > 1)) return null;
            const others = shuffle(sorted.slice(1)).slice(0, 4);
            if (others.length < 3) return null;
            const show = shuffle([best].concat(others));
            const opts = options4(best.name, others.slice(0, 3).map(o => o.name));
            if (!opts) return null;
            return {
                tag: 'EXPRESSION BY TISSUE',
                text: `Which tissue expresses ${gene} the most?`,
                options: opts,
                explain: `${best.name} lines sit at a median of ${best.med.toFixed(1)} on the log2 scale, about ${Math.round(Math.pow(2, best.med - sorted[1].med))} times the next tissue. Expression this tissue-specific usually marks a lineage factor.`,
                plot: (div) => Plotly.newPlot(div,
                    boxTraces(show.slice().sort((a, b) => a.med - b.med)),
                    plotLayout({
                        title: plotTitle(`${gene} expression by tissue`),
                        xaxis: {
                            title: { text: `${gene} mRNA, log2(TPM+1)`, font: { size: window.innerWidth <= 640 ? 11 : 12 }, standoff: 6 },
                            tickfont: { size: window.innerWidth <= 640 ? 10 : 11 }, zeroline: false
                        },
                        yaxis: { automargin: true, tickfont: { size: window.innerWidth <= 640 ? 9 : 10 } },
                        margin: { t: 40, r: 14, b: 44, l: 4 }
                    }), PLOT_CFG),
                usedGene: [gene],
                openInApp: {
                    label: `See ${gene} across tissues in the app`,
                    run: (a) => a.openGeneEffectModal(gene, 'tissue', { dataType: 'expr' })
                }
            };
        }));

    // F5. Mutation status strip plot, as the app's mutation inspect draws it:
    // one row per number of mutated copies, every cell line a dot.
    const STRIP_PAIRS = [
        { g: 'BRAF', m: 'BRAF' }, { g: 'KRAS', m: 'KRAS' }, { g: 'NRAS', m: 'NRAS' },
        { g: 'EGFR', m: 'EGFR' }, { g: 'PIK3CA', m: 'PIK3CA' }, { g: 'TP53', m: 'TP53' },
        { g: 'MDM2', m: 'TP53' }, { g: 'CTNNB1', m: 'CTNNB1' }
    ];
    FT('figStrip', 'MUTATION STATUS', () => hasPlotly() && !!A().mutations?.geneData,
        (c) => tryTimes(24, () => {
            const p = pick(STRIP_PAIRS);
            if (c.usedGene.has(p.g) || !A().mutations.geneData[p.m] || !D.row(p.g)) return null;
            const ids = D.ids(), row = D.row(p.g);
            const wt = [], m1 = [], m2 = [];
            for (let i = 0; i < ids.length; i++) {
                const v = D.ge(row, i);
                if (!isFinite(v)) continue;
                const lvl = D.hotspot(p.m, ids[i]);
                (lvl >= 2 ? m2 : lvl === 1 ? m1 : wt).push(v);
            }
            const mut = m1.concat(m2);
            if (mut.length < 12 || wt.length < 40) return null;
            const medWT = median(wt), medMut = median(mut);
            const d = medWT - medMut;
            const answer = Math.abs(d) < 0.15 ? `There is no clear difference`
                : d > 0 ? `Mutated lines depend more on ${p.g}` : `Wild-type lines depend more on ${p.g}`;
            const opts = options4(answer, shuffle([
                `Mutated lines depend more on ${p.g}`,
                `Wild-type lines depend more on ${p.g}`,
                'There is no clear difference',
                `Mutated lines grow faster`].filter(o => o !== answer)).slice(0, 3));
            if (!opts) return null;
            const jit = (base, n) => Array.from({ length: n }, () => base + (Math.random() - 0.5) * 0.5);
            return {
                tag: 'MUTATION STATUS',
                text: `Each dot is one cell line, split by ${p.m} mutation status. What does this plot show?`,
                options: opts,
                explain: `Median ${p.g} gene effect is ${medMut.toFixed(2)} in the ${num(mut.length)} ${p.m} mutated lines and ${medWT.toFixed(2)} in the ${num(wt.length)} wild-type lines.` +
                    (p.g === 'MDM2' ? ' Lines with working p53 need MDM2 to hold it down, so losing MDM2 hurts them and not the mutants.'
                        : p.g === p.m ? ' A mutant oncogene keeps the pathway on, and the cells then cannot do without it.' : ''),
                plot: (div) => {
                    const traces = [
                        { x: wt, y: jit(0, wt.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#888888', size: 5, opacity: 0.7 } },
                        { x: m1, y: jit(1, m1.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#3b82f6', size: 6, opacity: 0.8 } },
                        { x: m2, y: jit(2, m2.length), mode: 'markers', type: 'scatter', hoverinfo: 'skip', marker: { color: '#dc2626', size: 6, opacity: 0.85 } }
                    ];
                    return Plotly.newPlot(div, traces, plotLayout({
                        title: plotTitle(`${p.g} Gene Effect by ${p.m} mutation status`),
                        xaxis: zeroLineAxis(`${p.g} Gene Effect`),
                        yaxis: {
                            tickmode: 'array', tickvals: [0, 1, 2],
                            ticktext: [`WT (n=${wt.length})`, `1 mut (n=${m1.length})`, `2 mut (n=${m2.length})`],
                            range: [-0.6, 2.6], tickfont: { size: window.innerWidth <= 640 ? 9 : 10 }, automargin: true
                        },
                        margin: { t: 40, r: 14, b: 44, l: 4 }
                    }), PLOT_CFG);
                },
                usedGene: [p.g],
                openInApp: {
                    label: `See ${p.g} across tissues in the app`,
                    run: (a) => a.openGeneEffectModal(p.g, 'tissue')
                }
            };
        }));

    // F6. Wiki-style histogram: the cohort's distribution with this cell
    // line's own value as a red marker, exactly as the wiki panels look.
    const SIG_METRICS = [
        { key: 'Ploidy', label: 'Ploidy', dp: 2, why: 'Ploidy is the average number of chromosome copies. Around 2 is a normal set, and higher means the genome has been doubled.' },
        { key: 'Aneuploidy', label: 'Aneuploidy score', dp: 0, why: 'The aneuploidy score counts how many chromosome arms are gained or lost. Most cancer lines carry several.' },
        { key: 'LoHFraction', label: 'Loss of heterozygosity', dp: 2, why: 'Loss of heterozygosity means one parental copy is gone across that share of the genome, which is how many tumor suppressors are finished off.' }
    ];
    FT('figHist', 'WHERE DOES IT SIT', () => hasPlotly() && !!A().globalSignatures?.byCellLine,
        (c) => tryTimes(40, () => {
            const known = c.knownIdx && c.knownIdx.length ? c.knownIdx : c.idx;
            const i = pick(known);
            const ids = D.ids();
            const id = ids[i];
            if (c.usedCL.has(id)) return null;
            const name = D.name(id);
            const useDrug = Math.random() < 0.35 && knownCompounds().length >= 4;
            let vals = [], mine = null, label = '', why = '', dp = 2, title = '';
            if (useDrug) {
                const cands = knownCompounds().filter(cp => typeof cp.auc?.[id] === 'number');
                if (!cands.length) return null;
                const cp = pick(cands);
                const nice = drugName(cp.name);
                ids.forEach(x => { const v = cp.auc?.[x]; if (typeof v === 'number') vals.push(v); });
                mine = cp.auc[id];
                label = `${nice} AUC`;
                title = `${nice} response across the panel`;
                why = 'A low AUC means the drug killed the cells, and a high one means they carried on growing.';
                dp = 2;
            } else {
                const m = pick(SIG_METRICS);
                const sig = A().globalSignatures.byCellLine;
                ids.forEach(x => { const v = sig[x]?.[m.key]; if (typeof v === 'number') vals.push(v); });
                mine = typeof sig[id]?.[m.key] === 'number' ? sig[id][m.key] : null;
                label = m.label; why = m.why; dp = m.dp;
                title = `${m.label} across the panel`;
            }
            if (vals.length < 200 || mine == null) return null;
            const sorted = vals.slice().sort((a, b) => a - b);
            const pct = sorted.filter(v => v < mine).length / sorted.length * 100;
            // Only ask when the value sits well inside a band, never on a line.
            let answer = null;
            if (pct <= 20) answer = 'Among the lowest quarter';
            else if (pct >= 30 && pct <= 70) answer = 'Around the middle';
            else if (pct >= 80) answer = 'Among the highest quarter';
            if (!answer) return null;
            const opts = options4(answer, shuffle(['Among the lowest quarter', 'Around the middle',
                'Among the highest quarter', 'It was not measured'].filter(o => o !== answer)).slice(0, 3));
            if (!opts) return null;
            const lo = Math.min.apply(null, sorted), hi = Math.max.apply(null, sorted);
            return {
                tag: 'WHERE DOES IT SIT',
                text: `The bars are every cell line in the panel and the red line is ${name}. Where does ${name} sit for ${label}?`,
                options: opts,
                explain: `${name} sits at ${Number(mine).toFixed(dp)}, higher than ${Math.round(pct)} percent of the panel. ${why}`,
                plot: (div) => Plotly.newPlot(div, [{
                    type: 'histogram', x: vals,
                    marker: { color: '#9ca3af', line: { color: '#ffffff', width: 1 } },
                    xbins: { start: lo, end: hi, size: Math.max((hi - lo) / 30, 1e-6) },
                    hoverinfo: 'skip'
                }], plotLayout({
                    title: plotTitle(title, name + ' marked in red'),
                    xaxis: {
                        title: { text: label, font: { size: window.innerWidth <= 640 ? 11 : 12 }, standoff: 6 },
                        tickfont: { size: window.innerWidth <= 640 ? 10 : 11 },
                        showgrid: false, zeroline: false, showline: true, linecolor: '#d1d5db'
                    },
                    yaxis: { showgrid: false, zeroline: false, showticklabels: false, showline: false },
                    bargap: 0.15,
                    shapes: [{ type: 'line', x0: mine, x1: mine, y0: 0, y1: 1, yref: 'paper', line: { color: '#dc2626', width: 2 } }],
                    margin: { t: 46, r: 14, b: 46, l: 20 }
                }), PLOT_CFG),
                usedCL: [id],
                openInApp: { label: `Open the ${name} wiki`, run: (a) => a.openCellLineWiki(id) }
            };
        }));

    // F7. Hotspot frequency per tissue, as a bar chart. Only asked for genes
    // where the panel agrees with the textbook answer.
    FT('figMutFreq', 'MUTATION FREQUENCY',
        (c) => hasPlotly() && c.lineages.length >= 5 && !!A().mutations?.geneData,
        (c) => tryTimes(30, () => {
            const gene = pick(Object.keys(TEXTBOOK_HOTSPOT_TISSUE)
                .filter(g => A().mutations.geneData[g] && !c.usedGene.has(g)));
            if (!gene) return null;
            const want = TEXTBOOK_HOTSPOT_TISSUE[gene];
            const rows = c.lineages.map(l => {
                let mut = 0;
                for (const id of l.ids) if (D.hotspot(gene, id) >= 1) mut++;
                return { name: l.name, pct: mut / l.ids.length * 100, n: l.ids.length, mut };
            }).sort((a, b) => b.pct - a.pct);
            const top = rows[0];
            // The data has to agree with the textbook, or the question is
            // teaching a quirk of this panel rather than a fact.
            if (!top || top.name !== want || top.pct < 12) return null;
            const low = rows.slice(1).filter(r => r.pct <= top.pct / 2);
            if (low.length < 4) return null;
            const wrong = shuffle(low).slice(0, 4);
            const opts = options4(top.name, wrong.slice(0, 3).map(w => w.name));
            if (!opts) return null;
            const show = shuffle([top].concat(wrong)).sort((a, b) => a.pct - b.pct);
            return {
                tag: 'MUTATION FREQUENCY',
                text: `Which tissue has the highest share of ${gene} hotspot mutations?`,
                options: opts,
                explain: `${Math.round(top.pct)} of every 100 ${top.name.toLowerCase()} lines carry a hotspot mutation in ${gene} (${top.mut} of ${top.n}), which matches what is seen in patients.`,
                plot: (div) => Plotly.newPlot(div, [{
                    type: 'bar', orientation: 'h',
                    x: show.map(r => r.pct), y: show.map(r => r.name),
                    marker: { color: show.map(() => 'rgba(122, 185, 80, 0.85)') },
                    hoverinfo: 'skip'
                }], plotLayout({
                    title: plotTitle(`${gene} hotspot mutations`, 'percent of cell lines'),
                    xaxis: {
                        title: { text: 'Percent of lines mutated', font: { size: window.innerWidth <= 640 ? 11 : 12 }, standoff: 6 },
                        tickfont: { size: window.innerWidth <= 640 ? 10 : 11 }, zeroline: true, zerolinecolor: '#000'
                    },
                    yaxis: { automargin: true, tickfont: { size: window.innerWidth <= 640 ? 9 : 10 } },
                    margin: { t: 46, r: 14, b: 46, l: 4 }
                }), PLOT_CFG),
                usedGene: [gene],
                openInApp: {
                    label: `See ${gene} across tissues in the app`,
                    run: (a) => a.openGeneEffectModal(gene, 'tissue')
                }
            };
        }));

    // ------------------------------------------------------- bank assembly
    function buildContext(cohort) {
        const ids = D.ids();
        const counts = D.lineageCounts();
        const lineages = [...counts.entries()].filter(e => e[1] >= 10)
            .sort((a, b) => b[1] - a[1])
            .map(e => ({ name: e[0], count: e[1], ids: [] }));
        const byName = new Map(lineages.map(l => [l.name, l]));
        cohort.idx.forEach(i => {
            const l = byName.get(D.lineage(ids[i]));
            if (l) l.ids.push(ids[i]);
        });
        // Lineage-level questions need the whole panel behind them, so those
        // lists are filled from every line, not only the chosen tissue.
        const allByName = new Map(lineages.map(l => [l.name, { name: l.name, ids: [] }]));
        ids.forEach(id => { const l = allByName.get(D.lineage(id)); if (l) l.ids.push(id); });
        const lineageList = lineages.map(l => allByName.get(l.name)).filter(l => l.ids.length >= 10);

        const subtypes = new Set();
        cohort.idx.forEach(i => { const s = D.subtype(ids[i]); if (s) subtypes.add(s); });

        // The recognisable lines inside this cohort. Every question that names
        // a cell line draws from here first.
        const known = famousIds();
        const knownIdx = cohort.idx.filter(i => known.has(ids[i]));

        // Fusion questions only ask about the fusions people have heard of.
        const famousSet = new Set(FAMOUS_FUSIONS);
        const fusionLines = [], noFusionLines = [];
        knownIdx.forEach(i => {
            const id = ids[i];
            const calls = D.fusions(id).filter(f => famousSet.has(f.fusion));
            if (calls.length) fusionLines.push(id);
            else if (!D.fusions(id).length) noFusionLines.push(id);
        });

        const burden = [];
        knownIdx.forEach(i => {
            const n = D.damaging(ids[i]);
            if (typeof n === 'number' && n > 0) burden.push({ id: ids[i], n });
        });
        burden.sort((a, b) => b.n - a.n);

        return {
            key: cohort.key, idx: cohort.idx, knownIdx, lineages: lineageList, lineageCount: counts,
            subtypes, fusionLines, noFusionLines, burden,
            usedCL: new Set(), usedGene: new Set(), usedLineage: new Set()
        };
    }

    function commit(ctx, q) {
        (q.usedCL || []).forEach(x => ctx.usedCL.add(x));
        (q.usedGene || []).forEach(x => ctx.usedGene.add(x));
        (q.usedLineage || []).forEach(x => ctx.usedLineage.add(x));
    }

    // Roughly how hard each template is, so a game can start easy.
    const TEMPLATE_LEVEL = {
        figCorr: 1, figTissueDep: 1, figTissueExpr: 1, figMutFreq: 1, figHotspot: 2,
        figStrip: 2, figHist: 2,
        where: 1, howmany: 1, whoami: 2, dep: 2, expr: 2, hallmark: 2, pathway: 2,
        drug: 2, fusion: 2, burden: 3
    };

    // Templates are tried in turn, least used first. A template that cannot
    // find a fair question is set aside and tried again only if the rest run out.
    function drawFrom(list, ctx, n, cap) {
        const used = {}, dead = {};
        const out = [];
        const live = () => list.filter(t => !dead[t.id] && (used[t.id] || 0) < cap);
        let guard = 0;
        while (out.length < n && guard++ < 300) {
            let pool = live();
            if (!pool.length) { Object.keys(dead).forEach(k => delete dead[k]); pool = live(); }
            if (!pool.length) break;
            const lo = Math.min.apply(null, pool.map(t => used[t.id] || 0));
            const t = pick(pool.filter(x => (used[x.id] || 0) === lo));
            let q = null;
            try {
                if (t.needs(ctx)) q = t.build(ctx);
            } catch (e) { q = null; }
            if (!q) { dead[t.id] = true; continue; }
            q.template = t.id;
            q.level = q.level || TEMPLATE_LEVEL[t.id] || 2;
            used[t.id] = (used[t.id] || 0) + 1;
            commit(ctx, q);
            out.push(q);
        }
        return out;
    }

    // ------------------------------------------------- concept bank rotation
    const ASKED_KEY = 'correlateQuizAsked';
    const ASKED_MAX = 60;
    function loadAsked() {
        try {
            const raw = localStorage.getItem(ASKED_KEY);
            const arr = raw ? JSON.parse(raw) : [];
            return Array.isArray(arr) ? arr.filter(x => typeof x === 'string') : [];
        } catch (e) { return []; }
    }
    function rememberAsked(ids) {
        try {
            const merged = ids.concat(loadAsked()).slice(0, ASKED_MAX);
            localStorage.setItem(ASKED_KEY, JSON.stringify(merged));
        } catch (e) { }
    }

    function conceptToQuestion(entry) {
        const opts = shuffle(entry.o.map((label, i) => ({ label, correct: i === entry.a })));
        const id = entry.line ? lineId(entry.line) : null;
        return {
            tag: entry.cat, text: entry.q, options: opts, explain: entry.e,
            level: entry.level, conceptId: entry.id, template: 'concept',
            openInApp: id ? {
                label: 'Open the ' + entry.line + ' wiki',
                run: (a) => a.openCellLineWiki(id)
            } : null
        };
    }

    // Pick concept questions: never twice in one game, and a repeat player
    // gets the ones the last few games did not use. A cohort tissue pulls its
    // own questions forward without shutting the others out.
    function drawConcepts(ctx, n) {
        if (n <= 0) return [];
        const asked = new Set(loadAsked());
        const tissue = ctx.key !== 'all' ? ctx.key : null;
        const fresh = CONCEPTS.filter(e => !asked.has(e.id));
        const base = fresh.length >= n ? fresh : CONCEPTS;
        const onTopic = tissue ? shuffle(base.filter(e => e.tissue === tissue)) : [];
        const rest = shuffle(base.filter(e => onTopic.indexOf(e) < 0));
        // At most a third of the concept questions are tissue-specific, so a
        // Lung game still teaches the general material.
        const chosen = onTopic.slice(0, Math.max(1, Math.round(n / 3))).concat(rest).slice(0, n);
        rememberAsked(chosen.map(e => e.id));
        return chosen.map(conceptToQuestion);
    }

    // About half the game is a chart to read, a third is the concept bank and
    // the rest are the other data questions.
    function buildBank(ctx, n) {
        const wantFig = Math.round(n * 0.5);
        const wantConcept = Math.round(n * 0.3);
        const figs = drawFrom(FIG, ctx, wantFig, Math.max(2, Math.ceil(n / 7)));
        const concepts = drawConcepts(ctx, wantConcept);
        const rest = n - figs.length - concepts.length;
        const data = drawFrom(TEMPLATES, ctx, Math.max(rest, 0), Math.max(2, Math.ceil(n / 8)));
        let all = figs.concat(concepts, data);
        if (all.length < n) all = all.concat(drawFrom(FIG.concat(TEMPLATES), ctx, n - all.length, 3));
        // Shuffle first, then order by level, so the game ramps up without
        // always asking the same question in the same slot.
        return shuffle(all)
            .sort((a, b) => (a.level || 2) - (b.level || 2))
            .slice(0, n);
    }

    // ================================================================== UI
    function ensureDom() {
        // isConnected, not a plain null check: a host page that rewrites
        // document.body leaves S.root pointing at a detached node.
        if (S.root && S.root.isConnected) return S.root;
        if (!document.getElementById('cq-style')) {
            const st = document.createElement('style');
            st.id = 'cq-style';
            st.textContent = CSS;
            document.head.appendChild(st);
        }
        const root = document.createElement('div');
        root.id = 'cq-root';
        root.setAttribute('role', 'dialog');
        root.setAttribute('aria-label', 'Correlate Quest');
        root.innerHTML = '<button class="cq-corner" id="cq-sound" aria-label="Sound"></button>'
            + '<button class="cq-corner" id="cq-close" aria-label="Close">X</button>'
            + '<div class="cq-shell" id="cq-shell"></div>';
        document.body.appendChild(root);
        root.querySelector('#cq-close').addEventListener('click', () => api.close());
        const sb = root.querySelector('#cq-sound');
        sb.addEventListener('click', () => {
            S.sound = !S.sound;
            try { localStorage.setItem(SOUND_KEY, S.sound ? '1' : '0'); } catch (e) { }
            paintSound();
            if (S.sound) beep('blip');
        });
        S.root = root;
        try { S.sound = localStorage.getItem(SOUND_KEY) === '1'; } catch (e) { S.sound = false; }
        paintSound();
        if (!S.wired) {
            document.addEventListener('keydown', onKey, true);
            window.addEventListener('resize', onResize);
            S.wired = true;
        }
        return root;
    }
    function paintSound() {
        const b = S.root && S.root.querySelector('#cq-sound');
        if (b) {
            b.textContent = 'SFX';
            b.style.color = S.sound ? PAL.yellow : PAL.dim;
            b.title = S.sound ? 'Sound on' : 'Sound off';
        }
    }
    const shell = () => S.root.querySelector('#cq-shell');

    function onResize() {
        if (!S.opened || !S.figure || !S.figureEl) return;
        const w = Math.max(200, S.figureEl.clientWidth || 300);
        const cv = S.figureEl.querySelector('canvas');
        if (cv) { try { S.figure(cv, w); } catch (e) { } }
    }

    function onKey(e) {
        if (!S.opened) return;
        if (e.key === 'Escape') { e.stopPropagation(); api.close(); return; }
        if (S.screen !== 'play') return;
        if (!S.answered && /^[1-4]$/.test(e.key)) {
            const b = shell().querySelectorAll('.cq-ans')[Number(e.key) - 1];
            if (b) { e.preventDefault(); b.click(); }
            return;
        }
        if (S.answered && (e.key === 'Enter' || e.key === ' ')) {
            e.preventDefault();
            nextQuestion();
        }
    }

    function shake() {
        if (reduced()) return;
        const sh = shell();
        sh.classList.remove('cq-shake'); void sh.offsetWidth; sh.classList.add('cq-shake');
        S.root.classList.remove('cq-flash'); void S.root.offsetWidth; S.root.classList.add('cq-flash');
    }

    // ----------------------------------------------------------- screens
    function screenTitle() {
        S.screen = 'title';
        stopTimer();
        purgePlot();
        shell().classList.add('cq-center');
        shell().innerHTML =
            `<h1 class="cq-title">CORRELATE<br>QUEST</h1>
       <p class="cq-tag">A quiz built from ${num(D.ready() ? D.n() : 1208)} cancer cell lines</p>
       <p class="cq-blink">PRESS START</p>
       <button class="cq-btn cq-go" id="cq-start">START</button>
       <button class="cq-btn" id="cq-tourbtn">${TOUR.load() ? 'RESUME TOUR' : 'TOUR THE APP'}</button>
       <button class="cq-btn" id="cq-hs">HIGH SCORES</button>
       <button class="cq-btn cq-quiet" id="cq-how">HOW TO PLAY</button>
       <div id="cq-howbox"></div>`;
        shell().querySelector('#cq-start').onclick = () => { beep('blip'); screenSetup(); };
        shell().querySelector('#cq-tourbtn').onclick = () => { beep('blip'); TOUR.start(true); };
        shell().querySelector('#cq-hs').onclick = () => { beep('blip'); screenScores(null); };
        shell().querySelector('#cq-how').onclick = () => {
            const box = shell().querySelector('#cq-howbox');
            box.innerHTML = box.innerHTML ? '' :
                `<div class="cq-how">Ten questions, 20 seconds each, drawn from real cell line data.<br>
         A right answer scores 100, plus a bonus for speed and for a run of right answers.<br>
         Three lives. A wrong answer or a timeout costs one.</div>`;
        };
    }

    function screenSetup() {
        S.screen = 'setup';
        purgePlot();
        shell().classList.remove('cq-center');
        const counts = [...D.lineageCounts().entries()].filter(e => e[1] >= 10).sort((a, b) => b[1] - a[1]);
        let picked = { key: 'all', label: 'All cancers' };
        let nQ = 10;
        const draw = () => {
            shell().innerHTML =
                `<h2 class="cq-h2">PICK YOUR CELL LINES</h2>
         <div class="cq-list" id="cq-coh">
           <button class="cq-btn${picked.key === 'all' ? ' cq-sel' : ''}" data-k="all">ALL CANCERS (${num(D.n())})</button>
           ${counts.map(e => `<button class="cq-btn${picked.key === e[0] ? ' cq-sel' : ''}" data-k="${esc(e[0])}">${esc(e[0].toUpperCase())} (${num(e[1])})</button>`).join('')}
         </div>
         <h2 class="cq-h2">HOW MANY QUESTIONS</h2>
         <div class="cq-row">
           <button class="cq-btn${nQ === 10 ? ' cq-sel' : ''}" data-n="10">10</button>
           <button class="cq-btn${nQ === 20 ? ' cq-sel' : ''}" data-n="20">20</button>
         </div>
         <button class="cq-btn cq-go" id="cq-goBtn">GO</button>
         <button class="cq-btn cq-quiet" id="cq-back">BACK</button>`;
            shell().querySelectorAll('[data-k]').forEach(b => b.onclick = () => {
                const k = b.getAttribute('data-k');
                picked = k === 'all' ? { key: 'all', label: 'All cancers' } : { key: k, label: k };
                beep('blip'); draw();
            });
            shell().querySelectorAll('[data-n]').forEach(b => b.onclick = () => {
                nQ = Number(b.getAttribute('data-n')); beep('blip'); draw();
            });
            shell().querySelector('#cq-goBtn').onclick = () => startGame(picked, nQ);
            shell().querySelector('#cq-back').onclick = () => screenTitle();
        };
        draw();
    }

    function screenLoading(msg) {
        S.screen = 'loading';
        purgePlot();
        shell().classList.remove('cq-center');
        shell().innerHTML = `<div class="cq-spin"></div><p class="cq-tag">${esc(msg || 'LOADING...')}</p>`;
    }

    async function startGame(cohort, nQ) {
        screenLoading('LOADING...');
        const a = A();
        // Expression is a 51 MB lazy load. Wait for it, but never let a slow or
        // stalled load hold the game up: the expression questions simply drop.
        if (!a.expressionLoaded) {
            try {
                await Promise.race([
                    a.loadExpressionData(),
                    new Promise(r => setTimeout(r, 20000))
                ]);
            } catch (e) { /* expression questions drop out */ }
        }
        await new Promise(r => setTimeout(r, 30));
        const ids = D.ids();
        const idx = [];
        for (let i = 0; i < ids.length; i++) {
            if (cohort.key === 'all' || D.lineage(ids[i]) === cohort.key) idx.push(i);
        }
        S.cohort = { key: cohort.key, label: cohort.label, idx };
        S.nQ = nQ;
        const ctx = buildContext(S.cohort);
        S.bank = buildBank(ctx, nQ);
        if (!S.bank.length) {
            shell().innerHTML = '<p class="cq-tag">No questions could be built from this data. Try all cancers.</p>'
                + '<button class="cq-btn" id="cq-back">BACK</button>';
            shell().querySelector('#cq-back').onclick = () => screenSetup();
            return;
        }
        S.nQ = S.bank.length;
        S.qi = 0; S.score = 0; S.shown = 0; S.streak = 0; S.best = 0;
        S.lives = START_LIVES; S.correct = 0; S.lastEntry = null;
        renderQuestion();
    }

    function hudHtml() {
        const hearts = '*'.repeat(Math.max(S.lives, 0)) + '.'.repeat(Math.max(START_LIVES - S.lives, 0));
        return `<div class="cq-hud">
      <span>Q <b>${S.qi + 1}</b>/${S.nQ}</span>
      <span>SCORE <b id="cq-score">${num(S.shown)}</b></span>
      <span>RUN <b>${S.streak}</b></span>
      <span class="cq-hearts" aria-label="Lives">${hearts}</span>
    </div>`;
    }

    // A Plotly chart holds on to listeners and a WebGL-free but heavy DOM, so
    // it is purged before the card that holds it is thrown away.
    function purgePlot() {
        if (!S.plotDiv) return;
        try { if (typeof window.Plotly !== 'undefined') Plotly.purge(S.plotDiv); } catch (e) { }
        S.plotDiv = null;
    }

    function renderQuestion() {
        S.screen = 'play';
        shell().classList.remove('cq-center');
        purgePlot();
        S.answered = false; S.chosen = -1; S.gained = 0;
        S.shown = S.score;
        const q = S.bank[S.qi];
        S.figure = null; S.figureEl = null;
        shell().innerHTML = hudHtml()
            + `<div class="cq-bar"><i id="cq-timer"></i></div>`
            + `<div class="cq-cat">${esc(q.tag)}</div>`
            + (q.plot ? `<div class="cq-plot"><div id="cq-plot"></div></div>` : '')
            + `<p class="cq-q">${esc(q.text)}</p>`
            + (q.quote ? `<div class="cq-quote">${esc(q.quote)}</div>` : '')
            + `<div class="cq-fig" id="cq-figbox"></div>`
            + q.options.map((o, i) =>
                `<button class="cq-btn cq-ans" data-i="${i}"><span class="cq-key">${i + 1}</span><span>${esc(o.label)}</span></button>`).join('')
            + `<div id="cq-after"></div>`;
        if (q.figure && q.figureBefore) drawFigure(q.figure);
        shell().querySelectorAll('.cq-ans').forEach(b => {
            b.onclick = () => answer(Number(b.getAttribute('data-i')));
        });
        try { S.root.scrollTop = 0; } catch (e) { }
        // The clock waits for the chart: it is the question, and starting the
        // countdown against a blank box would be unfair.
        if (q.plot) {
            const div = shell().querySelector('#cq-plot');
            S.plotDiv = div;
            const go = () => { if (S.screen === 'play' && S.bank[S.qi] === q) startTimer(); };
            try {
                const r = q.plot(div);
                if (r && typeof r.then === 'function') r.then(go, go);
                else go();
            } catch (e) {
                const box = div && div.parentNode;
                if (box) box.remove();
                S.plotDiv = null;
                go();
            }
        } else {
            startTimer();
        }
    }

    function drawFigure(fn) {
        const box = shell().querySelector('#cq-figbox');
        if (!box) return;
        box.innerHTML = '<canvas></canvas>';
        const w = Math.max(200, box.clientWidth || 300);
        S.figure = fn; S.figureEl = box;
        try { fn(box.querySelector('canvas'), w); } catch (e) { box.innerHTML = ''; }
    }

    // ------------------------------------------------------------- timer
    function startTimer() {
        stopTimer();
        S.left = TIMER_MS;
        S.deadline = performance.now() + TIMER_MS;
        S.running = true;
        let lastTick = 6;
        const step = () => {
            if (!S.running) return;
            const left = Math.max(0, S.deadline - performance.now());
            S.left = left;
            const bar = shell().querySelector('#cq-timer');
            if (bar) {
                bar.style.width = (left / TIMER_MS * 100).toFixed(1) + '%';
                bar.classList.toggle('cq-low', left < 5000);
            }
            const secs = Math.ceil(left / 1000);
            if (secs <= 5 && secs !== lastTick) { lastTick = secs; beep('tick'); }
            if (left <= 0) { timeout(); return; }
            S.raf = requestAnimationFrame(step);
        };
        S.raf = requestAnimationFrame(step);
    }
    function stopTimer() {
        S.running = false;
        if (S.raf) cancelAnimationFrame(S.raf);
        S.raf = 0;
    }

    // The counter looks up its element every frame: the next question replaces
    // the whole HUD, and a held reference would leave a stale number on screen.
    function tickScore(target) {
        const step = () => {
            const el = shell().querySelector('#cq-score');
            if (!el) { S.shown = target; return; }
            if (S.shown >= target) { S.shown = target; el.textContent = num(target); return; }
            S.shown = Math.min(target, S.shown + Math.max(1, Math.ceil((target - S.shown) / 8)));
            el.textContent = num(S.shown);
            requestAnimationFrame(step);
        };
        step();
    }

    function timeout() { answer(-1); }

    function answer(i) {
        if (S.answered) return;
        S.answered = true;
        stopTimer();
        const q = S.bank[S.qi];
        const rightIdx = q.options.findIndex(o => o.correct);
        const ok = i === rightIdx;
        S.chosen = i;
        if (ok) {
            const timeBonus = Math.round(S.left / TIMER_MS * 100);
            const streakBonus = 25 * Math.min(S.streak, 5);
            S.gained = 100 + timeBonus + streakBonus;
            S.score += S.gained;
            S.streak++; S.correct++;
            S.best = Math.max(S.best, S.streak);
            beep('correct');
        } else {
            S.streak = 0;
            S.lives--;
            beep('wrong');
            shake();
        }
        shell().querySelectorAll('.cq-ans').forEach((b, k) => {
            b.setAttribute('disabled', 'disabled');
            b.onclick = null;
            if (k === rightIdx) b.classList.add('cq-right');
            else if (k === i) b.classList.add('cq-wrong');
        });
        const hearts = shell().querySelector('.cq-hearts');
        if (hearts) hearts.textContent = '*'.repeat(Math.max(S.lives, 0)) + '.'.repeat(Math.max(START_LIVES - S.lives, 0));
        tickScore(S.score);
        if (q.figure && !q.figureBefore) drawFigure(q.figure);
        const last = S.qi >= S.bank.length - 1 || S.lives <= 0;
        const verdict = ok
            ? `RIGHT +${num(S.gained)}`
            : (i === -1 ? 'OUT OF TIME' : 'WRONG');
        const after = shell().querySelector('#cq-after');
        after.innerHTML =
            `<div class="cq-exp" id="cq-exp">
        <div class="cq-verdict ${ok ? 'ok' : 'no'}">${verdict}</div>
        <p>${esc(q.explain)}</p>
        ${q.openInApp ? `<button class="cq-link" id="cq-open">${esc(q.openInApp.label)}</button>` : ''}
       </div>
       <button class="cq-btn cq-go" id="cq-next">${last ? 'SEE RESULT' : 'NEXT'}</button>`;
        const openBtn = after.querySelector('#cq-open');
        if (openBtn) openBtn.onclick = (ev) => {
            ev.stopPropagation();
            const a = A();
            api.close();
            try { q.openInApp.run(a); } catch (e) { console.warn('Correlate Quest could not open that view:', e); }
        };
        after.querySelector('#cq-exp').onclick = () => nextQuestion();
        after.querySelector('#cq-next').onclick = (ev) => { ev.stopPropagation(); nextQuestion(); };
        after.querySelector('#cq-next').scrollIntoView({ block: 'nearest' });
    }

    function nextQuestion() {
        if (!S.answered) return;
        if (S.lives <= 0 || S.qi >= S.bank.length - 1) { screenOver(); return; }
        S.qi++;
        renderQuestion();
    }

    // -------------------------------------------------------- game over
    function screenOver() {
        S.screen = 'over';
        purgePlot();
        shell().classList.remove('cq-center');
        stopTimer();
        beep('over');
        const isHigh = qualifies(S.score);
        shell().innerHTML =
            `<h2 class="cq-h2">GAME OVER</h2>
       <div class="cq-big">${num(S.score)}</div>
       <div class="cq-stat">
         RIGHT <span>${S.correct}</span> / ${S.qi + 1}<br>
         BEST RUN <span>${S.best}</span><br>
         LIVES LEFT <span>${Math.max(S.lives, 0)}</span>
       </div>
       <div id="cq-entry"></div>`;
        const box = shell().querySelector('#cq-entry');
        if (isHigh) initialsEntry(box);
        else {
            box.innerHTML = `<button class="cq-btn cq-go" id="cq-again">PLAY AGAIN</button>
        <button class="cq-btn" id="cq-copy">COPY RESULT</button>
        <button class="cq-btn cq-quiet" id="cq-scores">HIGH SCORES</button>
        <button class="cq-btn cq-quiet" id="cq-home">BACK TO TITLE</button>`;
            wireOverButtons(box);
        }
    }

    function initialsEntry(box) {
        const letters = ['A', 'A', 'A'];
        let slot = 0;
        const ABC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 ';
        const draw = () => {
            box.innerHTML = `<p class="cq-h2">NEW HIGH SCORE</p>
        <p class="cq-note">Enter your initials</p>
        <div class="cq-slots">${letters.map((l, i) =>
                `<div class="cq-slot">
             <button data-up="${i}">^</button>
             <div class="cq-ltr" data-slot="${i}" style="${i === slot ? 'border-color:' + PAL.yellow : ''}">${esc(l)}</div>
             <button data-dn="${i}">v</button>
           </div>`).join('')}</div>
        <button class="cq-btn cq-go" id="cq-save">SAVE</button>`;
            box.querySelectorAll('[data-up]').forEach(b => b.onclick = () => {
                const i = Number(b.getAttribute('data-up'));
                slot = i; letters[i] = ABC[(ABC.indexOf(letters[i]) + 1) % ABC.length]; beep('blip'); draw();
            });
            box.querySelectorAll('[data-dn]').forEach(b => b.onclick = () => {
                const i = Number(b.getAttribute('data-dn'));
                slot = i; letters[i] = ABC[(ABC.indexOf(letters[i]) - 1 + ABC.length) % ABC.length]; beep('blip'); draw();
            });
            box.querySelectorAll('[data-slot]').forEach(d => d.onclick = () => {
                slot = Number(d.getAttribute('data-slot')); draw();
            });
            box.querySelector('#cq-save').onclick = save;
        };
        const onType = (e) => {
            if (S.screen !== 'over') return;
            const k = e.key.toUpperCase();
            if (ABC.indexOf(k) >= 0 && k !== ' ') {
                letters[slot] = k; slot = Math.min(slot + 1, 2); draw(); e.preventDefault();
            } else if (e.key === 'ArrowRight') { slot = Math.min(slot + 1, 2); draw(); }
            else if (e.key === 'ArrowLeft') { slot = Math.max(slot - 1, 0); draw(); }
            else if (e.key === 'Enter') { save(); }
        };
        const save = () => {
            document.removeEventListener('keydown', onType);
            const entry = {
                initials: letters.join('').trim() || 'AAA',
                score: S.score,
                cohort: S.cohort.label,
                date: new Date().toISOString().slice(0, 10)
            };
            addScore(entry);
            S.lastEntry = entry;
            beep('correct');
            screenScores(entry);
        };
        document.addEventListener('keydown', onType);
        draw();
    }

    function wireOverButtons(box) {
        const again = box.querySelector('#cq-again');
        if (again) again.onclick = () => startGame(S.cohort, S.nQ);
        const home = box.querySelector('#cq-home');
        if (home) home.onclick = () => screenTitle();
        const sc = box.querySelector('#cq-scores');
        if (sc) sc.onclick = () => screenScores(S.lastEntry);
        const copy = box.querySelector('#cq-copy');
        if (copy) copy.onclick = () => copyResult(copy);
    }

    function copyResult(btn) {
        const url = location.origin + location.pathname + location.search;
        const txt = `I scored ${num(S.score)} in Correlate Quest (${S.cohort ? S.cohort.label : 'All cancers'}). Play at ${url}#quiz`;
        const done = () => { btn.textContent = 'COPIED!'; setTimeout(() => { btn.textContent = 'COPY RESULT'; }, 1600); };
        try {
            navigator.clipboard.writeText(txt).then(done, () => fallback());
        } catch (e) { fallback(); }
        function fallback() {
            const ta = document.createElement('textarea');
            ta.value = txt; ta.style.position = 'fixed'; ta.style.opacity = '0';
            document.body.appendChild(ta); ta.select();
            try { document.execCommand('copy'); done(); } catch (e2) { btn.textContent = 'COPY FAILED'; }
            document.body.removeChild(ta);
        }
    }

    function screenScores(highlight) {
        S.screen = 'scores';
        purgePlot();
        shell().classList.remove('cq-center');
        stopTimer();
        const list = loadScores();
        // The saved row comes back from storage as a new object, so the row
        // just played is matched on its fields, not on identity.
        const same = (e) => !!highlight && e.score === highlight.score
            && e.initials === highlight.initials && e.date === highlight.date;
        shell().innerHTML =
            `<h2 class="cq-h2">HIGH SCORES</h2>
       ${list.length ? `<table class="cq-table">${list.map((e, i) =>
                `<tr class="${same(e) ? 'cq-me' : ''}">
            <td>${i + 1}</td><td>${esc(e.initials || '???')}</td>
            <td class="cq-co">${esc(e.cohort || '')}<br>${esc(e.date || '')}</td>
            <td class="cq-r">${num(e.score)}</td></tr>`).join('')}</table>`
                : '<p class="cq-note">No scores yet. Play a round.</p>'}
       <p class="cq-note">High scores are stored on this device.</p>
       <button class="cq-btn cq-go" id="cq-again">PLAY AGAIN</button>
       <button class="cq-btn" id="cq-copy">COPY RESULT</button>
       <button class="cq-btn cq-quiet" id="cq-home">BACK TO TITLE</button>`;
        wireOverButtons(shell());
        if (!S.cohort) {
            const c = shell().querySelector('#cq-copy');
            if (c) c.remove();
            const ag = shell().querySelector('#cq-again');
            if (ag) ag.onclick = () => screenSetup();
        }
    }


    // ================================================================ tour
    // A guided walk through the app's real features. The quiz closes and a
    // compact bar takes over, watching the app's own state to know when a
    // mission is done. Nothing in app.js is patched or wrapped: the bar polls.
    const TOUR_KEY = 'correlateTourProgress';

    const el = (id) => document.getElementById(id);
    const shown = (id) => {
        const e = el(id);
        if (!e) return false;
        if (e.classList && e.classList.contains('active')) return true;
        const d = e.style.display;
        return d === 'flex' || d === 'block';
    };
    const val = (id) => (el(id)?.value || '');

    const MISSIONS = [
        {
            title: 'Open the Cell Line Browser',
            todo: 'Press Cell Line Browser in the Options card.',
            why: 'It is where you find the cell lines that fit a project.',
            hint: 'Top left of the page, next to Gene set analysis.',
            done: () => shown('cellLineBrowserModal')
        },
        {
            title: 'Narrow it to one tissue',
            todo: 'Choose a tissue in the Tissue dropdown.',
            why: 'Most questions in cancer biology are asked inside one cancer type.',
            hint: 'The dropdown at the top of the browser reads All tissues until you pick one.',
            done: () => shown('cellLineBrowserModal') && val('clbTissueFilter') !== ''
        },
        {
            title: 'Open a cell line page',
            todo: 'Click a cell line name in the list.',
            why: 'One page holds everything the app knows about that line.',
            hint: 'Any cell line name in the first column opens its page.',
            done: () => shown('clbWikiModal')
        },
        {
            title: 'Find the key genetic alterations',
            todo: 'Scroll the page down to Key genetic alterations.',
            why: 'The drivers are the first thing to check in a line you do not know.',
            hint: 'It sits about a third of the way down the page.',
            done: () => {
                if (!shown('clbWikiModal')) return false;
                const body = el('clbWikiBody');
                if (body) {
                    const hit = [...body.querySelectorAll('div,h3,h4,b,span')]
                        .find(n => /^Key genetic alterations/.test((n.textContent || '').trim()));
                    if (hit) {
                        const r = hit.getBoundingClientRect();
                        if (r.top < window.innerHeight * 0.9 && r.bottom > 0) return true;
                    }
                }
                // Some layouts scroll the modal itself rather than the body.
                for (const node of [el('clbWikiModal'), el('clbWikiBody')]) {
                    if (!node) continue;
                    const span = node.scrollHeight - node.clientHeight;
                    if (span > 40 && node.scrollTop > span * 0.3) return true;
                }
                return false;
            }
        },
        {
            title: 'Run a gene set analysis',
            todo: 'Close the wiki and the browser, put three or more genes in the gene box and press Run.',
            why: 'This is the core of the app: which genes rise and fall together across cell lines.',
            hint: 'Paste TP53, MDM2, CDKN1A, BAX, MDM4 into the gene box, then press Run.',
            arm: (ref) => { ref.results = A().results; },
            done: (ref) => {
                const r = A().results;
                return !!(r && r !== ref.results && r.correlations && r.correlations.length);
            }
        },
        {
            title: 'Lower the correlation cutoff',
            todo: 'Drag the Correlation Cutoff slider down, then press Run again.',
            why: 'The cutoff decides which links are strong enough to be drawn.',
            hint: 'The slider is in box 1, Parameters. Try 0.3, then press Run.',
            arm: (ref) => {
                ref.cutoff = parseFloat(val('correlationCutoff'));
                ref.results = A().results;
            },
            done: (ref) => {
                const now = parseFloat(val('correlationCutoff'));
                const r = A().results;
                return isFinite(now) && isFinite(ref.cutoff) && now < ref.cutoff
                    && !!r && r !== ref.results && !!r.correlations;
            }
        },
        {
            title: 'Open a correlation scatter',
            todo: 'Double-click a link in the network, or choose Other then Correlation.',
            why: 'The scatter shows the actual cell lines behind a link.',
            hint: 'In the network, double-click the line drawn between two genes.',
            done: () => shown('inspectModal')
        },
        {
            title: 'Color the scatter by a mutation',
            todo: 'In the scatter controls, pick a gene under Hotspot Overlay and set the mode to Color.',
            why: 'It shows whether a link only holds in the mutated cell lines.',
            hint: 'The Hotspot Overlay selectors sit in the scatter’s left-hand controls.',
            done: () => shown('inspectModal') && val('hotspotGene') !== '' && val('hotspotMode') === 'color'
        },
        {
            title: 'Open the gene effect view',
            todo: 'Choose Other, then Gene Effect, and pick a gene.',
            why: 'It shows one gene’s dependency across every tissue at once.',
            hint: 'Other is the third button in the Options card.',
            done: () => shown('geneEffectModal')
        },
        {
            title: 'Open the gene set heatmap',
            todo: 'Choose Other, then Gene set heatmap.',
            why: 'It puts many genes and many cell lines on one screen.',
            hint: 'Same Other menu, the item below Correlation.',
            done: () => shown('heatmapModal')
        },
        {
            title: 'Export a figure',
            desktopOnly: true,
            todo: 'Press Export image in any open view.',
            why: 'Figures leave the app ready to drop into a slide or a paper.',
            hint: 'Every popout has an Export image button in its top bar. Opening the dialog is enough.',
            done: () => shown('exportOptionsModal')
        },
        {
            title: 'Try a mutation analysis',
            todo: 'Choose Other, then Mutation analysis, pick a gene such as TP53 and press Run.',
            why: 'It ranks the genes whose dependency differs between mutated and wild-type lines.',
            hint: 'Pick TP53 under Hotspot Mutation, then press Run Mutation Analysis.',
            arm: (ref) => { ref.mut = A().mutationResults; },
            done: (ref) => !!A().mutationResults && A().mutationResults !== ref.mut
        }
    ];

    const isPhone = () => window.innerWidth <= 640;
    const tourMissions = () => MISSIONS.filter(m => !(m.desktopOnly && isPhone()));

    const TOUR = {
        idx: 0, points: 0, ref: {}, timer: 0, running: false, collapsed: false,
        hintShown: false, list: [], flash: false,

        load() {
            try {
                const raw = JSON.parse(localStorage.getItem(TOUR_KEY) || 'null');
                if (raw && typeof raw.idx === 'number') return raw;
            } catch (e) { }
            return null;
        },
        save() {
            try { localStorage.setItem(TOUR_KEY, JSON.stringify({ idx: this.idx, points: this.points })); }
            catch (e) { }
        },
        clear() { try { localStorage.removeItem(TOUR_KEY); } catch (e) { } },

        start(resume) {
            const saved = resume ? this.load() : null;
            this.list = tourMissions();
            this.idx = saved ? Math.min(saved.idx, this.list.length - 1) : 0;
            this.points = saved ? (saved.points || 0) : 0;
            this.running = true;
            this.collapsed = false;
            api.close();
            this.mount();
            this.enter();
            if (this.timer) clearInterval(this.timer);
            this.timer = setInterval(() => this.tick(), 500);
        },

        mount() {
            let bar = el('cq-tour');
            if (!bar) {
                bar = document.createElement('div');
                bar.id = 'cq-tour';
                bar.setAttribute('role', 'complementary');
                bar.setAttribute('aria-label', 'Correlate Quest tour');
                document.body.appendChild(bar);
            }
            bar.style.display = 'block';
            this.bar = bar;
            this.pad();
        },

        // The bar sits over the page, so the page needs room to scroll clear
        // of it. Without this the Run button hides underneath it on a phone.
        pad() {
            document.body.style.paddingBottom = isPhone() ? (this.collapsed ? '72px' : '34vh') : '';
        },

        enter() {
            const m = this.list[this.idx];
            this.hintShown = false;
            this.ref = {};
            if (m && m.arm) { try { m.arm(this.ref); } catch (e) { } }
            this.render();
        },

        render() {
            if (!this.bar) return;
            const m = this.list[this.idx];
            if (!m) return;
            const squares = this.list.map((_, i) =>
                `<i class="cq-sq${i < this.idx ? ' cq-sq-on' : ''}${i === this.idx ? ' cq-sq-now' : ''}"></i>`).join('');
            this.bar.className = this.collapsed ? 'cq-collapsed' : '';
            this.bar.innerHTML =
                `<div class="cq-tour-head">
           <span class="cq-tour-n">QUEST ${this.idx + 1}/${this.list.length}</span>
           <span class="cq-tour-pts">${num(this.points)}</span>
           <button class="cq-tour-chev" id="cq-tour-chev" aria-label="Collapse">${this.collapsed ? '^' : 'v'}</button>
         </div>
         <div class="cq-tour-body">
           <div class="cq-tour-title">${esc(m.title)}</div>
           <div class="cq-tour-do">${esc(m.todo)}</div>
           <div class="cq-tour-why">${esc(m.why)}</div>
           <div class="cq-tour-hint" id="cq-tour-hint">${this.hintShown ? esc(m.hint) : ''}</div>
           <div class="cq-tour-btns">
             <button class="cq-tour-b" id="cq-tour-hintbtn">HINT</button>
             <button class="cq-tour-b" id="cq-tour-skip">SKIP</button>
             <button class="cq-tour-b" id="cq-tour-exit">EXIT TOUR</button>
           </div>
           <div class="cq-tour-prog">${squares}</div>
         </div>`;
            el('cq-tour-chev').onclick = () => {
                this.collapsed = !this.collapsed; beep('blip'); this.pad(); this.render();
            };
            el('cq-tour-hintbtn').onclick = () => { this.hintShown = true; beep('blip'); this.render(); };
            el('cq-tour-skip').onclick = () => { beep('blip'); this.advance(0); };
            el('cq-tour-exit').onclick = () => this.exit(true);
        },

        tick() {
            if (!this.running || this.flash) return;
            const m = this.list[this.idx];
            if (!m) return;
            let ok = false;
            try { ok = !!m.done(this.ref); } catch (e) { ok = false; }
            if (ok) this.advance(100);
        },

        advance(points) {
            this.points += points;
            this.save();
            if (points > 0) {
                beep('fanfare');
                this.flash = true;
                const body = this.bar && this.bar.querySelector('.cq-tour-body');
                if (body) body.innerHTML = `<div class="cq-tour-done">DONE +${points}</div>`;
                const head = this.bar && this.bar.querySelector('.cq-tour-pts');
                if (head) head.textContent = num(this.points);
                setTimeout(() => { this.flash = false; this.step(); }, 1200);
            } else {
                this.step();
            }
        },

        step() {
            this.idx++;
            if (this.idx >= this.list.length) { this.finish(); return; }
            this.save();
            this.enter();
        },

        exit(keep) {
            this.running = false;
            if (this.timer) { clearInterval(this.timer); this.timer = 0; }
            if (this.bar) this.bar.style.display = 'none';
            document.body.style.paddingBottom = '';
            if (keep) this.save(); else this.clear();
        },

        finish() {
            const pts = this.points;
            this.exit(false);
            this.clear();
            ensureDom();
            S.opened = true;
            S.root.style.display = 'block';
            document.body.style.overflow = 'hidden';
            screenBadge(pts);
        }
    };

    function screenBadge(points) {
        S.screen = 'badge';
        purgePlot();
        shell().classList.add('cq-center');
        shell().innerHTML =
            `<h1 class="cq-title">CORRELATE<br>MASTER</h1>
       <p class="cq-tag">You have been everywhere in this app</p>
       <div class="cq-big">${num(points)}</div>
       <div class="cq-stat">MISSIONS <span>${tourMissions().length}</span></div>
       <button class="cq-btn cq-go" id="cq-badge-copy">COPY RESULT</button>
       <button class="cq-btn cq-quiet" id="cq-badge-home">BACK TO TITLE</button>`;
        shell().querySelector('#cq-badge-home').onclick = () => screenTitle();
        shell().querySelector('#cq-badge-copy').onclick = (ev) => {
            const btn = ev.currentTarget;
            const url = location.origin + location.pathname + location.search;
            const txt = `I finished the Correlate Quest tour with ${num(points)} points. Try it at ${url}#tour`;
            const done = () => { btn.textContent = 'COPIED!'; setTimeout(() => { btn.textContent = 'COPY RESULT'; }, 1600); };
            try { navigator.clipboard.writeText(txt).then(done, done); } catch (e) { done(); }
        };
    }

    // ================================================================ api
    const api = {
        open() {
            if (!D.ready()) {
                alert('The quiz needs the cell line data, which is still loading. Try again in a moment.');
                return;
            }
            ensureDom();
            S.opened = true;
            S.root.style.display = 'block';
            document.body.style.overflow = 'hidden';
            screenTitle();
        },
        close() {
            stopTimer();
            purgePlot();
            S.opened = false;
            if (S.root) S.root.style.display = 'none';
            document.body.style.overflow = '';
        },
        isOpen() { return S.opened; },
        // The guided tour, reachable from the title screen and the #tour route.
        tour() {
            if (!D.ready()) {
                alert('The tour needs the cell line data, which is still loading. Try again in a moment.');
                return;
            }
            ensureDom();
            TOUR.start(true);
        },
        // Read-only view of the questions in the current game, used by the
        // verification scripts. Never written to.
        peek() { return S.bank; }
    };

    window.CorrelateQuiz = api;
})();
