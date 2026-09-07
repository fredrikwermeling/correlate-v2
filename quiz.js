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

    // Genes worth naming in a question: the curated pathway panels, the
    // subtype hallmark lists and the hotspot mutation panel. Random genes out
    // of 18,000 make unreadable answer options.
    let NOTABLE = null;
    function notableGenes() {
        if (NOTABLE) return NOTABLE;
        const set = new Set();
        const pw = D.pathways();
        Object.keys(pw).forEach(k => (pw[k].genes || []).forEach(g => set.add(g)));
        const hl = D.hallmarks();
        Object.keys(hl).forEach(k => (hl[k].lookFor || []).forEach(g => set.add(g)));
        (A().mutations?.genes || Object.keys(A().mutations?.geneData || {})).forEach(g => set.add(g));
        SELECTIVE.forEach(g => set.add(g));
        const a = A();
        NOTABLE = [...set].filter(g => !a._isPolymorphicLocus?.(g) && a.geneIndex.has(g.toUpperCase()));
        return NOTABLE;
    }

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

    // Rank pairs with a cheap mean-imputed dot product first, then recompute
    // the honest pairwise-complete r for the handful that survive.
    function buildCorrelationPool() {
        const genes = notableGenes();
        const rows = [], names = [];
        for (const g of genes) {
            const row = D.row(g);
            if (!row) continue;
            let c = 0, s = 0, ss = 0;
            for (let i = 0; i < row.length; i++) {
                const v = D.ge(row, i);
                if (!isFinite(v)) continue;
                c++; s += v; ss += v * v;
            }
            if (c < 400) continue;
            const mu = s / c, sd = Math.sqrt(Math.max(ss / c - mu * mu, 0));
            if (sd < 0.15) continue;
            const z = new Float32Array(row.length);
            for (let i = 0; i < row.length; i++) {
                const v = D.ge(row, i);
                z[i] = isFinite(v) ? (v - mu) / sd : 0;
            }
            rows.push(z); names.push(g);
        }
        const nCL = D.n();
        const cand = [];
        for (let i = 0; i < rows.length; i++) {
            for (let j = i + 1; j < rows.length; j++) {
                let d = 0;
                const a = rows[i], b = rows[j];
                for (let k = 0; k < nCL; k++) d += a[k] * b[k];
                cand.push([d / nCL, i, j]);
            }
        }
        cand.sort((p, q) => p[0] - q[0]);
        const out = { pos: [], neg: [], zero: [] };
        const take = (list, bucket, want, test) => {
            for (const c of list) {
                if (bucket.length >= want) break;
                const g1 = names[c[1]], g2 = names[c[2]];
                const st = pearson(D.row(g1), D.row(g2));
                if (!st || !test(st.r)) continue;
                bucket.push({ g1, g2, r: st.r, n: st.n });
            }
        };
        take(cand.slice(0, 60), out.neg, 6, r => r <= -0.5);
        take(cand.slice(-60).reverse(), out.pos, 8, r => r >= 0.5);
        const mid = shuffle(cand.filter(c => Math.abs(c[0]) < 0.05)).slice(0, 60);
        take(mid, out.zero, 6, r => Math.abs(r) <= 0.1);
        return out;
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

    function scatterFigure(xs, ys, xLabel, yLabel) {
        return function (canvas, w) {
            const h = Math.min(Math.max(w * 0.78, 200), 300);
            const ctx = fitCanvas(canvas, w, h);
            ctx.fillStyle = PAL.panel; ctx.fillRect(0, 0, w, h);
            const L = 34, R = 10, T = 10, B = 34;
            let xlo = Infinity, xhi = -Infinity, ylo = Infinity, yhi = -Infinity;
            for (let i = 0; i < xs.length; i++) {
                xlo = Math.min(xlo, xs[i]); xhi = Math.max(xhi, xs[i]);
                ylo = Math.min(ylo, ys[i]); yhi = Math.max(yhi, ys[i]);
            }
            if (!(xhi > xlo)) { xhi = xlo + 1; }
            if (!(yhi > ylo)) { yhi = ylo + 1; }
            const px = v => L + (v - xlo) / (xhi - xlo) * (w - L - R);
            const py = v => (h - B) - (v - ylo) / (yhi - ylo) * (h - T - B);
            ctx.fillStyle = PAL.dim;
            ctx.fillRect(L - 3, T, 3, h - T - B + 3);
            ctx.fillRect(L - 3, h - B, w - L - R + 3, 3);
            ctx.fillStyle = PAL.blue;
            for (let i = 0; i < xs.length; i++) {
                ctx.fillRect(Math.round(px(xs[i])) - 2, Math.round(py(ys[i])) - 2, 4, 4);
            }
            ctx.fillStyle = PAL.yellow; ctx.font = pixFont(8); ctx.textBaseline = 'alphabetic';
            ctx.textAlign = 'center';
            ctx.fillText(clip(xLabel, Math.floor(w / 9)), (w + L) / 2, h - 8);
            ctx.save();
            ctx.translate(11, (h - B + T) / 2);
            ctx.rotate(-Math.PI / 2);
            ctx.fillText(clip(yLabel, Math.floor((h - T - B) / 9)), 0, 0);
            ctx.restore();
        };
    }

    // --------------------------------------------------------------- state
    const S = {
        opened: false, root: null, screen: 'title',
        cohort: null, nQ: 10,
        bank: [], qi: 0, score: 0, shown: 0, streak: 0, best: 0,
        lives: START_LIVES, correct: 0, answered: false, chosen: -1, gained: 0,
        deadline: 0, left: TIMER_MS, raf: 0, running: false,
        sound: false, actx: null, figure: null, figureEl: null, lastEntry: null, wired: false
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
                blip: [[520, 0]]
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
    function ctxCellLines(cohort) { return cohort.idx; }

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

    // 5. CORRELATION SIGN --------------------------------------------------
    T('corr', 'CORRELATION', (c) => !!(c.corr && (c.corr.pos.length || c.corr.neg.length || c.corr.zero.length)),
        (c) => tryTimes(30, () => {
            const buckets = [];
            if (c.corr.pos.length) buckets.push(['Positive', c.corr.pos]);
            if (c.corr.neg.length) buckets.push(['Negative', c.corr.neg]);
            if (c.corr.zero.length) buckets.push(['No correlation', c.corr.zero]);
            if (!buckets.length) return null;
            const [answer, list] = pick(buckets);
            const idx = ri(list.length);
            const p = list[idx];
            if (c.usedGene.has(p.g1) || c.usedGene.has(p.g2)) return null;
            list.splice(idx, 1);
            const r1 = D.row(p.g1), r2 = D.row(p.g2);
            const xs = [], ys = [];
            for (let i = 0; i < D.n(); i++) {
                const x = D.ge(r1, i), y = D.ge(r2, i);
                if (isFinite(x) && isFinite(y)) { xs.push(x); ys.push(y); }
            }
            if (xs.length < 200) return null;
            const opts = options4(answer, shuffle(['Positive', 'Negative', 'No correlation', 'Cannot tell']
                .filter(o => o !== answer)));
            if (!opts) return null;
            const word = answer === 'Positive'
                ? `lines that need ${p.g1} tend to need ${p.g2} as well`
                : answer === 'Negative'
                    ? `lines that need ${p.g1} tend not to need ${p.g2}`
                    : `knowing one score tells you nothing about the other`;
            return {
                tag: 'CORRELATION', text: `Each dot is one cell line. How do the gene effect scores of ${p.g1} and ${p.g2} relate?`,
                options: opts,
                explain: `The correlation is r = ${p.r.toFixed(2)} across ${num(p.n)} cell lines, so ${word}.`,
                figure: scatterFigure(xs, ys, p.g1, p.g2),
                figureBefore: true,
                usedGene: [p.g1, p.g2],
                openInApp: { label: 'See ' + p.g1 + ' across tissues', run: (a2) => a2.openGeneEffectModal(p.g1, 'tissue') }
            };
        }));

    // 6. MUTATION FREQUENCY ------------------------------------------------
    const HOTSPOT_GENES = ['BRAF', 'KRAS', 'NRAS', 'TP53', 'PIK3CA', 'EGFR', 'IDH1', 'CTNNB1', 'FLT3', 'JAK2', 'KIT'];
    T('mutfreq', 'MUTATIONS', (c) => c.lineages.length >= 4 && !!A().mutations?.geneData, (c) => tryTimes(60, () => {
        const gene = pick(HOTSPOT_GENES);
        if (c.usedGene.has(gene) || !A().mutations.geneData[gene]) return null;
        const rows = c.lineages.map(l => {
            let mut = 0;
            for (const id of l.ids) if (D.hotspot(gene, id) >= 1) mut++;
            return { name: l.name, frac: mut / l.ids.length, n: l.ids.length, mut };
        }).sort((a, b) => b.frac - a.frac);
        const top = rows[0];
        if (!top || top.frac < 0.12) return null;
        const low = rows.slice(1).filter(r => r.frac <= top.frac / 2);
        if (low.length < 3) return null;
        const wrong = shuffle(low).slice(0, 3);
        const opts = options4(top.name, wrong.map(w => w.name));
        if (!opts) return null;
        const bars = shuffle([top].concat(wrong)).map(r => ({
            label: r.name, value: Math.round(r.frac * 1000) / 10, hi: r.name === top.name
        }));
        return {
            tag: 'MUTATIONS', text: `In which tissue is ${gene} most often hotspot mutated?`,
            options: opts,
            explain: `${Math.round(top.frac * 100)} out of every 100 ${top.name.toLowerCase()} lines carry a ${gene} hotspot mutation (${top.mut} of ${top.n}).`,
            figure: barFigure(bars, { title: gene + ' hotspot, percent of lines', fmt: v => v.toFixed(0) + '%' }),
            figureBefore: false,
            usedGene: [gene]
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
        const good = (entry.lookFor || []).filter(g => !c.usedGene.has(g));
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
                if (c.usedGene.has(g)) continue;
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
        const gene = pick(mine.filter(g => !c.usedGene.has(g)) || []);
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
    T('drug', 'DRUG SCREEN', () => D.compounds().length >= 6, (c) => tryTimes(80, () => {
        const i = pick(ctxCellLines(c));
        const id = D.ids()[i];
        if (c.usedCL.has(id)) return null;
        const comps = D.compounds();
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
        const names = [hit.cp.name].concat(wrong.map(w => w.cp.name));
        if (!distinctNames(names)) return null;
        const opts = options4(hit.cp.name, wrong.map(w => w.cp.name));
        if (!opts) return null;
        const bars = shuffle([hit].concat(wrong)).map(x => ({
            label: x.cp.name, value: Math.round(x.v * 100) / 100, hi: x.cp.name === hit.cp.name
        }));
        const moa = hit.cp.moa || hit.cp.target || '';
        return {
            tag: 'DRUG SCREEN', text: `Which compound kills ${D.name(id)} best in the drug screen?`,
            options: opts,
            explain: `${hit.cp.name} leaves ${D.name(id)} at ${hit.v.toFixed(2)} on a scale where 1 means the cells are untouched.${moa ? ' It is ' + (/^[aeiou]/i.test(moa) ? 'an ' : 'a ') + moa.toLowerCase() + (hit.cp.target ? ', aimed at ' + hit.cp.target : '') + '.' : ''}`,
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
        const calls = D.fusions(hitId);
        if (!calls.length) return null;
        const call = calls[0];
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

        const fusionLines = [], noFusionLines = [];
        cohort.idx.forEach(i => {
            const id = ids[i];
            (D.fusions(id).length ? fusionLines : noFusionLines).push(id);
        });

        const burden = [];
        cohort.idx.forEach(i => {
            const n = D.damaging(ids[i]);
            if (typeof n === 'number' && n > 0) burden.push({ id: ids[i], n });
        });
        burden.sort((a, b) => b.n - a.n);

        return {
            key: cohort.key, idx: cohort.idx, lineages: lineageList, lineageCount: counts,
            subtypes, fusionLines, noFusionLines, burden,
            corr: null, usedCL: new Set(), usedGene: new Set(), usedLineage: new Set()
        };
    }

    function commit(ctx, q) {
        (q.usedCL || []).forEach(x => ctx.usedCL.add(x));
        (q.usedGene || []).forEach(x => ctx.usedGene.add(x));
        (q.usedLineage || []).forEach(x => ctx.usedLineage.add(x));
    }

    function buildBank(ctx, n) {
        const cap = Math.ceil(n / 10) * 3;
        const used = {}, dead = {};
        const out = [];
        const live = () => TEMPLATES.filter(t => !dead[t.id] && (used[t.id] || 0) < cap);
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
            used[t.id] = (used[t.id] || 0) + 1;
            commit(ctx, q);
            out.push(q);
        }
        return out;
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
        shell().classList.add('cq-center');
        shell().innerHTML =
            `<h1 class="cq-title">CORRELATE<br>QUEST</h1>
       <p class="cq-tag">A quiz built from ${num(D.ready() ? D.n() : 1208)} cancer cell lines</p>
       <p class="cq-blink">PRESS START</p>
       <button class="cq-btn cq-go" id="cq-start">START</button>
       <button class="cq-btn" id="cq-hs">HIGH SCORES</button>
       <button class="cq-btn cq-quiet" id="cq-how">HOW TO PLAY</button>
       <div id="cq-howbox"></div>`;
        shell().querySelector('#cq-start').onclick = () => { beep('blip'); screenSetup(); };
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
        try { ctx.corr = buildCorrelationPool(); } catch (e) { ctx.corr = { pos: [], neg: [], zero: [] }; }
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

    function renderQuestion() {
        S.screen = 'play';
        shell().classList.remove('cq-center');
        S.answered = false; S.chosen = -1; S.gained = 0;
        S.shown = S.score;
        const q = S.bank[S.qi];
        S.figure = null; S.figureEl = null;
        shell().innerHTML = hudHtml()
            + `<div class="cq-bar"><i id="cq-timer"></i></div>`
            + `<div class="cq-cat">${esc(q.tag)}</div>`
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
        startTimer();
        try { S.root.scrollTop = 0; } catch (e) { }
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
            S.opened = false;
            if (S.root) S.root.style.display = 'none';
            document.body.style.overflow = '';
        },
        isOpen() { return S.opened; }
    };

    window.CorrelateQuiz = api;
})();
