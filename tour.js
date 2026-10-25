// A tour of Correlate: fifteen pages in six sections, each with a real
// chart drawn from the loaded data and a button that opens the same view
// in the app, with the same filters and colouring. Charts are drawn on
// canvas here, sized to the card, so a phone gets the same picture as a
// desktop. Reads window.app only; the current page is kept in
// localStorage so the tour reopens where it was left.
(function () {
    'use strict';

    const A = () => window.app || null;
    const STEP_KEY = 'correlateTourStep';
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const phone = () => window.innerWidth <= 640;
    const gi = (g) => `<i>${esc(g)}</i>`;

    // ------------------------------------------------------------ data layer
    const D = {
        ready() {
            const a = A();
            return !!(a && a.metadata && a.metadata.cellLines && a.geneEffects && a.geneIndex);
        },
        ids() { return A().metadata.cellLines; },
        nGenes() { return A().nGenes || A().metadata.genes.length; },
        nExprLines() { return A().expressionMetadata?.nCellLines || 0; },
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
        // A gene's mRNA over the CRISPR cohort, in the same order as a GE row.
        exprRow(gene) {
            if (!A().expressionGeneIndex?.has(String(gene).toUpperCase())) return null;
            const n = D.ids().length, out = new Float64Array(n);
            for (let i = 0; i < n; i++) out[i] = D.expr(gene, i);
            return out;
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
        },
        async expression() {
            const a = A();
            if (a.expressionLoaded) return true;
            try { await a.loadExpressionData(); } catch (e) { }
            return !!a.expressionLoaded;
        }
    };
    const median = (arr) => {
        if (!arr.length) return NaN;
        const s = arr.slice().sort((a, b) => a - b);
        const m = s.length >> 1;
        return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
    };
    const quantile = (sorted, q) => {
        if (!sorted.length) return NaN;
        const pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
        return lo === hi ? sorted[lo] : sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
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
        const xs = [], ys = [], idx = [];
        for (let i = 0; i < D.ids().length; i++) {
            const x = D.ge(r1, i), y = D.ge(r2, i);
            if (isFinite(x) && isFinite(y)) { xs.push(x); ys.push(y); idx.push(i); }
        }
        return { xs, ys, idx };
    }
    // The same scatter split three ways by a gene's hotspot mutation.
    function scatterByMutation(g1, g2, mutGene) {
        const ids = D.ids(), rx = D.row(g1), ry = D.row(g2);
        const g = { wt: { xs: [], ys: [] }, m1: { xs: [], ys: [] }, m2: { xs: [], ys: [] } };
        for (let i = 0; i < ids.length; i++) {
            const x = D.ge(rx, i), y = D.ge(ry, i);
            if (!isFinite(x) || !isFinite(y)) continue;
            const lvl = D.hotspot(mutGene, ids[i]);
            const b = lvl >= 2 ? g.m2 : lvl === 1 ? g.m1 : g.wt;
            b.xs.push(x); b.ys.push(y);
        }
        return g;
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
    // Pairs above the cutoff among a gene set, on gene effect or on mRNA.
    function setEdges(genes, rows, cutoff) {
        const edges = [];
        for (let a = 0; a < genes.length; a++) for (let b = a + 1; b < genes.length; b++) {
            const st = pearson(rows[a], rows[b]);
            if (st && Math.abs(st.r) >= cutoff) edges.push({ a: genes[a], b: genes[b], r: st.r });
        }
        return edges;
    }
    const meanOf = (row) => { let n = 0, s = 0; for (let i = 0; i < D.ids().length; i++) { const v = D.ge(row, i); if (isFinite(v)) { n++; s += v; } } return n ? s / n : NaN; };

    // Welch's t-test of every gene between two groups of cell line indices,
    // the test the app's Inspect selection runs.
    function welchAll(selIdx, othIdx) {
        const a = A(), N = a.nCellLines, ge = a.geneEffects, out = [];
        for (let g = 0; g < a.nGenes; g++) {
            const off = g * N;
            let n1 = 0, s1 = 0, q1 = 0, n2 = 0, s2 = 0, q2 = 0;
            for (const i of selIdx) { const v = ge[off + i]; if (isFinite(v) && v > -990) { n1++; s1 += v; q1 += v * v; } }
            for (const i of othIdx) { const v = ge[off + i]; if (isFinite(v) && v > -990) { n2++; s2 += v; q2 += v * v; } }
            if (n1 < 3 || n2 < 3) continue;
            const m1 = s1 / n1, m2 = s2 / n2;
            const v1 = Math.max(0, (q1 - n1 * m1 * m1) / (n1 - 1)), v2 = Math.max(0, (q2 - n2 * m2 * m2) / (n2 - 1));
            const se = Math.sqrt(v1 / n1 + v2 / n2);
            if (!(se > 0)) continue;
            const t = (m1 - m2) / se;
            const df = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1));
            const p = a.tDistributionPValue(Math.abs(t), df);
            out.push({ gene: a.geneNames[g], d: m1 - m2, p });
        }
        return out;
    }

    // ------------------------------------------------------- canvas charts
    // The app's look: Arial, black zero lines, a light plot area, green
    // regression line, grey dots, blue and red for one and two mutated copies.
    const FONT = 'Arial, Helvetica, sans-serif';
    const C = { ink: '#374151', dim: '#6b7280', grid: '#e5e7eb', plotBg: '#fafafa', dot: '#9ca3af',
        green: '#6ba544', blue: '#3b82f6', red: '#dc2626', bar: 'rgba(122, 185, 80, 0.85)', amber: '#b45309' };
    function setupCanvas(div, cssH) {
        div.innerHTML = '';
        const W = Math.max(240, div.clientWidth || 360);
        const dpr = Math.min(2, window.devicePixelRatio || 1);
        const c = document.createElement('canvas');
        c.width = Math.round(W * dpr); c.height = Math.round(cssH * dpr);
        c.style.width = W + 'px'; c.style.height = cssH + 'px'; c.style.display = 'block';
        div.appendChild(c);
        const ctx = c.getContext('2d');
        ctx.scale(dpr, dpr);
        ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, W, cssH);
        return { c, ctx, W, H: cssH };
    }
    function niceTicks(lo, hi, want) {
        const span = hi - lo;
        if (!(span > 0)) return [lo];
        const raw = span / (want || 5);
        const p = Math.pow(10, Math.floor(Math.log10(raw)));
        const f = raw / p;
        const step = f < 1.5 ? p : f < 3 ? 2 * p : f < 7 ? 5 * p : 10 * p;
        const out = [];
        for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) out.push(+v.toFixed(10));
        return out;
    }
    const fmtTick = (v) => Number.isInteger(v) ? String(v) : String(+v.toFixed(2));
    function pad(vals) {
        let lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
        if (!(hi > lo)) { lo -= 1; hi += 1; }
        const m = (hi - lo) * 0.05;
        return [lo - m, hi + m];
    }
    function drawHeader(ctx, W, title, sub) {
        ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic'; ctx.fillStyle = C.ink;
        ctx.font = `bold ${phone() ? 12 : 14}px ${FONT}`;
        ctx.fillText(title, W / 2, 18);
        if (sub) { ctx.font = `10px ${FONT}`; ctx.fillStyle = C.dim; ctx.fillText(sub, W / 2, 33); }
        return sub ? 44 : 30;
    }
    // Axes, gridlines and labels for a numeric x axis; y is numeric unless
    // rows are given (a category axis with one label per row).
    function drawFrame(ctx, ar, xr, yr, o) {
        ctx.fillStyle = o.plain ? '#ffffff' : C.plotBg; ctx.fillRect(ar.x, ar.y, ar.w, ar.h);
        const sx = (v) => ar.x + (v - xr[0]) / (xr[1] - xr[0]) * ar.w;
        const sy = (v) => ar.y + ar.h - (v - yr[0]) / (yr[1] - yr[0]) * ar.h;
        const tf = phone() ? 10 : 11;
        ctx.font = `${tf}px ${FONT}`; ctx.fillStyle = C.ink; ctx.strokeStyle = C.grid; ctx.lineWidth = 1;
        for (const t of niceTicks(xr[0], xr[1], phone() ? 4 : 6)) {
            const x = Math.round(sx(t)) + 0.5;
            if (!o.plain) { ctx.beginPath(); ctx.moveTo(x, ar.y); ctx.lineTo(x, ar.y + ar.h); ctx.stroke(); }
            ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmtTick(t), x, ar.y + ar.h + 5);
        }
        if (o.rows) {
            const rh = ar.h / o.rows.length;
            ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
            o.rows.forEach((name, i) => ctx.fillText(name, ar.x - 6, ar.y + rh * (i + 0.5)));
        } else if (o.yTicks !== false) {
            for (const t of niceTicks(yr[0], yr[1], 5)) {
                const y = Math.round(sy(t)) + 0.5;
                ctx.beginPath(); ctx.moveTo(ar.x, y); ctx.lineTo(ar.x + ar.w, y); ctx.stroke();
                ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmtTick(t), ar.x - 6, y);
            }
        }
        ctx.strokeStyle = '#000'; ctx.lineWidth = 2;
        if (o.zeroX !== false && xr[0] < 0 && xr[1] > 0) { const x = Math.round(sx(0)); ctx.beginPath(); ctx.moveTo(x, ar.y); ctx.lineTo(x, ar.y + ar.h); ctx.stroke(); }
        if (o.zeroY && !o.rows && yr[0] < 0 && yr[1] > 0) { const y = Math.round(sy(0)); ctx.beginPath(); ctx.moveTo(ar.x, y); ctx.lineTo(ar.x + ar.w, y); ctx.stroke(); }
        if (o.baseline) { ctx.strokeStyle = '#d1d5db'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ar.x, ar.y + ar.h + 0.5); ctx.lineTo(ar.x + ar.w, ar.y + ar.h + 0.5); ctx.stroke(); }
        ctx.fillStyle = C.ink; ctx.font = `${phone() ? 11 : 12}px ${FONT}`;
        if (o.xLabel) { ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic'; ctx.fillText(o.xLabel, ar.x + ar.w / 2, ar.y + ar.h + 5 + tf + 16); }
        if (o.yLabel) {
            ctx.save(); ctx.translate(14, ar.y + ar.h / 2); ctx.rotate(-Math.PI / 2);
            ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(o.yLabel, 0, 0); ctx.restore();
        }
        return { sx, sy };
    }
    function drawLegend(ctx, x, y, items) {
        ctx.font = `${phone() ? 10 : 11}px ${FONT}`; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
        let cx = x;
        for (const it of items) {
            ctx.globalAlpha = it.hidden ? 0.35 : 1;
            ctx.fillStyle = it.color; ctx.beginPath(); ctx.arc(cx + 5, y, 4, 0, Math.PI * 2); ctx.fill();
            ctx.fillStyle = C.ink; ctx.fillText(it.name, cx + 14, y);
            cx += 14 + ctx.measureText(it.name).width + 16;
            ctx.globalAlpha = 1;
        }
    }
    // Colour between stops, t in [0,1].
    function lerpColor(stops, t) {
        t = Math.max(0, Math.min(1, t));
        let a = stops[0], b = stops[stops.length - 1];
        for (let i = 0; i < stops.length - 1; i++) if (t >= stops[i][0] && t <= stops[i + 1][0]) { a = stops[i]; b = stops[i + 1]; break; }
        const f = b[0] === a[0] ? 0 : (t - a[0]) / (b[0] - a[0]);
        const hex = (c) => [1, 3, 5].map(i => parseInt(c.slice(i, i + 2), 16));
        const ca = hex(a[1]), cb = hex(b[1]);
        return `rgb(${ca.map((v, i) => Math.round(v + (cb[i] - v) * f)).join(',')})`;
    }
    const GE_STOPS = [[0, '#e66101'], [0.5, '#f7f7f7'], [1, '#5e3c99']];
    const EXPR_STOPS = [[0, '#2166ac'], [0.5, '#f7f7f7'], [1, '#b2182b']];
    // A vertical colour key at the right edge, as the app's colorbar.
    function drawColorbar(ctx, x, y, h, stops, lo, hi, title) {
        const w = 12;
        for (let k = 0; k < h; k++) { ctx.fillStyle = lerpColor(stops, 1 - k / h); ctx.fillRect(x, y + k, w, 1.5); }
        ctx.strokeStyle = '#d1d5db'; ctx.lineWidth = 1; ctx.strokeRect(x + 0.5, y + 0.5, w, h);
        ctx.fillStyle = C.ink; ctx.font = `10px ${FONT}`; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        ctx.fillText(fmtTick(+hi.toFixed(1)), x + w + 4, y + 2);
        ctx.fillText(fmtTick(+lo.toFixed(1)), x + w + 4, y + h - 2);
        ctx.textBaseline = 'alphabetic';
        title.split('\n').forEach((t, i) => { ctx.font = i ? `10px ${FONT}` : `italic 10px ${FONT}`; ctx.fillText(t, x - 2, y - 16 + i * 12); });
    }

    // Scatter with an optional regression line and a legend for the groups.
    // A group may carry one colour per point (cols) with a colour key.
    function scatter(div, s) {
        const H = phone() ? 225 : 300;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const legendH = s.groups.length > 1 ? 22 : 0;
        const left = phone() ? 46 : 56;
        const right = s.colorbar ? 66 : 12;
        const ar = { x: left, y: top + 4, w: W - left - right, h: H - top - 4 - 46 - legendH };
        const shown = s.groups.filter(g => !g.hidden && g.xs.length);
        const allX = [].concat.apply([], s.groups.map(g => g.xs)), allY = [].concat.apply([], s.groups.map(g => g.ys));
        const xr = pad(allX), yr = pad(allY);
        const { sx, sy } = drawFrame(ctx, ar, xr, yr, { xLabel: s.xLabel, yLabel: s.yLabel, zeroY: true });
        for (const g of shown) {
            ctx.globalAlpha = g.alpha || 0.6;
            const r = g.size || 3;
            for (let i = 0; i < g.xs.length; i++) {
                ctx.fillStyle = g.cols ? g.cols[i] : g.color;
                ctx.beginPath(); ctx.arc(sx(g.xs[i]), sy(g.ys[i]), r, 0, Math.PI * 2); ctx.fill();
                if (g.cols) { ctx.strokeStyle = 'rgba(55,65,81,0.35)'; ctx.lineWidth = 0.5; ctx.stroke(); }
            }
        }
        ctx.globalAlpha = 1;
        if (s.regression && shown.length) {
            const xs = [].concat.apply([], shown.map(g => g.xs)), ys = [].concat.apply([], shown.map(g => g.ys));
            const n = xs.length; let a = 0, b = 0, aa = 0, ab = 0;
            for (let i = 0; i < n; i++) { a += xs[i]; b += ys[i]; aa += xs[i] * xs[i]; ab += xs[i] * ys[i]; }
            const den = n * aa - a * a;
            if (Math.abs(den) > 1e-9) {
                const slope = (n * ab - a * b) / den, ic = (b - slope * a) / n;
                ctx.save(); ctx.beginPath(); ctx.rect(ar.x, ar.y, ar.w, ar.h); ctx.clip();
                ctx.strokeStyle = C.green; ctx.lineWidth = 3; ctx.beginPath();
                ctx.moveTo(sx(xr[0]), sy(slope * xr[0] + ic)); ctx.lineTo(sx(xr[1]), sy(slope * xr[1] + ic)); ctx.stroke();
                ctx.restore();
            }
        }
        if (legendH) drawLegend(ctx, ar.x, H - 10, s.groups.map(g => ({ name: g.name, color: g.color, hidden: g.hidden })));
        if (s.colorbar) drawColorbar(ctx, ar.x + ar.w + 10, ar.y + 22, Math.min(120, ar.h - 30), s.colorbar.stops, s.colorbar.lo, s.colorbar.hi, s.colorbar.title);
    }

    // Horizontal box plots, one row per group, every point shown.
    function boxes(div, s) {
        const rows = s.rows;
        const rowH = phone() ? 22 : 30;
        const top = s.sub ? 44 : 30;
        const H = top + rows.length * rowH + 50;
        const { ctx, W } = setupCanvas(div, H);
        drawHeader(ctx, W, s.title, s.sub);
        ctx.font = `${phone() ? 10 : 11}px ${FONT}`;
        const labels = rows.map(r => `${r.name} (n=${r.vals.length})`);
        const lw = Math.max.apply(null, labels.map(l => ctx.measureText(l).width)) + 12;
        const ar = { x: lw, y: top, w: W - lw - 12, h: rows.length * rowH };
        const all = [].concat.apply([], rows.map(r => r.vals));
        const xr = pad(all);
        const { sx } = drawFrame(ctx, ar, xr, [0, 1], { xLabel: s.xLabel, rows: labels, zeroX: s.zero !== false });
        rows.forEach((r, i) => {
            const cy = ar.y + rowH * (i + 0.5);
            const sorted = r.vals.slice().sort((a, b) => a - b);
            const q1 = quantile(sorted, 0.25), q2 = quantile(sorted, 0.5), q3 = quantile(sorted, 0.75), iqr = q3 - q1;
            const lo = sorted.find(v => v >= q1 - 1.5 * iqr), hi = sorted.slice().reverse().find(v => v <= q3 + 1.5 * iqr);
            ctx.fillStyle = 'rgba(80,80,80,0.5)';
            for (const v of r.vals) { const jy = cy + (Math.random() - 0.5) * rowH * 0.5; ctx.beginPath(); ctx.arc(sx(v), jy, 2.2, 0, Math.PI * 2); ctx.fill(); }
            const bh = rowH * 0.55;
            ctx.strokeStyle = C.ink; ctx.lineWidth = 1.5;
            ctx.beginPath(); ctx.moveTo(sx(lo), cy); ctx.lineTo(sx(q1), cy); ctx.moveTo(sx(q3), cy); ctx.lineTo(sx(hi), cy); ctx.stroke();
            ctx.beginPath(); ctx.moveTo(sx(lo), cy - bh / 3); ctx.lineTo(sx(lo), cy + bh / 3); ctx.moveTo(sx(hi), cy - bh / 3); ctx.lineTo(sx(hi), cy + bh / 3); ctx.stroke();
            ctx.fillStyle = 'rgba(200,200,200,0.55)';
            ctx.fillRect(sx(q1), cy - bh / 2, Math.max(1, sx(q3) - sx(q1)), bh);
            ctx.strokeRect(sx(q1), cy - bh / 2, Math.max(1, sx(q3) - sx(q1)), bh);
            ctx.lineWidth = 2.5; ctx.beginPath(); ctx.moveTo(sx(q2), cy - bh / 2); ctx.lineTo(sx(q2), cy + bh / 2); ctx.stroke();
        });
    }

    // Dots in rows, jittered, one row per mutation state, median as a bar.
    function strip(div, s) {
        const rowH = phone() ? 34 : 46;
        const top = s.sub ? 44 : 30;
        const H = top + s.rows.length * rowH + 50;
        const { ctx, W } = setupCanvas(div, H);
        drawHeader(ctx, W, s.title, s.sub);
        ctx.font = `${phone() ? 10 : 11}px ${FONT}`;
        const labels = s.rows.map(r => `${r.name} (n=${r.vals.length})`);
        const lw = Math.max.apply(null, labels.map(l => ctx.measureText(l).width)) + 12;
        const ar = { x: lw, y: top, w: W - lw - 12, h: s.rows.length * rowH };
        const xr = pad([].concat.apply([], s.rows.map(r => r.vals)));
        const { sx } = drawFrame(ctx, ar, xr, [0, 1], { xLabel: s.xLabel, rows: labels });
        s.rows.forEach((r, i) => {
            const cy = ar.y + rowH * (i + 0.5);
            ctx.fillStyle = r.color; ctx.globalAlpha = 0.75;
            for (const v of r.vals) { ctx.beginPath(); ctx.arc(sx(v), cy + (Math.random() - 0.5) * rowH * 0.6, 2.6, 0, Math.PI * 2); ctx.fill(); }
            ctx.globalAlpha = 1;
            const m = median(r.vals);
            ctx.strokeStyle = r.color; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(sx(m), cy - rowH * 0.42); ctx.lineTo(sx(m), cy + rowH * 0.42); ctx.stroke();
        });
    }

    // Histogram of the panel. One value can be marked with a red line (a
    // cell line against the panel), or a subset drawn over it in red.
    function hist(div, s) {
        const H = phone() ? 195 : 250;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const legendH = s.over ? 20 : 0;
        const ar = { x: 24, y: top + 4, w: W - 36, h: H - top - 4 - 46 - legendH };
        const vals = s.vals.filter(v => isFinite(v));
        const lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
        const nb = 30, bw = (hi - lo) / nb || 1;
        const bin = (arr) => { const c = new Array(nb).fill(0); for (const v of arr) c[Math.min(nb - 1, Math.max(0, Math.floor((v - lo) / bw)))]++; return c; };
        const counts = bin(vals);
        const yr = [0, Math.max.apply(null, counts) * 1.05];
        const { sx, sy } = drawFrame(ctx, ar, [lo, hi], yr, { xLabel: s.xLabel, yTicks: false, zeroX: false, baseline: true, plain: true });
        const bars = (c, color) => { ctx.fillStyle = color; c.forEach((k, i) => { const x0 = sx(lo + i * bw), x1 = sx(lo + (i + 1) * bw); ctx.fillRect(x0 + 0.5, sy(k), Math.max(1, x1 - x0 - 1), ar.y + ar.h - sy(k)); }); };
        bars(counts, C.dot);
        if (s.over) bars(bin(s.over.vals.filter(v => isFinite(v))), s.over.color || C.red);
        if (typeof s.marker === 'number' && isFinite(s.marker)) {
            ctx.strokeStyle = C.red; ctx.lineWidth = 2.5; const x = sx(Math.min(hi, Math.max(lo, s.marker)));
            ctx.beginPath(); ctx.moveTo(x, ar.y); ctx.lineTo(x, ar.y + ar.h); ctx.stroke();
        }
        if (s.over) drawLegend(ctx, ar.x, H - 10, [{ name: s.allName || 'All cell lines', color: C.dot }, { name: s.over.name, color: s.over.color || C.red }]);
    }

    // A volcano as Inspect selection draws it: every gene grey, the ones
    // past both cutoffs in orange (more dependent in the selection) or
    // purple (less), the strongest few named.
    function volcano(div, s) {
        const H = phone() ? 240 : 320;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const left = phone() ? 46 : 56;
        const ar = { x: left, y: top + 4, w: W - left - 12, h: H - top - 4 - 46 };
        const pts = s.rows.map(r => ({ ...r, y: -Math.log10(Math.max(r.p, 1e-300)) }));
        const xr = pad(pts.map(p => p.d)), yr = [0, Math.max.apply(null, pts.map(p => p.y)) * 1.08];
        const { sx, sy } = drawFrame(ctx, ar, xr, yr, { xLabel: s.xLabel, yLabel: '-log10(p)', zeroY: false });
        const hit = (p) => Math.abs(p.d) >= s.dCut && p.p <= s.pCut;
        ctx.fillStyle = '#d1d5db'; ctx.globalAlpha = 0.55;
        for (const p of pts) if (!hit(p)) { ctx.beginPath(); ctx.arc(sx(p.d), sy(p.y), 1.8, 0, Math.PI * 2); ctx.fill(); }
        ctx.globalAlpha = 0.9;
        for (const p of pts) if (hit(p)) { ctx.fillStyle = p.d < 0 ? '#e66101' : '#5e3c99'; ctx.beginPath(); ctx.arc(sx(p.d), sy(p.y), 3.5, 0, Math.PI * 2); ctx.fill(); }
        ctx.globalAlpha = 1;
        const named = pts.filter(hit).sort((a, b) => a.p - b.p);
        const lab = [...named.filter(p => p.d < 0).slice(0, s.perSide), ...named.filter(p => p.d > 0).slice(0, s.perSide)];
        // Names pushed apart when two would overlap, each side on its own.
        ctx.font = `italic ${phone() ? 10 : 11}px ${FONT}`; ctx.fillStyle = C.ink; ctx.textBaseline = 'middle';
        for (const side of [-1, 1]) {
            const ps = lab.filter(p => Math.sign(p.d) === side).map(p => ({ p, x: sx(p.d), y: sy(p.y) - 6 })).sort((a, b) => a.y - b.y);
            for (let i = 1; i < ps.length; i++) if (ps[i].y - ps[i - 1].y < 12) ps[i].y = ps[i - 1].y + 12;
            for (const q of ps) {
                ctx.textAlign = side < 0 ? 'right' : 'left';
                ctx.fillText(q.p.gene, q.x + 8 * side, q.y);
                if (Math.abs(q.y - (sy(q.p.y) - 6)) > 1) { ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 0.8; ctx.beginPath(); ctx.moveTo(q.x + 4 * side, sy(q.p.y)); ctx.lineTo(q.x + 7 * side, q.y); ctx.stroke(); }
            }
        }
    }

    // Co-dependency of each member of a protein complex with one gene, as
    // signed bars: the Gene Browser's "does the complex act as one unit?".
    function rbars(div, s) {
        const rowH = phone() ? 22 : 26;
        const top = s.sub ? 44 : 30;
        const H = top + s.rows.length * rowH + 50;
        const { ctx, W } = setupCanvas(div, H);
        drawHeader(ctx, W, s.title, s.sub);
        ctx.font = `italic ${phone() ? 10 : 11}px ${FONT}`;
        const lw = Math.max.apply(null, s.rows.map(r => ctx.measureText(r.name).width)) + 14;
        const ar = { x: lw, y: top, w: W - lw - 48, h: s.rows.length * rowH };
        const { sx } = drawFrame(ctx, ar, [-1, 1], [0, 1], { xLabel: s.xLabel, rows: s.rows.map(r => r.name), zeroX: true });
        ctx.font = `italic ${phone() ? 10 : 11}px ${FONT}`; ctx.fillStyle = C.ink; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
        s.rows.forEach((r, i) => ctx.fillText(r.name, ar.x - 6, ar.y + rowH * (i + 0.5)));
        s.rows.forEach((r, i) => {
            const y = ar.y + rowH * i + rowH * 0.2, h = rowH * 0.6;
            if (!isFinite(r.r)) { ctx.fillStyle = C.dim; ctx.font = `10px ${FONT}`; ctx.textAlign = 'left'; ctx.fillText('no data', sx(0) + 6, y + h / 2); return; }
            const a = Math.min(1, Math.abs(r.r) / 0.6) * 0.75 + 0.15;
            ctx.fillStyle = r.r >= 0 ? `rgba(37,99,235,${a.toFixed(2)})` : `rgba(220,38,38,${a.toFixed(2)})`;
            const x0 = sx(Math.min(0, r.r)), x1 = sx(Math.max(0, r.r));
            ctx.fillRect(x0, y, Math.max(2, x1 - x0), h);
            ctx.fillStyle = C.ink; ctx.font = `bold 10px ${FONT}`; ctx.textAlign = r.r >= 0 ? 'left' : 'right';
            ctx.fillText(r.r.toFixed(2), r.r >= 0 ? x1 + 5 : x0 - 5, y + h / 2);
        });
    }

    // Horizontal bars, one per row.
    function hbars(div, s) {
        const rowH = phone() ? 17 : 24;
        const top = s.sub ? 44 : 30;
        const H = top + s.rows.length * rowH + 50;
        const { ctx, W } = setupCanvas(div, H);
        drawHeader(ctx, W, s.title, s.sub);
        ctx.font = `${phone() ? 10 : 11}px ${FONT}`;
        const lw = Math.max.apply(null, s.rows.map(r => ctx.measureText(r.name).width)) + 12;
        const ar = { x: lw, y: top, w: W - lw - 40, h: s.rows.length * rowH };
        const xr = [0, Math.max.apply(null, s.rows.map(r => r.value)) * 1.08];
        const { sx } = drawFrame(ctx, ar, xr, [0, 1], { xLabel: s.xLabel, rows: s.rows.map(r => r.name), zeroX: false });
        s.rows.forEach((r, i) => {
            const y = ar.y + rowH * i + rowH * 0.18;
            ctx.fillStyle = r.color || C.bar; ctx.fillRect(ar.x, y, sx(r.value) - ar.x, rowH * 0.64);
            ctx.fillStyle = C.ink; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
            ctx.fillText(r.text || num(r.value), sx(r.value) + 5, y + rowH * 0.32);
        });
    }

    // The gene set network as the app draws it: nodes on a ring, an edge
    // for every pair above the cutoff, blue positive and red negative,
    // wider the stronger. Nodes are green, or coloured by their mean gene
    // effect when asked, as Color by GE does.
    function network(div, s) {
        const H = phone() ? 240 : 340;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const keyW = s.colorNodes ? 70 : 0;
        const cx = (W - keyW) / 2, cy = top + (H - top - 44) / 2;
        const R = Math.min(W - keyW, H - top - 44) / 2 - (phone() ? 36 : 44);
        const nr = phone() ? 15 : 18;
        const pos = s.nodes.map((n, i) => { const a = -Math.PI / 2 + i * 2 * Math.PI / s.nodes.length; return { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a) }; });
        const idx = new Map(s.nodes.map((n, i) => [n.name, i]));
        for (const e of s.edges) {
            const a = pos[idx.get(e.a)], b = pos[idx.get(e.b)];
            if (!a || !b) continue;
            ctx.strokeStyle = e.r > 0 ? '#3182ce' : '#e53e3e'; ctx.lineWidth = 1 + 5 * Math.abs(e.r); ctx.globalAlpha = 0.85;
            ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
            ctx.globalAlpha = 1;
        }
        const vals = s.nodes.map(n => n.value).filter(isFinite);
        const lo = s.colorNodes ? Math.min.apply(null, vals) : 0, hi = s.colorNodes ? Math.max.apply(null, vals) : 0;
        s.nodes.forEach((n, i) => {
            const pnt = pos[i];
            ctx.fillStyle = s.colorNodes && isFinite(n.value) ? lerpColor(GE_STOPS, hi > lo ? (n.value - lo) / (hi - lo) : 0.5) : '#7ab950';
            ctx.strokeStyle = '#1f2937'; ctx.lineWidth = 1.5;
            ctx.beginPath(); ctx.arc(pnt.x, pnt.y, nr, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
            ctx.fillStyle = C.ink; ctx.font = `italic ${phone() ? 11 : 12}px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const outY = pnt.y + (pnt.y < cy - 1 ? -nr - 9 : nr + 9);
            ctx.fillText(n.name, pnt.x, Math.abs(pnt.y - cy) < 1 ? pnt.y + nr + 9 : outY);
        });
        if (s.colorNodes) drawColorbar(ctx, W - keyW + 18, top + 30, Math.min(120, H - top - 90), GE_STOPS, lo, hi, 'Mean\ngene effect');
        // Legend: the edges.
        const ly = H - 14;
        ctx.font = `${phone() ? 9 : 10}px ${FONT}`; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        let lx = 10;
        ctx.strokeStyle = '#3182ce'; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
        ctx.fillStyle = C.ink; ctx.fillText('r > 0', lx + 22, ly); lx += 22 + ctx.measureText('r > 0').width + 12;
        ctx.strokeStyle = '#e53e3e'; ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
        ctx.fillText('r < 0', lx + 22, ly); lx += 22 + ctx.measureText('r < 0').width + 16;
        ctx.fillText('wider = stronger', lx, ly);
        if (s.tag) {
            ctx.font = `bold 10px ${FONT}`; ctx.textAlign = 'right';
            ctx.fillStyle = s.tag.color; ctx.fillText(s.tag.text, W - 10, ly);
        }
    }

    // Colored grids come from the app's own drawer, the one the Matrix tab uses.
    function grid(div, spec) {
        const a = A();
        div.innerHTML = '';
        // Drawn into a child, because the drawer sizes its host to the grid.
        const inner = document.createElement('div');
        inner.style.margin = '0 auto';
        div.appendChild(inner);
        if (a && typeof a._drawCorrelationGrid === 'function') { a._drawCorrelationGrid(inner, spec); return; }
        div.innerHTML = '<div class="tour-loading">This chart is not available.</div>';
    }

    // ------------------------------------------------------------ actions
    // Open a pair's scatter with the hotspot overlay, and optionally the
    // hotspot filter, already set: the popout applies the preset itself
    // just before its first draw.
    function openPairWithHotspot(a, g1, g2, mutGene, filterLevel) {
        a._inspectPreset = { hotspotGene: mutGene, hotspotMode: 'color', filterGene: filterLevel ? mutGene : null, filterLevel: filterLevel || null };
        a.openInspectByGenes(g1, g2);
    }
    const P53_SET = ['TP53', 'MDM2', 'MDM4', 'CDKN1A', 'PPM1D', 'USP7'];
    const setText = () => P53_SET.filter(g => A().geneIndex.has(g)).join('\n');
    // Run the example set in the app, on gene effect or mRNA; optionally
    // land on one result tab or tick Color by GE once the results exist.
    function runExampleSet(o) {
        o = o || {};
        document.getElementById('modeGeneSetBtn')?.click();
        const basis = document.querySelector(`input[name="analysisBasis"][value="${o.basis || 'ge'}"]`);
        if (basis && !basis.checked) { basis.checked = true; basis.dispatchEvent(new Event('change', { bubbles: true })); }
        const ta = document.getElementById('geneTextarea');
        if (ta) { ta.value = setText(); ta.dispatchEvent(new Event('input', { bubbles: true })); }
        const before = A().results;
        setTimeout(() => document.getElementById('runAnalysis')?.click(), 150);
        if (!o.tab && !o.colorGE) return;
        let tries = 0;
        const arm = () => {
            if (A().results && A().results !== before && A().results.success) {
                if (o.tab) document.querySelector(`.nav-link[data-tab="${o.tab}"]`)?.click();
                if (o.colorGE) { const cb = document.getElementById('colorByGeneEffect'); if (cb && !cb.checked) { cb.checked = true; cb.dispatchEvent(new Event('change', { bubbles: true })); } }
                return;
            }
            if (tries++ < 80) setTimeout(arm, 250);
        };
        setTimeout(arm, 600);
    }
    // The Cell Line Browser filtered to a tissue and a hotspot gene.
    function openBrowserFiltered(a, tissue, hotspotGene) {
        a.openCellLineBrowser();
        setTimeout(() => {
            const set = (id, v) => { const el = document.getElementById(id); if (!el) return; el.value = v; el.dispatchEvent(new Event('change', { bubbles: true })); };
            if (tissue) set('clbTissueFilter', tissue);
            if (hotspotGene) set('clbHotspotFilter', hotspotGene);
        }, 400);
    }
    // Skin cell lines split by BRAF hotspot status, as CRISPR indices.
    function skinByBraf() {
        const ids = D.ids(), sel = [], oth = [];
        for (let i = 0; i < ids.length; i++) {
            if (D.lineage(ids[i]) !== 'Skin') continue;
            (D.hotspot('BRAF', ids[i]) >= 1 ? sel : oth).push(i);
        }
        return { sel, oth };
    }

    // ---------------------------------------------------------------- pages
    // Each page: a section, a title, a lead line (why it matters), a chart
    // drawn from the data, two or three plain sentences, and the action
    // that opens the same thing for real.
    const SECTIONS = ['The network', 'What a link is', 'Compare groups', 'One gene', 'Cell lines', 'Take it with you'];
    const PAGES = [
        // ---- 1. The network
        {
            section: 0, title: 'A gene set becomes a network',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4,
            lead: () => `Paste a gene set. Genes whose CRISPR knockout effects rise and fall together across ${num(D.ids().length)} cancer cell lines are linked.`,
            body: () => `<div class="tour-tiles">
                    <div class="tour-tile"><b>${num(D.ids().length)}</b><span>cell lines with CRISPR gene effect</span></div>
                    <div class="tour-tile"><b>${num(D.nExprLines() || 1718)}</b><span>cell lines with mRNA expression</span></div>
                    <div class="tour-tile"><b>${num(D.nGenes())}</b><span>genes, from DepMap and other resources</span></div>
                </div>
                <p>Here the p53 pathway: every link is a correlation at the default cutoff. Set the cutoff, drag the nodes, click a link for the plot behind it.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const rows = genes.map(g => D.row(g));
                network(div, { title: 'Gene set network, p53 pathway', sub: 'links at |r| of 0.5 and above, the default cutoff', nodes: genes.map(g => ({ name: g })), edges: setEdges(genes, rows, 0.5),
                    tag: { text: 'Built on gene effect (CRISPR)', color: '#3f6b25' } });
            },
            action: { label: 'Run this gene set', run: () => runExampleSet() }
        },
        {
            section: 0, title: 'Colour the network',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4,
            lead: () => 'Colour the nodes by gene effect, by mRNA level, or by your own numbers.',
            body: () => `<p>Paste a gene list <b>With Stats</b>, for example the log fold change and FDR from your own CRISPR screen or RNA-seq, and the network is coloured by them. Enrichr can highlight a pathway on top.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const rows = genes.map(g => D.row(g));
                network(div, { title: 'The same network, coloured by mean gene effect', sub: 'orange = the cell lines depend on it, purple = dispensable', colorNodes: true,
                    nodes: genes.map((g, i) => ({ name: g, value: meanOf(rows[i]) })), edges: setEdges(genes, rows, 0.5) });
            },
            action: { label: 'Run it with Color by GE', run: () => runExampleSet({ colorGE: true }) }
        },
        {
            section: 0, title: 'Gene effect or mRNA',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4 && !!A().loadExpressionData,
            lead: () => 'The same set, built on mRNA expression instead of CRISPR gene effect.',
            body: () => `<p>Every network, matrix and plot can be built on either layer, on the same cell lines, so the two are directly comparable. A tag on the network says which one you are looking at.</p>`,
            plot: async (div) => {
                div.innerHTML = '<div class="tour-loading">Loading expression data...</div>';
                if (!await D.expression()) { div.innerHTML = '<div class="tour-loading">Expression data is not available right now.</div>'; return; }
                const genes = P53_SET.filter(g => D.row(g) && D.exprRow(g));
                const rows = genes.map(g => D.exprRow(g));
                network(div, { title: 'p53 pathway, built on mRNA', sub: 'links at |r| of 0.5 and above', nodes: genes.map(g => ({ name: g })), edges: setEdges(genes, rows, 0.5),
                    tag: { text: 'Built on mRNA expression', color: '#6d28d9' } });
            },
            action: { label: 'Run this set on mRNA', run: () => runExampleSet({ basis: 'expr' }) }
        },
        // ---- 2. What a link is
        {
            section: 1, title: 'A link is a correlation',
            need: () => !!D.row('TP53') && !!D.row('MDM2'),
            lead: () => 'Behind every link is a pair of genes across all cell lines, one dot per cell line.',
            body: () => `<p>${gi('TP53')} and ${gi('MDM2')} correlate at r = ${(() => { const st = pearson(D.row('TP53'), D.row('MDM2')); return st ? st.r.toFixed(2) : '?'; })()}. Cell lines that need one need the other. Click any link in a network to open its plot; type any two genes to make your own.</p>`,
            plot: (div) => {
                const { xs, ys } = scatterXY('TP53', 'MDM2');
                scatter(div, { title: 'TP53 vs MDM2', sub: `n = ${num(xs.length)} cell lines`, xLabel: 'TP53 Gene Effect', yLabel: 'MDM2 Gene Effect',
                    groups: [{ xs, ys, color: C.dot, alpha: 0.6, size: 3 }], regression: true });
            },
            action: { label: 'Open TP53 vs MDM2', run: (a) => a.openInspectByGenes('TP53', 'MDM2') }
        },
        {
            section: 1, title: 'Colour by mutation, filter by cancer type',
            need: () => !!D.row('TP53') && !!D.row('MDM2') && !!A().mutations?.geneData?.TP53,
            lead: () => 'Colour the dots by a mutation. Wild-type TP53 cell lines depend on MDM2; mutated ones do not.',
            body: () => `<p>Colour by tissue, subtype, mutation, fusion or copy number. Filter to one cancer type, or to a combination of mutations, and the correlation is recomputed on those cell lines alone. A correlation that holds within one tissue is the one to trust.</p>`,
            plot: (div) => {
                const g = scatterByMutation('TP53', 'MDM2', 'TP53');
                const spec = { title: 'TP53 vs MDM2', sub: 'coloured by TP53 hotspot mutation', xLabel: 'TP53 Gene Effect', yLabel: 'MDM2 Gene Effect',
                    groups: [
                        { name: `WT (n=${g.wt.xs.length})`, xs: g.wt.xs, ys: g.wt.ys, color: C.dot, alpha: 0.6, size: 3 },
                        { name: `1 mut (n=${g.m1.xs.length})`, xs: g.m1.xs, ys: g.m1.ys, color: C.blue, alpha: 0.8, size: 3.2 },
                        { name: `2 mut (n=${g.m2.xs.length})`, xs: g.m2.xs, ys: g.m2.ys, color: C.red, alpha: 0.8, size: 3.2 }
                    ] };
                div._spec = spec;
                scatter(div, spec);
            },
            choices: [
                { label: 'Both groups', run: (div) => { div._spec.groups.forEach(g => g.hidden = false); scatter(div, div._spec); } },
                { label: 'Wild-type only', run: (div) => { div._spec.groups.forEach((g, i) => g.hidden = i !== 0); scatter(div, div._spec); } },
                { label: 'Mutated only', run: (div) => { div._spec.groups.forEach((g, i) => g.hidden = i === 0); scatter(div, div._spec); } }
            ],
            action: { label: 'Open it coloured by TP53', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', null) }
        },
        {
            section: 1, title: 'Which gene explains the dots?',
            need: () => !!D.row('SOX10') && !!D.row('MITF') && !!A().loadExpressionData,
            lead: () => 'Colour the plot by any gene. Or press Suggest, and every gene is ranked by how well it explains where the dots sit.',
            body: () => `<p>${gi('SOX10')} and ${gi('MITF')} are co-dependencies. Coloured by ${gi('TYR')} mRNA, the dependent cell lines light up red: they are the melanocyte lineage. Suggest found ${gi('TYR')} on its own, at the top of its list.</p>`,
            plot: async (div) => {
                div.innerHTML = '<div class="tour-loading">Loading expression data...</div>';
                if (!await D.expression()) { div.innerHTML = '<div class="tour-loading">Expression data is not available right now.</div>'; return; }
                const { xs, ys, idx } = scatterXY('SOX10', 'MITF');
                const v = idx.map(i => D.expr('TYR', i));
                const fin = v.filter(isFinite), lo = Math.min.apply(null, fin), hi = Math.max.apply(null, fin);
                const order = v.map((x, i) => i).sort((a, b) => (isFinite(v[a]) ? v[a] : -1e9) - (isFinite(v[b]) ? v[b] : -1e9));
                const cols = order.map(i => isFinite(v[i]) ? lerpColor(EXPR_STOPS, (v[i] - lo) / (hi - lo || 1)) : '#d1d5db');
                scatter(div, { title: 'SOX10 vs MITF, coloured by TYR mRNA', sub: `n = ${num(xs.length)} cell lines`, xLabel: 'SOX10 Gene Effect', yLabel: 'MITF Gene Effect',
                    groups: [{ xs: order.map(i => xs[i]), ys: order.map(i => ys[i]), cols, alpha: 0.85, size: 3.2 }],
                    colorbar: { stops: EXPR_STOPS, lo, hi, title: 'TYR\nmRNA' } });
            },
            action: {
                label: 'Open it coloured by TYR',
                run: async (a) => {
                    await D.expression();
                    a._inspectPreset = { colorByGene: 'TYR', colorByGeneType: 'expr' };
                    a.openInspectByGenes('SOX10', 'MITF');
                }
            }
        },
        // ---- 3. Compare groups
        {
            section: 2, title: 'Mutation analysis',
            need: () => !!D.row('BRAF') && !!A().mutations?.geneData?.BRAF,
            lead: () => 'Pick a cancer gene. Every gene is tested for a different dependency in mutated versus wild-type cell lines.',
            body: () => `<p>Here ${gi('BRAF')} by ${gi('BRAF')} mutation status: the mutated cell lines depend on it, the wild-type ones do not. The table ranks all ${num(D.nGenes())} genes, and Inspect opens any of them. Filter to one cancer type first to keep lineage out of the comparison.</p>`,
            plot: (div) => {
                const ids = D.ids(), row = D.row('BRAF');
                const wt = [], m1 = [], m2 = [];
                for (let i = 0; i < ids.length; i++) {
                    const v = D.ge(row, i); if (!isFinite(v)) continue;
                    const lvl = D.hotspot('BRAF', ids[i]);
                    (lvl >= 2 ? m2 : lvl === 1 ? m1 : wt).push(v);
                }
                strip(div, { title: 'BRAF Gene Effect by BRAF mutation status', xLabel: 'BRAF Gene Effect',
                    rows: [{ name: '2 mut', vals: m2, color: C.red }, { name: '1 mut', vals: m1, color: C.blue }, { name: 'WT', vals: wt, color: '#888888' }] });
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
            section: 2, title: 'Compare any two groups',
            need: () => !!D.row('BRAF') && !!A().mutations?.geneData?.BRAF && typeof A().tDistributionPValue === 'function' && skinByBraf().sel.length >= 5,
            lead: () => 'Any group of cell lines you can filter, against the rest: here BRAF-mutant melanoma against the other skin cell lines.',
            body: () => `<p>The volcano names what the mutated cell lines need more (${gi('BRAF')}, ${gi('MAP2K1')}) and less. The same comparison runs on mRNA and on the mutation spectrum, and any two gates drawn on a plot can be compared the same way.</p>`,
            plot: (div) => {
                const { sel, oth } = skinByBraf();
                const rows = welchAll(sel, oth);
                volcano(div, { title: `${sel.length} BRAF-mutant skin cell lines vs the other ${oth.length}`, sub: 'Welch t-test of gene effect, every gene', xLabel: 'Δ gene effect (BRAF-mutant minus the rest)', rows, dCut: 0.3, pCut: 0.05, perSide: 4 });
            },
            action: {
                label: 'Open this comparison',
                run: (a) => {
                    openBrowserFiltered(a, 'Skin', 'BRAF');
                    setTimeout(() => { try { a._resetGEInspectScope?.(); a._geInspectScope = 'lineage'; a._runBusy({ button: 'clbInspectGEBtn', label: 'Comparing' }, () => a.inspectSelectionGE()); } catch (e) { console.warn(e); } }, 1400);
                }
            }
        },
        // ---- 4. One gene
        {
            section: 3, title: 'The Gene Browser',
            need: () => !!D.row('KRAS'),
            lead: () => 'One gene, one page: dependency, expression, mutations, drugs, complexes, pathways, structure and sequence.',
            body: () => `<p>Each number is ranked against every other gene, so the page says in one line whether a gene is essential, selective or rarely needed, and whether that follows cancer type. The whole page exports as a PDF of several pages, or as one file for an AI.</p>`,
            plot: (div) => {
                const a = A(), row = D.row('KRAS');
                const vals = []; for (let i = 0; i < D.ids().length; i++) { const v = D.ge(row, i); if (isFinite(v)) vals.push(v); }
                const below = vals.filter(v => v < -0.5).length;
                let rank = '';
                try { const st = a._gbGeneStats?.('ge'); if (st) rank = `, more essential than ${Math.round((1 - a._gbRankIn(st.means, vals.reduce((s, v) => s + v, 0) / vals.length)) * 100)}% of genes`; } catch (e) { }
                const under = vals.filter(v => v < -0.5), over = vals.filter(v => v >= -0.5);
                hist(div, { title: 'KRAS gene effect across the panel', sub: `${Math.round(below / vals.length * 100)}% of cell lines below -0.5${rank}`, xLabel: 'KRAS Gene Effect', vals: over.concat(under), over: { vals: under, name: 'Dependent (below -0.5)', color: '#dc2626' }, allName: 'Not dependent' });
            },
            action: { label: 'Open KRAS in the Gene Browser', run: (a) => a.openGeneBrowser('KRAS') }
        },
        {
            section: 3, title: 'Does the complex act as one unit?',
            need: () => !!D.row('TSC1') && !!D.row('TSC2') && typeof A()._loadComplexes === 'function',
            lead: () => 'Every named protein complex the gene belongs to, with each member\'s co-dependency.',
            body: () => `<p>${gi('TSC2')} follows ${gi('TSC1')} almost perfectly across the cell lines; ${gi('TBC1D7')}, the third subunit, far less. The page does this for every complex in EBI Complex Portal and every Reactome pathway, so a member that does not follow stands out at once.</p>`,
            plot: async (div) => {
                div.innerHTML = '<div class="tour-loading">Loading protein complexes...</div>';
                const d = await A()._loadComplexes();
                const cpx = (d?.complexes || []).filter(c => c.genes.includes('TSC1'));
                const members = cpx.length ? cpx[0].genes : ['TSC1', 'TSC2', 'TBC1D7'];
                const me = D.row('TSC1');
                const rows = members.filter(g => g !== 'TSC1').map(g => ({ name: g, r: D.row(g) ? (pearson(me, D.row(g))?.r ?? NaN) : NaN })).sort((a, b) => (b.r || -9) - (a.r || -9));
                rbars(div, { title: cpx.length ? cpx[0].name : 'TSC1-TSC2 complex', sub: 'co-dependency of each member with TSC1 (r of gene effect)', xLabel: 'r with TSC1', rows });
            },
            action: { label: 'Open TSC1 in the Gene Browser', run: (a) => a.openGeneBrowser('TSC1') }
        },
        {
            section: 3, title: 'Co-expression, with neighbours flagged',
            need: () => !!A().loadExpressionData && !!D.row('ERBB2'),
            lead: () => 'The genes whose mRNA follows the gene\'s, against each other, and a warning when they sit next to it on the chromosome.',
            body: () => `<p>${gi('GRB7')} and ${gi('PGAP3')} correlate with ${gi('ERBB2')} because they share its amplicon on 17q12, not because they share its function. The Gene Browser marks such neighbours in amber, and the heatmap shows which genes form a real module.</p>`,
            plot: async (div) => {
                div.innerHTML = '<div class="tour-loading">Loading expression data...</div>';
                if (!await D.expression()) { div.innerHTML = '<div class="tour-loading">Expression data is not available right now.</div>'; return; }
                const genes = ['ERBB2', 'GRB7', 'PGAP3', 'NHSL3', 'DDR1', 'CTNND1', 'SSH3', 'KRT8'].filter(g => D.exprRow(g));
                const rows = genes.map(g => D.exprRow(g));
                const z = genes.map((g, i) => genes.map((h, j) => i === j ? 1 : (pearson(rows[i], rows[j])?.r ?? null)));
                grid(div, { title: 'Co-expression of ERBB2 and its partners', sub: 'mRNA r for every pair; GRB7 and PGAP3 are ERBB2\'s neighbours on 17q12', rowLabels: genes, colLabels: genes, z, zmin: -1, zmax: 1,
                    colorscale: EXPR_STOPS, showValues: true, colorbarTitle: 'r', cell: phone() ? 28 : 40, italicLabels: true });
            },
            action: { label: 'Open ERBB2 in the Gene Browser', run: (a) => a.openGeneBrowser('ERBB2') }
        },
        // ---- 5. Cell lines
        {
            section: 4, title: 'The Cell Line Browser',
            need: () => Array.isArray(A().drugResponse?.compounds) && A().drugResponse.compounds.some(c => /vemurafenib/i.test(c.name || '')) && !!A().mutations?.geneData?.BRAF,
            lead: () => 'Find cell lines by tissue, disease, mutation combinations, fusions, copy number, sex or drug response.',
            body: () => `<p>Here every cell line\'s response to vemurafenib, with the ${gi('BRAF')}-mutant ones in red. The browser lists what matches, sends the selection to a heatmap, a comparison or an export, and downloads every gene\'s expression, gene effect and mutations for each cell line.</p>`,
            plot: (div) => {
                const cp = A().drugResponse.compounds.find(c => /vemurafenib/i.test(c.name || ''));
                const all = [], mut = [];
                D.ids().forEach(id => { const v = cp.auc?.[id]; if (typeof v === 'number') { all.push(v); if (D.hotspot('BRAF', id) >= 1) mut.push(v); } });
                hist(div, { title: 'Vemurafenib response across the panel', sub: `${num(all.length)} cell lines tested in PRISM`, xLabel: 'AUC (0 = all cells killed, 1 = no effect)', vals: all, over: { vals: mut, name: `BRAF-mutant (n=${mut.length})` } });
            },
            action: { label: 'Open the browser on BRAF-mutant skin', run: (a) => openBrowserFiltered(a, 'Skin', 'BRAF') }
        },
        {
            section: 4, title: 'Every cell line has its own page',
            need: () => !!D.lineId('A375') && !!A().globalSignatures?.byCellLine,
            lead: () => 'Origin, drivers, dependencies, expression, copy number, drug response and authentication, gathered from DepMap, Cellosaurus and COSMIC.',
            body: () => `<p>Here A375\'s ploidy against the panel. Each page is a document of ten or more pages that exports as a PDF, as an image, or as one file for an AI.</p>`,
            plot: (div) => {
                const id = D.lineId('A375'), sig = A().globalSignatures.byCellLine;
                const vals = []; D.ids().forEach(x => { const v = sig[x]?.Ploidy; if (typeof v === 'number') vals.push(v); });
                hist(div, { title: 'Ploidy across the panel', sub: 'A375 marked in red', xLabel: 'Ploidy', vals, marker: sig[id]?.Ploidy });
            },
            action: { label: 'Open the A375 page', run: (a) => a.openCellLineWiki(D.lineId('A375')) }
        },
        // ---- 6. Take it with you
        {
            section: 5, title: 'Take the results with you',
            lead: () => 'Every view exports the same way: a figure for a slide, a table for a spreadsheet, a methods paragraph, or a file to discuss with an AI.',
            body: () => `<div class="tour-cards">
                    <div class="tour-card"><b>Export image</b><span>PNG, SVG, TIFF, PDF or PowerPoint at print size. Every image carries its settings: open it again in Correlate and the analysis comes back.</span></div>
                    <div class="tour-card"><b>Export for AI</b><span>One file with the data behind the view, written for a language model. Upload it and ask: which of these dependencies has a drug? Which cell line should I use?</span></div>
                    <div class="tour-card"><b>Methods</b><span>How every number was calculated, as a paragraph for a methods section and as a plain walk-through.</span></div>
                    <div class="tour-card"><b>Tables</b><span>CSV of every result, and of every gene in a cell line or every cell line for a gene.</span></div>
                </div>`,
            actions: [
                { label: 'Run the p53 example, then try the exports', run: () => runExampleSet() }
            ]
        },
        {
            section: 5, title: 'Where to go next',
            lead: () => 'More tools, one line each.',
            body: () => `<ul class="tour-more">
                    <li><b>Heatmap</b>: a gene set across the cell lines, clustered and grouped.</li>
                    <li><b>Compilation</b>: a grid of plots, one per gene, sharing one set of filters.</li>
                    <li><b>Threshold scan</b>: how many cell lines pass a cutoff, and which.</li>
                    <li><b>Gates</b>: draw a box on a plot and carry those cell lines into the next analysis.</li>
                    <li><b>PCA and UMAP</b>: the cell lines laid out by a gene class, with the genes that drive the axes.</li>
                    <li><b>Green Listed</b>: design the CRISPR library for the genes you found.</li>
                </ul>`,
            actions: [
                { label: 'Run the p53 example', run: () => runExampleSet() },
                { label: 'Open the Gene Browser', run: (a) => a.openGeneBrowser('KRAS') },
                { label: 'Open the Cell Line Browser', run: (a) => a.openCellLineBrowser() },
                { label: 'The statistics in depth', run: () => { const m = document.getElementById('infographicModal'); if (m) m.style.display = 'flex'; } }
            ]
        }
    ];

    // ------------------------------------------------------------------ UI
    const CSS = `
#tourModal .modal { max-width: 860px; }
#tourModal .modal-header { background: linear-gradient(90deg, #f3faee, #fff); }
#tourModal .modal-body { overflow-y: auto; padding: 12px 24px 16px; }
#tourModal .tour-sections { display: flex; gap: 4px; flex-wrap: wrap; margin: 0 0 12px; }
#tourModal .tour-sec { flex: 1 1 0; min-width: 0; border: none; border-radius: 4px; padding: 5px 4px; font-size: 11px; font-weight: 600; color: #6b7280; background: #f3f4f6; cursor: pointer; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
#tourModal .tour-sec.done { background: #e4f2d9; color: #3f6b25; }
#tourModal .tour-sec.on { background: var(--green-600); color: #fff; }
#tourModal .tour-sec:hover:not(.on) { background: #e5e7eb; }
#tourModal .tour-sec.done:hover:not(.on) { background: #d4ebc4; }
#tourModal .tour-lead { font-size: 16px; line-height: 1.45; font-weight: 600; color: #1f2937; margin: 0 0 10px; }
#tourModal .tour-text { font-size: 13px; line-height: 1.6; color: #374151; }
#tourModal .tour-text p { margin: 0 0 10px; }
#tourModal .tour-text ul { margin: 0 0 10px 18px; padding: 0; }
#tourModal .tour-text li { margin-bottom: 4px; }
#tourModal .tour-chart { min-height: 120px; margin: 4px 0 12px; border: 1px solid var(--gray-200); border-radius: 6px; overflow-x: auto; overflow-y: hidden; padding: 6px 0; background: #fff; }
#tourModal .tour-loading { padding: 40px 12px; text-align: center; color: #6b7280; font-size: 12px; }
#tourModal .tour-choices { display: flex; flex-wrap: nowrap; gap: 6px; margin: -4px 0 10px; }
#tourModal .tour-choices:empty { display: none; }
#tourModal .tour-choice { flex: 1 1 0; min-width: 0; min-height: 34px; font-size: 12px; padding: 4px 4px; white-space: nowrap; }
#tourModal .tour-choice.on { background: var(--green-600); color: #fff; border-color: var(--green-600); }
#tourModal .tour-actions { display: flex; flex-wrap: wrap; gap: 8px; margin: 4px 0 6px; align-items: center; }
#tourModal .tour-actions .btn { min-height: 40px; }
#tourModal .tour-actions .btn.tour-primary { background: var(--green-600); color: #fff; border-color: var(--green-600); font-weight: 600; padding-left: 16px; padding-right: 16px; }
#tourModal .tour-actions .btn.tour-primary:hover { background: var(--green-700); }
#tourModal .tour-tiles { display: flex; gap: 10px; margin: 0 0 12px; }
#tourModal .tour-tile { flex: 1 1 0; background: #f3faee; border: 1px solid #cfe6bb; border-radius: 6px; padding: 10px 12px; }
#tourModal .tour-tile b { display: block; font-size: 22px; line-height: 1.1; color: #3f6b25; }
#tourModal .tour-tile span { display: block; font-size: 11px; color: #4b5563; margin-top: 4px; line-height: 1.3; }
#tourModal .tour-cards { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin: 0 0 12px; }
#tourModal .tour-card { background: #fff; border: 1px solid var(--gray-200); border-left: 3px solid var(--green-500); border-radius: 6px; padding: 10px 12px; }
#tourModal .tour-card b { display: block; font-size: 13px; color: #1f2937; margin-bottom: 3px; }
#tourModal .tour-card span { font-size: 12px; color: #4b5563; line-height: 1.5; }
#tourModal .tour-more { list-style: none; margin: 0 0 12px !important; padding: 0; display: grid; grid-template-columns: 1fr 1fr; gap: 4px 16px; }
#tourModal .tour-more li { font-size: 12.5px; }
#tourModal .modal-footer { display: flex; align-items: center; justify-content: space-between; gap: 10px; padding: 10px 24px; border-top: 1px solid var(--gray-200); flex-shrink: 0; }
#tourModal .tour-progress { flex: 1; height: 6px; background: #e5e7eb; border-radius: 3px; overflow: hidden; }
#tourModal .tour-progress > div { height: 100%; background: var(--green-600); border-radius: 3px; transition: width .25s; }
#tourModal .tour-step { font-size: 11px; color: #6b7280; white-space: nowrap; }
#tourModal .tour-nav { min-width: 84px; min-height: 40px; }
@media (max-width: 640px) {
  #tourModal .modal-body { padding: 10px 14px 12px; }
  #tourModal .tour-sections { gap: 3px; }
  #tourModal .tour-sec { font-size: 10px; padding: 4px 2px; }
  #tourModal .tour-lead { font-size: 15px; }
  #tourModal .tour-text { font-size: 14px; line-height: 1.45; }
  #tourModal .tour-text p { margin: 0 0 6px; }
  #tourModal .tour-chart { margin: 2px 0 8px; }
  #tourModal .tour-tiles { gap: 6px; }
  #tourModal .tour-tile { padding: 8px; }
  #tourModal .tour-tile b { font-size: 18px; }
  #tourModal .tour-cards, #tourModal .tour-more { grid-template-columns: 1fr; }
  #tourModal .tour-nav { min-width: 72px; min-height: 44px; font-size: 14px; }
  #tourModal .tour-actions .btn { flex: 1 1 100%; min-height: 44px; font-size: 14px; }
  #tourModal .tour-choice { font-size: 13px; min-height: 40px; }
}`;

    let pages = [], step = 0, built = false, resizeTimer = null;

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
                    <div class="tour-sections" id="tourSections"></div>
                    <div class="tour-lead" id="tourLead"></div>
                    <div class="tour-chart" id="tourChart" style="display:none;"></div>
                    <div class="tour-choices" id="tourChoices"></div>
                    <div class="tour-text" id="tourText"></div>
                    <div class="tour-actions" id="tourActions"></div>
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-outline btn-sm tour-nav" id="tourBack">Back</button>
                    <div class="tour-progress"><div id="tourBar"></div></div>
                    <span class="tour-step" id="tourStep"></span>
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
        // A chart is drawn at the card's width, so a turned phone gets a
        // fresh drawing rather than a stretched one.
        window.addEventListener('resize', () => {
            if (!isOpen()) return;
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(drawChart, 200);
        });
    }

    const isOpen = () => document.getElementById('tourModal')?.classList.contains('active');
    function saveStep(i) { try { localStorage.setItem(STEP_KEY, String(i)); } catch (e) { } }
    function loadStep() { try { return parseInt(localStorage.getItem(STEP_KEY), 10) || 0; } catch (e) { return 0; } }

    function actionButton(act, primary) {
        const b = document.createElement('button');
        b.type = 'button';
        b.className = 'btn btn-sm ' + (primary ? 'tour-primary' : 'btn-outline');
        b.textContent = act.label;
        b.addEventListener('click', () => {
            close();
            try { const r = act.run(A()); if (r && r.catch) r.catch(e => console.warn('Tour action failed:', e)); } catch (e) { console.warn('Tour action failed:', e); }
        });
        return b;
    }

    function drawChart() {
        const p = pages[step];
        const chart = document.getElementById('tourChart');
        if (!p || !p.plot) { chart.style.display = 'none'; chart.innerHTML = ''; return; }
        chart.style.display = 'block';
        const mine = step;
        Promise.resolve().then(() => { if (mine === step) return p.plot(chart); })
            .catch(e => { console.warn('Tour chart failed:', e); chart.innerHTML = '<div class="tour-loading">This chart could not be drawn.</div>'; });
    }

    function show(i) {
        step = Math.max(0, Math.min(pages.length - 1, i));
        saveStep(step);
        const p = pages[step];
        document.getElementById('tourTitle').textContent = p.title;
        document.getElementById('tourLead').innerHTML = p.lead ? p.lead() : '';
        document.getElementById('tourText').innerHTML = p.body();
        // Section pills: the current section, those passed, and a jump to each.
        const secs = document.getElementById('tourSections');
        secs.innerHTML = '';
        const firstOf = (k) => pages.findIndex(pg => pg.section === k);
        SECTIONS.forEach((name, k) => {
            const first = firstOf(k);
            if (first < 0) return;
            const b = document.createElement('button');
            b.type = 'button';
            b.className = 'tour-sec' + (p.section === k ? ' on' : first < step ? ' done' : '');
            b.textContent = `${k + 1}. ${name}`;
            b.title = name;
            b.addEventListener('click', () => show(first));
            secs.appendChild(b);
        });
        const ch = document.getElementById('tourChoices');
        ch.innerHTML = '';
        if (p.choices) {
            p.choices.forEach((c, k) => {
                const b = document.createElement('button');
                b.type = 'button';
                b.className = 'btn btn-outline btn-sm tour-choice' + (k === 0 ? ' on' : '');
                b.textContent = c.label;
                b.addEventListener('click', () => {
                    ch.querySelectorAll('.tour-choice').forEach(x => x.classList.remove('on'));
                    b.classList.add('on');
                    const div = document.getElementById('tourChart');
                    if (div && div._spec) { try { c.run(div); } catch (e) { } }
                });
                ch.appendChild(b);
            });
        }
        const acts = document.getElementById('tourActions');
        acts.innerHTML = '';
        (p.actions || (p.action ? [p.action] : [])).forEach((a, k) => acts.appendChild(actionButton(a, k === 0)));
        document.getElementById('tourBack').style.visibility = step === 0 ? 'hidden' : 'visible';
        document.getElementById('tourNext').textContent = step === pages.length - 1 ? 'Finish' : 'Next';
        document.getElementById('tourBar').style.width = `${Math.round((step + 1) / pages.length * 100)}%`;
        document.getElementById('tourStep').textContent = `${step + 1} of ${pages.length}`;
        document.querySelector('#tourModal .modal-body').scrollTop = 0;
        drawChart();
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
        document.getElementById('tourModal')?.classList.remove('active');
    }

    window.CorrelateTour = { open, close, isOpen };
})();
