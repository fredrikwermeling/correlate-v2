// A guided tour of Correlate: one page per feature, each with a real chart
// drawn from the loaded data and a button that opens that view for real.
// Charts are drawn on canvas here, sized exactly to the card, so a phone
// gets the same picture as a desktop. Reads window.app only. No network
// calls; the current page is kept in localStorage so the tour reopens
// where it was left.
(function () {
    'use strict';

    const A = () => window.app || null;
    const STEP_KEY = 'correlateTourStep';
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const phone = () => window.innerWidth <= 640;

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
        const xs = [], ys = [];
        for (let i = 0; i < D.ids().length; i++) {
            const x = D.ge(r1, i), y = D.ge(r2, i);
            if (isFinite(x) && isFinite(y)) { xs.push(x); ys.push(y); }
        }
        return { xs, ys };
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

    // ------------------------------------------------------- canvas charts
    // The app's look: Arial, black zero lines, a light plot area, green
    // regression line, grey dots, blue and red for one and two mutated copies.
    const FONT = 'Arial, Helvetica, sans-serif';
    const C = { ink: '#374151', dim: '#6b7280', grid: '#e5e7eb', plotBg: '#fafafa', dot: '#9ca3af',
        green: '#6ba544', blue: '#3b82f6', red: '#dc2626', bar: 'rgba(122, 185, 80, 0.85)' };

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

    // Scatter with an optional regression line and a legend for the groups.
    function scatter(div, s) {
        const H = phone() ? 250 : 300;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const legendH = s.groups.length > 1 ? 22 : 0;
        const left = phone() ? 46 : 56;
        const ar = { x: left, y: top + 4, w: W - left - 12, h: H - top - 4 - 46 - legendH };
        const shown = s.groups.filter(g => !g.hidden && g.xs.length);
        const allX = [].concat.apply([], s.groups.map(g => g.xs)), allY = [].concat.apply([], s.groups.map(g => g.ys));
        const xr = pad(allX), yr = pad(allY);
        const { sx, sy } = drawFrame(ctx, ar, xr, yr, { xLabel: s.xLabel, yLabel: s.yLabel, zeroY: true });
        for (const g of shown) {
            ctx.fillStyle = g.color; ctx.globalAlpha = g.alpha || 0.6;
            const r = g.size || 3;
            for (let i = 0; i < g.xs.length; i++) { ctx.beginPath(); ctx.arc(sx(g.xs[i]), sy(g.ys[i]), r, 0, Math.PI * 2); ctx.fill(); }
        }
        ctx.globalAlpha = 1;
        if (s.regression && shown.length) {
            const xs = [].concat.apply([], shown.map(g => g.xs)), ys = [].concat.apply([], shown.map(g => g.ys));
            const n = xs.length; let a = 0, b = 0, aa = 0, ab = 0;
            for (let i = 0; i < n; i++) { a += xs[i]; b += ys[i]; aa += xs[i] * xs[i]; ab += xs[i] * ys[i]; }
            const den = n * aa - a * a;
            if (Math.abs(den) > 1e-9) {
                const slope = (n * ab - a * b) / den, ic = (b - slope * a) / n;
                ctx.strokeStyle = C.green; ctx.lineWidth = 3; ctx.beginPath();
                ctx.moveTo(sx(xr[0]), sy(slope * xr[0] + ic)); ctx.lineTo(sx(xr[1]), sy(slope * xr[1] + ic)); ctx.stroke();
            }
        }
        if (legendH) drawLegend(ctx, ar.x, H - 10, s.groups.map(g => ({ name: g.name, color: g.color, hidden: g.hidden })));
    }

    // Horizontal box plots, one row per group, every point shown.
    function boxes(div, s) {
        const rows = s.rows;
        const rowH = phone() ? 26 : 30;
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
        const rowH = phone() ? 40 : 46;
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

    // Histogram of the panel with one line marked in red, as the cell line pages draw it.
    function hist(div, s) {
        const H = phone() ? 220 : 250;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const ar = { x: 24, y: top + 4, w: W - 36, h: H - top - 4 - 46 };
        const vals = s.vals.filter(v => isFinite(v));
        const lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
        const nb = 30, bw = (hi - lo) / nb || 1;
        const counts = new Array(nb).fill(0);
        for (const v of vals) counts[Math.min(nb - 1, Math.floor((v - lo) / bw))]++;
        const yr = [0, Math.max.apply(null, counts) * 1.05];
        const { sx, sy } = drawFrame(ctx, ar, [lo, hi], yr, { xLabel: s.xLabel, yTicks: false, zeroX: false, baseline: true, plain: true });
        ctx.fillStyle = C.dot;
        counts.forEach((c, i) => { const x0 = sx(lo + i * bw), x1 = sx(lo + (i + 1) * bw); ctx.fillRect(x0 + 0.5, sy(c), Math.max(1, x1 - x0 - 1), ar.y + ar.h - sy(c)); });
        if (typeof s.marker === 'number' && isFinite(s.marker)) {
            ctx.strokeStyle = C.red; ctx.lineWidth = 2.5; const x = sx(Math.min(hi, Math.max(lo, s.marker)));
            ctx.beginPath(); ctx.moveTo(x, ar.y); ctx.lineTo(x, ar.y + ar.h); ctx.stroke();
        }
    }

    // Horizontal bars, one per row.
    function hbars(div, s) {
        const rowH = phone() ? 22 : 24;
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
            ctx.fillStyle = C.bar; ctx.fillRect(ar.x, y, sx(r.value) - ar.x, rowH * 0.64);
            ctx.fillStyle = C.ink; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
            ctx.fillText(num(r.value), sx(r.value) + 5, y + rowH * 0.32);
        });
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

    // The gene set network as the app draws it: nodes on a ring, an edge for
    // every pair above the cutoff, blue positive and red negative, wider the
    // stronger, and nodes colored by their mean gene effect.
    function network(div, s) {
        const H = phone() ? 300 : 340;
        const { ctx, W } = setupCanvas(div, H);
        const top = drawHeader(ctx, W, s.title, s.sub);
        const cx = W / 2, cy = top + (H - top - 44) / 2;
        const R = Math.min(W, H - top - 44) / 2 - (phone() ? 36 : 44);
        const nr = phone() ? 15 : 18;
        const pos = s.nodes.map((n, i) => { const a = -Math.PI / 2 + i * 2 * Math.PI / s.nodes.length; return { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a) }; });
        const idx = new Map(s.nodes.map((n, i) => [n.name, i]));
        ctx.font = `9px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        for (const e of s.edges) {
            const a = pos[idx.get(e.a)], b = pos[idx.get(e.b)];
            if (!a || !b) continue;
            ctx.strokeStyle = e.r > 0 ? '#3182ce' : '#e53e3e'; ctx.lineWidth = 1 + 5 * Math.abs(e.r); ctx.globalAlpha = 0.85;
            ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
            ctx.globalAlpha = 1;
        }
        s.nodes.forEach((n, i) => {
            const pnt = pos[i];
            ctx.fillStyle = '#7ab950';
            ctx.strokeStyle = '#1f2937'; ctx.lineWidth = 1.5;
            ctx.beginPath(); ctx.arc(pnt.x, pnt.y, nr, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
            ctx.fillStyle = C.ink; ctx.font = `italic ${phone() ? 11 : 12}px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const outY = pnt.y + (pnt.y < cy - 1 ? -nr - 9 : nr + 9);
            ctx.fillText(n.name, pnt.x, Math.abs(pnt.y - cy) < 1 ? pnt.y + nr + 9 : outY);
        });
        // Legend: edges and the node color scale.
        const ly = H - 14;
        ctx.font = `${phone() ? 9 : 10}px ${FONT}`; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        let lx = 10;
        ctx.strokeStyle = '#3182ce'; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
        ctx.fillStyle = C.ink; ctx.fillText('r > 0', lx + 22, ly); lx += 22 + ctx.measureText('r > 0').width + 12;
        ctx.strokeStyle = '#e53e3e'; ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
        ctx.fillText('r < 0', lx + 22, ly); lx += 22 + ctx.measureText('r < 0').width + 16;
        ctx.fillText('wider = stronger', lx, ly);
    }

    // Colored grids come from the app's own drawer, the one the Matrix tab uses.
    function grid(div, spec) {
        const a = A();
        div.innerHTML = '';
        // Drawn into a child, because the drawer sizes its host to the grid.
        const inner = document.createElement('div');
        div.appendChild(inner);
        if (a && typeof a._drawCorrelationGrid === 'function') { a._drawCorrelationGrid(inner, spec); return; }
        div.innerHTML = '<div class="tour-loading">This chart is not available.</div>';
    }

    // Open a pair's scatter with the hotspot overlay, and optionally the
    // hotspot filter, already set: the popout applies the preset itself
    // just before its first draw.
    function openPairWithHotspot(a, g1, g2, mutGene, filterLevel) {
        a._inspectPreset = { hotspotGene: mutGene, hotspotMode: 'color', filterGene: filterLevel ? mutGene : null, filterLevel: filterLevel || null };
        a.openInspectByGenes(g1, g2);
    }

    const P53_SET = ['TP53', 'MDM2', 'MDM4', 'CDKN1A', 'PPM1D', 'USP7'];
    const setText = () => P53_SET.filter(g => A().geneIndex.has(g)).join('\n');
    // Run the example set in the app; optionally land on one result tab once
    // the results exist.
    function runExampleSet(tab) {
        document.getElementById('modeGeneSetBtn')?.click();
        const ta = document.getElementById('geneTextarea');
        if (ta) { ta.value = setText(); ta.dispatchEvent(new Event('input', { bubbles: true })); }
        const before = A().results;
        setTimeout(() => document.getElementById('runAnalysis')?.click(), 150);
        if (!tab) return;
        let tries = 0;
        const arm = () => {
            if (A().results && A().results !== before && A().results.success) { document.querySelector(`.nav-link[data-tab="${tab}"]`)?.click(); return; }
            if (tries++ < 60) setTimeout(arm, 250);
        };
        setTimeout(arm, 600);
    }

    // ---------------------------------------------------------------- pages
    // Each page: a title, a few plain sentences, an optional chart drawn from
    // the data, and the action that opens the same thing for real.
    const PAGES = [
        {
            title: 'Welcome',
            body: () => `<p>Correlate compiles data on ${num(D.ids().length)} human cancer cell lines from DepMap and other resources.</p>
                <p>From DepMap, the Cancer Dependency Map: CRISPR knockout screens showing which genes each cell line depends on, plus mutations, copy number, mRNA levels and drug response. From other resources: Cellosaurus, Oncotree, curated driver gene and fusion lists, published breast cancer subtypes, and a retroelement signal from public RNA-seq.</p>
                <p>Each page of this tour shows one feature with a real chart and a button that opens it in the app.</p>`,
        },
        {
            title: 'Gene set analysis',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4,
            body: () => `<p>Gene set analysis is the central tool of Correlate. Paste a set of genes and the app correlates their dependency scores across all cell lines. Pairs above the cutoff become links, and the set becomes a network.</p>
                <p>Here, the p53 pathway. Blue links are positive correlations, red negative. Nodes can be colored by gene effect, or by your own statistics with the "With Stats" input.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const edges = [];
                for (let a = 0; a < genes.length; a++) for (let b = a + 1; b < genes.length; b++) {
                    const st = pearson(D.row(genes[a]), D.row(genes[b]));
                    if (st && Math.abs(st.r) >= 0.5) edges.push({ a: genes[a], b: genes[b], r: st.r });
                }
                const nodes = genes.map(g => {
                    const row = D.row(g); let n = 0, sum = 0;
                    for (let i = 0; i < D.ids().length; i++) { const v = D.ge(row, i); if (isFinite(v)) { n++; sum += v; } }
                    return { name: g, value: n ? sum / n : 0 };
                });
                network(div, { title: 'Gene set network, p53 pathway', sub: 'links at |r| of 0.5 and above, the default cutoff', nodes, edges });
            },
            action: {
                label: 'Run this gene set in the app',
                run: () => runExampleSet(null)
            }
        },
        {
            title: 'What one link means',
            need: () => !!D.row('TP53') && !!D.row('MDM2'),
            body: () => `<p>One link is a scatter: each dot is a cell line, placed by its score for each gene. Here r = ${(() => { const st = pearson(D.row('TP53'), D.row('MDM2')); return st ? st.r.toFixed(2) : '?'; })()} across ${num(scatterXY('TP53', 'MDM2').xs.length)} cell lines.</p>
                <p>Cell lines with working p53 need MDM2 to keep it in check; cell lines that have lost p53 do not. The correlation reveals that relationship.</p>`,
            plot: (div) => {
                const { xs, ys } = scatterXY('TP53', 'MDM2');
                scatter(div, { title: 'TP53 vs MDM2', sub: `n = ${num(xs.length)} cell lines`, xLabel: 'TP53 Gene Effect', yLabel: 'MDM2 Gene Effect',
                    groups: [{ xs, ys, color: C.dot, alpha: 0.6, size: 3 }], regression: true });
            },
            action: { label: 'Open TP53 vs MDM2 in the app', run: (a) => a.openInspectByGenes('TP53', 'MDM2') }
        },
        {
            title: 'The same pair, split by TP53 mutation',
            need: () => !!D.row('TP53') && !!D.row('MDM2') && !!A().mutations?.geneData?.TP53,
            body: () => `<p>The same scatter, colored by TP53 mutation: grey wild-type, blue one mutated copy, red both. The wild-type cell lines are the ones that depend on MDM2.</p>
                <p>Any scatter can be colored this way with the Hotspot overlay, or limited to one group with the Hotspot filter.</p>`,
            plot: (div) => {
                const g = scatterByMutation('TP53', 'MDM2', 'TP53');
                const spec = { title: 'TP53 vs MDM2', sub: 'colored by TP53 hotspot mutation', xLabel: 'TP53 Gene Effect', yLabel: 'MDM2 Gene Effect',
                    groups: [
                        { name: `WT (n=${g.wt.xs.length})`, xs: g.wt.xs, ys: g.wt.ys, color: C.dot, alpha: 0.6, size: 3 },
                        { name: `1 mut (n=${g.m1.xs.length})`, xs: g.m1.xs, ys: g.m1.ys, color: C.blue, alpha: 0.8, size: 3.2 },
                        { name: `2 mut (n=${g.m2.xs.length})`, xs: g.m2.xs, ys: g.m2.ys, color: C.red, alpha: 0.8, size: 3.2 }
                    ] };
                div._spec = spec;
                scatter(div, spec);
            },
            // Which group to show, chosen by the reader.
            choices: [
                { label: 'Both groups', run: (div) => { div._spec.groups.forEach(g => g.hidden = false); scatter(div, div._spec); } },
                { label: 'Wild-type only', run: (div) => { div._spec.groups.forEach((g, i) => g.hidden = i !== 0); scatter(div, div._spec); } },
                { label: 'Mutated only', run: (div) => { div._spec.groups.forEach((g, i) => g.hidden = i === 0); scatter(div, div._spec); } }
            ],
            actions: [
                { label: 'Open with the TP53 overlay', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', null) },
                { label: 'Open only the TP53 wild-type lines', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', '0') },
                { label: 'Open only the TP53 mutated lines', run: (a) => openPairWithHotspot(a, 'TP53', 'MDM2', 'TP53', '1+2') }
            ]
        },
        {
            title: 'The correlation matrix',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4,
            body: () => `<p>The Matrix tab shows every pair in the set, above or below the cutoff, as one grid: red positive, blue negative, the number is r.</p>
                <p>Click a cell to open its scatter; blank the pairs below the cutoff; export as image or CSV.</p>`,
            plot: (div) => {
                const genes = P53_SET.filter(g => D.row(g));
                const z = genes.map(g1 => genes.map(g2 => g1 === g2 ? 1 : (pearson(D.row(g1), D.row(g2))?.r ?? null)));
                grid(div, { title: 'Correlation of gene effects, p53 pathway', rowLabels: genes, colLabels: genes, z, zmin: -1, zmax: 1,
                    colorscale: [[0, '#2166ac'], [0.5, '#f7f7f7'], [1, '#b2182b']], showValues: true, colorbarTitle: 'r' });
            },
            action: {
                label: 'Run this set and open the Matrix tab',
                run: () => runExampleSet('matrix')
            }
        },
        {
            title: 'The gene effect score',
            need: () => !!D.row('SOX10'),
            body: () => `<p>A gene effect score says how much a cell line needs a gene: 0 means no effect, about -1 is a typical essential gene, more negative is a stronger need.</p>
                <p>SOX10 across tissues: skin cell lines sit far to the left, the others near zero.</p>`,
            plot: (div) => {
                const groups = tissueGroups('SOX10', false, 14).sort((a, b) => a.med - b.med);
                const keep = groups.slice(0, 3).concat(groups.slice(-4));
                boxes(div, { title: 'SOX10 Gene Effect by tissue', xLabel: 'SOX10 Gene Effect', rows: keep });
            },
            action: { label: 'Open SOX10 in the Gene Effect view', run: (a) => a.openGeneEffectModal('SOX10', 'tissue') }
        },
        {
            title: 'Coloring a scatter by a mutation',
            need: () => !!D.row('BRAF') && !!D.row('MAPK1') && !!A().mutations?.geneData?.BRAF,
            body: () => `<p>BRAF against MAPK1, colored by BRAF mutation. The mutated cell lines sit low on both axes: they depend on BRAF and on the kinase below it.</p>
                <p>A dependency that follows a mutation is what a targeted drug is built on.</p>`,
            plot: (div) => {
                const g = scatterByMutation('BRAF', 'MAPK1', 'BRAF');
                scatter(div, { title: 'BRAF vs MAPK1', sub: 'colored by BRAF hotspot mutation', xLabel: 'BRAF Gene Effect', yLabel: 'MAPK1 Gene Effect',
                    groups: [
                        { name: `WT (n=${g.wt.xs.length})`, xs: g.wt.xs, ys: g.wt.ys, color: C.dot, alpha: 0.6, size: 3 },
                        { name: `1 mut (n=${g.m1.xs.length})`, xs: g.m1.xs, ys: g.m1.ys, color: C.blue, alpha: 0.8, size: 3.2 },
                        { name: `2 mut (n=${g.m2.xs.length})`, xs: g.m2.xs, ys: g.m2.ys, color: C.red, alpha: 0.8, size: 3.2 }
                    ] });
            },
            action: { label: 'Open this scatter with the overlay', run: (a) => openPairWithHotspot(a, 'BRAF', 'MAPK1', 'BRAF', null) }
        },
        {
            title: 'Mutation analysis',
            need: () => !!D.row('BRAF') && !!A().mutations?.geneData?.BRAF,
            body: () => `<p>Mutation analysis splits the panel by one mutation and ranks the genes whose dependency differs most between the two groups.</p>
                <p>The simplest case: BRAF's own score by BRAF mutation status. The thick mark is the median.</p>`,
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
            title: 'The Cell Line Browser',
            body: () => `<p>The browser finds the cell lines that fit a project: filter by tissue, subtype, disease, sex, mutation, fusion or copy number, tick the ones you want, and send them on to a heatmap, an export or a comparison.</p>
                <p>The chart shows the panel by tissue.</p>`,
            plot: (div) => {
                const counts = new Map();
                D.ids().forEach(id => { const l = D.lineage(id); if (l) counts.set(l, (counts.get(l) || 0) + 1); });
                const rows = [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 12).map(e => ({ name: e[0], value: e[1] }));
                hbars(div, { title: 'Cell lines per tissue', sub: 'the twelve largest groups', xLabel: 'Cell lines', rows });
            },
            action: { label: 'Open the Cell Line Browser', run: (a) => a.openCellLineBrowser() }
        },
        {
            title: 'A cell line\'s page',
            need: () => !!D.lineId('A375') && !!A().globalSignatures?.byCellLine,
            body: () => `<p>Every cell line has a page: origin, drivers, copy number, fusions, dependencies, expression, drug response and how to authenticate a stock.</p>
                <p>Charts on the page place the cell line among all others. Here: ploidy, with A375 in red.</p>`,
            plot: (div) => {
                const id = D.lineId('A375'), sig = A().globalSignatures.byCellLine;
                const vals = []; D.ids().forEach(x => { const v = sig[x]?.Ploidy; if (typeof v === 'number') vals.push(v); });
                hist(div, { title: 'Ploidy across the panel', sub: 'A375 marked in red', xLabel: 'Ploidy', vals, marker: sig[id]?.Ploidy });
            },
            action: { label: 'Open the A375 page', run: (a) => a.openCellLineWiki(D.lineId('A375')) }
        },
        {
            title: 'Expression',
            need: () => !!A().loadExpressionData,
            body: () => `<p>The app also carries mRNA levels, as log2(TPM+1). Every gene view and scatter can switch from gene effect to mRNA.</p>
                <p>MITF, the melanocyte lineage factor, is expressed far above the rest in skin cell lines.</p>`,
            plot: async (div) => {
                if (!A().expressionLoaded) {
                    div.innerHTML = '<div class="tour-loading">Loading expression data...</div>';
                    try { await A().loadExpressionData(); } catch (e) { }
                    if (!A().expressionLoaded) { div.innerHTML = '<div class="tour-loading">Expression data is not available right now.</div>'; return; }
                }
                const groups = tissueGroups('MITF', true, 14).sort((a, b) => b.med - a.med);
                const keep = groups.slice(0, 4).concat(groups.slice(-3));
                boxes(div, { title: 'MITF expression by tissue', xLabel: 'MITF mRNA, log2(TPM+1)', rows: keep, zero: false });
            },
            action: { label: 'Open MITF expression by tissue', run: (a) => a.openGeneEffectModal('MITF', 'tissue', { dataType: 'expr' }) }
        },
        {
            title: 'The gene set heatmap',
            need: () => P53_SET.filter(g => D.row(g)).length >= 4 && typeof A()._hmOpenModal === 'function',
            body: () => `<p>The heatmap shows many genes across many cell lines at once. Each row is a gene, each column a cell line, colored by that cell line's score, z-scored per gene.</p>
                <p>Here the p53 pathway across skin cell lines. Orange is a dependency, purple means the knockout helped growth.</p>`,
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
                const use = cols.slice(0, phone() ? 30 : 45);
                const z = stats.map(st => use.map(i => { const v = D.ge(st.row, i); return isFinite(v) ? (v - st.m) / st.sd : null; }));
                grid(div, { title: 'p53 pathway genes across skin lines', sub: 'gene effect, z-scored per gene', rowLabels: genes,
                    colLabels: use.map(i => D.name(ids[i])), showColLabels: false, xLabel: `${use.length} skin cell lines`, z, zmin: -2.5, zmax: 2.5,
                    colorscale: [[0, '#e66101'], [0.5, '#f7f7f7'], [1, '#5e3c99']], showValues: false, colorbarTitle: 'z' });
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
            body: () => `<p>Each cell line's page carries its PRISM drug screen: AUC from 0 (all cells killed) to 1 (no effect).</p>
                <p>Vemurafenib, the BRAF V600E inhibitor, across the panel, with A375 in red. A375 carries BRAF V600E and sits among the most sensitive.</p>`,
            plot: (div) => {
                const id = D.lineId('A375');
                const cp = A().drugResponse.compounds.find(c => /vemurafenib/i.test(c.name || ''));
                const vals = []; D.ids().forEach(x => { const v = cp.auc?.[x]; if (typeof v === 'number') vals.push(v); });
                hist(div, { title: 'Vemurafenib response across the panel', sub: 'A375 marked in red', xLabel: 'AUC (0 = all cells killed, 1 = no effect)', vals, marker: cp.auc?.[id] });
            },
            action: { label: 'Open the A375 page', run: (a) => a.openCellLineWiki(D.lineId('A375')) }
        },
        {
            title: 'Where to go next',
            body: () => `<p>Three ways to start: run the p53 example and swap in your own genes; open the Cell Line Browser and filter to your tissue; open a correlation for two genes you know.</p>
                <p>Every popout has Export buttons for figures and data. The last button opens the longer explanation of the statistics.</p>`,
            actions: [
                { label: 'Run the p53 example', run: () => runExampleSet(null) },
                { label: 'Open the Cell Line Browser', run: (a) => a.openCellLineBrowser() },
                { label: 'Open TP53 vs MDM2', run: (a) => a.openInspectByGenes('TP53', 'MDM2') },
                { label: 'The statistics in more depth', run: () => { const m = document.getElementById('infographicModal'); if (m) m.style.display = 'flex'; } }
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
#tourModal .tour-chart { min-height: 120px; margin: 4px 0 12px; border: 1px solid var(--gray-200); border-radius: 6px; overflow-x: auto; overflow-y: hidden; padding: 6px 0; background: #fff; }
#tourModal .tour-loading { padding: 40px 12px; text-align: center; color: #6b7280; font-size: 12px; }
#tourModal .tour-choices { display: flex; flex-wrap: nowrap; gap: 6px; margin: -4px 0 10px; }
#tourModal .tour-choices:empty { display: none; }
#tourModal .tour-choice { flex: 1 1 0; min-width: 0; min-height: 34px; font-size: 12px; padding: 4px 4px; white-space: nowrap; }
#tourModal .tour-choice.on { background: var(--green-600); color: #fff; border-color: var(--green-600); }
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
                    <div class="tour-text" id="tourText"></div>
                    <div class="tour-chart" id="tourChart" style="display:none;"></div>
                    <div class="tour-choices" id="tourChoices"></div>
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
        document.getElementById('tourTitle').innerHTML = `Tour: ${esc(p.title)} <span class="tour-step">${step + 1} of ${pages.length}</span>`;
        document.getElementById('tourText').innerHTML = p.body();
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
