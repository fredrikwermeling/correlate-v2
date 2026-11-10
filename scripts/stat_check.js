// Checks the app's own statistics against reference values from scipy.
// Run: node scripts/stat_check.js   (exits 1 on a mismatch)
// The methods are lifted out of app.js as text, so this tests the shipped code.
const fs = require('fs');
const path = require('path');
const src = fs.readFileSync(path.join(__dirname, '..', 'app.js'), 'utf8');

function grab(name) {
    const i = src.indexOf('\n    ' + name + '(');
    if (i < 0) throw new Error('method not found: ' + name);
    let depth = 0, j = src.indexOf('{', i);
    for (; ; j++) {
        if (src[j] === '{') depth++;
        if (src[j] === '}') { depth--; if (!depth) break; }
    }
    return src.slice(i + 5, j + 1);
}

const names = ['_chi2P1', 'normalUpperTail', '_slopeDiffP', 'tDistributionPValue', 'incompleteBeta',
    'betaCF', 'logGamma', '_pearsonP', '_welchFromSums', '_mannWhitneyP', '_normalCdf'];
const T = eval('(class {' + names.map(grab).join('\n') + '})');
const t = new T();

const pts = (rows) => rows.map(([x, y]) => ({ x, y }));
const g1 = pts([[1, 1.1], [2, 1.9], [3, 3.2], [4, 3.9], [5, 5.1], [6, 6.2]]);
const g2 = pts([[1, 0.9], [2, 1.2], [3, 1.4], [4, 2.1], [5, 2.0], [6, 2.9], [7, 3.1]]);
const a = [1.2, 3.4, 2.2, 5.1, 4.4, 6.0, 2.9], b = [0.5, 1.1, 2.0, 1.7, 0.9, 2.4, 1.3, 0.2];
const sums = (v) => [v.reduce((s, x) => s + x, 0), v.reduce((s, x) => s + x * x, 0), v.length];

// [label, value from the app, scipy reference, relative tolerance]
const cases = [
    ['chi2 sf(5.98, df 1)', t._chi2P1(5.98), 0.01446900220244473, 1e-6],
    ['chi2 sf(12.87, df 1)', t._chi2P1(12.87), 0.00033389160933658526, 1e-6],
    ['chi2 sf(0.5, df 1)', t._chi2P1(0.5), 0.47950012218695337, 1e-6],
    ['Pearson p, r 0.9, n 5', Number(t._pearsonP(0.9, 5)), 0.03738607346849863, 2e-3],
    ['Pearson p, r -0.35, n 40', Number(t._pearsonP(-0.35, 40)), 0.026832078464118343, 2e-3],
    ['slope difference p', t._slopeDiffP(g1, g2), 5.491526426653061e-07, 1e-6],
    ['Welch t p', t._welchFromSums(...sums(a), ...sums(b)).p, 0.009242968731417984, 1e-6],
    ['Mann-Whitney p', t._mannWhitneyP(a, b).p, 0.009217797347054925, 1e-4],
];

let bad = 0;
for (const [label, got, want, tol] of cases) {
    const ok = isFinite(got) && Math.abs(got - want) <= tol * Math.abs(want);
    if (!ok) bad++;
    console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: ${got} (scipy ${want})`);
}
if (isFinite(t._slopeDiffP(g1.slice(0, 2), g2))) { bad++; console.log('FAIL slope difference with 2 points should be NaN'); }
process.exit(bad ? 1 : 0);
