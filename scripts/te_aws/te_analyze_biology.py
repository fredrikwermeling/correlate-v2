#!/usr/bin/env python3
"""The question the whole exercise was for: does retroelement transcription
predict anything?

Three tests, in increasing order of interest:
  1. Sanity     - do the counts still behave (controls, spread)?
  2. Viral mimicry - do de-repressed lines carry an interferon signature?
  3. Dependency - do they lean harder on the brakes (ADAR1, TREX1, SAMHD1,
                  RNASEH2), which is the therapeutically interesting claim?

Usage: python3 te_analyze_biology.py counts_all.json [/path/to/web_data]
"""
import json, sys, gzip, array, math
from collections import defaultdict

WEB = sys.argv[2] if len(sys.argv) > 2 else '/Users/fredrikwermeling/Documents/correlate_v2/web_data'
COUNTS = sys.argv[1] if len(sys.argv) > 1 else 'counts_all.json'

# ---------- helpers ---------------------------------------------------------
def pearson(a, b):
    n = len(a)
    if n < 8: return float('nan'), float('nan'), n
    ma, mb = sum(a)/n, sum(b)/n
    va = sum((x-ma)**2 for x in a); vb = sum((y-mb)**2 for y in b)
    if va <= 0 or vb <= 0: return float('nan'), float('nan'), n
    r = sum((x-ma)*(y-mb) for x, y in zip(a, b)) / math.sqrt(va*vb)
    r = max(-0.999999, min(0.999999, r))
    # two-sided p from the t approximation; adequate at these n
    t = abs(r) * math.sqrt((n-2)/(1-r*r))
    # Abramowitz-Stegun normal tail, t is near-normal for n in the hundreds
    z = t
    p = 2 * (1 - (1 - 0.5*math.exp(-0.717*z - 0.416*z*z))) if z < 6 else 1e-9
    p = min(1.0, max(1e-300, p))
    return r, p, n

def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0]*len(v); i = 0
        while i < len(order):
            j = i
            while j+1 < len(order) and v[order[j+1]] == v[order[i]]: j += 1
            avg = (i+j)/2 + 1
            for k in range(i, j+1): r[order[k]] = avg
            i = j+1
        return r
    return pearson(rank(a), rank(b))

def load_matrix(meta_path, bin_path):
    meta = json.load(open(meta_path))
    a = array.array('h'); a.frombytes(gzip.open(bin_path, 'rb').read())
    return meta, a

# ---------- 1. retroelement scores per cell line -----------------------------
res = [r for r in json.load(open(COUNTS)) if r.get('total_mapped') and not r.get('error')]
panel = {e['id']: e for e in json.load(open('te_panel.json'))}
FAM = defaultdict(list)
for i, e in panel.items(): FAM[e.get('role') == 'element' and e['family'] or e['role']].append(i)

te = {}
for r in res:
    m = r['total_mapped']/1e6
    cpm = lambda ids: sum((r['counts_unique'].get(i) or 0) for i in ids)/m
    el = [i for i, e in panel.items() if e['role'] == 'element']
    vals = [(r['counts_unique'].get(i) or 0)/m for i in el]
    te[r['model_id']] = {
        'name': r.get('name'), 'lineage': r.get('lineage'),
        'total': sum(vals), 'n_active': sum(1 for v in vals if v > 0.5),
        'L1': cpm(FAM['L1']), 'HERVK': cpm(FAM['HERVK']), 'SVA': cpm(FAM['SVA']),
        'passenger': cpm([i for i, e in panel.items() if e['role'] == 'control_passenger']),
        'housekeeping': cpm([i for i, e in panel.items() if e['role'] == 'control_housekeeping']),
    }
print(f"cell lines with retroelement data: {len(te)}")

# ---------- 2. interferon score, same formula the app uses -------------------
emeta, earr = load_matrix(f'{WEB}/expression_metadata.json', f'{WEB}/expression.bin.gz')
enC, eSF, eNA = emeta['nCellLines'], emeta['scaleFactor'], emeta['naValue']
egi = {g.upper(): i for i, g in enumerate(emeta['genes'])}
ISG = ['ISG15','IFIT1','IFIT2','IFIT3','MX1','MX2','OAS1','OAS2','OAS3','OASL','IFI6','IFI27',
       'IFI44','IFI44L','RSAD2','STAT1','STAT2','IRF7','IRF9','BST2','XAF1','HERC5','USP18',
       'RIGI','IFIH1','CMPK2','EPSTI1','LY6E','SAMD9','SAMD9L','IFITM1','IFITM3','PARP9','DTX3L']
zrows = []
for g in ISG:
    gi = egi.get(g)
    if gi is None: continue
    row = [earr[gi*enC + i] for i in range(enC)]
    ok = [v/eSF for v in row if v != eNA]
    if len(ok) < 10: continue
    mu = sum(ok)/len(ok); sd = math.sqrt(sum((v-mu)**2 for v in ok)/(len(ok)-1))
    if sd <= 0: continue
    zrows.append((gi, mu, sd))
ifn = {}
for i, cl in enumerate(emeta['cellLines']):
    zs = []
    for gi, mu, sd in zrows:
        v = earr[gi*enC + i]
        if v != eNA: zs.append((v/eSF - mu)/sd)
    if len(zs) >= 0.6*len(zrows): ifn[cl] = sum(zs)/len(zs)
print(f"IFN score computed for {len(ifn)} lines from {len(zrows)} ISGs")

# ---------- 3. gene effect ---------------------------------------------------
gmeta, garr = load_matrix(f'{WEB}/metadata.json', f'{WEB}/geneEffects.bin.gz')
gnC, gSF, gNA = gmeta['nCellLines'], gmeta['scaleFactor'], gmeta['naValue']
ggi = {g.upper(): i for i, g in enumerate(gmeta['genes'])}
gcl = {cl: i for i, cl in enumerate(gmeta['cellLines'])}
def dep(gene, cl):
    gi, ci = ggi.get(gene.upper()), gcl.get(cl)
    if gi is None or ci is None: return None
    v = garr[gi*gnC + ci]
    return None if v == gNA else v/gSF

# ---------- report -----------------------------------------------------------
print("\n=== 1. CONTROLS (should stay low) ===")
ids = [c for c in te]
for lab, key in [('vs genic-L1 passenger', 'passenger'), ('vs housekeeping', 'housekeeping')]:
    a = [te[c]['total'] for c in ids]; b = [te[c][key] for c in ids]
    r, p, n = pearson(a, b)
    print(f"  element signal {lab:24s} r = {r:+.2f}  (n={n})")

print("\n=== 2. VIRAL MIMICRY: retroelements vs interferon ===")
shared = [c for c in ids if c in ifn]
print(f"  cell lines with both: {len(shared)}")
for lab, key in [('all elements', 'total'), ('LINE-1', 'L1'), ('HERV-K', 'HERVK'), ('SVA', 'SVA'), ('elements switched on', 'n_active')]:
    a = [te[c][key] for c in shared]; b = [ifn[c] for c in shared]
    r, p, n = pearson(a, b); rs, ps, _ = spearman(a, b)
    print(f"  {lab:22s} Pearson r = {r:+.3f} (p = {p:.2g})   Spearman rho = {rs:+.3f}")

print("\n=== 3. DEPENDENCY: do de-repressed lines lean on the brakes? ===")
GENES = ['ADAR','TREX1','SAMHD1','RNASEH2A','RNASEH2B','RNASEH2C','CGAS','STING1','IFIH1',
         'RIGI','MAVS','EIF2AK2','TBK1','STAT1','JAK1','SETDB1','TRIM28','MPHOSPH8','TASOR','MORC2']
rows = []
for g in GENES:
    pairs = [(te[c]['total'], dep(g, c)) for c in ids if dep(g, c) is not None]
    if len(pairs) < 30: continue
    a = [x for x, _ in pairs]; b = [y for _, y in pairs]
    r, p, n = pearson(a, b)
    rows.append((g, r, p, n))
rows.sort(key=lambda x: x[1])
print(f"  {'gene':10s} {'r':>7s} {'p':>10s} {'n':>5s}   negative r = more dependent when elements are high")
for g, r, p, n in rows:
    star = ' *' if p < 0.05/len(rows) else ''
    print(f"  {g:10s} {r:+7.3f} {p:10.2g} {n:5d}{star}")
print(f"\n  * survives Bonferroni across {len(rows)} genes")

print("\n=== top 10 lines by retroelement signal ===")
for c in sorted(ids, key=lambda c: -te[c]['total'])[:10]:
    i = ifn.get(c)
    print(f"  {str(te[c]['name'])[:18]:18s} {str(te[c]['lineage'])[:22]:22s} total {te[c]['total']:7.1f} CPM  "
          f"active {te[c]['n_active']:3d}  IFN {('%+.2f' % i) if i is not None else '   n/a'}")
