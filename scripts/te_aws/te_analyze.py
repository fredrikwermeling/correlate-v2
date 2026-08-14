#!/usr/bin/env python3
"""Validate the pilot: do the numbers behave like biology, or like artefact?

Three checks the data has to survive before any of it reaches the app:
  1. Library health   - housekeeping counts sane, comparable across lines.
  2. Passenger control- intergenic elements must NOT simply track genic L1s.
                        If they do, we are measuring transcription in general.
  3. Dynamic range    - do lines actually differ, and is the difference
                        carried by many elements or by one outlier locus?
"""
import json, sys, statistics as st

res = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'te_counts_pilot.json'))
panel = {e['id']: e for e in json.load(open('te_panel.json'))}
res = [r for r in res if not r.get('error') and r.get('total_mapped')]
print(f"lines with data: {len(res)}\n")

def cpm(r, ids):
    m = r['total_mapped'] / 1e6
    v = [r['counts_unique'].get(i) for i in ids]
    v = [x for x in v if x is not None]
    return [x / m for x in v]

E  = [i for i, e in panel.items() if e['role'] == 'element']
L1 = [i for i, e in panel.items() if e['role']=='element' and e['family']=='L1']
HK = [i for i, e in panel.items() if e['role'] == 'control_housekeeping']
PS = [i for i, e in panel.items() if e['role'] == 'control_passenger']

rows = []
for r in res:
    el, hk, ps = cpm(r, E), cpm(r, HK), cpm(r, PS)
    l1 = cpm(r, L1)
    rows.append({'name': r['name'], 'lineage': r.get('lineage'),
                 'L1TD1': r.get('L1TD1_expr'),
                 'mapped_M': r['total_mapped']/1e6,
                 'hk': sum(hk), 'elem_sum': sum(el), 'l1_sum': sum(l1),
                 'passenger': sum(ps),
                 'n_active': sum(1 for x in el if x > 0.5)})   # >0.5 CPM = on

print("=== 1. LIBRARY HEALTH (housekeeping CPM) ===")
hks = [x['hk'] for x in rows]
print(f"  median {st.median(hks):.0f}  min {min(hks):.0f}  max {max(hks):.0f}"
      f"  -> {'OK, comparable' if max(hks) < 20*min(hks) else 'WARNING: uneven'}")

print("\n=== 2. PASSENGER CONTROL ===")
def pearson(a, b):
    n=len(a); ma, mb = sum(a)/n, sum(b)/n
    num = sum((x-ma)*(y-mb) for x,y in zip(a,b))
    da = (sum((x-ma)**2 for x in a))**.5; db=(sum((y-mb)**2 for y in b))**.5
    return num/(da*db) if da and db else float('nan')
r_pass = pearson([x['elem_sum'] for x in rows], [x['passenger'] for x in rows])
r_hk   = pearson([x['elem_sum'] for x in rows], [x['hk'] for x in rows])
print(f"  intergenic elements vs genic-L1 passenger : r = {r_pass:+.2f}")
print(f"  intergenic elements vs housekeeping       : r = {r_hk:+.2f}")
print("  -> " + ("CLEAN: element signal is not just bulk transcription"
                 if abs(r_pass) < 0.7 and abs(r_hk) < 0.7 else
                 "SUSPECT: tracks general transcription, needs a stricter filter"))

print("\n=== 3. DYNAMIC RANGE (elements switched on, >0.5 CPM) ===")
rows.sort(key=lambda x: -x['n_active'])
print(f"  {'line':16s} {'lineage':22s} {'active':>7s} {'elemCPM':>9s} {'L1TD1':>6s}")
for x in rows[:8]:
    print(f"  {x['name']:16s} {str(x['lineage'])[:22]:22s} {x['n_active']:7d} {x['elem_sum']:9.1f} {x['L1TD1']:6.2f}")
print("   ...")
for x in rows[-5:]:
    print(f"  {x['name']:16s} {str(x['lineage'])[:22]:22s} {x['n_active']:7d} {x['elem_sum']:9.1f} {x['L1TD1']:6.2f}")
na = [x['n_active'] for x in rows]
print(f"\n  active elements: median {st.median(na):.0f}, range {min(na)}-{max(na)}"
      f"  -> {'usable spread' if max(na) > 2*max(1,min(na)) else 'TOO FLAT to be a browser column'}")
json.dump(rows, open('pilot_summary.json','w'), indent=1)
