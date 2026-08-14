#!/usr/bin/env python3
"""Does the streaming method reproduce the index-based pilot?

Same cell lines, same panel, two different access paths. If per-locus counts
agree, the cloud run is trustworthy; if they do not, nothing downstream is.
"""
import json, sys
new = {r['model_id']: r for r in json.load(open(sys.argv[1])) if not r.get('error')}
old = {r['model_id']: r for r in json.load(open(sys.argv[2])) if not r.get('error')}
shared = sorted(set(new) & set(old))
if not shared:
    print("no overlapping cell lines to compare"); sys.exit(1)
print(f"comparing {len(shared)} cell lines\n")
allok = True
for m in shared:
    a, b = new[m], old[m]
    ids = set(a['counts_unique']) | set(b['counts_unique'])
    diffs = [(i, a['counts_unique'].get(i, 0), b['counts_unique'].get(i, 0)) for i in ids]
    nz = [(i, x, y) for i, x, y in diffs if x or y]
    exact = sum(1 for _, x, y in nz if x == y)
    close = sum(1 for _, x, y in nz if x != y and abs(x-y) <= max(2, 0.05*max(x, y)))
    bad = [(i, x, y) for i, x, y in nz if x != y and abs(x-y) > max(2, 0.05*max(x, y))]
    rd = abs(a['total_mapped'] - b['total_mapped']) / max(1, b['total_mapped'])
    ok = not bad and rd < 0.15
    allok &= ok
    print(f"  {a.get('name','?'):16s} loci nonzero={len(nz):4d} exact={exact:4d} close={close:3d} "
          f"off={len(bad):3d} libsize_diff={rd*100:4.1f}%  {'OK' if ok else 'CHECK'}")
    for i, x, y in bad[:3]:
        print(f"      {i}: stream={x} indexed={y}")
print("\n" + ("VERIFIED - streaming reproduces the pilot" if allok else
      "MISMATCH - do not run the full job until this is understood"))
