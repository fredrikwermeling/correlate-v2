#!/usr/bin/env python3
"""Summarize a scan: which cell lines carry which virus, and how that agrees
with the Cellosaurus records the app already has (web_data/virus_status.json).

A call needs reads AND breadth: at least 10 reads per million unmapped reads
scanned, and at least 20% of the genome's 50 bins touched. Breadth is what
separates an infection from a vector fragment or a shared low-complexity
stretch.

Usage:  python3 analyze.py counts_all.json
"""
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_RPM, MIN_BREADTH, MIN_READS = 10, 0.2, 5


def main():
    res = [r for r in json.load(open(sys.argv[1])) if not r.get('error')]
    vs = json.load(open(os.path.join(HERE, '..', '..', 'web_data', 'virus_status.json')))['byCellLine']
    known = {cl: {a['agent'].replace('-', '') for a in recs} for cl, recs in vs.items()}
    calls = defaultdict(list)
    for r in res:
        n = max(1, r.get('unmapped_scanned', 0))
        for v, k in r.get('hits', {}).items():
            rpm, b = k * 1e6 / n, r['breadth'].get(v, 0)
            if k >= MIN_READS and rpm >= MIN_RPM and b >= MIN_BREADTH:
                calls[v].append((r['name'], r['lineage'], k, round(rpm, 1), b, r['model_id']))
    print(f'{len(res)} cell lines scanned\n')
    for v, rows in sorted(calls.items(), key=lambda x: -len(x[1])):
        print(f'== {v}: {len(rows)} cell lines')
        for name, lin, k, rpm, b, cid in sorted(rows, key=lambda x: -x[3]):
            rec = 'Cellosaurus agrees' if any(v.startswith(a[:3]) for a in known.get(cid, ())) else ''
            print(f'   {name:14s} {lin[:22]:22s} reads={k:<7} per M unmapped={rpm:<8} breadth={b}  {rec}')
    print('\nCellosaurus-recorded agents not found:')
    called = {(c[5], v) for v, rows in calls.items() for c in rows}
    for r in res:
        for a in known.get(r['model_id'], ()):
            if not any(cid == r['model_id'] and v.startswith(a[:3]) for cid, v in called):
                print(f"   {r['name']:14s} {a}")


if __name__ == '__main__':
    main()
