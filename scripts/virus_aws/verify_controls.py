#!/usr/bin/env python3
"""Check the five control cell lines before trusting a run.

Cellosaurus records Raji and NAMALWA as EBV-transformed, PLC/PRF/5 as
carrying hepatitis B and MJ as HTLV-1 positive; RT4 has no record and is the
negative control for those four viruses.

Usage:  python3 verify_controls.py counts_verify.json
"""
import json
import sys

EXPECT = {'Raji': 'EBV', 'NAMALWA': 'EBV', 'PLC/PRF/5': 'HBV', 'MJ': 'HTLV1'}
MIN_READS, MIN_BREADTH = 20, 0.2


def main():
    res = {r['name']: r for r in json.load(open(sys.argv[1]))}
    ok = True
    for name, virus in EXPECT.items():
        r = res.get(name, {})
        n, b = r.get('hits', {}).get(virus, 0), r.get('breadth', {}).get(virus, 0)
        good = n >= MIN_READS and b >= MIN_BREADTH
        ok &= good
        print(f"{'ok  ' if good else 'FAIL'} {name:10s} {virus:6s} reads={n} breadth={b}")
    rt4 = res.get('RT4', {})
    for virus in set(EXPECT.values()):
        n = rt4.get('hits', {}).get(virus, 0)
        good = n < MIN_READS
        ok &= good
        print(f"{'ok  ' if good else 'FAIL'} RT4        {virus:6s} reads={n} (negative control)")
    print('VERIFIED' if ok else 'MISMATCH: do not start the full run')


if __name__ == '__main__':
    main()
