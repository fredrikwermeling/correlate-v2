#!/usr/bin/env python3
"""Rebuild Green Listed's globalSignatures.json from DepMap's
OmicsGlobalSignatures.csv.

Per-cell-line genome-wide features computed once for the panel:
  Ploidy — average chromosome dosage (from PureCN); ~2 in a diploid line,
           ~4 in a whole-genome-doubled line.
  WGD    — whole-genome-doubling flag (boolean).
  CIN    — chromosomal instability score.
  LoHFraction — loss-of-heterozygosity fraction.
  MSIScore — MSIsensor2 score (>= 20 = MSI-high).
  Aneuploidy — Ben-David 2021 aneuploidy score.

We filter rows to default-entry only (one row per ModelID) and emit a
byCellLine dict that the cnService loads on CN-mode activation. Used to
rescale relative CN → ≈ N actual copies via measured ploidy.

Output: greenlistedv2/globalSignatures.json
"""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, 'OmicsGlobalSignatures.csv')
OUT  = '/Users/fredrikwermeling/Documents/greenlistedv2/globalSignatures.json'
FIELDS = ['MSIScore', 'LoHFraction', 'WGD', 'CIN', 'Ploidy', 'Aneuploidy']

def parse_float(s):
    if s is None or s == '' or s.upper() == 'NA': return None
    try: return float(s)
    except ValueError: return None

def main():
    if not os.path.exists(SRC):
        print(f'OmicsGlobalSignatures.csv not found at {SRC}', file=sys.stderr); sys.exit(1)
    by_cl = {}
    n_seen = n_default = 0
    with open(SRC, newline='') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            n_seen += 1
            if row.get('IsDefaultEntryForModel') != 'Yes': continue
            n_default += 1
            mid = row['ModelID']
            entry = {}
            for fld in FIELDS:
                v = parse_float(row.get(fld, ''))
                if v is None: continue
                if fld == 'WGD':       entry[fld] = bool(int(v))
                elif fld == 'Aneuploidy': entry[fld] = int(v)
                else: entry[fld] = round(v, 3)
            if entry: by_cl[mid] = entry
    out = {
        '_description': 'DepMap OmicsGlobalSignatures — per-cell-line genome-wide features. Used by the CN service to rescale DepMap relative copy number into estimated actual copy counts (via measured Ploidy). Filtered to default-entry rows (one per ModelID).',
        'schemaVersion': 1,
        'fields': FIELDS,
        'byCellLine': by_cl
    }
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    size_kb = os.path.getsize(OUT) / 1024
    n_with_ploidy = sum(1 for v in by_cl.values() if 'Ploidy' in v)
    print(f'  scanned {n_seen} rows, kept {n_default} default-entry rows')
    print(f'  cell lines emitted: {len(by_cl)}, with measured ploidy: {n_with_ploidy}')
    # Coverage check vs CN matrix.
    cn_path = '/Users/fredrikwermeling/Documents/greenlistedv2/cn_metadata.json'
    if os.path.exists(cn_path):
        cn = json.load(open(cn_path))
        cn_set = set(cn['cellLines'])
        covered = sum(1 for id in cn_set if id in by_cl and 'Ploidy' in by_cl[id])
        print(f'  CN coverage: {covered} / {len(cn_set)} lines now have measured ploidy')
    print(f'  wrote {OUT} ({size_kb:.0f} KB)')

if __name__ == '__main__':
    main()
