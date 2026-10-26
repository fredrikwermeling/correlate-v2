#!/usr/bin/env python3
"""Download the gene annotation of each reference genome from NCBI into
refs/features.json: per reference, its genes (CDS, plus annotated RNAs such
as EBV's EBERs) with name, product and 1-based segments on the genome.

Run once after fetch_refs.py:  python3 fetch_features.py
"""
import json
import os
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=ft&retmode=text'
KEEP = {'CDS', 'ncRNA', 'misc_RNA', 'rRNA'}


def parse(txt):
    """NCBI five-column feature table -> list of features."""
    feats, cur, gene_names = [], None, {}
    for line in txt.splitlines():
        if not line or line.startswith('>'):
            continue
        cols = line.split('\t')
        if cols[0] and cols[1]:                      # a new feature, or its first segment
            try:
                a, b = int(cols[0].strip('<>^')), int(cols[1].strip('<>^'))
            except ValueError:
                continue
            if len(cols) > 2 and cols[2]:
                cur = {'type': cols[2], 'seg': [[a, b]], 'q': {}}
                feats.append(cur)
            elif cur:
                cur['seg'].append([a, b])            # a further segment (spliced)
        elif len(cols) >= 5 and cur:                 # qualifier
            cur['q'].setdefault(cols[3], cols[4])
    for f in feats:                                  # gene names by span, for CDS that lack one
        if f['type'] == 'gene' and f['q'].get('gene'):
            gene_names[(min(map(min, f['seg'])), max(map(max, f['seg'])))] = f['q']['gene']
    out = []
    for f in feats:
        if f['type'] not in KEEP:
            continue
        lo, hi = min(map(min, f['seg'])), max(map(max, f['seg']))
        name = f['q'].get('gene') or gene_names.get((lo, hi)) or f['q'].get('product') or f['q'].get('locus_tag') or f['type']
        out.append({'name': name, 'type': f['type'], 'product': f['q'].get('product', ''),
                    'strand': -1 if f['seg'][0][0] > f['seg'][0][1] else 1,
                    'seg': [[min(s), max(s)] for s in f['seg']]})
    return out


def main():
    refs = json.load(open(os.path.join(HERE, 'references.json')))['references']
    out = {}
    for r in refs:
        txt = urllib.request.urlopen(EFETCH.format(r['acc']), timeout=120).read().decode()
        out[r['key']] = parse(txt)
        print(f"{r['key']:6s} {len(out[r['key']])} features")
        time.sleep(0.4)
    os.makedirs(os.path.join(HERE, 'refs'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'refs', 'features.json'), 'w'))


if __name__ == '__main__':
    main()
