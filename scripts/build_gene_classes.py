#!/usr/bin/env python3
"""
Build web_data/gene_classes.json from Green Listed's curated gene lists
(~/Documents/greenlistedv2/geneSets.json, built and audited there by
tools/build_gene_sets.py), so both apps offer the same classes.

    python3 scripts/build_gene_classes.py [path/to/geneSets.json]

Only the human lists are used. Symbols are mapped onto the names the app's
matrices use through web_data/synonyms.json; genes the data lacks are dropped.
Oncogenes and tumour suppressors are split here from OncoKB's own gene list
(its Gene Type column; a gene typed as both goes in both), placed right after
Green Listed's Cancer genes, which come from the same list.
"""
import json
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web_data')
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/Documents/greenlistedv2/geneSets.json')
ONCOKB = 'https://www.oncokb.org/api/v1/utils/cancerGeneList.txt'


def oncokb_roles():
    req = urllib.request.Request(ONCOKB, headers={'User-Agent': 'correlate-gene-classes/1.0'})
    rows = urllib.request.urlopen(req, timeout=120).read().decode('utf-8', 'replace').splitlines()
    head = rows[0].split('\t')
    si, ti = head.index('Hugo Symbol'), head.index('Gene Type')
    onc, tsg = [], []
    for r in rows[1:]:
        c = r.split('\t')
        if len(c) <= ti:
            continue
        t = c[ti].strip().upper()
        if 'ONCOGENE' in t:
            onc.append(c[si])
        if 'TSG' in t:
            tsg.append(c[si])
    return [
        {'key': 'oncogene', 'label': 'Oncogenes', 'description': 'Cancer genes OncoKB types as oncogenes (including those it also types as tumor suppressors).', 'source': 'OncoKB Cancer Gene List, Gene Type', 'genes': onc},
        {'key': 'tsg', 'label': 'Tumor suppressors', 'description': 'Cancer genes OncoKB types as tumor suppressors (including those it also types as oncogenes).', 'source': 'OncoKB Cancer Gene List, Gene Type', 'genes': tsg},
    ]


def main():
    gl = json.load(open(SRC))
    syn = json.load(open(os.path.join(WEB, 'synonyms.json')))
    data = set()
    for f in ('metadata.json', 'expression_metadata.json'):
        data |= {g.upper() for g in json.load(open(os.path.join(WEB, f)))['genes']}
    out = {'_doc': "Curated gene classes from Green Listed's geneSets.json (human), mapped to the app's gene names. Built by scripts/build_gene_classes.py; do not edit by hand.",
           'classes': {}}
    sets = []
    for s in gl['species']['Human']:
        sets.append(s)
        if s['key'] == 'cancer':
            sets += oncokb_roles()
    for s in sets:
        genes = set()
        for g in s['genes']:
            up = g.upper()
            t = up if up in data else (syn.get(up) or {}).get('d', '').upper()
            if t in data:
                genes.add(t)
        out['classes'][s['key']] = {'label': s['label'], 'description': s.get('description', ''),
                                    'source': s.get('source', ''), 'n': len(genes), 'genes': sorted(genes)}
        print(f"{s['key']:12s} {len(s['genes']):5d} listed, {len(genes):5d} in the data")
    with open(os.path.join(WEB, 'gene_classes.json'), 'w') as f:
        json.dump(out, f, separators=(',', ':'))


if __name__ == '__main__':
    main()
