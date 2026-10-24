#!/usr/bin/env python3
"""
Build web_data/synonyms.json: every human gene name a user may type, mapped
to the name this app's data uses.

Sources (download fresh, pass the paths):
  NCBI  https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
  HGNC  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

    python3 scripts/build_synonyms.py Homo_sapiens.gene_info.gz hgnc_complete_set.txt

Keys are upper-cased: approved symbols, previous symbols, aliases (NCBI and
HGNC) and exact full names. Rules that keep it honest:
  - a name that is another gene's approved symbol is never an alias of this one;
  - a gene name the data itself uses always means that gene;
  - a key claimed by several genes goes to the one claiming it most strongly
    (symbol > previous symbol > alias / full name) and is dropped on a tie.
The target "d" is the name the data uses for the gene (its current or an older
symbol, whichever the matrices carry), else the current symbol.

Output: {KEY: {"d": TARGET, "r": "l"}}
"""
import csv
import gzip
import json
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web_data')
ncbi_path, hgnc_path = sys.argv[1], sys.argv[2]

data_genes = set()
for f in ('metadata.json', 'expression_metadata.json', 'cn_metadata.json'):
    data_genes |= {g.upper() for g in json.load(open(os.path.join(WEB, f)))['genes']}

split = lambda s, sep: [x.strip() for x in (s or '').split(sep) if x.strip() and x.strip() != '-']

# One record per gene, keyed by HGNC id when there is one, else NCBI GeneID.
genes = {}
hgnc_by_id = {}
with open(hgnc_path, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row.get('status') != 'Approved':
            continue
        hgnc_by_id[row['hgnc_id']] = row
with gzip.open(ncbi_path, 'rt', newline='') as f:
    rd = csv.reader(f, delimiter='\t')
    next(rd)
    for r in rd:
        if r[0] != '9606':
            continue
        xrefs = split(r[5], '|')
        hid = next((x[len('HGNC:'):] for x in xrefs if x.startswith('HGNC:HGNC:')), None)
        key = hid or 'NCBI:' + r[1]
        genes[key] = {'symbol': r[2], 'aliases': split(r[4], '|'), 'names': [r[8]] if r[8] != '-' else [],
                      'prev': [], 'hgnc': hgnc_by_id.get(hid)}
for hid, row in hgnc_by_id.items():
    g = genes.setdefault(hid, {'symbol': row['symbol'], 'aliases': [], 'names': [], 'prev': [], 'hgnc': row})
    g['hgnc'] = row
for g in genes.values():
    row = g['hgnc']
    if row:
        g['symbol'] = row['symbol']
        g['prev'] = split(row.get('prev_symbol', '').strip('"'), '|')
        g['aliases'] += split(row.get('alias_symbol', '').strip('"'), '|')
        g['names'] += [row['name']] + split(row.get('prev_name', '').strip('"'), '|')

approved = {g['symbol'].upper() for g in genes.values()}
# Keep HGNC-named genes and anything the data carries; NCBI's other records
# are mostly LOC / unnamed loci nobody types.
genes = {k: g for k, g in genes.items() if g['hgnc'] or g['symbol'].upper() in data_genes}

claims = defaultdict(list)  # key -> [(strength, target)]
for g in genes.values():
    sym = g['symbol'].upper()
    target = next((x for x in [sym] + [p.upper() for p in g['prev']] if x in data_genes), sym)
    for k, strength in ([(sym, 3)] + [(p, 2) for p in g['prev']]
                        + [(a, 1) for a in g['aliases']] + [(n, 1) for n in g['names']]):
        k = k.upper().strip()
        if not k:
            continue
        if strength < 3 and k in approved:
            continue
        claims[k].append((strength, target))

out = {}
dropped = 0
for k, c in claims.items():
    if k in data_genes:
        out[k] = {'d': k, 'r': 'l'}
        continue
    best = max(s for s, _ in c)
    targets = {t for s, t in c if s == best}
    if len(targets) == 1:
        out[k] = {'d': targets.pop(), 'r': 'l'}
    else:
        dropped += 1
for k in data_genes:
    out.setdefault(k, {'d': k, 'r': 'l'})

with open(os.path.join(WEB, 'synonyms.json'), 'w') as f:
    json.dump(dict(sorted(out.items())), f, separators=(',', ':'))
print(f'{len(genes)} genes, {len(out)} names, {dropped} ambiguous names dropped')
