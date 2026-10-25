#!/usr/bin/env python3
"""
Build web_data/complexes.json: named protein complexes (EBI Complex Portal)
and named Reactome pathways, each with its members mapped to the app's gene
names, for the Gene Browser's "Protein complexes and pathways" section.

    python3 scripts/build_complexes.py [workdir]

Downloads (into workdir, default a temp folder):
  Complex Portal human complexes  ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv
  UniProt accession -> gene       rest.uniprot.org (human, all entries)
  Reactome pathways               reactome.org/download/current/ReactomePathways.gmt.zip

Complex members are the expanded participant list (sub-complexes resolved to
proteins); non-protein participants (CHEBI, RNA) are dropped. Reactome keeps
pathways of 5 to 150 genes: bigger ones are broad parents ("Signal
Transduction") whose co-membership says little.
"""
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web_data')
CPX_URL = 'https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv'
UP_URL = 'https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606&fields=accession,gene_primary&format=tsv'
RX_URL = 'https://reactome.org/download/current/ReactomePathways.gmt.zip'
RX_MIN, RX_MAX = 5, 150
EVIDENCE = {'ECO:0000353': 'experimental'}  # everything else is curator inference (orthology, paralogy, background knowledge)


def fetch(url, path):
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={'User-Agent': 'correlate-complexes/1.0'})
        with urllib.request.urlopen(req, timeout=600) as r, open(path, 'wb') as f:
            f.write(r.read())
    return path


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp()
    syn = json.load(open(os.path.join(WEB, 'synonyms.json')))
    data = set()
    for f in ('metadata.json', 'expression_metadata.json'):
        data |= {g.upper() for g in json.load(open(os.path.join(WEB, f)))['genes']}

    def to_app(sym):
        up = sym.upper()
        if up in data:
            return up
        t = (syn.get(up) or {}).get('d', '').upper()
        return t if t in data else None

    up = {}
    for line in open(fetch(UP_URL, os.path.join(work, 'up.tsv'))).read().splitlines()[1:]:
        c = line.split('\t')
        if len(c) >= 2 and c[1].strip():
            up[c[0]] = c[1].split(';')[0].strip()

    complexes = []
    rows = open(fetch(CPX_URL, os.path.join(work, '9606.tsv')), encoding='utf-8').read().splitlines()
    head = rows[0].lstrip('#').split('\t')
    ix = {h: i for i, h in enumerate(head)}
    for line in rows[1:]:
        c = line.split('\t')
        members = []
        for part in c[ix['Expanded participant list']].split('|'):
            acc = re.sub(r'\(.*\)$', '', part).strip()
            acc = re.sub(r'[-_].*$', '', acc)  # isoform / chain -> entry
            g = up.get(acc)
            g = to_app(g) if g else None
            if g and g not in members:
                members.append(g)
        if len(members) < 2:
            continue
        dis = []
        for d in c[ix['Disease']].split('|'):
            m = re.match(r'\s*([^\[:]+?)\s*(\[|:)', d)
            if m and len(m.group(1)) < 70:
                dis.append(m.group(1))
        aliases = [a for a in c[ix['Aliases for complex']].split('|') if a and a != '-'][:3]
        desc = c[ix['Description']].strip()
        desc = '' if desc == '-' else re.sub(r'\s+', ' ', desc)
        go = [m.group(1) for m in re.finditer(r'GO:\d+\(([^)]+)\)', c[ix['Go Annotations']])]
        complexes.append({
            'id': c[ix['Complex ac']], 'name': c[ix['Recommended name']], 'aka': aliases,
            'desc': desc, 'go': go[:6], 'disease': dis[:4],
            'ev': EVIDENCE.get(c[ix['Evidence Code']].split('(')[0], 'inferred'),
            'genes': members,
        })
    complexes.sort(key=lambda x: x['id'])

    pathways = []
    with zipfile.ZipFile(fetch(RX_URL, os.path.join(work, 'rx.zip'))) as zf:
        gmt = zf.read([n for n in zf.namelist() if n.endswith('.gmt')][0]).decode('utf-8')
    for line in gmt.splitlines():
        c = line.split('\t')
        genes = []
        for g in c[2:]:
            a = to_app(g)
            if a and a not in genes:
                genes.append(a)
        if RX_MIN <= len(genes) <= RX_MAX:
            pathways.append({'id': c[1], 'name': c[0], 'genes': genes})

    out = {
        '_doc': 'Named protein complexes (EBI Complex Portal, human) and Reactome pathways of '
                f'{RX_MIN}-{RX_MAX} genes, members mapped to the app gene names. Built by scripts/build_complexes.py.',
        'built': time.strftime('%Y-%m-%d'),
        'complexes': complexes, 'pathways': pathways,
    }
    path = os.path.join(WEB, 'complexes.json')
    with open(path, 'w') as f:
        json.dump(out, f, separators=(',', ':'), ensure_ascii=False)
    print(f'{len(complexes)} complexes, {len(pathways)} pathways, {os.path.getsize(path) / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
