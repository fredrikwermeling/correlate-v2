#!/usr/bin/env python3
"""
Cell-line usage and availability -> web_data/cellline_usage.json

Usage:  python3 scripts/build_cellline_usage.py <cellosaurus.txt> [Model26Q1.csv] [--offline]
        cellosaurus.txt from https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt

For every DepMap model with an RRID:
  - distributors: cell-line collections that sell it now (Cellosaurus DR lines,
    minus entries flagged "Discontinued" in CC lines)
  - pubs: Europe PMC hit count for papers that write "<name> cell(s)" for any
    of the line's name variants. Counts are cached in scripts/.cache so
    reruns only query new lines.

Name variants: the Cellosaurus ID plus synonyms that collapse to the same
alphanumeric form, a shorter prefix of it (LNCaP for "LNCaP clone FGC") or a
suffix with letters and digits (H1581 for NCI-H1581). Variant suffixes (MCF7/WT) and
unrelated aliases are dropped.
"""
from concurrent.futures import ThreadPoolExecutor
import csv, json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
CELLO = ARGS[0] if ARGS else '/tmp/cellosaurus.txt'
MODEL = ARGS[1] if len(ARGS) > 1 else os.path.join(ROOT, 'Model26Q1.csv')
OUT = os.path.join(ROOT, 'web_data', 'cellline_usage.json')
CACHE = os.path.join(HERE, '.cache', 'europepmc_counts.json')

# Public cell-line collections and vendors that sell the line itself.
COLLECTIONS = {
    'ATCC': 'ATCC', 'DSMZ': 'DSMZ', 'ECACC': 'ECACC', 'JCRB': 'JCRB', 'RCB': 'RIKEN',
    'KCLB': 'KCLB', 'CLS': 'CLS', 'BCRJ': 'BCRJ', 'BCRC': 'BCRC', 'TKG': 'TKG',
    'Coriell': 'Coriell', 'NCI-DTP': 'NCI-DTP', 'IZSLER': 'IZSLER', 'CCTCC': 'CCTCC',
    'Ximbio': 'Ximbio', 'Kerafast': 'Kerafast', 'AddexBio': 'AddexBio', 'Lonza': 'Lonza',
}


# Names whose "<name> cells" search mostly finds something else: ordinary
# words or abbreviations (WILL-1 reads as "will 1 cells", SEM as scanning EM,
# HEL / HLF as human embryonic / lung fibroblasts), or another line (PK-1 as
# LLC-PK1). Checked by hand against the joined spelling; left uncounted.
AMBIGUOUS = {'WILL-1', 'TEN', 'SEM', 'HEL', 'HLF', 'GCT', 'CMS', 'DEL', 'MON', 'VAL',
             'BEN', 'SCH', 'SKN', 'CII', 'HTK-', 'MTA', 'PK-1', 'DOK', 'NMB', 'HEY'}
# Cellosaurus adds a bracketed note to a name another line also uses, so its
# count would mix the two; C32 is kept, the melanoma line is the usual one.
SHARED_OK = {'C32 [Human melanoma]'}


def strip(s):
    return re.sub(r'[^A-Z0-9]', '', s.upper())


def parse_cellosaurus(path, wanted):
    rec = None
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            tag, val = line[:2], line[5:].rstrip('\n')
            if tag == 'ID':
                rec = {'id': val, 'sy': [], 'dr': [], 'disc': set(), 'ac': None}
            elif rec is None:
                continue
            elif tag == 'AC':
                rec['ac'] = val.strip()
            elif tag == 'SY':
                rec['sy'] += [s.strip() for s in val.split(';') if s.strip()]
            elif tag == 'DR':
                db, _, acc = val.partition(';')
                rec['dr'].append((db.strip(), acc.strip()))
            elif tag == 'CC' and val.startswith('Discontinued:'):
                m = re.match(r'Discontinued: ([^;]+); ([^;.]+)', val)
                if m:
                    rec['disc'].add((m.group(1).strip(), m.group(2).strip()))
            elif line.startswith('//'):
                if rec['ac'] in wanted:
                    yield rec
                rec = None


def name_variants(main, synonyms):
    m = strip(main)
    keep = {main}
    base = re.match(r'^(.*?)[ ,]+clone\b', main, re.I)
    if base:
        keep.add(base.group(1))
    for s in synonyms:
        t = strip(s)
        if not t or '/' in s:
            continue
        if t == m or (len(t) >= 4 and m.startswith(t)) or \
           (len(t) >= 4 and m.endswith(t) and re.search(r'\d', t) and re.search(r'[A-Z]', t)):
            keep.add(s)
    # Europe PMC splits on punctuation: "MCF-7" also finds "MCF 7"; the joined
    # form needs its own term.
    out = {}
    for k in keep:
        k = re.sub(r'\s*\[.*?\]', '', k)          # Cellosaurus disambiguators
        k = re.sub(r'[^A-Za-z0-9.\-_ ]+', ' ', k).strip()
        for v in (k, re.sub(r'[\s\-.:_]', '', k)):
            if len(strip(v)) >= 3:
                # Search is case-insensitive and splits on punctuation.
                out.setdefault(re.sub(r'[\s\-.:_]+', ' ', v).upper().strip(), v)
    return sorted(out.values(), key=len)


def epmc_count(query, tries=4):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode(
        {'query': query, 'format': 'json', 'pageSize': 1, 'resultType': 'idlist'})
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.load(r)
            if 'hitCount' in d:
                return d['hitCount']
            print('  error', d.get('errMsg'), query[:80])
            return None
        except Exception as e:
            print('  retry', i + 1, e)
            time.sleep(3 * (i + 1))
    return None


def build_query(variants):
    # Name must sit next to "cell(s)" / "cell line" somewhere, which drops gene,
    # enzyme and catalog-number homonyms (PC3 the convertase, HEL the antigen).
    if not variants:
        return None
    q = ' OR '.join('"%s cells" OR "%s cell"' % (v, v) for v in variants)
    if len(q) > 1400:   # API limit is 1500 characters
        q = ' OR '.join('"%s cells"' % v for v in variants)[:1400].rsplit(' OR ', 1)[0]
    return q


def main():
    rows = list(csv.DictReader(open(MODEL, encoding='utf-8')))
    by_rrid = {}
    for r in rows:
        if r.get('RRID'):
            by_rrid.setdefault(r['RRID'], []).append(r['ModelID'])
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)

    recs = list(parse_cellosaurus(CELLO, set(by_rrid)))
    print(len(recs), 'Cellosaurus records for', len(by_rrid), 'RRIDs')
    queries = {rec['ac']: build_query(name_variants(rec['id'], rec['sy'])) for rec in recs}
    todo = sorted({q for q in queries.values() if q and cache.get(q) is None})
    if '--offline' in sys.argv:   # rebuild the JSON from the cache only
        todo = []
    print(len(todo), 'Europe PMC queries to run')
    # Each full-text phrase search takes 10-60 s server side; a few in parallel.
    with ThreadPoolExecutor(6) as pool:
        for n, (q, c) in enumerate(zip(todo, pool.map(epmc_count, todo))):
            cache[q] = c
            if n % 20 == 0:
                json.dump(cache, open(CACHE, 'w'))
                print(n, c, q[:60], flush=True)

    out = {}
    for rec in recs:
        dist = {}
        for db, acc in rec['dr']:
            if db in COLLECTIONS and (db, acc) not in rec['disc']:
                dist.setdefault(COLLECTIONS[db], acc)
        unclear = rec['id'] in AMBIGUOUS or ('[' in rec['id'] and rec['id'] not in SHARED_OK)
        pubs = None if unclear else cache.get(queries[rec['ac']])
        for mid in by_rrid[rec['ac']]:
            out[mid] = {'pubs': pubs, 'dist': dist}
            if unclear:
                out[mid]['amb'] = 1
    json.dump(cache, open(CACHE, 'w'))

    json.dump({'source': 'Europe PMC hit counts and Cellosaurus distributors',
               'built': time.strftime('%Y-%m-%d'), 'lines': out},
              open(OUT, 'w'), separators=(',', ':'))
    print('wrote', OUT, len(out))


if __name__ == '__main__':
    main()
