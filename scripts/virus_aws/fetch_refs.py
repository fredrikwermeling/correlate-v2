#!/usr/bin/env python3
"""Download the reference genomes in references.json from NCBI into refs/.

Run once, on the laptop or the instance:  python3 fetch_refs.py
Skips files already present. Checks each header names the expected accession.
"""
import json
import os
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'


def main():
    refs = json.load(open(os.path.join(HERE, 'references.json')))['references']
    os.makedirs(os.path.join(HERE, 'refs'), exist_ok=True)
    for r in refs:
        path = os.path.join(HERE, 'refs', f"{r['key']}.fa")
        if os.path.exists(path) and os.path.getsize(path) > 1000:
            continue
        txt = urllib.request.urlopen(EFETCH.format(r['acc']), timeout=120).read().decode()
        if not txt.startswith('>' + r['acc']):
            raise SystemExit(f"{r['key']}: unexpected header {txt[:80]!r}")
        open(path, 'w').write(txt)
        print(f"{r['key']:6s} {r['acc']}  {len(txt):,} bytes")
        time.sleep(0.4)  # NCBI asks for at most 3 requests a second without a key
    print('references ready')


if __name__ == '__main__':
    main()
