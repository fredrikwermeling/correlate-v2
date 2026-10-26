#!/usr/bin/env python3
"""Turn a scan result into web_data/virus_rna.json for the app.

Usage:  python3 build_web_data.py counts_sample.json [counts_extra.json ...]

Per cell line it keeps the reads scanned and every virus with 3 or more
reads touching 3 or more of the 50 genome bins (a one-bin hit is a shared
or conserved stretch). A virus counts as DETECTED with at least 5 reads, at
least 10 per million unmapped and breadth 0.2 or more (the rule in
analyze.py); below that it is a trace. A detection is graded by reads per
million unmapped: strong (1,000+), moderate (300+) or low. The cut-offs
come from the first full screen: real infections sit far above 1,000, and a
continuous haze of 10 to 250 per million follows the strongest producers. A real infection
is usually strong (22Rv1 XMRV ~350,000); low levels can be carry-over from
another library on the same sequencing run.
"""
import json
import os
import sys
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_READS, MIN_RPM, MIN_BREADTH, KEEP = 5, 10, 0.2, 3


def main():
    res = [r for src in sys.argv[1:] for r in json.load(open(src)) if not r.get('error')]
    refs = json.load(open(os.path.join(HERE, 'references.json')))['references']
    mode = res[0].get('mode', 'sample') if res else 'sample'
    lines = {}
    for r in res:
        n = r.get('unmapped_scanned', 0)
        v = {}
        for key, k in r.get('hits', {}).items():
            rpm = k * 1e6 / max(1, n)
            b = r['breadth'].get(key, 0)
            det = k >= MIN_READS and rpm >= MIN_RPM and b >= MIN_BREADTH
            if not det and (k < KEEP or b < 0.06):
                continue
            v[key] = {'reads': k, 'rpm': round(rpm, 1), 'breadth': b, 'detected': det}
            if det:
                v[key]['level'] = 'strong' if rpm >= 1000 else 'moderate' if rpm >= 300 else 'low'
        lines[r['model_id']] = {'scanned': n, 'v': v}
    out = {
        '_description': 'Viral RNA in CCLE RNA-seq: reads that did not map to the human genome, matched by 25-mers to 18 viral (and one mycoplasma) genomes. scanned = unmapped reads examined; per virus: reads, rpm (reads per million unmapped reads scanned), breadth (share of 50 genome bins covered), detected (reads >= 5, rpm >= 10, breadth >= 0.2), level (strong: rpm >= 1000, moderate: >= 300, low: below; low can be carry-over from another library on the same run). Undetected viruses are listed only with 3+ reads over 3+ bins.',
        'mode': mode,
        'modeNote': 'A sample: the last 25 MB of each file, about 330,000 unmapped reads per cell line. Strongly expressed viruses are found; a weak signal can be missed.' if mode == 'sample' else 'Every unmapped read in the file.',
        'source': 'Computed for Correlate from the public DepMap CCLE RNA-seq BAMs, hg19 and hg38 (s3://depmap-omics-ccle); references from NCBI RefSeq.',
        'built': date.today().isoformat(),
        'thresholds': {'reads': MIN_READS, 'rpm': MIN_RPM, 'breadth': MIN_BREADTH},
        'viruses': {r['key']: {'label': r['label'], 'acc': r['acc'], **({'note': r['note']} if r.get('note') else {})} for r in refs},
        'byCellLine': lines,
    }
    dst = os.path.join(HERE, '..', '..', 'web_data', 'virus_rna.json')
    json.dump(out, open(dst, 'w'), separators=(',', ':'))
    det = sum(1 for l in lines.values() if any(x['detected'] for x in l['v'].values()))
    print(f'{len(lines)} cell lines, {det} with a detected virus -> {os.path.normpath(dst)}')


if __name__ == '__main__':
    main()
