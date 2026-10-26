#!/usr/bin/env python3
"""Turn a scan result into web_data/virus_rna.json for the app.

Usage:  python3 build_web_data.py counts_sample.json [counts_extra.json ...] [--genes counts_positive.json]

--genes: a rescan of the positive cell lines that recorded where on each
genome the reads fall (virus_scan.py 'pos'). Reads are then assigned to the
annotated genes (refs/features.json) and each gene gets an approximate TPM:
its share of the unmapped reads sampled, times the unmapped share of the
file (bytes), times 1e6, scaled by 2,000 / gene length (2 kb as the average
expressed transcript length). Approximate, within about twofold.

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
    args = sys.argv[1:]
    gsrc = None
    if '--genes' in args:
        i = args.index('--genes'); gsrc = args[i + 1]; args = args[:i] + args[i + 2:]
    res = [r for src in args for r in json.load(open(src)) if not r.get('error')]
    feats = json.load(open(os.path.join(HERE, 'refs', 'features.json')))
    pos_by = {r['model_id']: r for r in json.load(open(gsrc)) if not r.get('error')} if gsrc else {}
    seen = {}  # virus -> feature indices that received reads anywhere
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
            g = pos_by.get(r['model_id'])
            if det and g and g.get('pos', {}).get(key) and g.get('file_bytes') and g.get('unmapped_start') is not None:
                pb, ns = g['posbin'], max(1, g['unmapped_scanned'])
                ushare = (g['file_bytes'] - g['unmapped_start']) / g['file_bytes']
                frac = lambda reads: reads / ns * ushare
                fl = feats.get(key, [])
                per = {}
                for b, c in g['pos'][key]:
                    x = b * pb + pb // 2 + 1                       # 1-based bin centre
                    for fi, f in enumerate(fl):
                        if any(a <= x <= e for a, e in f['seg']):
                            per[fi] = per.get(fi, 0) + c
                genes = []
                for fi, c in sorted(per.items(), key=lambda t: -t[1])[:15]:
                    L = sum(e - a + 1 for a, e in fl[fi]['seg'])
                    genes.append([fi, c, round(frac(c) * 1e6 * 2000 / max(L, 200), 1)])
                    seen.setdefault(key, set()).add(fi)
                glen = max(1, (len(open(os.path.join(HERE, 'refs', f'{key}.fa')).read().split('\n', 1)[1].replace('\n', ''))))
                nb = 100
                cov = [0] * nb
                for b, c in g['pos'][key]:
                    cov[min(nb - 1, (b * pb) * nb // glen)] += c
                v[key].update({'genes': genes, 'cov': cov, 'rpmAll': round(frac(g['hits'].get(key, k)) * 1e6, 1)})
        lines[r['model_id']] = {'scanned': n, 'v': v}
    out = {
        '_description': 'Viral RNA in CCLE RNA-seq: reads that did not map to the human genome, matched by 25-mers to 18 viral (and one mycoplasma) genomes. scanned = unmapped reads examined; per virus: reads, rpm (reads per million unmapped reads scanned), breadth (share of 50 genome bins covered), detected (reads >= 5, rpm >= 10, breadth >= 0.2), level (strong: rpm >= 1000, moderate: >= 300, low: below; low can be carry-over from another library on the same run). Undetected viruses are listed only with 3+ reads over 3+ bins.',
        'mode': mode,
        'modeNote': 'A sample: the last 25 MB of each file, about 330,000 unmapped reads per cell line. Strongly expressed viruses are found; a weak signal can be missed.' if mode == 'sample' else 'Every unmapped read in the file.',
        'source': 'Computed for Correlate from the public DepMap CCLE RNA-seq BAMs, hg19 and hg38 (s3://depmap-omics-ccle); references from NCBI RefSeq.',
        'built': date.today().isoformat(),
        'thresholds': {'reads': MIN_READS, 'rpm': MIN_RPM, 'breadth': MIN_BREADTH},
        'viruses': {r['key']: {'label': r['label'], 'acc': r['acc'], **({'note': r['note']} if r.get('note') else {})} for r in refs},
        # Annotated genes that received reads somewhere, indexed as in 'genes'.
        'features': {k: {str(fi): {'name': feats[k][fi]['name'], 'product': feats[k][fi]['product'], 'seg': feats[k][fi]['seg'], 'strand': feats[k][fi]['strand']} for fi in sorted(ix)} for k, ix in seen.items()},
        'genomeLength': {r['key']: len(open(os.path.join(HERE, 'refs', f"{r['key']}.fa")).read().split('\n', 1)[1].replace('\n', '')) for r in refs},
        'tpmNote': 'genes: [feature index, reads, approximate TPM]. TPM is estimated from the share of all reads (unmapped share of the file in bytes times the share of sampled unmapped reads) and a 2 kb average transcript length; approximate, within about twofold. rpmAll: virus reads per million of all reads. cov: reads in 100 equal stretches of the genome.',
        'byCellLine': lines,
    }
    dst = os.path.join(HERE, '..', '..', 'web_data', 'virus_rna.json')
    json.dump(out, open(dst, 'w'), separators=(',', ':'))
    det = sum(1 for l in lines.values() if any(x['detected'] for x in l['v'].values()))
    print(f'{len(lines)} cell lines, {det} with a detected virus -> {os.path.normpath(dst)}')


if __name__ == '__main__':
    main()
