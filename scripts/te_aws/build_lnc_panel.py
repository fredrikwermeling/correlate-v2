#!/usr/bin/env python3
"""Build the non-coding panel: lncRNAs the app cannot see today.

DepMap's expression matrix is protein-coding only, so every lncRNA is invisible
to Correlate - including XIST, whose absence is why 506 of 1,208 cell lines
carry an "unknown" sex call. This panel is measured by the same streaming
quantifier as the retroelement panel.

Curated, not comprehensive, on purpose: lncRNA annotation is unstable between
databases, most lncRNAs are near-zero in any given line, and CCLE libraries are
poly-A selected. A few dozen abundant, well-annotated genes are trustworthy; a
genome-wide dump would be mostly noise.

Intervals are MERGED EXONS, not gene spans, so intronic and overlapping
transcription is excluded. A spliced read crossing two exons is counted in
both, which inflates a gene's total by a constant factor; since every
comparison the browser makes is between cell lines on the same gene, that
cancels.

Usage: python3 build_lnc_panel.py refGene.txt.gz lnc_panel.json
"""
import gzip, json, sys
from collections import defaultdict

GENES = {
    # gene: why it is in the panel
    'XIST':      'X-inactivation; absence from the protein-coding matrix is why 506 lines have an unknown sex call',
    'TSIX':      'antisense to XIST; with unstranded libraries it cannot be separated from XIST, kept to make that explicit',
    'TERC':      'telomerase RNA component; the app already carries TERT, ATRX, DAXX, POT1',
    'MALAT1':    'abundant nuclear lncRNA, metastasis-associated',
    'NEAT1':     'paraspeckle scaffold; induced by innate immune signalling',
    'NORAD':     'genomic stability, PUMILIO sequestration',
    'H19':       'imprinted, IGF2 locus',
    'MEG3':      'imprinted tumour suppressor, DLK1-DIO3 locus',
    'MEG8':      'DLK1-DIO3 locus',
    'KCNQ1OT1':  'imprinted, 11p15 locus',
    'HOTAIR':    'HOX locus, PRC2-associated',
    'PVT1':      'MYC locus; co-amplified with MYC',
    'CDKN2B-AS1':'ANRIL, 9p21 locus the wiki already discusses',
    'GAS5':      'growth arrest, glucocorticoid receptor decoy',
    'CYTOR':     'LINC00152, proliferation',
    'SNHG1':     'snoRNA host gene',
    'SNHG5':     'snoRNA host gene',
    'SNHG16':    'snoRNA host gene',
    'DANCR':     'differentiation-associated',
    'CCAT1':     '8q24 MYC enhancer region',
    'CCAT2':     '8q24 MYC enhancer region',
    'UCA1':      'bladder cancer-associated',
    'TUG1':      'taurine-upregulated, PRC2-associated',
    'FIRRE':     'nuclear architecture, X-linked',
    'FENDRR':    'FOXF1 adjacent',
    'RMRP':      'RNase MRP RNA; mutated in cartilage-hair hypoplasia',
    'RPPH1':     'RNase P RNA',
    'VTRNA1-1':  'vault RNA',
    'MIR155HG':  'miR-155 host; immune activation',
    'MIR17HG':   'miR-17-92 host; oncogenic cluster',
    'MIR31HG':   'miR-31 host; 9p21',
    'LINC00473': 'CREB-regulated',
    'HOTTIP':    'HOXA locus',
    'HOXA-AS2':  'HOXA locus',
    'DLEU1':     '13q14 deletion region, CLL',
    'DLEU2':     '13q14, miR-15a/16-1 host',
    'SAMMSON':   'melanoma-specific dependency',
    'OIP5-AS1':  'cyrano, widely expressed',
    'ZFAS1':     'snoRNA host, widely expressed',
    'LUCAT1':    'smoke and oxidative stress induced',
    'LINC01133': 'TGF-beta associated',
    'CRNDE':     'colorectal neoplasia differentially expressed',
    'HULC':      'liver cancer-associated',
    'MIAT':      'myocardial infarction-associated, neural',
    'NRIR':      'negative regulator of the interferon response',
    'EPCAM-DT':  'EPCAM divergent transcript',
}

def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else 'refGene.txt.gz'
    out = sys.argv[2] if len(sys.argv) > 2 else 'lnc_panel.json'
    exons = defaultdict(list)      # gene -> [(chrom,strand,s,e)]
    coding = defaultdict(list)     # chrom -> [(s,e,name)] protein-coding spans
    with gzip.open(ref, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            acc, chrom, strand, name2 = f[1], f[2], f[3], f[12]
            if '_' in chrom: continue
            if acc.startswith('NM_'):
                coding[chrom].append((int(f[4]), int(f[5]), name2))
            if name2 not in GENES: continue
            starts = [int(x) for x in f[9].rstrip(',').split(',') if x]
            ends   = [int(x) for x in f[10].rstrip(',').split(',') if x]
            for s, e in zip(starts, ends):
                exons[name2].append((chrom, strand, s, e))

    panel, report = [], []
    for gene in sorted(exons):
        rows = exons[gene]
        chrom = rows[0][0]; strand = rows[0][1]
        iv = sorted((s, e) for c, st, s, e in rows if c == chrom)
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
            else: merged.append([s, e])
        # Does any merged exon sit inside a protein-coding gene? If so its
        # counts are contaminated by that gene and must be read with care.
        ov = set()
        for s, e in merged:
            for cs, ce, cname in coding.get(chrom, []):
                if cs <= e and ce >= s: ov.add(cname)
        bp = sum(e - s for s, e in merged)
        for i, (s, e) in enumerate(merged, 1):
            panel.append({'chrom': chrom.replace('chr', ''), 'start': s, 'end': e,
                          'name': gene, 'family': 'lncRNA', 'strand': strand,
                          'id': f'LNC_{gene}_e{i:02d}', 'role': 'lncrna_exon',
                          'gene': gene,
                          'overlaps_coding': sorted(ov) or None})
        report.append((gene, chrom, len(merged), bp, sorted(ov)))

    json.dump(panel, open(out, 'w'), indent=0)
    print(f"{len(panel)} intervals over {len(report)} genes -> {out}\n")
    print(f"{'gene':14s} {'chrom':7s} {'exons':>5s} {'bp':>7s}  overlaps protein-coding")
    for g, c, n, bp, ov in report:
        flag = ', '.join(ov[:3]) + ('...' if len(ov) > 3 else '') if ov else '-'
        print(f"{g:14s} {c:7s} {n:5d} {bp:7d}  {flag}")
    miss = [g for g in GENES if g not in exons]
    if miss: print("\nNOT FOUND in this annotation:", ', '.join(miss))

if __name__ == '__main__':
    main()
