"""Build web_data/cytobands.json: GRCh38 G-banding for the Gene Browser ideogram.

Source: UCSC hg38 cytoBand.txt.gz
(https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz).
Usage: python3 scripts/build_cytobands.py path/to/cytoBand.txt.gz
Each band is [start_kb, end_kb, name, stain]; stain keeps UCSC's gieStain
(gneg, gpos25..gpos100, acen, gvar, stalk).
"""
import gzip, json, re, sys

src = sys.argv[1]
out = {}
with gzip.open(src, 'rt') as f:
    for line in f:
        chrom, start, end, name, stain = line.rstrip('\n').split('\t')
        m = re.fullmatch(r'chr(\d+|X|Y)', chrom)
        if not m:
            continue
        out.setdefault(m.group(1), []).append([round(int(start) / 1000), round(int(end) / 1000), name, stain])
for bands in out.values():
    bands.sort(key=lambda b: b[0])
doc = {
    '_doc': 'GRCh38 cytogenetic bands (UCSC hg38 cytoBand). Per chromosome: [start_kb, end_kb, band, gieStain]. Built by scripts/build_cytobands.py.',
    'assembly': 'GRCh38',
    'bands': out,
}
with open('web_data/cytobands.json', 'w') as f:
    json.dump(doc, f, separators=(',', ':'))
print(sum(len(v) for v in out.values()), 'bands on', len(out), 'chromosomes')
