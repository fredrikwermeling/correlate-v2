#!/usr/bin/env python3
"""
Build web_data/common_essentials.json from DepMap's CRISPRInferredCommonEssentials.csv
(one "SYMBOL (entrez)" per line under an "Essentials" header). Sorted symbol list.

    python3 scripts/build_common_essentials.py [26Q1/CRISPRInferredCommonEssentials.csv]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, '26Q1', 'CRISPRInferredCommonEssentials.csv')
lines = [l.strip() for l in open(src).read().splitlines()[1:] if l.strip()]
genes = sorted({l.split(' (')[0].strip() for l in lines})
with open(os.path.join(ROOT, 'web_data', 'common_essentials.json'), 'w') as f:
    json.dump(genes, f)
print(f'{len(genes)} common essentials from {src}')
