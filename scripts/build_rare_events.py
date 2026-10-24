#!/usr/bin/env python3
"""
Build web_data/rare_events.json: genes with a damaging mutation or a fusion in
only 1-2 cell lines, which damaging_mutations.json and translocations.json leave
out (they keep genes with at least 3, the minimum for a mutated-vs-rest group).
Same rules as those builders: app cohort only, default profile per model for
mutations, every fusion row for fusions.

    python3 scripts/build_rare_events.py

Output: {"_doc", "damaging": {GENE: [ACH...]}, "fusion": {GENE: {ACH: [partners]}}}
"""
import json
import os
import re
from collections import defaultdict

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
meta = json.load(open(os.path.join(ROOT, 'web_data', 'metadata.json')))
cohort = set(meta['cellLines'])
# Only genes the app can show a page for (screened or expressed).
known = set(meta['genes']) | set(json.load(open(os.path.join(ROOT, 'web_data', 'expression_metadata.json')))['genes'])

dm = pd.read_csv(os.path.join(ROOT, '26Q1', 'OmicsSomaticMutationsMatrixDamaging.csv'), low_memory=False)
dm = dm[(dm['IsDefaultEntryForModel'] == 'Yes') & dm['ModelID'].isin(cohort)].drop_duplicates('ModelID')
gene_cols = [c for c in dm.columns if re.match(r'^.+?\s*\(\d+\)$', c)]
vals = dm[gene_cols].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy() > 0
models = dm['ModelID'].to_numpy()
damaging = {}
for j, col in enumerate(gene_cols):
    hit = vals[:, j]
    k = int(hit.sum())
    if 1 <= k < 3:
        damaging.setdefault(re.match(r'^(.+?)\s*\(\d+\)$', col).group(1).strip(), []).extend(models[hit].tolist())

fu = pd.read_csv(os.path.join(ROOT, '26Q1', 'OmicsFusionFiltered.csv'), usecols=['ModelID', 'Gene1', 'Gene2'])
fu = fu[fu['ModelID'].isin(cohort)]
partners = defaultdict(lambda: defaultdict(set))
for m, a, b in fu.itertuples(index=False):
    a = str(a).split(' (')[0].strip(); b = str(b).split(' (')[0].strip()
    if not a or not b:
        continue
    partners[a][m].add(b); partners[b][m].add(a)
fusion = {g: {m: sorted(p) for m, p in d.items()} for g, d in partners.items() if 1 <= len(d) < 3}

out = {
    '_doc': 'Genes with a damaging mutation or a fusion in only 1-2 cell lines of the app cohort (DepMap 26Q1), which damaging_mutations.json and translocations.json leave out because they keep genes with at least 3. Mutations: default profile per model.',
    'damaging': {g: v for g, v in sorted(damaging.items()) if g in known},
    'fusion': {g: v for g, v in sorted(fusion.items()) if g in known},
}
with open(os.path.join(ROOT, 'web_data', 'rare_events.json'), 'w') as f:
    json.dump(out, f, separators=(',', ':'))
print(f"{len(damaging)} rare damaging genes, {len(fusion)} rare fusion genes")
