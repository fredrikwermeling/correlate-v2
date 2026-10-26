"""
Whole-transcriptome identity check -> web_data/transcriptome_identity.json

For every CRISPR-cohort cell line with mRNA data: its five most similar cell
lines by whole-transcriptome profile, and whether any of them shares its
tissue (Oncotree lineage). Cell lines from the same patient (DepMap
PatientID) are never counted as neighbors, so a parental line and its
derivatives cannot vouch for each other.

Profile: the 5,000 most variable genes (log2 TPM+1), z-scored per gene;
similarity is the Pearson r between two cell lines' profiles. Tested
against 1,000 / 3,000 / 10,000 genes and Spearman: all within 1 point.

Status:
  match      2 or more of the 5 nearest share its tissue
  mixed      exactly 1 does, and 3 or more share one other tissue
  other      none does
  unchecked  tissue Unknown / Other, or fewer than 5 cell lines with mRNA

Run: python3 scripts/build_transcriptome_identity.py
"""
import csv
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from matrix_io import read_float  # noqa: E402

W = os.path.join(HERE, '..', 'web_data')
N_GENES, K = 5000, 5


def main():
    em = json.load(open(os.path.join(W, 'expression_metadata.json')))
    md = json.load(open(os.path.join(W, 'metadata.json')))
    cm = json.load(open(os.path.join(W, 'cellLineMetadata.json')))
    X = read_float(os.path.join(W, 'expression.bin.gz'), em)
    ecol = {c: i for i, c in enumerate(em['cellLines'])}
    ids = [c for c in md['cellLines'] if c in ecol]
    lin = [cm['lineage'].get(c) or 'Unknown' for c in ids]
    patient = {r['ModelID']: r['PatientID'] for r in csv.DictReader(open(os.path.join(HERE, 'Model_25Q3.csv')))}

    A = np.nan_to_num(X[:, [ecol[c] for c in ids]].T.astype(np.float64))
    top = np.argsort(-A.var(0))[:N_GENES]
    B = A[:, top]
    B = (B - B.mean(0)) / (B.std(0) + 1e-9)
    B -= B.mean(1, keepdims=True)
    B /= np.linalg.norm(B, axis=1, keepdims=True)
    S = B @ B.T
    P = np.array([patient.get(c) or c for c in ids])
    S[P[:, None] == P[None, :]] = -9  # self and same-patient lines
    nn = np.argsort(-S, 1)[:, :K]

    count = Counter(lin)
    out, tally = {}, Counter()
    for i, c in enumerate(ids):
        near = [ids[j] for j in nn[i]]
        tis = Counter(lin[j] for j in nn[i])
        own = tis.get(lin[i], 0)
        like, likeN = next(((t, n) for t, n in tis.most_common() if t != lin[i]), (None, 0))
        if lin[i] in ('Unknown', 'Other') or count[lin[i]] < 5:
            s = 'unchecked'
        elif own >= 2:
            s = 'match'
        elif own == 1:
            s = 'mixed' if likeN >= 3 else 'match'
        else:
            s = 'other'
        tally[s] += 1
        out[c] = {'s': s, 'nn': [[near[k], round(float(S[i, nn[i][k]]), 3)] for k in range(K)]}
        if s in ('mixed', 'other'):
            out[c]['like'] = like

    checked = sum(v for k, v in tally.items() if k != 'unchecked')
    res = {
        '_description': 'Whole-transcriptome identity check: for each CRISPR-cohort cell line with mRNA data, its five most similar cell lines (Pearson r over the 5,000 most variable genes, z-scored; cell lines from the same patient excluded) and whether they share its tissue. s: match (2+ of 5 share its tissue), mixed (1 does and 3+ share another tissue), other (none does), unchecked (tissue Unknown/Other or under 5 cell lines). like: the other tissue it most resembles.',
        'source': 'DepMap 26Q1 OmicsExpressionTPMLogp1; patient IDs from the DepMap Model table',
        'genes': N_GENES, 'k': K,
        'counts': dict(tally),
        'cohortOwnTissueRate': round((tally['match'] + tally['mixed']) / checked, 4),
        'byCellLine': out,
    }
    json.dump(res, open(os.path.join(W, 'transcriptome_identity.json'), 'w'), separators=(',', ':'))
    print(dict(tally), res['cohortOwnTissueRate'])


if __name__ == '__main__':
    main()
