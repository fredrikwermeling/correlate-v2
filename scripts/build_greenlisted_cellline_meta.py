#!/usr/bin/env python3
"""Rebuild Green Listed's cellLineMetadata.json from DepMap Model.csv.

The previous snapshot covered 1186 cell lines, but the current CN matrix
(OmicsCNGene 24Q4) carries 1929 lines, leaving 743 with no name / sex /
cancer-type annotation. They surfaced in the CN picker only as bare
ACH-XXXXXX DepMap IDs. DepMap publishes a complete annotation table in
Model.csv (one row per ModelID); reading the latest release and emitting
the slim shape Green Listed expects fills the gap.

Output: greenlistedv2/cellLineMetadata.json
Schema (one list + parallel dicts keyed by ACH-ModelID):
    cellLines, cellLineName, strippedCellLineName, lineage, primaryDisease,
    subtype, sex, age, ageCategory, patientRace, primaryOrMetastasis,
    sampleCollectionSite, growthPattern, patientTumorGrade,
    depmapModelType, oncotreeSubtype, oncotreeCode, patientSubtypeFeatures,
    rrid, engineeredModel, culturedResistanceDrug

`sexByExpression` and `strProfile` from the old V2 file are omitted —
they come from separate DepMap files (Y-chromosome expression / STR
profile tables) and Green Listed doesn't read them.
"""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_CSV = os.path.join(ROOT, 'Model26Q1.csv')
OUT = '/Users/fredrikwermeling/Documents/greenlistedv2/cellLineMetadata.json'

# Model.csv column → one or more output dict keys. A single source column
# can mirror into more than one destination (subtype + oncotreeSubtype are
# duplicates kept for backwards-compat with the old V2-shaped metadata).
COLS = [
    ('CellLineName',           ['cellLineName']),
    ('StrippedCellLineName',   ['strippedCellLineName']),
    ('OncotreeLineage',        ['lineage']),
    ('OncotreePrimaryDisease', ['primaryDisease']),
    ('OncotreeSubtype',        ['subtype', 'oncotreeSubtype']),
    ('Sex',                    ['sex']),
    ('Age',                    ['age']),
    ('AgeCategory',            ['ageCategory']),
    ('PatientRace',            ['patientRace']),
    ('PrimaryOrMetastasis',    ['primaryOrMetastasis']),
    ('SampleCollectionSite',   ['sampleCollectionSite']),
    ('GrowthPattern',          ['growthPattern']),
    ('PatientTumorGrade',      ['patientTumorGrade']),
    ('DepmapModelType',        ['depmapModelType']),
    ('OncotreeCode',           ['oncotreeCode']),
    ('PatientSubtypeFeatures', ['patientSubtypeFeatures']),
    ('RRID',                   ['rrid']),
    ('EngineeredModel',        ['engineeredModel']),
    ('CulturedResistanceDrug', ['culturedResistanceDrug']),
]

def main():
    if not os.path.exists(MODEL_CSV):
        print(f'Model file not found at {MODEL_CSV}', file=sys.stderr); sys.exit(1)
    print(f'reading {MODEL_CSV} ...')
    with open(MODEL_CSV, newline='') as f:
        rdr = csv.DictReader(f)
        header_set = set(rdr.fieldnames)
        missing = [src for src, _ in COLS if src not in header_set]
        if missing:
            print(f'  warning: Model.csv missing columns: {missing}', file=sys.stderr)
        out = {'cellLines': []}
        for _, dsts in COLS:
            for d in dsts: out[d] = {}
        for row in rdr:
            mid = row.get('ModelID', '').strip()
            if not mid: continue
            out['cellLines'].append(mid)
            for src, dsts in COLS:
                val = (row.get(src) or '').strip()
                if val:
                    for d in dsts: out[d][mid] = val
    n = len(out['cellLines'])
    n_named = len(out['cellLineName'])
    print(f'  models: {n}, with name: {n_named}')

    # Quick coverage check against the CN matrix.
    cn_path = '/Users/fredrikwermeling/Documents/greenlistedv2/cn_metadata.json'
    if os.path.exists(cn_path):
        cn = json.load(open(cn_path))
        cn_set = set(cn['cellLines'])
        named_set = set(out['cellLineName'].keys())
        covered = cn_set & named_set
        print(f'  CN coverage: {len(covered)} / {len(cn_set)} lines now have a name')
        missing_cn = cn_set - named_set
        if missing_cn:
            print(f'  still missing: {sorted(missing_cn)[:10]}{" ..." if len(missing_cn) > 10 else ""}')

    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    size_kb = os.path.getsize(OUT) / 1024
    print(f'  wrote {OUT} ({size_kb:.0f} KB)')

if __name__ == '__main__':
    main()
