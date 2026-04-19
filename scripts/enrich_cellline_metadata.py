#!/usr/bin/env python3
"""
Enrich web_data/cellLineMetadata.json with additional fields from DepMap 25Q3
Model.csv — everything the Cell Line Wiki modal needs to render.

Adds (keyed by ACH id):
  age, ageCategory, patientRace, primaryOrMetastasis, sampleCollectionSite,
  growthPattern, patientTumorGrade, depmapModelType, oncotreeSubtype,
  oncotreeCode, patientSubtypeFeatures, RRID, engineeredModel,
  culturedResistanceDrug

Idempotent. Run once; re-run whenever Model.csv is refreshed.
"""

import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_CSV = os.path.join(HERE, "Model_25Q3.csv")
METADATA_JSON = os.path.join(HERE, "..", "web_data", "cellLineMetadata.json")

# DepMap column → JSON key (camelCase)
FIELD_MAP = {
    "Age": "age",
    "AgeCategory": "ageCategory",
    "PatientRace": "patientRace",
    "PrimaryOrMetastasis": "primaryOrMetastasis",
    "SampleCollectionSite": "sampleCollectionSite",
    "GrowthPattern": "growthPattern",
    "PatientTumorGrade": "patientTumorGrade",
    "DepmapModelType": "depmapModelType",
    "OncotreeSubtype": "oncotreeSubtype",
    "OncotreeCode": "oncotreeCode",
    "PatientSubtypeFeatures": "patientSubtypeFeatures",
    "RRID": "rrid",
    "EngineeredModel": "engineeredModel",
    "CulturedResistanceDrug": "culturedResistanceDrug",
}


def clean(v):
    if v is None:
        return ""
    s = str(v).strip()
    return s if s and s.lower() != "nan" else ""


def main():
    # Load Model.csv
    rows_by_id = {}
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            rows_by_id[row["ModelID"]] = row
    print(f"Loaded {len(rows_by_id)} rows from {MODEL_CSV}")

    # Load metadata
    with open(METADATA_JSON) as f:
        meta = json.load(f)
    n_cells = len(meta["cellLines"])
    print(f"Enriching {n_cells} cell lines in {METADATA_JSON}")

    # For every target field, build a per-cell-line map
    for depmap_col, json_key in FIELD_MAP.items():
        field_map = {}
        for cl in meta["cellLines"]:
            row = rows_by_id.get(cl)
            if not row:
                continue
            v = clean(row.get(depmap_col, ""))
            if v:
                field_map[cl] = v
        meta[json_key] = field_map
        non_empty = len(field_map)
        print(f"  {json_key:24s}: {non_empty}/{n_cells} cells have value")

    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)
    print(f"\nDone. {METADATA_JSON} now includes the above {len(FIELD_MAP)} enrichment fields.")


if __name__ == "__main__":
    main()
