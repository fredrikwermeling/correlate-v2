#!/usr/bin/env python3
"""
Build CoExpress's web_data from Correlate V2's 26Q1 build.

CoExpress is the V1 app with EXPRESSION as its primary matrix: the file the
app loads as geneEffects.bin.gz/metadata.json holds the 26Q1 log2(TPM+1)
matrix (19,215 genes x 1,719 default-entry models), so every feature
(correlations, network, clusters, mutation analysis) runs on expression.

- metadata.json          expression genes/cellLines (+genesFull from the CSV header)
- geneEffects.bin.gz     byte-copy of expression.bin.gz (int16, scale 1800)
- expression.bin.gz/_metadata  copied too (keeps the secondary expression layer working)
- cellLineMetadata.json  base + enrichment for the 1,719 cohort from Model26Q1.csv,
                         sexByExpression/strProfile carried from V2 where known
- everything not cohort-scoped is copied from V2's web_data verbatim; the
  cohort-scoped layers (translocations, damaging, validated fusions,
  functional loss, subtypes, signatures) are rebuilt by the copied scripts
  run from coexpress-app afterwards.
"""

import csv
import gzip
import json
import os
import re
import shutil

V2 = "/Users/fredrikwermeling/Documents/correlate_v2"
CO = "/Users/fredrikwermeling/Documents/coexpress-app"
Q = os.path.join(V2, "26Q1")
W2 = os.path.join(V2, "web_data")
WC = os.path.join(CO, "web_data")

FIELD_MAP = {
    "Age": "age", "AgeCategory": "ageCategory", "PatientRace": "patientRace",
    "PrimaryOrMetastasis": "primaryOrMetastasis", "SampleCollectionSite": "sampleCollectionSite",
    "GrowthPattern": "growthPattern", "PatientTumorGrade": "patientTumorGrade",
    "DepmapModelType": "depmapModelType", "OncotreeSubtype": "oncotreeSubtype",
    "OncotreeCode": "oncotreeCode", "PatientSubtypeFeatures": "patientSubtypeFeatures",
    "RRID": "rrid", "EngineeredModel": "engineeredModel", "CulturedResistanceDrug": "culturedResistanceDrug",
}

os.makedirs(WC, exist_ok=True)

# 1. metadata.json from the expression layer, with genesFull from the CSV header
em = json.load(open(os.path.join(W2, "expression_metadata.json")))
with open(os.path.join(Q, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")) as f:
    header = f.readline().rstrip("\n").split(",")
genes_full_by_name = {}
for col in header[6:]:
    m = re.match(r"^(.+?)\s*\(\d+\)$", col)
    genes_full_by_name[(m.group(1).strip() if m else col.strip())] = col
genes_full = [genes_full_by_name.get(g, g) for g in em["genes"]]
meta = {
    "nGenes": em["nGenes"], "nCellLines": em["nCellLines"],
    "scaleFactor": em["scaleFactor"], "naValue": em["naValue"],
    "genes": em["genes"], "genesFull": genes_full, "cellLines": em["cellLines"],
}
json.dump(meta, open(os.path.join(WC, "metadata.json"), "w"))
print(f"metadata.json: {meta['nGenes']} genes x {meta['nCellLines']} lines, scale {meta['scaleFactor']}")

# 2. primary matrix = expression bytes
shutil.copy2(os.path.join(W2, "expression.bin.gz"), os.path.join(WC, "geneEffects.bin.gz"))
print("geneEffects.bin.gz <- expression.bin.gz")

# 3. secondary expression layer kept identical
shutil.copy2(os.path.join(W2, "expression.bin.gz"), os.path.join(WC, "expression.bin.gz"))
shutil.copy2(os.path.join(W2, "expression_metadata.json"), os.path.join(WC, "expression_metadata.json"))

# 4. cellLineMetadata for the 1,719 cohort
model = {}
with open(os.path.join(V2, "Model26Q1.csv")) as f:
    for row in csv.DictReader(f):
        if row.get("ModelID"):
            model[row["ModelID"]] = row
old = json.load(open(os.path.join(W2, "cellLineMetadata.json")))
cohort = meta["cellLines"]
out = {"cellLines": cohort, "cellLineName": {}, "strippedCellLineName": {},
       "lineage": {}, "primaryDisease": {}, "subtype": {}}
for k in FIELD_MAP.values():
    out[k] = {}
missing = 0
for cl in cohort:
    m = model.get(cl)
    if not m:
        missing += 1
        out["cellLineName"][cl] = cl
        out["strippedCellLineName"][cl] = cl
        out["lineage"][cl] = "Unknown"
        out["primaryDisease"][cl] = "Unknown"
        out["subtype"][cl] = ""
        continue
    out["cellLineName"][cl] = m.get("CellLineName", "")
    out["strippedCellLineName"][cl] = m.get("StrippedCellLineName", "")
    out["lineage"][cl] = m.get("OncotreeLineage", "")
    out["primaryDisease"][cl] = m.get("OncotreePrimaryDisease", "")
    out["subtype"][cl] = m.get("OncotreeSubtype", "")
    out["sex"] = out.get("sex", {})
    out["sex"][cl] = m.get("Sex", "")
    for src, dst in FIELD_MAP.items():
        v = (m.get(src) or "").strip()
        if v:
            out[dst][cl] = v
cs = set(cohort)
for field in ("sexByExpression", "strProfile"):
    if field in old:
        out[field] = {cl: v for cl, v in old[field].items() if cl in cs}
json.dump(out, open(os.path.join(WC, "cellLineMetadata.json"), "w"))
print(f"cellLineMetadata.json: {len(cohort)} lines ({missing} not in Model26Q1)")

# 5. verbatim copies of the non-cohort-scoped layers
COPY = ["mutations.json", "growth_rate.json", "drug_response.json", "cn.bin.gz",
        "cn_metadata.json", "gene_locations.json", "synonyms.json", "orthologs.json",
        "corum_partners.json", "reactome_partners.json", "lehmann_tnbc.json",
        "problematic_lines.json", "hla_cn.json", "common_essentials.json",
        "clinical_fusions.json", "clinical_cn.json", "curated_fusions.json",
        "coexpress_logo.png"]
for f in COPY:
    src = os.path.join(W2, f)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(WC, f))
    elif not os.path.exists(os.path.join(WC, f)):
        print(f"  NOTE: {f} not in V2 web_data and not already present")
print("copied shared layers")
print("DONE (cohort-scoped rebuilds run separately from coexpress-app/scripts)")
