#!/usr/bin/env python3
"""
Rebuild Correlate V2's core matrices from the DepMap 26Q1 release.

Inputs (in ../26Q1/):
  CRISPRGeneEffect.csv                                  rows=ModelID, cols="GENE (Entrez)"
  OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv    rows=profiles, 6 meta cols, filter IsDefaultEntryForModel

Outputs (in ../web_data/):
  geneEffects.bin.gz + metadata.json            int16, scaleFactor 5000, NA -32768, gene-major
  expression.bin.gz + expression_metadata.json  int16, scaleFactor 1800, NA -32768, gene-major
  cellLineMetadata.json                         BASE fields for the gene-effect cohort, from ../Model26Q1.csv.
                                                sexByExpression + strProfile are restored from the previous
                                                file for carried-over lines (XIST is absent from the
                                                protein-coding expression file, so the imputation cannot be
                                                re-run without the AllGenes download; STR needs Cellosaurus).

Cohort rule (same as 25Q3): every row of CRISPRGeneEffect.csv, sorted by ModelID.
"""

import csv
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
Q = os.path.join(ROOT, "26Q1")
WEB = os.path.join(ROOT, "web_data")

GE_CSV = os.path.join(Q, "CRISPRGeneEffect.csv")
EXPR_CSV = os.path.join(Q, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
MODEL_CSV = os.path.join(ROOT, "Model26Q1.csv")

GE_SCALE = 5000
EXPR_SCALE = 1800
NA = -32768


def parse_gene(col):
    m = re.match(r"^(.+?)\s*\(\d+\)$", col)
    return m.group(1).strip() if m else col.strip()


def quantize(mat, scale):
    out = np.full(mat.shape, NA, dtype=np.int16)
    valid = ~np.isnan(mat)
    out[valid] = np.clip(np.round(mat[valid] * scale).astype(np.int32), -32767, 32767).astype(np.int16)
    return out


def write_bin(path, int16_matrix):
    with gzip.open(path, "wb") as f:
        f.write(int16_matrix.flatten().tobytes())
    print(f"  wrote {path} ({os.path.getsize(path)/1e6:.1f} MB)")


def build_gene_effect():
    print("=== Gene effect (CRISPRGeneEffect.csv) ===")
    with open(GE_CSV) as f:
        reader = csv.reader(f)
        header = next(reader)
        genes_full = header[1:]
        genes = [parse_gene(c) for c in genes_full]
        rows = {}
        for row in reader:
            model = row[0].strip()
            vals = np.array([float(v) if v not in ("", "NA") else np.nan for v in row[1:]], dtype=np.float32)
            rows[model] = vals
    cell_lines = sorted(rows)
    print(f"  {len(genes)} genes x {len(cell_lines)} cell lines")
    mat = np.vstack([rows[cl] for cl in cell_lines])            # [cell][gene]
    ge = quantize(mat, GE_SCALE).T.copy()                        # gene-major
    write_bin(os.path.join(WEB, "geneEffects.bin.gz"), ge)
    meta = {
        "nGenes": len(genes),
        "nCellLines": len(cell_lines),
        "scaleFactor": GE_SCALE,
        "naValue": NA,
        "genes": genes,
        "genesFull": genes_full,
        "cellLines": cell_lines,
    }
    with open(os.path.join(WEB, "metadata.json"), "w") as f:
        json.dump(meta, f)
    print(f"  min {np.nanmin(mat):.3f}  max {np.nanmax(mat):.3f}  NaN {np.isnan(mat).sum()}")
    return cell_lines


def build_expression():
    print("=== Expression (OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv) ===")
    with open(EXPR_CSV) as f:
        reader = csv.reader(f)
        header = next(reader)
        model_col = header.index("ModelID")
        default_col = header.index("IsDefaultEntryForModel")
        gene_start = 6
        genes_full = header[gene_start:]
        genes = [parse_gene(c) for c in genes_full]
        rows = {}
        for row in reader:
            if row[default_col].strip() != "Yes":
                continue
            model = row[model_col].strip()
            vals = np.array([float(v) if v not in ("", "NA") else np.nan for v in row[gene_start:]], dtype=np.float32)
            rows[model] = vals
    cell_lines = sorted(rows)
    print(f"  {len(genes)} genes x {len(cell_lines)} cell lines")
    mat = np.vstack([rows[cl] for cl in cell_lines])
    ex = quantize(mat, EXPR_SCALE).T.copy()
    write_bin(os.path.join(WEB, "expression.bin.gz"), ex)
    meta = {
        "nGenes": len(genes),
        "nCellLines": len(cell_lines),
        "scaleFactor": EXPR_SCALE,
        "naValue": NA,
        "genes": genes,
        "cellLines": cell_lines,
    }
    with open(os.path.join(WEB, "expression_metadata.json"), "w") as f:
        json.dump(meta, f)
    print(f"  min {np.nanmin(mat):.3f}  max {np.nanmax(mat):.3f}")


def build_cell_line_metadata(cohort):
    print("=== cellLineMetadata.json (base fields) ===")
    old_path = os.path.join(WEB, "cellLineMetadata.json")
    with open(old_path) as f:
        old = json.load(f)

    model = {}
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            mid = row.get("ModelID", "")
            if mid:
                model[mid] = row

    out = {
        "cellLines": cohort,
        "cellLineName": {},
        "strippedCellLineName": {},
        "lineage": {},
        "primaryDisease": {},
        "subtype": {},
    }
    missing = []
    for cl in cohort:
        m = model.get(cl)
        if not m:
            missing.append(cl)
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

    # Carry forward what cannot be recomputed from these downloads.
    for field in ("sexByExpression", "strProfile"):
        if field in old:
            kept = {cl: v for cl, v in old[field].items() if cl in set(cohort)}
            out[field] = kept
            print(f"  carried forward {field} for {len(kept)} lines")

    with open(old_path, "w") as f:
        json.dump(out, f)
    if missing:
        print(f"  WARNING: {len(missing)} cohort lines missing from Model26Q1.csv: {missing[:5]}")
    print(f"  {len(cohort)} cell lines written")


if __name__ == "__main__":
    cohort = build_gene_effect()
    build_expression()
    build_cell_line_metadata(cohort)
    print("=== DONE ===")
