#!/usr/bin/env python3
"""
Impute sex from Y-chromosome markers + XIST expression.

Merges DepMap Sex annotation with expression-based inference into
web_data/cellLineMetadata.json.

Input:
  /Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanAllGenes.csv
  scripts/Model_25Q3.csv

Output: web_data/cellLineMetadata.json gets two independent fields keyed by ACH-id:
  - sex:             DepMap annotation ('Male' | 'Female' | 'Unknown')
  - sexByExpression: expression-based call ('male' | 'female' | 'unknown') — computed
                     for every cell line, independent of annotation.

Classifier (two-rule, independent thresholds):
  - Y_mean > Y_THR           -> 'male'    (Y presence is unambiguous; takes precedence)
  - XIST  > X_THR            -> 'female'  (only if Y is below threshold)
  - otherwise                -> 'unknown' (Y-loss males AND XIST-silenced females live here)

Thresholds Y_THR=1.0, X_THR=1.0 chosen to separate the two populations robustly in
both directions. Note: XIST silencing is common in breast and other cancers, so a
large minority of annotated females land in 'unknown' — this is biology, not a bug.

Script is idempotent.
"""

import csv
import json
import os
import re
import sys
from collections import Counter

import numpy as np

EXPRESSION_CSV = "/Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanAllGenes.csv"
HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_CSV = os.path.join(HERE, "Model_25Q3.csv")
METADATA_JSON = os.path.join(HERE, "..", "web_data", "cellLineMetadata.json")

Y_MARKERS = ["RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y"]
XIST = "XIST"

# Classifier thresholds (log-TPM units). See tuning in commit history.
Y_THR = 1.0
X_THR = 1.0


def parse_gene_name(col_header):
    m = re.match(r"^(.+?)\s*\(\d+\)$", col_header)
    return m.group(1).strip() if m else col_header.strip()


def main():
    print(f"Reading header from {EXPRESSION_CSV}")
    with open(EXPRESSION_CSV) as f:
        header = next(csv.reader(f))

    model_id_col = header.index("ModelID")
    is_default_col = header.index("IsDefaultEntryForModel")

    target_genes = Y_MARKERS + [XIST]
    gene_col_idx = {}
    for i, h in enumerate(header[6:], start=6):
        name = parse_gene_name(h)
        if name in target_genes:
            gene_col_idx[name] = i
    missing = [g for g in target_genes if g not in gene_col_idx]
    if missing:
        sys.exit(f"ERROR: markers missing from expression file: {missing}")
    print(f"  Located all {len(target_genes)} markers")

    print("Streaming rows (filtering IsDefaultEntryForModel=Yes)")
    cell_vals = {}
    with open(EXPRESSION_CSV) as f:
        reader = csv.reader(f)
        next(reader)
        for row_num, row in enumerate(reader):
            if row_num % 500 == 0 and row_num > 0:
                print(f"  Row {row_num}")
            if row[is_default_col].strip() != "Yes":
                continue
            cl = row[model_id_col].strip()
            vals = {}
            for g, idx in gene_col_idx.items():
                raw = row[idx].strip() if idx < len(row) else ""
                try:
                    vals[g] = float(raw) if raw else np.nan
                except ValueError:
                    vals[g] = np.nan
            cell_vals[cl] = vals
    print(f"  Loaded markers for {len(cell_vals)} cell lines")

    # Per cell line: mean(Y) and XIST
    scores = {}
    for cl, v in cell_vals.items():
        y_vals = [v[g] for g in Y_MARKERS if not np.isnan(v.get(g, np.nan))]
        x = v.get(XIST, np.nan)
        y_mean = float(np.mean(y_vals)) if y_vals else np.nan
        scores[cl] = (y_mean, float(x) if not np.isnan(x) else np.nan)

    print(f"Reading Sex labels from {MODEL_CSV}")
    depmap_sex = {}
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            s = row["Sex"].strip()
            depmap_sex[row["ModelID"]] = s if s else "Unknown"

    # Classifier: Y_THR takes precedence over X_THR so XIST-positive males
    # (e.g. Klinefelter) still get called male.
    def classify(y, x):
        if np.isnan(y) or np.isnan(x):
            return "unknown"
        if y > Y_THR:
            return "male"
        if x > X_THR:
            return "female"
        return "unknown"

    print(f"  Thresholds: Y_mean > {Y_THR} -> male, XIST > {X_THR} (Y low) -> female")

    # Validation on known-label cell lines
    male_y, male_x, female_y, female_x = [], [], [], []
    for cl, (y, x) in scores.items():
        if np.isnan(y) or np.isnan(x):
            continue
        label = depmap_sex.get(cl, "Unknown")
        if label == "Male":
            male_y.append(y); male_x.append(x)
        elif label == "Female":
            female_y.append(y); female_x.append(x)

    m_calls = Counter(classify(y, x) for y, x in zip(male_y, male_x))
    f_calls = Counter(classify(y, x) for y, x in zip(female_y, female_x))
    print("  Validation:")
    print(
        f"    Of {len(male_y)} known males:  "
        f"{m_calls.get('male', 0)} male, "
        f"{m_calls.get('female', 0)} female ({100 * m_calls.get('female', 0) / len(male_y):.2f}% FP), "
        f"{m_calls.get('unknown', 0)} unknown"
    )
    print(
        f"    Of {len(female_y)} known females: "
        f"{f_calls.get('male', 0)} male ({100 * f_calls.get('male', 0) / len(female_y):.2f}% FP), "
        f"{f_calls.get('female', 0)} female, "
        f"{f_calls.get('unknown', 0)} unknown"
    )

    by_expression = {cl: classify(y, x) for cl, (y, x) in scores.items()}

    print(f"Updating {METADATA_JSON}")
    with open(METADATA_JSON) as f:
        meta = json.load(f)

    sex_map = {}
    exp_map = {}
    for cl in meta["cellLines"]:
        sex_map[cl] = depmap_sex.get(cl, "Unknown")
        exp_map[cl] = by_expression.get(cl, "unknown")

    meta["sex"] = sex_map
    meta["sexByExpression"] = exp_map
    # Remove any previous field name
    meta.pop("sexImputed", None)

    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print("\nFinal distribution in cellLineMetadata.json:")
    print("  sex (annotation):    ", dict(Counter(sex_map.values())))
    print("  sexByExpression:     ", dict(Counter(exp_map.values())))

    # Agreement crosstab
    ct = Counter((sex_map[cl], exp_map[cl]) for cl in meta["cellLines"])
    print("\nAgreement crosstab (annotation × expression):")
    for (a, e), n in sorted(ct.items()):
        print(f"    {a:8s} × {e:8s}: {n}")


if __name__ == "__main__":
    main()
