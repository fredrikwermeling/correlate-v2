#!/usr/bin/env python3
"""
Impute sex from Y-chromosome markers + XIST expression.

Merges DepMap Sex annotation with expression-based inference into
web_data/cellLineMetadata.json.

Input:
  /Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanAllGenes.csv
  scripts/Model_25Q3.csv

Output: web_data/cellLineMetadata.json gets two new fields keyed by ACH-id:
  - sex:        DepMap annotation ('Male' | 'Female' | 'Unknown')
  - sexImputed: expression-based call ('likely_male' | 'likely_female' | 'unknown')

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

    # Score per cell line: mean(Y) - XIST
    #   Male: Y high, XIST low -> score high
    #   Female: Y low, XIST high -> score low
    #   Y-loss male or low-confidence: near zero
    scores = {}
    for cl, v in cell_vals.items():
        y_vals = [v[g] for g in Y_MARKERS if not np.isnan(v.get(g, np.nan))]
        x = v.get(XIST, np.nan)
        y_mean = float(np.mean(y_vals)) if y_vals else np.nan
        if np.isnan(y_mean) or np.isnan(x):
            scores[cl] = (np.nan, np.nan, np.nan)
        else:
            scores[cl] = (y_mean, float(x), y_mean - float(x))

    print(f"Reading Sex labels from {MODEL_CSV}")
    depmap_sex = {}
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            s = row["Sex"].strip()
            depmap_sex[row["ModelID"]] = s if s else "Unknown"

    male_s, female_s = [], []
    for cl, (_, _, s) in scores.items():
        if np.isnan(s):
            continue
        label = depmap_sex.get(cl, "Unknown")
        if label == "Male":
            male_s.append(s)
        elif label == "Female":
            female_s.append(s)
    male_s = np.array(male_s)
    female_s = np.array(female_s)
    print(f"  Training scores: {len(male_s)} males, {len(female_s)} females")
    print(
        f"    Male   median={np.median(male_s):+.3f}, 1st pct={np.percentile(male_s, 1):+.3f}"
    )
    print(
        f"    Female median={np.median(female_s):+.3f}, 99th pct={np.percentile(female_s, 99):+.3f}"
    )

    # Conservative thresholds: <1% false-positive rate on each side
    thr_male = float(np.percentile(female_s, 99))
    thr_female = float(np.percentile(male_s, 1))
    print(
        f"  Thresholds: >{thr_male:+.3f} -> likely_male, <{thr_female:+.3f} -> likely_female"
    )
    if thr_male <= thr_female:
        print("  WARNING: thresholds overlap — the two populations are not separable.")

    def classify(s):
        if np.isnan(s):
            return "unknown"
        if s > thr_male:
            return "likely_male"
        if s < thr_female:
            return "likely_female"
        return "unknown"

    # Validation on known-label cell lines
    male_calls = Counter(classify(s) for s in male_s)
    female_calls = Counter(classify(s) for s in female_s)
    print("  Validation:")
    print(
        f"    Of {len(male_s)} known males:  "
        f"{male_calls.get('likely_male', 0)} likely_male, "
        f"{male_calls.get('likely_female', 0)} likely_female, "
        f"{male_calls.get('unknown', 0)} unknown"
    )
    print(
        f"    Of {len(female_s)} known females: "
        f"{female_calls.get('likely_male', 0)} likely_male, "
        f"{female_calls.get('likely_female', 0)} likely_female, "
        f"{female_calls.get('unknown', 0)} unknown"
    )

    imputed = {cl: classify(s) for cl, (_, _, s) in scores.items()}

    print(f"Updating {METADATA_JSON}")
    with open(METADATA_JSON) as f:
        meta = json.load(f)

    sex_map = {}
    imputed_map = {}
    for cl in meta["cellLines"]:
        sex_map[cl] = depmap_sex.get(cl, "Unknown")
        imputed_map[cl] = imputed.get(cl, "unknown")

    meta["sex"] = sex_map
    meta["sexImputed"] = imputed_map

    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print("\nFinal distribution in cellLineMetadata.json:")
    print("  sex (DepMap):", dict(Counter(sex_map.values())))
    print("  sexImputed:  ", dict(Counter(imputed_map.values())))

    unknown_cls = [cl for cl, s in sex_map.items() if s == "Unknown"]
    unk_imp = Counter(imputed_map[cl] for cl in unknown_cls)
    print(f"\nFor {len(unknown_cls)} DepMap-Unknown cell lines:")
    for k, v in unk_imp.items():
        print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
