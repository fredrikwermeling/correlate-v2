#!/usr/bin/env python3
"""
Build web_data/mutations.json from DepMap's OmicsSomaticMutationsMatrixHotspot.csv.

Same recipe as the original build (kept in coexpress/preprocess_coexpress.py):
genes with >= 3 mutated models, top 49 by mutated-model count, values capped at 2.
Counts run over every model in the matrix; the app scopes to its cohort at use time.

Usage: python3 build_hotspot_mutations.py <OmicsSomaticMutationsMatrixHotspot.csv>
"""

import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "web_data", "mutations.json")


def parse_gene(col):
    m = re.match(r"^(.+?)\s*\(\d+\)$", col)
    return m.group(1).strip() if m else col.strip()


def main():
    csv_path = sys.argv[1]
    with open(csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        model_col = header.index("ModelID")
        default_col = header.index("IsDefaultEntryForModel") if "IsDefaultEntryForModel" in header else None
        gene_cols = [(i, c) for i, c in enumerate(header) if re.match(r"^.+\s*\(\d+\)$", c)]
        data = {parse_gene(c): {"column": c, "mutations": {}, "counts": {"0": 0, "1": 0, "2": 0}}
                for _, c in gene_cols}
        for row in reader:
            # 26Q1: one row per profile; keep each model's default entry only.
            if default_col is not None and row[default_col].strip() != "Yes":
                continue
            model = row[model_col]
            for i, c in gene_cols:
                gene = parse_gene(c)
                v = row[i].strip() if i < len(row) else "0"
                try:
                    vi = min(int(float(v)) if v else 0, 2)
                except ValueError:
                    vi = 0
                data[gene]["counts"][str(vi)] += 1
                if vi > 0:
                    data[gene]["mutations"][model] = vi

    kept = {g: d for g, d in data.items() if len(d["mutations"]) >= 3}
    for g, d in kept.items():
        d["total_mutated"] = len(d["mutations"])
    top = sorted(kept, key=lambda g: kept[g]["total_mutated"], reverse=True)[:49]
    out = {
        "genes": top,
        "geneCounts": {g: kept[g]["total_mutated"] for g in top},
        "geneData": {g: kept[g] for g in top},
    }
    with open(OUT, "w") as f:
        json.dump(out, f)
    print(f"{len(top)} hotspot genes -> {OUT}")
    print("top 5:", [(g, out['geneCounts'][g]) for g in top[:5]])


if __name__ == "__main__":
    main()
