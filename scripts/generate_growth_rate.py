#!/usr/bin/env python3
"""
Generate growth_rate.json for Correlate V2.
Averages per-screen CRISPR-inferred growth rates per ModelID.

Input: CRISPRInferredModelGrowthRate.csv (DepMap 25Q3)
Output: web_data/growth_rate.json
"""

import csv
import json
import os
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV = os.path.join(SCRIPT_DIR, "CRISPRInferredModelGrowthRate.csv")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "..", "web_data")

def main():
    print("Reading CRISPRInferredModelGrowthRate.csv...")
    vals = defaultdict(list)
    with open(INPUT_CSV, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            model_id = row[0].split('_')[0]  # ACH-000001
            for v in row[1:]:
                v = v.strip()
                if v:
                    vals[model_id].append(float(v))

    # Average per model
    growth_rates = {}
    for model_id, vlist in vals.items():
        if vlist:
            growth_rates[model_id] = round(sum(vlist) / len(vlist), 6)

    print(f"  {len(growth_rates)} cell lines with growth rate data")

    # Stats
    all_vals = sorted(growth_rates.values())
    print(f"  Range: {all_vals[0]:.4f} to {all_vals[-1]:.4f}")
    print(f"  Median: {all_vals[len(all_vals)//2]:.4f}")

    # Check overlap with gene effect metadata
    meta_path = os.path.join(OUTPUT_DIR, "metadata.json")
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
        ge_cells = set(meta['cellLines'])
        overlap = ge_cells & set(growth_rates.keys())
        print(f"  Overlap with gene effect cell lines: {len(overlap)}/{len(ge_cells)} ({len(overlap)/len(ge_cells)*100:.1f}%)")

    output_path = os.path.join(OUTPUT_DIR, "growth_rate.json")
    with open(output_path, 'w') as f:
        json.dump(growth_rates, f)

    size_kb = os.path.getsize(output_path) / 1024
    print(f"\n=== DONE ===")
    print(f"  growth_rate.json: {size_kb:.1f} KB, {len(growth_rates)} cell lines")

if __name__ == '__main__':
    main()
