#!/usr/bin/env python3
"""
Regenerate expression.bin.gz and expression_metadata.json for Correlate V2.
Uses scale factor 1800 (max ~18.2) instead of 3000 (max ~10.9) to avoid clipping.

Input: DepMap OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
Output: web_data/expression.bin.gz, web_data/expression_metadata.json
"""

import csv
import json
import gzip
import re
import os
import numpy as np

# Paths
EXPRESSION_CSV = "/Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
OUTPUT_DIR = "/Users/fredrikwermeling/Documents/correlate_v2/web_data"

SCALE_FACTOR = 1800  # Max = 32767/1800 = 18.2, covers max observed 17.36
NA_VALUE = -32768

def parse_gene_name(col_header):
    match = re.match(r'^(.+?)\s*\(\d+\)$', col_header)
    return match.group(1).strip() if match else col_header.strip()

def main():
    print("Reading expression CSV...")
    with open(EXPRESSION_CSV, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

    meta_col_names = header[:6]
    gene_columns = header[6:]
    gene_names = [parse_gene_name(col) for col in gene_columns]
    n_genes = len(gene_names)
    print(f"  Found {n_genes} genes")

    model_id_col = header.index('ModelID')
    is_default_col = header.index('IsDefaultEntryForModel')

    print("Reading expression values (filtering to IsDefaultEntryForModel=Yes)...")
    cell_lines = []
    expression_data = []

    with open(EXPRESSION_CSV, 'r') as f:
        reader = csv.reader(f)
        next(reader)
        for row_num, row in enumerate(reader):
            if row_num % 200 == 0:
                print(f"  Processing row {row_num}...")

            if row[is_default_col].strip() != 'Yes':
                continue

            model_id = row[model_id_col].strip()
            cell_lines.append(model_id)

            values = []
            for i in range(6, len(row)):
                try:
                    values.append(float(row[i].strip()))
                except (ValueError, IndexError):
                    values.append(float('nan'))
            while len(values) < n_genes:
                values.append(float('nan'))

            expression_data.append(values)

    n_cell_lines = len(cell_lines)
    print(f"  Loaded {n_cell_lines} cell lines x {n_genes} genes")

    # Convert to Int16 with new scale factor
    print(f"Converting to Int16 binary matrix (scale={SCALE_FACTOR}, max={32767/SCALE_FACTOR:.2f})...")
    data_array = np.array(expression_data, dtype=np.float32)

    # Check range
    valid = data_array[~np.isnan(data_array)]
    print(f"  Value range: {valid.min():.4f} to {valid.max():.4f}")
    clipped = np.sum(valid > 32767 / SCALE_FACTOR)
    print(f"  Values that would be clipped: {clipped}")

    int16_data = np.full(data_array.shape, NA_VALUE, dtype=np.int16)
    valid_mask = ~np.isnan(data_array)
    int16_data[valid_mask] = np.clip(
        np.round(data_array[valid_mask] * SCALE_FACTOR).astype(np.int32),
        -32767, 32767
    ).astype(np.int16)

    # Transpose to gene-major order
    int16_transposed = int16_data.T
    flat_data = int16_transposed.flatten()
    raw_bytes = flat_data.tobytes()

    output_bin = os.path.join(OUTPUT_DIR, 'expression.bin.gz')
    print(f"  Writing {output_bin}...")
    with gzip.open(output_bin, 'wb') as f:
        f.write(raw_bytes)

    uncompressed_mb = len(raw_bytes) / (1024 * 1024)
    compressed_mb = os.path.getsize(output_bin) / (1024 * 1024)
    print(f"  Uncompressed: {uncompressed_mb:.1f} MB, Compressed: {compressed_mb:.1f} MB")

    # Write metadata
    print("Creating expression_metadata.json...")
    metadata = {
        'nGenes': n_genes,
        'nCellLines': n_cell_lines,
        'scaleFactor': SCALE_FACTOR,
        'naValue': NA_VALUE,
        'genes': gene_names,
        'cellLines': cell_lines
    }

    with open(os.path.join(OUTPUT_DIR, 'expression_metadata.json'), 'w') as f:
        json.dump(metadata, f)

    print(f"\n=== DONE ===")
    print(f"  expression.bin.gz: {compressed_mb:.1f} MB")
    print(f"  Scale factor: {SCALE_FACTOR} (max value: {32767/SCALE_FACTOR:.2f})")
    print(f"  {n_genes} genes x {n_cell_lines} cell lines")

if __name__ == '__main__':
    main()
