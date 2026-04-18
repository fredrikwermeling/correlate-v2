#!/usr/bin/env python3
"""
Regenerate geneEffects.bin.gz with a corrected scale factor.

The original scaleFactor=10000 causes Int16 overflow for gene effect values
below -3.2767 (Int16 min usable = -32767, divided by 10000 = -3.2767).
DepMap data goes as low as -5.78, so ~6000 values were being clipped.

Fix: use scaleFactor=5000 → range ±6.5534, precision 0.0002 (sufficient).
"""

import csv
import json
import gzip
import os
import numpy as np

# Paths
CRISPR_CSV = "/Users/fredrikwermeling/Documents/CRISPR screen script depmap/crispr screen script depmap conections/CRISPRGeneEffect.csv"
OUTPUT_DIR = "/Users/fredrikwermeling/Documents/correlate_v2/web_data"
V1_OUTPUT_DIR = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app/web_data"

NEW_SCALE_FACTOR = 5000
NA_VALUE = -32768

def main():
    # Load existing metadata to preserve structure
    meta_path = os.path.join(OUTPUT_DIR, 'metadata.json')
    with open(meta_path, 'r') as f:
        metadata = json.load(f)

    existing_genes = metadata['genes']
    existing_cell_lines = metadata['cellLines']
    n_genes = metadata['nGenes']
    n_cell_lines = metadata['nCellLines']
    old_sf = metadata['scaleFactor']

    print(f"Existing: {n_genes} genes x {n_cell_lines} cell lines, scaleFactor={old_sf}")
    print(f"New scaleFactor: {NEW_SCALE_FACTOR}")
    print(f"  Old range: ±{32767/old_sf:.4f}")
    print(f"  New range: ±{32767/NEW_SCALE_FACTOR:.4f}")

    # Read CSV
    print(f"\nReading {CRISPR_CSV}...")
    with open(CRISPR_CSV, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

    # Parse gene names from header (format: "GENE (EntrezID)")
    import re
    csv_genes = {}
    for i, col in enumerate(header[1:], 1):
        match = re.match(r'^(.+?)\s*\(\d+\)$', col)
        gene_name = match.group(1).strip() if match else col.strip()
        csv_genes[gene_name] = i

    # Build cell line -> row mapping
    print("Reading cell line data...")
    csv_data = {}  # model_id -> {gene_name: value}
    with open(CRISPR_CSV, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row_num, row in enumerate(reader):
            model_id = row[0].strip()
            csv_data[model_id] = row

    print(f"  CSV has {len(csv_data)} cell lines, {len(csv_genes)} genes")

    # Build the matrix in the same order as existing metadata
    print("Building matrix in existing gene/cell line order...")
    data_matrix = np.full((n_genes, n_cell_lines), np.nan, dtype=np.float32)

    genes_found = 0
    cl_found = 0
    for gi, gene in enumerate(existing_genes):
        if gene in csv_genes:
            genes_found += 1
            col_idx = csv_genes[gene]
            for ci, cl in enumerate(existing_cell_lines):
                if cl in csv_data:
                    row = csv_data[cl]
                    if col_idx < len(row):
                        val = row[col_idx].strip()
                        if val and val != 'NA':
                            try:
                                data_matrix[gi, ci] = float(val)
                            except ValueError:
                                pass

    for cl in existing_cell_lines:
        if cl in csv_data:
            cl_found += 1

    print(f"  Matched {genes_found}/{n_genes} genes, {cl_found}/{n_cell_lines} cell lines")

    # Stats
    valid = ~np.isnan(data_matrix)
    print(f"\n  Total values: {valid.sum()}")
    print(f"  NaN values: {(~valid).sum()}")
    print(f"  Min: {np.nanmin(data_matrix):.4f}")
    print(f"  Max: {np.nanmax(data_matrix):.4f}")
    below_old = np.sum(data_matrix[valid] < -32767/old_sf)
    below_new = np.sum(data_matrix[valid] < -32767/NEW_SCALE_FACTOR)
    print(f"  Values that were clipped (old SF={old_sf}): {below_old}")
    print(f"  Values that would be clipped (new SF={NEW_SCALE_FACTOR}): {below_new}")

    # Quantize to Int16
    print(f"\nQuantizing with scaleFactor={NEW_SCALE_FACTOR}...")
    int16_data = np.full(data_matrix.shape, NA_VALUE, dtype=np.int16)
    int16_data[valid] = np.clip(
        np.round(data_matrix[valid] * NEW_SCALE_FACTOR).astype(np.int32),
        -32767, 32767
    ).astype(np.int16)

    # Write binary
    flat = int16_data.flatten()
    raw_bytes = flat.tobytes()

    for out_dir, label in [(OUTPUT_DIR, "V2"), (V1_OUTPUT_DIR, "V1")]:
        out_path = os.path.join(out_dir, 'geneEffects.bin.gz')
        print(f"\nWriting {label}: {out_path}")
        with gzip.open(out_path, 'wb') as f:
            f.write(raw_bytes)

        size_mb = os.path.getsize(out_path) / (1024 * 1024)
        print(f"  Size: {size_mb:.1f} MB")

        # Update metadata
        meta_file = os.path.join(out_dir, 'metadata.json')
        with open(meta_file, 'r') as f:
            meta = json.load(f)
        meta['scaleFactor'] = NEW_SCALE_FACTOR
        with open(meta_file, 'w') as f:
            json.dump(meta, f)
        print(f"  Updated scaleFactor in metadata.json: {old_sf} → {NEW_SCALE_FACTOR}")

    print("\n=== DONE ===")

if __name__ == '__main__':
    main()
