#!/usr/bin/env python3
"""
Extract HLA / antigen-processing-machinery copy number from DepMap's
OmicsCNGeneWGS.csv into a compact web_data/hla_cn.json.

Used by Correlate's class-I antigen-presentation inference. The full WGS
CN matrix is ~400 MB; we only need ~12 genes to call likely class-I loss
in combination with B2M damaging mutations and HLA expression. Output is
~50 KB.

Genes extracted:
  Class-I HLA:    HLA-A, HLA-B, HLA-C
  β2-microglobulin: B2M
  Antigen processing: TAP1, TAP2, TAPBP, ERAP1, ERAP2, NLRC5
  Class-II HLA (for completeness): HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DQA1

Format detection: DepMap's WGS CN matrix is log2(CN/2) ratio. Threshold
guidance for downstream calls:
  log2 ratio < -0.4 : suggestive of one-copy loss
  log2 ratio < -0.7 : strong evidence
  log2 ratio < -1.0 : near-homozygous deletion

Input
  Path to OmicsCNGeneWGS.csv from DepMap 25Q3 (~400 MB).
  Download: https://depmap.org/portal/data_page/?tab=currentRelease
            (file: OmicsCNGeneWGS.csv)

Usage
  python scripts/process_hla_cn.py /path/to/OmicsCNGeneWGS.csv
"""

import csv
import json
import os
import re
import sys


HLA_PANEL = [
    # Class-I
    'HLA-A', 'HLA-B', 'HLA-C',
    # β2-microglobulin
    'B2M',
    # Antigen-processing machinery
    'TAP1', 'TAP2', 'TAPBP', 'ERAP1', 'ERAP2', 'NLRC5',
    # Class-II (for completeness — not used by class-I inference)
    'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1',
]


def parse_gene_column_name(col):
    """DepMap CN columns are typically 'GENE (entrez_id)' — strip the suffix."""
    m = re.match(r'^([A-Za-z0-9.\-]+)(?:\s*\(\d+\))?\s*$', col)
    return m.group(1).upper() if m else col.strip().upper()


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_hla_cn.py /path/to/OmicsCNGeneWGS.csv")
    src = sys.argv[1]
    if not os.path.exists(src):
        sys.exit(f"File not found: {src}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")
    out_path = os.path.join(web, "hla_cn.json")

    # Load valid cell lines so we only emit data for lines we have GE for.
    with open(os.path.join(web, "metadata.json")) as f:
        valid_cell_lines = set(json.load(f)["cellLines"])
    print(f"Filtering to {len(valid_cell_lines)} cell lines from metadata.json")

    panel_set = {g.upper() for g in HLA_PANEL}
    by_cl = {}
    n_rows = 0
    sample_values = []

    print(f"Reading {src} ...")
    with open(src, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        # First column is ModelID; remaining are gene columns.
        gene_cols = header[1:]
        gene_lookup = {}
        for i, raw in enumerate(gene_cols):
            g = parse_gene_column_name(raw)
            if g in panel_set:
                gene_lookup[i] = g
        print(f"  matched {len(gene_lookup)} of {len(panel_set)} panel genes in CSV header")
        if len(gene_lookup) < len(panel_set):
            missing = panel_set - set(gene_lookup.values())
            print(f"  missing from CSV: {sorted(missing)}")

        for row in reader:
            n_rows += 1
            cl = row[0].strip()
            if cl not in valid_cell_lines:
                continue
            entry = {}
            for i, g in gene_lookup.items():
                v = row[i + 1].strip() if i + 1 < len(row) else ""
                if not v or v.lower() in {"nan", "na", "none", ""}:
                    continue
                try:
                    fv = float(v)
                except ValueError:
                    continue
                entry[g] = round(fv, 3)
                if len(sample_values) < 200:
                    sample_values.append(fv)
            if entry:
                by_cl[cl] = entry

    # Detect format: log2-ratio (most values |<2|, near 0 = diploid) vs
    # absolute (mostly 0–6, near 2 = diploid).
    if sample_values:
        n_neg = sum(1 for v in sample_values if v < 0)
        n_le2 = sum(1 for v in sample_values if v <= 2)
        looks_log2 = (n_neg / len(sample_values) > 0.1) and (n_le2 / len(sample_values) > 0.9)
        cn_format = "log2_ratio" if looks_log2 else "absolute"
    else:
        cn_format = "unknown"
    print(f"  detected CN format: {cn_format}")

    output = {
        "_description": (
            "Per-cell-line copy number values for HLA / antigen-processing-"
            "machinery genes, extracted from DepMap OmicsCNGeneWGS.csv. Used by "
            "Correlate's class-I antigen-presentation inference."
        ),
        "schemaVersion": "1.0",
        "source": "DepMap 25Q3 OmicsCNGeneWGS.csv",
        "cnFormat": cn_format,
        "panel": HLA_PANEL,
        "thresholds": {
            "log2_ratio": {
                "suggestive_loss": -0.4,
                "strong_loss": -0.7,
                "deep_deletion": -1.0
            },
            "note": "If cnFormat is 'absolute', subtract 2 and divide by 2 to get approximate log2-ratio."
        },
        "byCellLine": by_cl,
    }
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten {out_path} ({size_kb:.0f} KB) — {len(by_cl)} cell lines")
    # Summary
    n_with_loss = 0
    for cl, e in by_cl.items():
        class1 = [e.get('HLA-A'), e.get('HLA-B'), e.get('HLA-C')]
        valid = [v for v in class1 if v is not None]
        if not valid:
            continue
        mean_class1 = sum(valid) / len(valid)
        if cn_format == "log2_ratio" and mean_class1 < -0.4:
            n_with_loss += 1
        elif cn_format == "absolute" and mean_class1 < 1.6:
            n_with_loss += 1
    print(f"  {n_with_loss} cell lines show suggestive class-I CN loss (mean HLA-A/B/C threshold).")


if __name__ == "__main__":
    main()
