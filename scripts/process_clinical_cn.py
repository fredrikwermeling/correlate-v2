#!/usr/bin/env python3
"""
Build web_data/clinical_cn.json from DepMap's OmicsCNGeneWGS.csv plus the
curated panels in scripts/clinical_cn_anchor.json. Per cell line, calls
clinically actionable focal AMPLIFICATIONS and EXTENDED DELETIONS using
relative-CN thresholds:

    amp        relative CN ≥ 3.0    (≥6 actual copies; high-level focal amp)
    strong_amp relative CN ≥ 5.0    (focal high-level)
    del        relative CN ≤ 0.5    (single-copy loss)
    deep_del   relative CN ≤ 0.3    (near-homozygous, biological LoF)

The 8 TSGs in DepMap's inferred-LoF panel (RB1, TP53, PTEN, NF1, CDKN2A,
VHL, MTAP, APC) are deliberately excluded from extended_deletions —
inferred.functionalLoss already covers them with integrated CN+mut+expr
evidence, no need to duplicate.

Input
  Path to OmicsCNGeneWGS.csv (~400 MB, DepMap 25Q3).

Usage
  python scripts/process_clinical_cn.py /path/to/OmicsCNGeneWGS.csv
"""

import csv
import json
import os
import re
import sys
from collections import defaultdict


def parse_gene_column_name(col):
    m = re.match(r'^([A-Za-z0-9.\-]+)(?:\s*\(\d+\))?\s*$', col)
    return m.group(1).upper() if m else col.strip().upper()


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_clinical_cn.py /path/to/OmicsCNGeneWGS.csv")
    src = sys.argv[1]
    if not os.path.exists(src):
        sys.exit(f"File not found: {src}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")

    anchor_path = os.path.join(script_dir, "clinical_cn_anchor.json")
    with open(anchor_path) as f:
        anchor = json.load(f)
    amp_panel = {e["gene"].upper(): e["context"] for e in anchor["amplifications"]}
    del_panel = {e["gene"].upper(): e["context"] for e in anchor["extended_deletions"]}
    panel_set = set(amp_panel) | set(del_panel)
    th = anchor["thresholds"]

    with open(os.path.join(web, "metadata.json")) as f:
        valid_cell_lines = set(json.load(f)["cellLines"])

    by_cl = {}  # cl -> { 'amplifications': [...], 'deletions': [...] }
    n_kept = 0

    print(f"Reading {src} ...")
    with open(src, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        try:
            model_idx = header.index("ModelID")
        except ValueError:
            sys.exit("CSV is missing 'ModelID' column.")
        try:
            default_idx = header.index("IsDefaultEntryForModel")
        except ValueError:
            default_idx = None
        first_gene_col = max([
            i for i, c in enumerate(header)
            if c in {"ModelID", "SequencingID", "IsDefaultEntryForModel",
                     "ModelConditionID", "IsDefaultEntryForMC"} or c == ""
        ]) + 1
        gene_lookup = {}
        for i in range(first_gene_col, len(header)):
            g = parse_gene_column_name(header[i])
            if g in panel_set:
                gene_lookup[i] = g
        print(f"  matched {len(gene_lookup)} of {len(panel_set)} panel genes in CSV header")
        if len(gene_lookup) < len(panel_set):
            missing = panel_set - set(gene_lookup.values())
            print(f"  missing from CSV: {sorted(missing)}")

        for row in reader:
            if default_idx is not None and len(row) > default_idx:
                if (row[default_idx] or "").strip().lower() != "yes":
                    continue
            cl = row[model_idx].strip() if len(row) > model_idx else ""
            if not cl or cl not in valid_cell_lines:
                continue
            n_kept += 1
            amps = []
            dels = []
            for i, g in gene_lookup.items():
                v = row[i].strip() if i < len(row) else ""
                if not v or v.lower() in {"nan", "na", "none", ""}:
                    continue
                try:
                    fv = float(v)
                except ValueError:
                    continue
                if g in amp_panel:
                    if fv >= th["strong_amp"]:
                        amps.append({"gene": g, "cn": round(fv, 2), "tier": "strong_amp", "context": amp_panel[g]})
                    elif fv >= th["amp"]:
                        amps.append({"gene": g, "cn": round(fv, 2), "tier": "amp", "context": amp_panel[g]})
                if g in del_panel:
                    if fv <= th["deep_del"]:
                        dels.append({"gene": g, "cn": round(fv, 2), "tier": "deep_del", "context": del_panel[g]})
                    elif fv <= th["del"]:
                        dels.append({"gene": g, "cn": round(fv, 2), "tier": "del", "context": del_panel[g]})
            if amps or dels:
                amps.sort(key=lambda x: -x["cn"])
                dels.sort(key=lambda x: x["cn"])
                entry = {}
                if amps: entry["amplifications"] = amps
                if dels: entry["deletions"] = dels
                by_cl[cl] = entry

    print(f"  {n_kept} default-entry cell lines processed; {len(by_cl)} have at least one CN call")

    # Counts per gene for sanity check
    amp_counts = defaultdict(int)
    del_counts = defaultdict(int)
    for entry in by_cl.values():
        for a in entry.get("amplifications", []):
            amp_counts[a["gene"]] += 1
        for d in entry.get("deletions", []):
            del_counts[d["gene"]] += 1

    print("\nTop 10 amplifications by cell-line count:")
    for g, n in sorted(amp_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {g:8s} n={n}")
    print("\nTop 10 extended deletions by cell-line count:")
    for g, n in sorted(del_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {g:8s} n={n}")

    output = {
        "_description": (
            "Clinically actionable focal CN events per cell line: amplifications "
            "(oncogene high-level CN gain) and extended deletions (TSGs not in "
            "DepMap's inferred-LoF panel). Source: OmicsCNGeneWGS.csv; thresholds "
            "in 'thresholds' field. The existing 8 TSGs in inferred.functionalLoss "
            "are deliberately excluded from this list to avoid duplication."
        ),
        "schemaVersion": "1.0",
        "thresholds": th,
        "amplificationPanel": [
            {"gene": g, "context": amp_panel[g]} for g in sorted(amp_panel.keys())
        ],
        "deletionPanel": [
            {"gene": g, "context": del_panel[g]} for g in sorted(del_panel.keys())
        ],
        "byCellLine": by_cl,
    }
    out_path = os.path.join(web, "clinical_cn.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
