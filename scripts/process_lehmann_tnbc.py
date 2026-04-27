#!/usr/bin/env python3
"""
Map Lehmann TNBC subtype classifications (cell-line names) to DepMap ACH
IDs and emit web_data/lehmann_tnbc.json.

Source: scripts/lehmann_tnbc_anchor.json (curated from Lehmann et al.,
JCI 2011 PMID 21633166 + 2016 reclassification PMID 27310713).

Usage
  python scripts/process_lehmann_tnbc.py
"""

import json
import os
import re
import sys


def normalize_name(s):
    """Strip non-alphanumeric and uppercase — matches DepMap's
    strippedCellLineName convention (BT-549 → 'BT549',
    MDA-MB-231 → 'MDAMB231', SUM-149PT → 'SUM149PT')."""
    return re.sub(r"[^A-Za-z0-9]", "", s or "").upper()


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")
    anchor_path = os.path.join(script_dir, "lehmann_tnbc_anchor.json")
    metadata_path = os.path.join(web, "cellLineMetadata.json")

    with open(anchor_path) as f:
        anchor = json.load(f)
    with open(metadata_path) as f:
        meta = json.load(f)

    # Build a lookup: stripped-uppercase name → ACH ID
    name_to_ach = {}
    for ach, name in (meta.get("strippedCellLineName") or {}).items():
        n = normalize_name(name)
        if n:
            name_to_ach[n] = ach
    # Some lines also have alt names in cellLineName — fall back to those.
    for ach, name in (meta.get("cellLineName") or {}).items():
        n = normalize_name(name)
        if n and n not in name_to_ach:
            name_to_ach[n] = ach

    by_cl = {}
    matched = 0
    unmatched = []
    for entry in anchor["cellLines"]:
        norm = normalize_name(entry["name"])
        ach = name_to_ach.get(norm)
        if not ach:
            unmatched.append(entry["name"])
            continue
        matched += 1
        by_cl[ach] = {
            "name": entry["name"],
            "tnbcType6": entry.get("tnbcType6"),
            "tnbcType4": entry.get("tnbcType4"),
        }

    print(f"Matched {matched} / {len(anchor['cellLines'])} Lehmann cell lines to DepMap ACH IDs")
    if unmatched:
        print(f"  unmatched: {unmatched}")

    output = {
        "_description": (
            "Lehmann TNBC molecular-subtype classifications mapped to DepMap "
            "ACH IDs. Two label systems carried per cell line: tnbcType6 (the "
            "original 2011 six-class scheme: BL1, BL2, IM, M, MSL, LAR) and "
            "tnbcType4 (the 2016 reclassification that collapsed IM and MSL "
            "to immune / stromal contamination signal — IM lines have null "
            "tnbcType4)."
        ),
        "schemaVersion": "1.0",
        "source_2011": anchor.get("source_2011"),
        "source_2016": anchor.get("source_2016"),
        "subtypes_2011": anchor.get("subtypes_2011"),
        "subtypes_2016": anchor.get("subtypes_2016"),
        "byCellLine": by_cl,
    }
    out_path = os.path.join(web, "lehmann_tnbc.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"Written {out_path} ({size_kb:.1f} KB) — {len(by_cl)} cell lines")


if __name__ == "__main__":
    main()
