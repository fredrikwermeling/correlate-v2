#!/usr/bin/env python3
"""
Generate web_data/clinical_fusions.json from a curated anchor list and orthogonal
evidence (lineage match + partner expression z-score + partner dependency z-score).

Input:
  - <OmicsFusionFiltered.csv> from DepMap (same file used by process_translocations.py)
  - scripts/clinical_fusions_anchor.json (curated driver list)
  - web_data/cellLineMetadata.json (lineage)
  - web_data/expression_metadata.json + web_data/expression.bin.gz (int16 / scaleFactor)
  - web_data/metadata.json + web_data/geneEffects.bin.gz (int16 / scaleFactor)

Output:
  web_data/clinical_fusions.json
    {
      "_description": ...,
      "schemaVersion": "1.0",
      "anchors": [...]                  // curated list passed through
      "fusions": [...],                 // canonical fusion names emitted at least once
      "fusionData": {                   // per-fusion: cell lines + tier + evidence
        "BCR-ABL1": {
          "cellLines": {
            "ACH-000076": {
              "tier": "high",
              "atypicalLineage": false,
              "lineage": "Myeloid",
              "exprZ": 0.42,
              "depZ": -2.1
            }, ...
          },
          "expectedLineages": ["Myeloid","Lymphoid"],
          "validatePartner": "ABL1",
          "diseaseContext": "CML, Ph+ B-ALL"
        }
      },
      "byCellLine": {                   // per-cell-line index for the cell line browser
        "ACH-000076": [
          { "fusion":"BCR-ABL1","tier":"high","atypicalLineage":false }
        ]
      }
    }

Tier rules (lineage is a modifier, not a gate):
  curated + expression match + dependency match     -> high
  curated + (expression OR dependency) match        -> medium
  curated + neither support but lineage match       -> medium
  curated + neither support and lineage mismatch    -> low
  atypicalLineage flag set whenever lineage doesn't match expected list.

A signal "matches" when its z-score crosses a threshold:
  expression: exprZ >= +1.5  (partner overexpressed vs same-lineage negatives)
  dependency: depZ  <= -1.0  (partner more essential than baseline)

Usage:
  python scripts/process_clinical_fusions.py <path_to_OmicsFusionFiltered.csv>
"""

import csv
import gzip
import json
import math
import os
import struct
import sys
from collections import defaultdict


EXPR_Z_THRESHOLD = 1.5     # partner overexpression
DEP_Z_THRESHOLD = -1.0     # partner more essential than baseline (lower GE = more essential)


def load_int16_matrix(bin_path, n_rows, n_cols, scale_factor, na_value):
    """Read int16 little-endian matrix, return list-of-lists in float space.
    Layout: row-major, gene-major (row = gene, col = cell line)."""
    with gzip.open(bin_path, "rb") as f:
        raw = f.read()
    expected = n_rows * n_cols * 2
    if len(raw) != expected:
        sys.exit(f"Binary size mismatch for {bin_path}: got {len(raw)}, expected {expected}")
    arr = struct.unpack(f"<{n_rows * n_cols}h", raw)
    out = []
    for r in range(n_rows):
        row = []
        base = r * n_cols
        for c in range(n_cols):
            v = arr[base + c]
            row.append(None if v == na_value else v / scale_factor)
        out.append(row)
    return out


def zscore_for_value(value, others):
    """Z-score of value vs a list of float values (None entries skipped)."""
    pool = [x for x in others if x is not None and not math.isnan(x)]
    if value is None or math.isnan(value) or len(pool) < 5:
        return None
    mean = sum(pool) / len(pool)
    var = sum((x - mean) ** 2 for x in pool) / len(pool)
    sd = math.sqrt(var)
    if sd <= 0:
        return None
    return (value - mean) / sd


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_clinical_fusions.py <OmicsFusionFiltered.csv>")
    fusion_csv = sys.argv[1]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")
    anchor_path = os.path.join(script_dir, "clinical_fusions_anchor.json")

    # ---------- load anchor list ----------
    with open(anchor_path) as f:
        anchor = json.load(f)
    anchors = list(anchor["fusions"])
    kmt2a_partners = anchor.get("kmt2aPartners", [])

    # ---------- load metadata ----------
    with open(os.path.join(web, "cellLineMetadata.json")) as f:
        cl_meta = json.load(f)
    lineage_map = cl_meta.get("lineage", {})
    name_map = cl_meta.get("strippedCellLineName", {})
    disease_map = cl_meta.get("primaryDisease", {})

    with open(os.path.join(web, "metadata.json")) as f:
        ge_meta = json.load(f)
    ge_genes = ge_meta["genes"]
    ge_cls = ge_meta["cellLines"]
    ge_n_genes = ge_meta["nGenes"]
    ge_n_cls = ge_meta["nCellLines"]
    ge_scale = ge_meta["scaleFactor"]
    ge_na = ge_meta["naValue"]
    ge_gene_idx = {g: i for i, g in enumerate(ge_genes)}
    ge_cl_idx = {cl: i for i, cl in enumerate(ge_cls)}

    with open(os.path.join(web, "expression_metadata.json")) as f:
        ex_meta = json.load(f)
    ex_genes = ex_meta["genes"]
    ex_cls = ex_meta["cellLines"]
    ex_n_genes = ex_meta["nGenes"]
    ex_n_cls = ex_meta["nCellLines"]
    ex_scale = ex_meta["scaleFactor"]
    ex_na = ex_meta["naValue"]
    ex_gene_idx = {g: i for i, g in enumerate(ex_genes)}
    ex_cl_idx = {cl: i for i, cl in enumerate(ex_cls)}

    valid_cell_lines = set(ge_cls)

    # ---------- helpers ----------
    def get_ge_row(gene):
        idx = ge_gene_idx.get(gene)
        if idx is None:
            return None
        return ge_data[idx]

    def get_ex_row(gene):
        idx = ex_gene_idx.get(gene)
        if idx is None:
            return None
        return ex_data[idx]

    print(f"Loading expression matrix ({ex_n_genes} genes x {ex_n_cls} cell lines)...")
    ex_data = load_int16_matrix(
        os.path.join(web, "expression.bin.gz"), ex_n_genes, ex_n_cls, ex_scale, ex_na
    )
    print(f"Loading gene effect matrix ({ge_n_genes} genes x {ge_n_cls} cell lines)...")
    ge_data = load_int16_matrix(
        os.path.join(web, "geneEffects.bin.gz"), ge_n_genes, ge_n_cls, ge_scale, ge_na
    )

    # ---------- parse raw fusion calls ----------
    # Detect supplementary (Arriba-annotated) vs standard format. The
    # supplementary file carries `confidence` (high/medium/low), `type`
    # (translocation/inversion/deletion/read-through/...) and `reading_frame`
    # — we use these to filter out the noise the user described (passenger
    # junctions in chromothriptic regions, read-through transcripts, etc.).
    print(f"Parsing {fusion_csv} ...")
    fusion_pairs = defaultdict(set)  # (g5,g3) -> set(cellLineIDs)
    n_seen = n_kept = 0
    drop_low_conf = drop_readthrough = 0

    with open(fusion_csv, newline="") as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        model_col = next(
            (c for c in cols if c in ("ModelID", "DepMap_ID", "CellLineID", "ModelId")),
            cols[0],
        )
        gene1_col = gene2_col = fusion_name_col = None
        confidence_col = type_col = None
        for c in cols:
            cl_norm = c.lower().replace("_", "").replace(" ", "")
            if cl_norm.startswith("gene1") or cl_norm in ("leftgene", "genea", "leftgenename"):
                gene1_col = c
            elif cl_norm.startswith("gene2") or cl_norm in ("rightgene", "geneb", "rightgenename"):
                gene2_col = c
            elif cl_norm in ("fusionname", "fusion", "name", "canonicalfusionname"):
                fusion_name_col = c
            elif cl_norm == "confidence":
                confidence_col = c
            elif cl_norm == "type":
                type_col = c

        is_supplementary = confidence_col is not None
        print(
            "  format: "
            + ("supplementary (Arriba-annotated)" if is_supplementary else "standard")
        )

        for row in reader:
            n_seen += 1
            cl = row[model_col].strip()
            if cl not in valid_cell_lines:
                continue

            # Quality gates only meaningful in supplementary
            if is_supplementary:
                conf = (row.get(confidence_col) or "").strip().lower()
                if conf and conf != "high":
                    drop_low_conf += 1
                    continue
                fusion_type = (row.get(type_col) or "").strip().lower() if type_col else ""
                # Read-through transcripts are adjacent-gene co-transcription,
                # not real DNA rearrangements. Skip them for clinical calls.
                if "read-through" in fusion_type or "readthrough" in fusion_type:
                    drop_readthrough += 1
                    continue

            if gene1_col and gene2_col:
                g1 = row[gene1_col].split(" (")[0].strip()
                g2 = row[gene2_col].split(" (")[0].strip()
            elif fusion_name_col:
                fn = row[fusion_name_col].strip()
                sep = "--" if "--" in fn else "-"
                parts = fn.split(sep, 1)
                if len(parts) != 2:
                    continue
                g1, g2 = parts[0].strip(), parts[1].strip()
            else:
                continue
            if not g1 or not g2:
                continue
            fusion_pairs[(g1, g2)].add(cl)
            n_kept += 1

    print(
        f"  {n_seen} rows seen, {n_kept} kept after quality filter "
        f"(dropped {drop_low_conf} non-high-confidence, {drop_readthrough} read-through)"
    )
    print(f"  {len(fusion_pairs)} distinct ordered fusion pairs")

    # ---------- evaluate each anchor ----------
    fusion_data = {}
    by_cell_line = defaultdict(list)

    def evaluate(fusion_name, gene5, gene3, validate_partner, validate_expr,
                 validate_dep, expected_lineages, disease_context):
        # Match either orientation
        cls = set(fusion_pairs.get((gene5, gene3), set())) | set(fusion_pairs.get((gene3, gene5), set()))
        if not cls:
            return
        # Build pools for z-score calculation
        partner_ex_row = get_ex_row(validate_partner) if validate_expr else None
        partner_ge_row = get_ge_row(validate_partner) if validate_dep else None

        per_cell = {}
        for cl in sorted(cls):
            actual_lineage = lineage_map.get(cl, "")
            atypical = bool(expected_lineages) and (actual_lineage not in expected_lineages)

            expr_z = None
            if partner_ex_row is not None:
                idx = ex_cl_idx.get(cl)
                if idx is not None:
                    val = partner_ex_row[idx]
                    # pool: lineage-matched negatives if available, else all negatives
                    pool_cls = [c for c in ex_cls if c not in cls and lineage_map.get(c) == actual_lineage]
                    if len(pool_cls) < 5:
                        pool_cls = [c for c in ex_cls if c not in cls]
                    others = [partner_ex_row[ex_cl_idx[c]] for c in pool_cls]
                    expr_z = zscore_for_value(val, others)

            dep_z = None
            if partner_ge_row is not None:
                idx = ge_cl_idx.get(cl)
                if idx is not None:
                    val = partner_ge_row[idx]
                    pool_cls = [c for c in ge_cls if c not in cls and lineage_map.get(c) == actual_lineage]
                    if len(pool_cls) < 5:
                        pool_cls = [c for c in ge_cls if c not in cls]
                    others = [partner_ge_row[ge_cl_idx[c]] for c in pool_cls]
                    dep_z = zscore_for_value(val, others)

            expr_match = expr_z is not None and expr_z >= EXPR_Z_THRESHOLD
            dep_match = dep_z is not None and dep_z <= DEP_Z_THRESHOLD

            # Tier rules — only count signals the curator marked as informative
            # (validateExpression / validateDependency). Lineage is a modifier,
            # not a gate: strong orthogonal evidence keeps the call even when
            # the cell line's tissue is atypical for the canonical disease.
            n_enabled = int(validate_expr) + int(validate_dep)
            n_passed = (int(expr_match) if validate_expr else 0) + (int(dep_match) if validate_dep else 0)

            if n_enabled == 2:
                if n_passed == 2:
                    tier = "high"
                elif n_passed == 1:
                    tier = "medium"
                else:
                    tier = "medium" if not atypical else "low"
            elif n_enabled == 1:
                if n_passed == 1:
                    tier = "high"
                else:
                    tier = "medium" if not atypical else "low"
            else:  # neither signal informative — call rests on curated + lineage only
                tier = "medium" if not atypical else "low"

            per_cell[cl] = {
                "tier": tier,
                "atypicalLineage": atypical,
                "lineage": actual_lineage,
                "exprZ": None if expr_z is None else round(expr_z, 2),
                "depZ": None if dep_z is None else round(dep_z, 2),
            }
            by_cell_line[cl].append({
                "fusion": fusion_name,
                "tier": tier,
                "atypicalLineage": atypical,
            })

        fusion_data[fusion_name] = {
            "cellLines": per_cell,
            "expectedLineages": expected_lineages,
            "validatePartner": validate_partner,
            "diseaseContext": disease_context,
        }

    for a in anchors:
        evaluate(
            fusion_name=a["fusion"], gene5=a["gene5"], gene3=a["gene3"],
            validate_partner=a["validatePartner"],
            validate_expr=a.get("validateExpression", False),
            validate_dep=a.get("validateDependency", False),
            expected_lineages=a.get("expectedLineages", []),
            disease_context=a.get("diseaseContext", ""),
        )

    # KMT2A — variable 3' partner, KMT2A is always 5'. Validation = KMT2A dependency.
    for k in kmt2a_partners:
        evaluate(
            fusion_name=k["fusion"], gene5="KMT2A", gene3=k["partner"],
            validate_partner="KMT2A",
            validate_expr=False, validate_dep=True,
            expected_lineages=k.get("expectedLineages", []),
            disease_context=k.get("diseaseContext", ""),
        )

    # ---------- write output ----------
    fusions_emitted = sorted(fusion_data.keys())
    output = {
        "_description": (
            "Clinical driver fusions called by process_clinical_fusions.py. Each call has "
            "a confidence tier ('high', 'medium', 'low') and an atypicalLineage flag. "
            "Lineage is a modifier, not a gate — strong orthogonal evidence keeps the call "
            "even if the cell line's tissue doesn't match the canonical disease context."
        ),
        "schemaVersion": "1.0",
        "thresholds": {
            "expressionZ": EXPR_Z_THRESHOLD,
            "dependencyZ": DEP_Z_THRESHOLD,
        },
        "anchors": anchors + [
            {
                "fusion": k["fusion"], "gene5": "KMT2A", "gene3": k["partner"],
                "validatePartner": "KMT2A", "validateExpression": False, "validateDependency": True,
                "expectedLineages": k.get("expectedLineages", []),
                "diseaseContext": k.get("diseaseContext", ""),
            }
            for k in kmt2a_partners
        ],
        "fusions": fusions_emitted,
        "fusionData": fusion_data,
        "byCellLine": dict(by_cell_line),
    }

    out_path = os.path.join(web, "clinical_fusions.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    # ---------- summary ----------
    n_calls = sum(len(v["cellLines"]) for v in fusion_data.values())
    tier_counts = defaultdict(int)
    atypical_count = 0
    for fd in fusion_data.values():
        for cl_info in fd["cellLines"].values():
            tier_counts[cl_info["tier"]] += 1
            if cl_info["atypicalLineage"]:
                atypical_count += 1

    print(f"\nWritten {out_path}")
    print(f"  {len(fusions_emitted)} fusions emitted across {n_calls} cell-line calls")
    print(f"  tier breakdown: high={tier_counts['high']}, medium={tier_counts['medium']}, low={tier_counts['low']}")
    print(f"  atypical-lineage calls: {atypical_count}")
    print(f"\nTop 10 fusions by call count:")
    for fname in sorted(fusion_data.keys(), key=lambda f: -len(fusion_data[f]["cellLines"]))[:10]:
        n = len(fusion_data[fname]["cellLines"])
        tiers = defaultdict(int)
        for c in fusion_data[fname]["cellLines"].values():
            tiers[c["tier"]] += 1
        print(f"  {fname}: n={n}  high={tiers['high']} medium={tiers['medium']} low={tiers['low']}")


if __name__ == "__main__":
    main()
