#!/usr/bin/env python3
"""
Build web_data/validated_fusions.json — the curated/validated fusion set used by
the mutational-analysis "Fusion (validated)" axis.

Source: web_data/clinical_fusions.json (curated driver fusions, validated per cell
line by tissue match + partner expression z-score + partner CRISPR dependency
z-score; each call carries a tier in {high, medium, low}). This is the same data
the Cell Line Browser shows under "Clinically relevant fusions".

We keep only HIGH + MEDIUM tier calls (drop the handful of low-tier ones), index
each call under both participating genes (split the fusion name on '-'), and emit
a structure that MIRRORS translocations.json exactly so the existing fusion-axis
code consumes it unchanged:

{
  "genes": [...],                       # alphabetical (UI sorts by count itself)
  "geneCounts": {gene: N_lines},
  "geneData": {
    GENE: {
      "translocations": {cellLine: n_partners, ...},  # validated fusions involving GENE
      "partners":       {cellLine: [partnerGene, ...]},
      "counts": {"0": N_wt, "1": N, "2": N},
      "total_translocated": N
    }
  }
}

Unlike process_translocations.py there is NO minimum-cell-line filter — every
validated fusion is curated and kept, however rare.

Usage:
    python build_validated_fusions.py
"""

import json
import os
from collections import defaultdict

ALLOWED_TIERS = {"high", "medium"}


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    metadata_file = os.path.join(project_dir, "web_data", "metadata.json")
    clinical_file = os.path.join(project_dir, "web_data", "clinical_fusions.json")
    output_file = os.path.join(project_dir, "web_data", "validated_fusions.json")

    with open(metadata_file) as f:
        cohort = set(json.load(f)["cellLines"])
    print(f"Loaded {len(cohort)} cell lines from metadata.json")

    with open(clinical_file) as f:
        clinical = json.load(f)
    fusion_data = clinical["fusionData"]
    print(f"{len(fusion_data)} curated fusions in clinical_fusions.json")

    # gene -> cell line -> set of partner genes (from high/medium-tier calls only)
    gene_partners = defaultdict(lambda: defaultdict(set))
    n_calls = 0
    for fusion_name, fd in fusion_data.items():
        genes = [g.strip() for g in fusion_name.split("-") if g.strip()]
        if len(genes) < 2:
            print(f"  WARNING: could not split fusion name into genes: {fusion_name!r}")
        for cell_line, info in fd.get("cellLines", {}).items():
            if info.get("tier") not in ALLOWED_TIERS:
                continue
            if cell_line not in cohort:
                continue
            n_calls += 1
            for g in genes:
                partners = set(genes) - {g}
                gene_partners[g][cell_line].update(partners)

    print(f"{n_calls} high/medium-tier per-line calls kept")

    genes_sorted = sorted(gene_partners.keys())
    gene_counts = {}
    gene_data = {}
    for gene in genes_sorted:
        cell_data = gene_partners[gene]
        total = len(cell_data)
        gene_counts[gene] = total
        translocations = {}
        partners = {}
        counts = defaultdict(int)
        for cell_line, partner_set in cell_data.items():
            n = len(partner_set)
            translocations[cell_line] = n
            partners[cell_line] = sorted(partner_set)
            counts[str(min(n, 2))] += 1
        gene_data[gene] = {
            "translocations": translocations,
            "partners": partners,
            "counts": {
                "0": len(cohort) - total,
                "1": counts.get("1", 0),
                "2": counts.get("2", 0),
            },
            "total_translocated": total,
        }

    output = {"genes": genes_sorted, "geneCounts": gene_counts, "geneData": gene_data}
    with open(output_file, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"Written {output_file}")
    print(f"  {len(genes_sorted)} genes")
    print(f"  Top genes by line count: {sorted(gene_counts.items(), key=lambda x: -x[1])[:12]}")


if __name__ == "__main__":
    main()
