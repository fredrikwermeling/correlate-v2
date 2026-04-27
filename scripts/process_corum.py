#!/usr/bin/env python3
"""
Build web_data/corum_partners.json from a CORUM human-complexes table.

CORUM (https://mips.helmholtz-muenchen.de/corum/) is a manually curated
catalog of mammalian protein complexes. The Correlate AI export uses it to
populate the focal-gene pathway-partner expression whitelist for arbitrary
focal genes — without it, only the ~9 hand-curated complexes (NEDD8,
Proteasome, Hippo, MYC, TP53, BRCA, mTOR, BCL2, splicing) plus the wiki
cancer pathways have partner expression carried through the variance filter.
With CORUM, ~3,000 human genes pick up their complex co-members
automatically (SMARCA4 → BAF subunits, MCM4 → MCM2-7 helicase,
RUVBL1 → INO80/SRCAP, etc.).

Input
  Path to a CORUM human-complexes file. Two accepted forms:
    - .txt / .tsv  : tab-separated, downloaded directly from CORUM
    - .zip         : the downloadable zip from CORUM containing the TSV
  Download from: https://mips.helmholtz-muenchen.de/corum/download
  Pick the "Human" complexes file.

Output
  web_data/corum_partners.json
    {
      "_description": "...",
      "schemaVersion": "1.0",
      "source": "CORUM <release>",
      "nGenes": <count>,
      "nComplexes": <count>,
      "partners": {
        "SMARCA4": ["ARID1A","ARID1B","SMARCB1","SMARCC1", ...],
        "NEDD8":   ["UBA3","NAE1","CUL1","CUL2","RBX1", ...],
        ...
      }
    }

Usage
  python scripts/process_corum.py /path/to/humanComplexes.txt
  python scripts/process_corum.py /path/to/corum_humanComplexes.zip
"""

import csv
import io
import json
import os
import re
import sys
import zipfile
from collections import defaultdict


# Possible column-name variants for the gene-symbol member list — CORUM
# changes the header subtly between releases.
GENE_COL_CANDIDATES = [
    "subunits(Gene name)",
    "subunits(Gene names)",
    "subunits Gene name",
    "subunits.Gene.name.",
    "subunits.Gene.name",
    "gene_name",
    "GeneName",
    "subunits(GeneName)",
    "subunits_gene_name",
]

ORGANISM_COL_CANDIDATES = ["Organism", "organism"]


def open_tsv(path):
    """Return an iterable of rows (dicts) regardless of zip / plain TSV."""
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            # Find the first .txt or .tsv inside.
            inner = [n for n in zf.namelist() if n.lower().endswith((".txt", ".tsv"))]
            if not inner:
                sys.exit(f"No .txt or .tsv inside {path}")
            with zf.open(inner[0]) as f:
                text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
                reader = csv.DictReader(text, delimiter="\t")
                rows = list(reader)
                cols = reader.fieldnames
        return cols, rows
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        cols = reader.fieldnames
    return cols, rows


def detect(cols, candidates):
    norm = {c.strip().lower(): c for c in cols}
    for cand in candidates:
        key = cand.strip().lower()
        if key in norm:
            return norm[key]
    return None


def parse_gene_list(value):
    """CORUM separates gene symbols with semicolons. Strip whitespace and
    drop empties / numeric-only entries (which sometimes leak through when
    the source field accidentally has Entrez IDs)."""
    if not value:
        return []
    parts = re.split(r"[;,]\s*", value)
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        # Skip numeric-only entries (Entrez ID leakage)
        if p.isdigit():
            continue
        # Skip unmapped placeholders
        if p in {"-", "NA", "None", "?"}:
            continue
        out.append(p.upper())
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_corum.py <humanComplexes.txt | corum.zip>")
    src = sys.argv[1]
    if not os.path.exists(src):
        sys.exit(f"File not found: {src}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")

    print(f"Reading {src} ...")
    cols, rows = open_tsv(src)
    if not cols:
        sys.exit("Could not parse columns from the file.")
    print(f"  columns: {cols[:8]}{' ...' if len(cols) > 8 else ''}")

    gene_col = detect(cols, GENE_COL_CANDIDATES)
    if not gene_col:
        sys.exit(
            "Could not find a gene-symbol member column. Tried: "
            + ", ".join(GENE_COL_CANDIDATES)
            + f". Got: {cols}"
        )
    organism_col = detect(cols, ORGANISM_COL_CANDIDATES)
    print(f"  using gene column: {gene_col!r}")
    if organism_col:
        print(f"  using organism column: {organism_col!r}")

    # Build gene → set of co-complex members.
    gene_partners = defaultdict(set)
    n_complexes = 0
    n_skipped_organism = 0
    for row in rows:
        if organism_col:
            org = (row.get(organism_col) or "").strip().lower()
            # Skip if this is a multi-organism file and the row isn't human.
            # CORUM's "Human" download is already filtered, but be safe.
            if org and "human" not in org and "homo sapiens" not in org:
                n_skipped_organism += 1
                continue
        members = parse_gene_list(row.get(gene_col))
        if len(members) < 2:
            continue
        n_complexes += 1
        for g in members:
            for m in members:
                if m != g:
                    gene_partners[g].add(m)

    if n_skipped_organism:
        print(f"  skipped {n_skipped_organism} non-human rows")
    print(f"  {n_complexes} complexes ingested, {len(gene_partners)} genes with at least one partner")

    # Sort partners deterministically and convert to lists.
    partners = {g: sorted(s) for g, s in sorted(gene_partners.items())}

    # Sanity-check key cancer / cell-bio genes.
    sanity = ["NEDD8", "SMARCA4", "MCM4", "RUVBL1", "BRCA1", "CDK4", "TP53", "MYC"]
    print("\nSanity check (sample focal genes → partner count):")
    for g in sanity:
        n = len(partners.get(g, []))
        sample = ", ".join(partners.get(g, [])[:6])
        print(f"  {g:10s} n={n:3d}  ({sample}{'...' if n > 6 else ''})")

    output = {
        "_description": (
            "Per-gene complex co-members from CORUM, used by the Correlate AI "
            "export to carry pathway-partner expression through the variance "
            "filter for arbitrary focal genes."
        ),
        "schemaVersion": "1.0",
        "source": f"CORUM ({os.path.basename(src)})",
        "nGenes": len(partners),
        "nComplexes": n_complexes,
        "partners": partners,
    }
    out_path = os.path.join(web, "corum_partners.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
