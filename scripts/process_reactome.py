#!/usr/bin/env python3
"""
Build web_data/reactome_partners.json by fetching Reactome's pathway GMT
and computing per-gene pathway co-members.

Reactome (https://reactome.org/) is a curated database of human pathways
covering signalling cascades, metabolism, transcription, gene expression,
etc. — broader than CORUM (which is specifically physical complexes).

Used by the Correlate AI export as a second tier in the focal-gene partner
lookup, layered after CORUM. Catches pathway/cascade relationships CORUM
misses — e.g. IL4R → JAK1, STAT6 (signalling cascade, not a stable complex).

Filtering:
  - Pathway size 2–150 genes. Smaller is degenerate; larger pathways are
    typically broad parents like "Signal Transduction" with ~2,600 members,
    where being in the same pathway tells you almost nothing about
    functional relatedness.
  - Per-gene partner cap (default 250) — for very promiscuous genes, keep
    the partners shared across the most pathways (most-co-occurring), drop
    the long tail.

Source: https://reactome.org/download/current/ReactomePathways.gmt.zip

Output: web_data/reactome_partners.json
"""

import io
import json
import os
import ssl
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict, Counter


GMT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
USER_AGENT = "correlate-v2/process_reactome.py"

# Reactome's CA chain isn't always in the default Python bundle.
_SSL_CTX = ssl._create_unverified_context()

MIN_PATHWAY_SIZE = 5
MAX_PATHWAY_SIZE = 100
MAX_PARTNERS_PER_GENE = 120
# Tightened after a first run produced 11 MB of mostly noise: hub
# genes (TP53, MYC, BRCA1, JAK1) all maxed out at 250 partners,
# pulling in genes that share a broad parent pathway (e.g. "Cytokine
# Signaling in Immune system" with ~80 members) without functional
# meaning. Specific pathways live in the 5-50 range; everything bigger
# is typically a parent grouping that conveys little partner signal.


def fetch_gmt():
    """Return list-of-lists: each pathway = [name, reactome_id, ...gene_symbols]."""
    print(f"Fetching {GMT_URL} ...")
    req = urllib.request.Request(GMT_URL, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60, context=_SSL_CTX) as r:
        raw = r.read()
    print(f"  downloaded {len(raw)} bytes")
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        inner = [n for n in zf.namelist() if n.lower().endswith(".gmt")]
        if not inner:
            sys.exit("No .gmt file inside the zip.")
        with zf.open(inner[0]) as f:
            text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
            pathways = []
            for line in text:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    pathways.append(parts)
    return pathways


def from_gmt_file(path):
    print(f"Reading {path} ...")
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            inner = [n for n in zf.namelist() if n.lower().endswith(".gmt")]
            if not inner:
                sys.exit("No .gmt file inside the zip.")
            with zf.open(inner[0]) as f:
                text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
                pathways = [line.rstrip("\n").split("\t") for line in text if line.strip()]
        return [p for p in pathways if len(p) >= 3]
    with open(path, encoding="utf-8", errors="replace") as f:
        pathways = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    return [p for p in pathways if len(p) >= 3]


def main():
    if len(sys.argv) >= 2 and sys.argv[1] not in ("--help", "-h"):
        all_pathways = from_gmt_file(sys.argv[1])
    else:
        all_pathways = fetch_gmt()
    print(f"  {len(all_pathways)} pathways total")

    # Filter pathway size and collect per-gene partner counters.
    kept = []
    skipped_too_small = 0
    skipped_too_big = 0
    for parts in all_pathways:
        name, rid, *genes = parts
        # GMT may include duplicates / empty entries — clean them.
        genes = [g.strip().upper() for g in genes if g.strip()]
        genes = list(dict.fromkeys(genes))  # dedup preserve order
        if len(genes) < MIN_PATHWAY_SIZE:
            skipped_too_small += 1
            continue
        if len(genes) > MAX_PATHWAY_SIZE:
            skipped_too_big += 1
            continue
        kept.append((name, rid, genes))
    print(
        f"  kept {len(kept)} pathways "
        f"(dropped {skipped_too_small} <{MIN_PATHWAY_SIZE} genes, {skipped_too_big} >{MAX_PATHWAY_SIZE} genes)"
    )

    # gene → Counter(partner → number of shared pathways). Tracking the
    # co-occurrence count lets us cap per-gene partners at the most-shared
    # ones rather than dropping arbitrarily.
    co_counts = defaultdict(Counter)
    for _name, _rid, genes in kept:
        for g in genes:
            for m in genes:
                if m != g:
                    co_counts[g][m] += 1

    # Cap per-gene partners — keep the top MAX_PARTNERS_PER_GENE most-shared.
    partners = {}
    over_cap = 0
    for g, ctr in co_counts.items():
        if len(ctr) > MAX_PARTNERS_PER_GENE:
            over_cap += 1
            top = ctr.most_common(MAX_PARTNERS_PER_GENE)
            partners[g] = sorted(p for p, _ in top)
        else:
            partners[g] = sorted(ctr.keys())

    # Sort gene keys deterministically.
    partners = {g: partners[g] for g in sorted(partners.keys())}

    print(f"  {len(partners)} genes have at least one pathway partner")
    print(f"  {over_cap} genes hit the {MAX_PARTNERS_PER_GENE}-partner cap")

    sanity = ["IL4R", "JAK1", "STAT6", "NEDD8", "SMARCA4", "MCM4", "RUVBL1",
              "BRCA1", "TP53", "MYC", "EGFR", "KRAS"]
    print("\nSanity check (focal gene → partner count):")
    for g in sanity:
        n = len(partners.get(g, []))
        sample = ", ".join(partners.get(g, [])[:6])
        print(f"  {g:10s} n={n:3d}  ({sample}{'...' if n > 6 else ''})")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")
    out_path = os.path.join(web, "reactome_partners.json")
    output = {
        "_description": (
            "Per-gene pathway co-members from Reactome. Used by the Correlate "
            "AI export as a second tier (after CORUM) in the focal-gene "
            "partner lookup, to catch pathway / signalling-cascade "
            "relationships CORUM doesn't cover (e.g. IL4R → JAK / STAT)."
        ),
        "schemaVersion": "1.0",
        "source": f"Reactome GMT (https://reactome.org/, {time.strftime('%Y-%m-%d')})",
        "filter": {
            "minPathwaySize": MIN_PATHWAY_SIZE,
            "maxPathwaySize": MAX_PATHWAY_SIZE,
            "maxPartnersPerGene": MAX_PARTNERS_PER_GENE,
            "rationale": "Pathways larger than 150 genes (e.g. 'Signal Transduction', 'Metabolism') are broad parents whose co-membership conveys little functional information. Cap on per-gene partners keeps the JSON tight and prefers the most-co-occurring partners for promiscuous hub genes.",
        },
        "nGenes": len(partners),
        "nPathways": len(kept),
        "partners": partners,
    }
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
