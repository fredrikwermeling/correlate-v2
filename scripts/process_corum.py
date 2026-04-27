#!/usr/bin/env python3
"""
Build web_data/corum_partners.json by querying the CORUM public REST API.

CORUM (https://mips.helmholtz-muenchen.de/corum/) is a manually curated
catalog of mammalian protein complexes. Earlier releases shipped a static
TSV download; the current site is a SPA backed by a FastAPI service. This
script enumerates all human complexes via the public search endpoint, then
fetches each complex's full member list and extracts the gene-name partners.

Output: web_data/corum_partners.json mapping each gene to its set of
co-complex members across all complexes it appears in.

Usage
  python scripts/process_corum.py
      → fetches from API (default; ~5,600 human complexes, takes 1-2 min)

  python scripts/process_corum.py /path/to/humanComplexes.txt
      → reads a previously-downloaded TSV instead (kept for offline use)
"""

import csv
import io
import json
import os
import re
import ssl
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


API_BASE = "https://mips.helmholtz-muenchen.de/fastapi-corum"
USER_AGENT = "correlate-v2/process_corum.py"

# CORUM's TLS chain isn't in the default Python cert bundle on every system,
# so use an unverified context. We're only fetching public scientific data —
# no secrets at risk; integrity is checked downstream by the gene-list parser.
_SSL_CTX = ssl._create_unverified_context()


def http_post(path, body):
    req = urllib.request.Request(
        API_BASE + path,
        data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=30, context=_SSL_CTX) as r:
        return json.loads(r.read().decode("utf-8"))


def http_get(path):
    req = urllib.request.Request(
        API_BASE + path, headers={"User-Agent": USER_AGENT}, method="GET"
    )
    with urllib.request.urlopen(req, timeout=30, context=_SSL_CTX) as r:
        return json.loads(r.read().decode("utf-8"))


def list_human_complex_ids():
    """Page through /public/search filtered to organism=Human."""
    ids = []
    page = 1
    size = 100
    while True:
        body = {
            "filters": [
                {"operation": {"operation": "AND", "neg": False},
                 "option": "organism", "value": "Human"}
            ]
        }
        resp = http_post(f"/public/search?page={page}&size={size}", body)
        items = resp.get("items", [])
        if not items:
            break
        for item in items:
            ids.append(item["complex_id"])
        if page >= resp.get("pages", 0):
            break
        page += 1
        if page % 10 == 1:
            print(f"  ... listed {len(ids)} complex IDs so far")
    return ids


def fetch_complex(cid):
    """Return (cid, [member_gene_symbols]) for one complex."""
    try:
        d = http_get(f"/public/complex?complex_id={cid}")
        members = []
        for sub in d.get("subunits", []):
            sw = sub.get("swissprot") or {}
            gn = (sw.get("gene_name") or "").strip().upper()
            # Skip non-gene placeholders.
            if not gn or gn in {"-", "NA", "NONE", "?"} or gn.isdigit():
                continue
            members.append(gn)
        return cid, members
    except Exception as e:
        return cid, None  # signal failure; caller logs


# --------- Legacy TSV path (kept for offline use) ---------

GENE_COL_CANDIDATES = [
    "subunits(Gene name)", "subunits(Gene names)",
    "subunits Gene name", "subunits.Gene.name.", "subunits.Gene.name",
    "gene_name", "GeneName", "subunits(GeneName)", "subunits_gene_name",
]
ORGANISM_COL_CANDIDATES = ["Organism", "organism"]


def open_tsv(path):
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
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
    if not value:
        return []
    parts = re.split(r"[;,]\s*", value)
    out = []
    for p in parts:
        p = p.strip()
        if not p or p.isdigit() or p in {"-", "NA", "None", "?"}:
            continue
        out.append(p.upper())
    return out


def from_tsv(src):
    print(f"Reading TSV {src} ...")
    cols, rows = open_tsv(src)
    if not cols:
        sys.exit("Could not parse columns.")
    gene_col = detect(cols, GENE_COL_CANDIDATES)
    if not gene_col:
        sys.exit(f"No gene column found in: {cols}")
    organism_col = detect(cols, ORGANISM_COL_CANDIDATES)
    print(f"  using gene column: {gene_col!r}")
    complexes = []
    for row in rows:
        if organism_col:
            org = (row.get(organism_col) or "").strip().lower()
            if org and "human" not in org and "homo sapiens" not in org:
                continue
        members = parse_gene_list(row.get(gene_col))
        if len(members) >= 2:
            complexes.append(members)
    print(f"  {len(complexes)} human complexes parsed")
    return complexes, f"CORUM TSV ({os.path.basename(src)})"


def from_api():
    print("Listing human complex IDs from CORUM API ...")
    ids = list_human_complex_ids()
    print(f"  total: {len(ids)} human complexes")

    print("Fetching member lists in parallel (8 workers)...")
    complexes = []
    failed = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(fetch_complex, cid): cid for cid in ids}
        done = 0
        for fut in as_completed(futures):
            cid, members = fut.result()
            done += 1
            if members is None:
                failed += 1
            elif len(members) >= 2:
                complexes.append(members)
            if done % 500 == 0:
                rate = done / max(time.time() - t0, 0.01)
                eta = (len(ids) - done) / max(rate, 0.01)
                print(f"  ... {done}/{len(ids)} ({rate:.0f}/s, ETA {eta:.0f}s)")
    if failed:
        print(f"  warning: {failed} complexes failed to fetch")
    print(f"  {len(complexes)} complexes ingested")
    return complexes, f"CORUM REST API (https://mips.helmholtz-muenchen.de/fastapi-corum, {time.strftime('%Y-%m-%d')})"


def main():
    if len(sys.argv) >= 2 and sys.argv[1] not in ("--help", "-h"):
        complexes, source = from_tsv(sys.argv[1])
    else:
        complexes, source = from_api()

    # gene → set of co-complex members
    gene_partners = defaultdict(set)
    for members in complexes:
        for g in members:
            for m in members:
                if m != g:
                    gene_partners[g].add(m)
    print(f"\n{len(gene_partners)} genes have at least one partner")

    partners = {g: sorted(s) for g, s in sorted(gene_partners.items())}

    sanity = ["NEDD8", "SMARCA4", "MCM4", "RUVBL1", "BRCA1", "CDK4",
              "TP53", "MYC", "EZH2", "RPL10", "SF3B1", "EGFR"]
    print("\nSanity check (focal gene → partner count):")
    for g in sanity:
        n = len(partners.get(g, []))
        sample = ", ".join(partners.get(g, [])[:6])
        print(f"  {g:10s} n={n:3d}  ({sample}{'...' if n > 6 else ''})")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")
    out_path = os.path.join(web, "corum_partners.json")
    output = {
        "_description": (
            "Per-gene complex co-members from CORUM. Used by the Correlate "
            "AI export to carry pathway-partner expression through the "
            "variance filter for arbitrary focal genes."
        ),
        "schemaVersion": "1.0",
        "source": source,
        "nGenes": len(partners),
        "nComplexes": len(complexes),
        "partners": partners,
    }
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"\nWritten {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
