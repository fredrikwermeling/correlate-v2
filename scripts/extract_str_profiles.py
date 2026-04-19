#!/usr/bin/env python3
"""
Extract STR (Short Tandem Repeat) profiles from Cellosaurus for every cell line
in cellLineMetadata.json, keyed by RRID.

Input:  cellosaurus.txt from https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
        (user downloads offline; script expects path)
Output: augments web_data/cellLineMetadata.json with a new `strProfile` map:
            rrid → { marker: "allele1,allele2", ... }

STR markers of interest (the core ANSI/ATCC panel + extensions commonly seen
on commercial authentication reports):
  Amelogenin, CSF1PO, D3S1358, D5S818, D7S820, D8S1179, D13S317, D16S539,
  D18S51, D21S11, D19S433, D2S1338, FGA, TH01, TPOX, vWA, Penta D, Penta E.
"""

import json
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
METADATA_JSON = os.path.join(HERE, "..", "web_data", "cellLineMetadata.json")
CELLOSAURUS = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cellosaurus.txt"

STR_MARKERS = {
    "Amelogenin", "CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179",
    "D13S317", "D16S539", "D18S51", "D21S11", "D19S433", "D2S1338",
    "FGA", "TH01", "TPOX", "vWA", "Penta D", "Penta E",
    "D1S1656", "D2S441", "D10S1248", "D12S391", "D22S1045", "D6S1043"
}


def parse_cellosaurus(path):
    """Iterate Cellosaurus .txt, yield dict per cell line with {rrid, name, str}."""
    rrid = None
    name = None
    # marker → Counter of profiles (because multiple labs may report)
    marker_profiles = {}

    def flush():
        if rrid and marker_profiles:
            resolved = {}
            for m, c in marker_profiles.items():
                # Pick the most frequently reported profile; tie → longest (more alleles)
                top = c.most_common()
                if not top:
                    continue
                best = max(top, key=lambda x: (x[1], len(x[0])))
                resolved[m] = best[0]
            if resolved:
                return {"rrid": rrid, "name": name, "str": resolved}
        return None

    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("ID   "):
                name = line[5:].strip()
            elif line.startswith("AC   "):
                rrid = line[5:].strip()
            elif line.startswith("ST   ") and ":" in line:
                body = line[5:].strip()
                if body.startswith("Source(s):"):
                    continue
                # Strip source annotation: "D13S317: 12,13.3 (ATCC=CCL-2; ...)"
                m = re.match(r"^([^:]+):\s*([^()]+?)(\s*\(.*\))?$", body)
                if not m:
                    continue
                marker = m.group(1).strip()
                alleles = m.group(2).strip()
                if marker not in STR_MARKERS:
                    continue
                marker_profiles.setdefault(marker, Counter())[alleles] += 1
            elif line == "//":
                rec = flush()
                if rec:
                    yield rec
                rrid = name = None
                marker_profiles = {}


def main():
    print(f"Parsing {CELLOSAURUS}...")
    by_rrid = {}
    count = 0
    for rec in parse_cellosaurus(CELLOSAURUS):
        by_rrid[rec["rrid"]] = rec["str"]
        count += 1
    print(f"  Parsed STR profiles for {len(by_rrid)} cell lines (of {count} scanned)")

    print(f"Loading {METADATA_JSON}...")
    with open(METADATA_JSON) as f:
        meta = json.load(f)

    profiles = {}
    matched = 0
    for cl in meta["cellLines"]:
        rrid = meta.get("rrid", {}).get(cl)
        if not rrid:
            continue
        prof = by_rrid.get(rrid)
        if prof:
            profiles[cl] = prof
            matched += 1

    print(f"  Matched STR for {matched}/{len(meta['cellLines'])} DepMap cell lines")

    meta["strProfile"] = profiles
    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print(f"Wrote strProfile to {METADATA_JSON}")


if __name__ == "__main__":
    main()
