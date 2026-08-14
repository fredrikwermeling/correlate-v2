#!/usr/bin/env python3
"""
Derive the V1 build (correlate.cmm.se) from this V2 repo.

V1 = V2 minus the blocks marked `V2 ONLY` that Fredrik has decided are
extended-build features, minus the version badge / changelog, minus the
"extended version" tagline, and named plain "Correlate".

app.js is copied verbatim: every removable feature degrades safely when its
markup is absent (the UMAP wiring is gated on #clbUmapSection, the analysis
basis falls back to gene effect when the radios are missing).

Usage:  python3 scripts/build_v1.py [--check]
        --check only reports what would be removed.
"""

import os
import re
import shutil
import sys

V2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V1 = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app"

# id -> human name, for the contiguous <div id=...> blocks that come out
REMOVE_DIVS = [
    ("clbUmapSection", "PCA / UMAP dimensionality reduction"),
    ("changelogModal", "changelog modal"),
    ("basisParams", "expression basis toggle"),
]

ASSETS = ["network_example.svg", "scatter_example.svg", "tsc_pathway.svg",
          "favicon-32x32.png", "favicon.ico",
          # UI images live in web_data but are build assets, not release data,
          # so they ride along with every rebuild rather than waiting for a
          # DepMap release copy.
          "web_data/clb_banner.png"]


def cut_div(lines, div_id):
    """Remove <div id="div_id"> ... </div> plus any immediately preceding
    HTML comment block, by tracking div depth."""
    start = next((i for i, l in enumerate(lines) if f'id="{div_id}"' in l), None)
    if start is None:
        return lines, 0
    # absorb a preceding comment that introduces the block
    first = start
    j = start - 1
    while j >= 0 and lines[j].strip() == "":
        j -= 1
    if j >= 0 and "-->" in lines[j]:
        k = j
        while k >= 0 and "<!--" not in lines[k]:
            k -= 1
        if k >= 0:
            first = k
    depth = 0
    end = None
    for i in range(start, len(lines)):
        depth += lines[i].count("<div") - lines[i].count("</div")
        if depth <= 0 and i > start:
            end = i
            break
    if end is None:
        raise SystemExit(f"could not find the end of #{div_id}")
    removed = end - first + 1
    return lines[:first] + lines[end + 1:], removed


def main():
    check = "--check" in sys.argv
    src = open(os.path.join(V2, "index.html")).readlines()
    report = []
    for div_id, name in REMOVE_DIVS:
        src, n = cut_div(src, div_id)
        report.append(f"  {name} (#{div_id}): {n} lines")
    html = "".join(src)

    badge = re.search(r'\s*<a href="#" id="versionBadge".*?</a>', html, re.S)
    if badge:
        html = html.replace(badge.group(0), "")
        report.append("  version badge")
    tag = re.search(r'\s*<div class="logo-tagline"[^>]*>— extended version —</div>', html)
    if tag:
        html = html.replace(tag.group(0), "")
        report.append("  'extended version' tagline")

    html = html.replace("<title>Correlate V2</title>", "<title>Correlate</title>")
    html = html.replace("<p>Correlate V2 | Part of the", "<p>Correlate | Part of the")

    print("Removing for V1:")
    print("\n".join(report))
    if check:
        print("\n--check: nothing written")
        return

    open(os.path.join(V1, "index.html"), "w").write(html)
    shutil.copy2(os.path.join(V2, "app.js"), os.path.join(V1, "app.js"))
    for a in ASSETS:
        s = os.path.join(V2, a)
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(V1, a))
    print(f"\nWrote index.html + app.js + {len(ASSETS)} assets to:\n  {V1}")
    print("Data files (web_data/) are copied separately, only when the release changes.")


if __name__ == "__main__":
    main()
