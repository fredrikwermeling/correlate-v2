#!/usr/bin/env python3
"""
Build web_data/derivative_lines.json: cell lines that are a treated or
engineered version of another line in the same panel.

The panel carries a handful of matched pairs, a parental line and a version of
it selected under a drug, and nothing in the app said so. They are exactly the
comparison a resistance experiment wants, and exactly the pair that will look
like a striking "difference between two cell lines" if you do not know they
are the same line twice.

The pairing is read off the names, which follow DepMap's convention of the
parental name plus what was done and a suffix (A-375_DAB_R, HCC-827-GR5,
A549_CRAF_KD). A suffix is only accepted when every token in it is recognised,
so SK-N-SH is not filed as a derivative of SKN and HT-29 is not filed as a
derivative of HT.

Two kinds of change are distinguished, because they are not the same
experiment: a line GROWN IN a drug until it survived it, and a line with a
GENE KNOCKED DOWN in it. A549_CRAF_KD is the second kind, and describing it as
a drug exposure (which an earlier version did, by keeping drugs and genes in
one dictionary) is simply wrong.

    python3 scripts/build_derivative_lines.py

Output: web_data/derivative_lines.json
    { "_description": ..., "built": "YYYY-MM-DD",
      "pairs": [ {parent, parentName, derivative, derivativeName,
                  agents: [...], relationship: "drug-resistant" | "engineered",
                  note} ],
      "byCellLine": { "ACH-...": {role: "parent"|"derivative", partner: "ACH-..."} } }
"""

import json
import os
import re
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
WEB = os.path.join(os.path.dirname(HERE), "web_data")

# Drug tokens that can appear in a derivative's suffix. These are compounds
# the line was grown in.
AGENTS = {
    "DAB": "dabrafenib (BRAF inhibitor)",
    "TRAM": "trametinib (MEK inhibitor)",
    "VEM": "vemurafenib (BRAF inhibitor)",
    "ROX": "roxadustat",
    "SCH772984": "SCH772984 (ERK inhibitor)",
    "GR": "gefitinib (EGFR inhibitor), selected for resistance",
}
# Gene tokens. A gene in the suffix is the thing that was knocked down, NOT
# something the cells were treated with, and keeping the two in one dictionary
# is how A549_CRAF_KD came out described as a drug exposure.
GENES = {
    "CRAF": "RAF1 (CRAF)",
}
# Suffix tokens that are neither, marking what was done.
MARKERS = {
    "R": "selected for resistance",
    "KD": "knocked down",
}


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", (s or "").upper())


def main():
    meta = json.load(open(os.path.join(WEB, "cellLineMetadata.json")))
    names = meta["cellLineName"]
    lines = meta["cellLines"]

    by_norm = {}
    for cl in lines:
        by_norm.setdefault(norm(names.get(cl, "")), cl)

    pairs = []
    for cl in lines:
        raw = (names.get(cl, "") or "").upper()
        # Longest parent wins: HT-144_DAB_R is a derivative of HT-144, not of HT.
        for i in range(len(raw) - 1, 0, -1):
            if raw[i] not in "_- ":
                continue
            parent = by_norm.get(norm(raw[:i]))
            suffix = raw[i + 1:]
            if not parent or parent == cl or not suffix:
                continue
            tokens = [t for t in re.split(r"[_\- ]+", suffix) if t]
            agents, genes, markers, ok = [], [], [], True
            for t in tokens:
                base = re.sub(r"\d+$", "", t) or t
                if t in AGENTS or base in AGENTS:
                    label = AGENTS.get(t) or AGENTS[base]
                    # A token like GR5 carries both the drug and the fact that
                    # the line was selected under it.
                    if ", selected for resistance" in label:
                        label = label.replace(", selected for resistance", "")
                        markers.append("selected for resistance")
                    agents.append(label)
                elif t in GENES or base in GENES:
                    genes.append(GENES.get(t) or GENES[base])
                elif t in MARKERS or base in MARKERS:
                    markers.append(MARKERS.get(t) or MARKERS[base])
                else:
                    ok = False
                    break
            if not ok or not markers or not (agents or genes):
                break  # longest parent matched but the suffix says nothing known
            resistant = "selected for resistance" in markers
            knockdown = "knocked down" in markers
            if knockdown and genes:
                relationship = "gene-knockdown"
                change = "had " + " and ".join(genes) + " knocked down"
                changeShort = " and ".join(genes) + " knockdown"
            elif resistant and agents:
                relationship = "drug-resistant"
                change = "was grown in " + " and ".join(agents) + " until it survived it"
                changeShort = "resistant to " + " and ".join(agents)
            else:
                relationship = "derived"
                change = "was derived under " + " and ".join(agents + genes)
                changeShort = " and ".join(agents + genes)
            pairs.append({
                "parent": parent,
                "parentName": names[parent],
                "derivative": cl,
                "derivativeName": names[cl],
                "relationship": relationship,
                "agents": agents,          # drugs the line was grown in
                "genesTargeted": genes,    # genes knocked down in it
                "change": change,
                "changeShort": changeShort,
                "note": (
                    f"{names[cl]} is {names[parent]} that {change}. The same line "
                    "twice, not two lines: a difference between them is what that "
                    "change did, and anything they share is the background they "
                    "both came from."
                ),
            })
            break

    # A parent can have several derivatives (A-375 has three), so the partner
    # is a list; keeping one silently hid the others.
    by_cl = {}
    for p in pairs:
        by_cl.setdefault(p["parent"], {"role": "parent", "partners": []})["partners"].append(p["derivative"])
        by_cl.setdefault(p["derivative"], {"role": "derivative", "partners": []})["partners"].append(p["parent"])

    payload = {
        "_description": (
            "Cell lines that are a treated or engineered version of another line "
            "in the same panel, paired with that parental line. Read off the "
            "DepMap names, which follow the convention of the parental name plus "
            "an agent and a suffix; a suffix is only accepted when every token in "
            "it names a recognised agent or marker, so unrelated lines that merely "
            "share a name prefix are not paired."
        ),
        "built": date.today().isoformat(),
        "agents": AGENTS,
        "genes": GENES,
        "pairs": pairs,
        "byCellLine": by_cl,
    }
    dest = os.path.join(WEB, "derivative_lines.json")
    json.dump(payload, open(dest, "w"))
    print(f"{len(pairs)} pairs")
    for p in pairs:
        print(f"  {p['parentName']:12s} -> {p['derivativeName']:22s} [{p['relationship']}] {p['changeShort']}")
    print("->", dest)


if __name__ == "__main__":
    main()
