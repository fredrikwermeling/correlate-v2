"""
Culture medium per cell line -> web_data/cellLineMetadata.json["culture"].

Why: dependencies on nutrient-handling genes (iron uptake, purine / folate /
one-carbon, serine-glycine, cystine, pyruvate, glutamine) can reflect what is
in the medium rather than the biology of the line. RPMI, for example, has no
iron salts and no hypoxanthine, while DMEM:F12 has both.

Per line (every line in cellLineMetadata["cellLines"] gets an entry, so a
blank is always explained):
  medium       raw formulation string, e.g. "RPMI + 10% FBS + 2mM Glutamine"
  base         the base medium as written, e.g. "RPMI", "DMEM:F12", "AlphaMEM"
  family       RPMI | DMEM | DMEM-F12 | F12 | IMDM | MEM | McCoy | L-15 | other | unknown
  serum        { type, pct } when a serum is named, e.g. { "type": "FBS", "pct": 10 }
  serumFree    true when DepMap flags the medium serum-free
  supplements  the remaining " + " parts, as written
  formulationId  DepMap's MF-xxx-xxx id
  source       "screen" | "model_default" | "missing"
  why          present when source is "missing"

Sources, in order of preference:
  1. The condition the CRISPR screen was run in (ScreenSequenceMap ->
     ModelConditionID -> ModelCondition). Needs 26Q1/ScreenSequenceMap.csv
     and 26Q1/ModelCondition.csv. When they are present this script prints
     their headers and stops, so the join is written against the real
     column names rather than assumed ones.
  2. The model's onboarded (default) medium from Model26Q1.csv.
     NOTE: in this release the Model.csv columns read shifted: the MF-xxx id
     sits under "OnboardedMedia" and the formulation text under
     "FormulationID". Both are detected by content below, not by header.

Run: python3 scripts/build_culture_media.py   (idempotent)
"""
import csv
import json
import os
import re
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_CSV = os.path.join(ROOT, "Model26Q1.csv")
SCREEN_MAP = os.path.join(ROOT, "26Q1", "ScreenSequenceMap.csv")
CONDITION = os.path.join(ROOT, "26Q1", "ModelCondition.csv")
METADATA_JSON = os.path.join(ROOT, "web_data", "cellLineMetadata.json")

MF_ID = re.compile(r"^MF-\d{3}-\d{3}$")


def family_of(base):
    b = base.strip().lower().replace("’", "'")
    if not b:
        return "unknown"
    # Mixtures of two bases other than DMEM:F12 are their own thing.
    if b in ("dmem:f12", "dmem/f12", "dmem-f12", "dmem:f-12"):
        return "DMEM-F12"
    if ":" in b or "/" in b:
        return "other"
    if b.startswith("rpmi"):
        return "RPMI"
    if b.startswith("dmem"):
        return "DMEM"
    if b in ("f12", "f-12", "f-12k", "f12k", "ham's f12"):
        return "F12"
    if b.startswith("imdm"):
        return "IMDM"
    if b in ("mem", "emem", "alphamem", "alpha-mem", "alpha mem", "bme"):
        return "MEM"
    if b.startswith("mccoy"):
        return "McCoy"
    if b.startswith("l-15") or b.startswith("l15"):
        return "L-15"
    return "other"


def parse_medium(text):
    parts = [p.strip() for p in re.split(r"\s*\+\s*", text) if p.strip()]
    base = parts[0] if parts else ""
    # "DMEM:KSFM (2:1)" -> base "DMEM:KSFM (2:1)", family from the part before the ratio
    base_core = re.sub(r"\s*\(.*?\)\s*$", "", base)
    serum = None
    supplements = []
    for p in parts[1:]:
        m = re.match(r"^(\d+(?:\.\d+)?)\s*%\s*(FBS|FCS|horse serum|calf serum|serum)$", p, re.I)
        if m and serum is None:
            kind = m.group(2)
            serum = {"type": "FBS" if kind.upper() in ("FBS", "FCS") else kind.lower(),
                     "pct": float(m.group(1)) if "." in m.group(1) else int(m.group(1))}
        else:
            supplements.append(p)
    # A serum written inside the base part ("DMEM + 10%FBS" without a space) is
    # caught by the loop; a base-only string has neither.
    return base, family_of(base_core), serum, supplements


def model_default():
    out = {}
    with open(MODEL_CSV, newline="") as f:
        for r in csv.DictReader(f):
            cells = [r.get("OnboardedMedia", "") or "", r.get("FormulationID", "") or ""]
            mf = next((c for c in cells if MF_ID.match(c.strip())), "")
            text = next((c for c in cells if c.strip() and not MF_ID.match(c.strip())), "")
            out[r["ModelID"]] = {
                "text": text.strip(),
                "mf": mf.strip(),
                "serumFree": (r.get("SerumFreeMedia", "") or "").strip().lower() == "true",
            }
    return out


def main():
    if os.path.exists(SCREEN_MAP) and os.path.exists(CONDITION):
        for p in (SCREEN_MAP, CONDITION):
            with open(p, newline="") as f:
                print(os.path.basename(p), "columns:", next(csv.reader(f)))
        sys.exit("Screen-condition files found. Write the ScreenSequenceMap -> ModelCondition join "
                 "against the columns printed above, then re-run.")

    with open(METADATA_JSON) as f:
        meta = json.load(f)
    models = model_default()

    culture = {}
    for cl in meta["cellLines"]:
        m = models.get(cl)
        if not m or not m["text"]:
            culture[cl] = {
                "family": "unknown", "source": "missing",
                "why": "no medium recorded for this model in DepMap Model.csv"
                       if m else "model not found in DepMap Model.csv",
            }
            continue
        base, fam, serum, supp = parse_medium(m["text"])
        d = {"medium": m["text"], "base": base, "family": fam, "source": "model_default"}
        if serum:
            d["serum"] = serum
        if m["serumFree"]:
            d["serumFree"] = True
        if supp:
            d["supplements"] = supp
        if m["mf"]:
            d["formulationId"] = m["mf"]
        culture[cl] = d

    meta["culture"] = culture
    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print("culture:", len(culture), "lines")
    print("  source:", dict(Counter(d["source"] for d in culture.values())))
    print("  family:", dict(Counter(d["family"] for d in culture.values()).most_common()))
    lin = meta.get("lineage", {})
    blood = Counter(culture[cl]["family"] for cl in culture if lin.get(cl) in ("Lymphoid", "Myeloid"))
    print("  family in Lymphoid + Myeloid:", dict(blood.most_common()))


if __name__ == "__main__":
    main()
