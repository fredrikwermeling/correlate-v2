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
  why          present when source is not "screen", saying why
  screenMedia  every distinct medium among the line's screens, when more than one
  modelDefault the model's default medium, when the screen used a different one

Sources, in order of preference:
  1. The condition the CRISPR screen was run in: ScreenSequenceMap rows that
     feed the combined gene effect (PassesQC True, ExcludeFromCRISPRCombined
     False, ScreenType 2DS) -> ModelConditionID -> ModelCondition. Needs
     26Q1/ScreenSequenceMap.csv and 26Q1/ModelCondition.csv. ModelCondition
     has the same shift as Model.csv (MF id under GrowthMedia, text under
     FormulationID). Sanger screen conditions usually record no medium; a
     line whose screens record none falls back to 2. A line screened in two
     different media keeps both in screenMedia, the Broad one first.
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


def screen_media():
    """Model -> list of (text, mf, serumFree) for the media its combined screens ran in,
    Broad conditions first. Empty list when the screens record no medium."""
    if not (os.path.exists(SCREEN_MAP) and os.path.exists(CONDITION)):
        return None
    cond = {}
    with open(CONDITION, newline="") as f:
        for r in csv.DictReader(f):
            cells = [r.get("GrowthMedia", "") or "", r.get("FormulationID", "") or ""]
            cond[r["ModelConditionID"]] = {
                "text": next((c.strip() for c in cells if c.strip() and not MF_ID.match(c.strip())), ""),
                "mf": next((c.strip() for c in cells if MF_ID.match(c.strip())), ""),
                "serumFree": (r.get("SerumFreeMedia", "") or "").strip().lower() == "true",
                "source": r.get("DataSource", ""),
            }
    per = {}
    with open(SCREEN_MAP, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("PassesQC") != "True" or r.get("ExcludeFromCRISPRCombined") != "False" or r.get("ScreenType") != "2DS":
                continue
            per.setdefault(r["ModelID"], set()).add(r["ModelConditionID"])
    out = {}
    for model, cids in per.items():
        seen, media = set(), []
        for c in sorted((cond[c] for c in cids if c in cond), key=lambda c: c["source"] != "BROAD"):
            if c["text"] and c["text"] not in seen:
                seen.add(c["text"])
                media.append(c)
        out[model] = media
    return out


def entry_from(text, mf, serum_free, source):
    base, fam, serum, supp = parse_medium(text)
    d = {"medium": text, "base": base, "family": fam, "source": source}
    if serum:
        d["serum"] = serum
    if serum_free:
        d["serumFree"] = True
    if supp:
        d["supplements"] = supp
    if mf:
        d["formulationId"] = mf
    return d


def main():

    with open(METADATA_JSON) as f:
        meta = json.load(f)
    models = model_default()
    screens = screen_media()
    if screens is None:
        print("Screen-condition files not found in 26Q1/, using the model default medium only.")

    culture = {}
    for cl in meta["cellLines"]:
        m = models.get(cl)
        sm = (screens or {}).get(cl, [])
        if sm:
            first = sm[0]
            d = entry_from(first["text"], first["mf"], first["serumFree"], "screen")
            if len(sm) > 1:
                d["screenMedia"] = [c["text"] for c in sm]
            if m and m["text"] and m["text"] != first["text"]:
                d["modelDefault"] = m["text"]
            culture[cl] = d
            continue
        screen_why = ("the conditions of this line's screens record no medium (typical of Sanger screens)"
                      if screens is not None else "screen-condition tables not available")
        if m and m["text"]:
            d = entry_from(m["text"], m["mf"], m["serumFree"], "model_default")
            d["why"] = screen_why + ", so this is the medium DepMap lists for the model"
            culture[cl] = d
        else:
            culture[cl] = {
                "family": "unknown", "source": "missing",
                "why": screen_why + ", and " + ("no medium is recorded for this model in DepMap Model.csv"
                                                if m else "the model is not in DepMap Model.csv"),
            }

    meta["culture"] = culture
    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print("culture:", len(culture), "lines")
    print("  source:", dict(Counter(d["source"] for d in culture.values())))
    print("  screened in more than one medium:", sum(1 for d in culture.values() if d.get("screenMedia")))
    print("  screen medium differs from model default:", sum(1 for d in culture.values() if d.get("modelDefault")))
    print("  family:", dict(Counter(d["family"] for d in culture.values()).most_common()))
    lin = meta.get("lineage", {})
    blood = Counter(culture[cl]["family"] for cl in culture if lin.get(cl) in ("Lymphoid", "Myeloid"))
    print("  family in Lymphoid + Myeloid:", dict(blood.most_common()))


if __name__ == "__main__":
    main()
