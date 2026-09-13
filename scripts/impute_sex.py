#!/usr/bin/env python3
"""
Impute sex from Y-chromosome markers + XIST expression.
NOTE: XIST is a lncRNA, so it is absent from the protein-coding
expression file the app ships. This reads DepMap's AllGenes
expression export; keep it in step with the release when possible
(shared values move by <0.001 log-TPM between releases, so an
older AllGenes file still classifies correctly).

Merges DepMap Sex annotation with expression-based inference into
web_data/cellLineMetadata.json.

Input:
  /Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanAllGenes.csv
  scripts/Model_25Q3.csv

Output: web_data/cellLineMetadata.json gets two independent fields keyed by ACH-id:
  - sex:             DepMap annotation ('Male' | 'Female' | 'Unknown')
  - sexByExpression: expression-based call ('male' | 'female' | 'unknown') — computed
                     for every cell line, independent of annotation.
  - sexChromosomes:  per line with expression data only (absent = never measured):
                       y      mean log-TPM of the six Y-linked markers
                       xist   XIST log-TPM
                       xcn    median relative copy number over non-PAR chrX genes
                              (from web_data/cn.bin.gz, 1.0 = the line's modal baseline;
                              one X in a diploid line reads ~0.5), omitted when no CN
                       status one of:
                         y_present      Y-linked genes expressed
                         y_loss         annotated male, Y-linked genes silent: FUNCTIONAL loss
                                        of Y (an expression call, not a DNA one)
                         xist_present   no Y, XIST expressed (inactive X present)
                         xi_lost        annotated female, XIST off, chrX CN < 0.75: the
                                        inactive X was lost, one X left
                         xist_silenced  annotated female, XIST off, chrX CN >= 0.75: two X
                                        copies but no XIST (Xi erosion, or Xa duplicated)
                         both_low       no Y, no XIST, and annotation or CN cannot split it

Classifier (two-rule, independent thresholds):
  - Y_mean > Y_THR           -> 'male'    (Y presence is unambiguous; takes precedence)
  - XIST  > X_THR            -> 'female'  (only if Y is below threshold)
  - otherwise                -> 'unknown' (Y-loss males AND XIST-silenced females live here)

Thresholds Y_THR=1.0, X_THR=1.0 chosen to separate the two populations robustly in
both directions. Note: XIST silencing is common in breast and other cancers, so a
large minority of annotated females land in 'unknown' — this is biology, not a bug.

Script is idempotent.
"""

import csv
import json
import os
import re
import sys
from collections import Counter

import numpy as np

EXPRESSION_CSV = "/Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanAllGenes.csv"
HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_CSV = os.path.join(HERE, "..", "Model26Q1.csv")
METADATA_JSON = os.path.join(HERE, "..", "web_data", "cellLineMetadata.json")

Y_MARKERS = ["RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y"]
XIST = "XIST"

# Classifier thresholds (log-TPM units). See tuning in commit history.
Y_THR = 1.0
X_THR = 1.0
# chrX relative CN below this = one X copy against a diploid baseline (0.5 is
# the ideal; 0.75 is the midpoint to two copies).
XCN_THR = 0.75
# GRCh38 PAR1 ends at 2.78 Mb, PAR2 starts at 155.7 Mb; both are present on Y
# and are excluded from the chrX measurement.
PAR1_END = 2_800_000
PAR2_START = 155_700_000
CN_BIN = os.path.join(HERE, "..", "web_data", "cn.bin.gz")
CN_META = os.path.join(HERE, "..", "web_data", "cn_metadata.json")
GENE_LOC = os.path.join(HERE, "..", "web_data", "gene_locations.json")


def chrx_median_cn():
    """Per-line median relative CN over non-PAR chrX genes, from the shipped matrix."""
    import gzip
    cm = json.load(open(CN_META))
    loc = json.load(open(GENE_LOC))["genes"]
    cn = np.frombuffer(gzip.open(CN_BIN).read(), dtype=np.int16)
    cn = cn.reshape(cm["nGenes"], cm["nCellLines"]).astype(float)
    cn[cn == cm["naValue"]] = np.nan
    cn /= cm["scaleFactor"]
    rows = [i for i, g in enumerate(cm["genes"])
            if loc.get(g, {}).get("chr") == "X"
            and PAR1_END < (loc[g].get("start") or 0) < PAR2_START]
    sub = cn[rows, :]
    out = {}
    with np.errstate(all="ignore"):
        med = np.nanmedian(sub, axis=0)
    for ci, cl in enumerate(cm["cellLines"]):
        if not np.isnan(med[ci]):
            out[cl] = float(med[ci])
    print(f"  chrX CN: {len(rows)} non-PAR genes, {len(out)} lines")
    return out


def chromosome_status(y, x, xcn, annotation):
    if np.isnan(y) or np.isnan(x):
        return None
    if y > Y_THR:
        return "y_present"
    if x > X_THR:
        return "xist_present"
    if annotation == "Male":
        return "y_loss"
    if annotation == "Female" and xcn is not None:
        return "xi_lost" if xcn < XCN_THR else "xist_silenced"
    return "both_low"


def parse_gene_name(col_header):
    m = re.match(r"^(.+?)\s*\(\d+\)$", col_header)
    return m.group(1).strip() if m else col_header.strip()


def main():
    print(f"Reading header from {EXPRESSION_CSV}")
    with open(EXPRESSION_CSV) as f:
        header = next(csv.reader(f))

    model_id_col = header.index("ModelID")
    is_default_col = header.index("IsDefaultEntryForModel")

    target_genes = Y_MARKERS + [XIST]
    gene_col_idx = {}
    for i, h in enumerate(header[6:], start=6):
        name = parse_gene_name(h)
        if name in target_genes:
            gene_col_idx[name] = i
    missing = [g for g in target_genes if g not in gene_col_idx]
    if missing:
        sys.exit(f"ERROR: markers missing from expression file: {missing}")
    print(f"  Located all {len(target_genes)} markers")

    print("Streaming rows (filtering IsDefaultEntryForModel=Yes)")
    cell_vals = {}
    with open(EXPRESSION_CSV) as f:
        reader = csv.reader(f)
        next(reader)
        for row_num, row in enumerate(reader):
            if row_num % 500 == 0 and row_num > 0:
                print(f"  Row {row_num}")
            if row[is_default_col].strip() != "Yes":
                continue
            cl = row[model_id_col].strip()
            vals = {}
            for g, idx in gene_col_idx.items():
                raw = row[idx].strip() if idx < len(row) else ""
                try:
                    vals[g] = float(raw) if raw else np.nan
                except ValueError:
                    vals[g] = np.nan
            cell_vals[cl] = vals
    print(f"  Loaded markers for {len(cell_vals)} cell lines")

    # Per cell line: mean(Y) and XIST
    scores = {}
    for cl, v in cell_vals.items():
        y_vals = [v[g] for g in Y_MARKERS if not np.isnan(v.get(g, np.nan))]
        x = v.get(XIST, np.nan)
        y_mean = float(np.mean(y_vals)) if y_vals else np.nan
        scores[cl] = (y_mean, float(x) if not np.isnan(x) else np.nan)

    print(f"Reading Sex labels from {MODEL_CSV}")
    depmap_sex = {}
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            s = row["Sex"].strip()
            depmap_sex[row["ModelID"]] = s if s else "Unknown"

    # Classifier: Y_THR takes precedence over X_THR so XIST-positive males
    # (e.g. Klinefelter) still get called male.
    def classify(y, x):
        if np.isnan(y) or np.isnan(x):
            return "unknown"
        if y > Y_THR:
            return "male"
        if x > X_THR:
            return "female"
        return "unknown"

    print(f"  Thresholds: Y_mean > {Y_THR} -> male, XIST > {X_THR} (Y low) -> female")

    # Validation on known-label cell lines
    male_y, male_x, female_y, female_x = [], [], [], []
    for cl, (y, x) in scores.items():
        if np.isnan(y) or np.isnan(x):
            continue
        label = depmap_sex.get(cl, "Unknown")
        if label == "Male":
            male_y.append(y); male_x.append(x)
        elif label == "Female":
            female_y.append(y); female_x.append(x)

    m_calls = Counter(classify(y, x) for y, x in zip(male_y, male_x))
    f_calls = Counter(classify(y, x) for y, x in zip(female_y, female_x))
    print("  Validation:")
    print(
        f"    Of {len(male_y)} known males:  "
        f"{m_calls.get('male', 0)} male, "
        f"{m_calls.get('female', 0)} female ({100 * m_calls.get('female', 0) / len(male_y):.2f}% FP), "
        f"{m_calls.get('unknown', 0)} unknown"
    )
    print(
        f"    Of {len(female_y)} known females: "
        f"{f_calls.get('male', 0)} male ({100 * f_calls.get('male', 0) / len(female_y):.2f}% FP), "
        f"{f_calls.get('female', 0)} female, "
        f"{f_calls.get('unknown', 0)} unknown"
    )

    by_expression = {cl: classify(y, x) for cl, (y, x) in scores.items()}

    print(f"Updating {METADATA_JSON}")
    with open(METADATA_JSON) as f:
        meta = json.load(f)

    sex_map = {}
    exp_map = {}
    for cl in meta["cellLines"]:
        sex_map[cl] = depmap_sex.get(cl, "Unknown")
        exp_map[cl] = by_expression.get(cl, "unknown")

    meta["sex"] = sex_map
    meta["sexByExpression"] = exp_map

    xcn_map = chrx_median_cn()
    chrom = {}
    for cl in meta["cellLines"]:
        if cl not in scores:
            continue
        y, x = scores[cl]
        if np.isnan(y) or np.isnan(x):
            continue
        xcn = xcn_map.get(cl)
        d = {"y": round(y, 3), "xist": round(x, 3),
             "status": chromosome_status(y, x, xcn, sex_map[cl])}
        if xcn is not None:
            d["xcn"] = round(xcn, 3)
        chrom[cl] = d
    meta["sexChromosomes"] = chrom
    # Remove any previous field name
    meta.pop("sexImputed", None)

    with open(METADATA_JSON, "w") as f:
        json.dump(meta, f)

    print("\nFinal distribution in cellLineMetadata.json:")
    print("  sex (annotation):    ", dict(Counter(sex_map.values())))
    print("  sexByExpression:     ", dict(Counter(exp_map.values())))
    print("  sexChromosomes:      ", len(chrom), "lines,", dict(Counter(d["status"] for d in chrom.values())))

    # Agreement crosstab
    ct = Counter((sex_map[cl], exp_map[cl]) for cl in meta["cellLines"])
    print("\nAgreement crosstab (annotation × expression):")
    for (a, e), n in sorted(ct.items()):
        print(f"    {a:8s} × {e:8s}: {n}")


if __name__ == "__main__":
    main()
