"""
Data-driven check on DepMap's TP53 loss-of-function call ->
web_data/cellLineMetadata.json["p53Function"].

DepMap's call (OmicsInferredMolecularSubtypes TP53_LoF, shown in the app as
"functional loss of TP53") is genomic: a deep deletion, a loss-of-function
mutation at high allele fraction, or no expression. It can be wrong about
whether p53 still WORKS. Two readouts in the same data measure that directly:

  1. TP53 gene effect. Where p53 is active, knocking it out lets cells outgrow
     the pool (positive gene effect); where it is already lost, it does nothing.
  2. p53 target-gene score: mean z-score (across the expression cohort) of
     canonical direct p53 targets (TARGETS below).

A logistic model trained on DepMap's own calls turns the two into P(loss).
Probabilities are out-of-fold (10-fold, unshuffled, so deterministic): no line
is judged by a model that saw its own label. Lines without expression data use
a gene-effect-only model.

Per line: { pLoss, targetScore?, geneEffect, call: "loss" | "intact",
            verdict?: "active_despite_loss" | "low_without_call", readouts }
  active_despite_loss  DepMap calls loss but pLoss < 0.1
  low_without_call     DepMap calls intact but pLoss > 0.9. Weaker: p53 can be
                       held down without a TP53 lesion (MDM2 / MDM4, viral
                       proteins, low basal signalling), so this says the pathway
                       looks quiet, not that TP53 is lost.

Only TP53 is checked. For the other genes in DepMap's call (PTEN, RB1, CDKN2A,
MTAP, NF1, APC, VHL) the gene's own knockout effect does not separate lost from
intact lines well enough (it flags 786-O, a textbook VHL-null line, as intact),
and there is no equivalent target signature, so no verdict is attempted.

Run: python3 scripts/build_p53_function.py   (idempotent; needs scikit-learn)
"""
import json
import os
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, cross_val_predict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from matrix_io import read_float  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, "web_data")
TARGETS = ["CDKN1A", "MDM2", "EDA2R", "ZMAT3", "RPS27L", "SESN1", "SESN2", "TNFRSF10B", "FDXR",
           "AEN", "DDB2", "RRM2B", "SPATA18", "TP53I3", "BAX", "PHLDA3", "TRIAP1", "GDF15"]


def main():
    m = json.load(open(os.path.join(WEB, "metadata.json")))
    em = json.load(open(os.path.join(WEB, "expression_metadata.json")))
    ge = read_float(os.path.join(WEB, "geneEffects.bin.gz"), m)
    ex = read_float(os.path.join(WEB, "expression.bin.gz"), em)
    loss_set = set(json.load(open(os.path.join(WEB, "functional_loss.json")))["geneData"]["TP53"]["mutations"])
    meta_path = os.path.join(WEB, "cellLineMetadata.json")
    meta = json.load(open(meta_path))

    cls = m["cellLines"]
    eidx = {c: i for i, c in enumerate(em["cellLines"])}
    targets = [g for g in TARGETS if g in em["genes"]]
    z = []
    for g in targets:
        v = ex[em["genes"].index(g)]
        z.append((v - np.nanmean(v)) / np.nanstd(v))
    score_all = np.nanmean(np.array(z), axis=0)
    score = np.array([score_all[eidx[c]] if c in eidx else np.nan for c in cls])
    tp53 = ge[m["genes"].index("TP53")]
    y = np.array([c in loss_set for c in cls]).astype(int)

    cv = KFold(n_splits=10, shuffle=False)
    p = np.full(len(cls), np.nan)
    both = ~np.isnan(tp53) & ~np.isnan(score)
    p[both] = cross_val_predict(LogisticRegression(), np.column_stack([tp53, score])[both], y[both],
                                cv=cv, method="predict_proba")[:, 1]
    only_ge = ~np.isnan(tp53) & np.isnan(score)
    ge_ok = ~np.isnan(tp53)
    if only_ge.any():
        # Fit the gene-effect-only model on every line with a gene effect, then
        # score the ones that lack expression (out-of-fold for them as well).
        p_ge = cross_val_predict(LogisticRegression(), tp53[ge_ok].reshape(-1, 1), y[ge_ok],
                                 cv=cv, method="predict_proba")[:, 1]
        full = np.full(len(cls), np.nan)
        full[ge_ok] = p_ge
        p[only_ge] = full[only_ge]

    out = {}
    for i, cl in enumerate(cls):
        if np.isnan(p[i]):
            continue
        d = {"pLoss": round(float(p[i]), 3), "geneEffect": round(float(tp53[i]), 3),
             "call": "loss" if y[i] else "intact",
             "readouts": "gene effect + p53 targets" if both[i] else "gene effect only"}
        if not np.isnan(score[i]):
            d["targetScore"] = round(float(score[i]), 2)
        if y[i] and p[i] < 0.1:
            d["verdict"] = "active_despite_loss"
        elif not y[i] and p[i] > 0.9:
            d["verdict"] = "low_without_call"
        out[cl] = d
    meta["p53Function"] = out
    json.dump(meta, open(meta_path, "w"))

    name = meta.get("cellLineName", {})
    act = sorted([c for c, d in out.items() if d.get("verdict") == "active_despite_loss"], key=lambda c: out[c]["pLoss"])
    low = [c for c, d in out.items() if d.get("verdict") == "low_without_call"]
    print(f"p53Function: {len(out)} lines ({int(both.sum())} with both readouts, {int(only_ge.sum())} gene effect only)")
    print(f"  active despite loss call ({len(act)}):", ", ".join(f"{name.get(c, c)} P={out[c]['pLoss']}" for c in act))
    print(f"  low without a call: {len(low)}")


if __name__ == "__main__":
    main()
