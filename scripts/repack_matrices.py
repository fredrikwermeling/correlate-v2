"""
Repack web_data's int16 matrices for shipping: requantise to 3 decimals
(scaleFactor 1000, the precision every export already rounds to) and store the
low and high bytes as separate planes. Gene effect 37 -> 27 MB, expression
52 -> 45 MB, copy number 62 -> 48 MB, which is most of the first-load wait on
a phone. Idempotent: a file whose metadata already says byteSplit is skipped.

Run after any script that rewrites a matrix (build_26q1_core.py,
process_cn_matrix.py):  python3 scripts/repack_matrices.py
"""
import json
import os
import numpy as np
from matrix_io import read_float, write_bytesplit

WEB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "web_data")
TARGET_SCALE = 1000
PAIRS = [("geneEffects.bin.gz", "metadata.json"),
         ("expression.bin.gz", "expression_metadata.json"),
         ("cn.bin.gz", "cn_metadata.json")]

for bin_name, meta_name in PAIRS:
    bp, mp = os.path.join(WEB, bin_name), os.path.join(WEB, meta_name)
    meta = json.load(open(mp))
    if meta.get("byteSplit"):
        print(f"{bin_name}: already packed, skipped")
        continue
    before = os.path.getsize(bp)
    v = read_float(bp, meta)
    na = np.isnan(v)
    q = np.clip(np.round(np.nan_to_num(v) * TARGET_SCALE), -32767, 32767).astype(np.int16)
    q[na] = meta["naValue"]
    write_bytesplit(bp, q)
    meta["scaleFactor"] = TARGET_SCALE
    meta["byteSplit"] = True
    json.dump(meta, open(mp, "w"))
    print(f"{bin_name}: {before/1e6:.1f} -> {os.path.getsize(bp)/1e6:.1f} MB")
