"""
Read / repack the gzipped int16 matrices in web_data/ (geneEffects, expression, cn).

Two on-disk layouts, told apart by the metadata's "byteSplit" flag:
  plain      little-endian int16 pairs, gene-major (what the build scripts write)
  byteSplit  every low byte, then every high byte (what ships: gzip packs it
             about a quarter smaller); written by repack_matrices.py

Always read through read_int16 so a script works on either layout.
"""
import gzip
import numpy as np


def read_int16(bin_path, meta):
    raw = gzip.open(bin_path, "rb").read()
    if meta.get("byteSplit"):
        n = len(raw) // 2
        lo = np.frombuffer(raw[:n], dtype=np.uint8).astype(np.uint16)
        hi = np.frombuffer(raw[n:], dtype=np.uint8).astype(np.uint16)
        arr = ((hi << 8) | lo).view(np.int16)
    else:
        arr = np.frombuffer(raw, dtype="<i2")
    return arr[: meta["nGenes"] * meta["nCellLines"]].reshape(meta["nGenes"], meta["nCellLines"])


def read_float(bin_path, meta):
    a = read_int16(bin_path, meta)
    out = a.astype(np.float32)
    out[a == meta["naValue"]] = np.nan
    out /= meta["scaleFactor"]
    return out


def write_bytesplit(bin_path, int16_matrix):
    b = np.ascontiguousarray(int16_matrix.astype("<i2")).reshape(-1).view(np.uint8).reshape(-1, 2)
    with gzip.open(bin_path, "wb", compresslevel=9) as f:
        f.write(b[:, 0].tobytes())
        f.write(b[:, 1].tobytes())
