#!/usr/bin/env python3
"""Process DepMap OmicsCNGeneWGS.csv into a compact binary matrix.

Mirrors the geneEffects pipeline: gene-major int16 array with a
scaleFactor for fixed-point encoding and -32768 reserved for NaN.

Output:
  web_data/cn.bin.gz          — gzipped int16 matrix, nGenes × nCellLines
  web_data/cn_metadata.json   — gene + cell-line lists, scaleFactor, naValue

The cohort is whatever cell lines DepMap publishes CN for (filtered to
default-entry rows). The app's renderCellLineList uses metadata.cellLines
(GE cohort) as the master list; lines without CN data show "—" in the
sort-value column.

Scale factor: 3000 → range ±10.92 in CN units. CN values are relative
to diploid (1.0); typical range is 0 to ~6 with rare extreme amps
that may clip near 11. That's acceptable — anything ≥ 6 already
qualifies as "strong amplification".
"""
import csv, gzip, json, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_PATH = os.path.join(ROOT, 'OmicsCNGeneWGS.csv')
OUT_BIN = os.path.join(ROOT, 'web_data', 'cn.bin.gz')
OUT_META = os.path.join(ROOT, 'web_data', 'cn_metadata.json')
SCALE_FACTOR = 3000
NA_VALUE = -32768

# Metadata columns at the start of each row — skip when collecting
# gene values.
META_COLS = {'SequencingID', 'ModelID', 'IsDefaultEntryForModel',
             'ModelConditionID', 'IsDefaultEntryForMC'}
# First unnamed column is the row index.

def main():
    print(f'reading {IN_PATH} ...')
    if not os.path.exists(IN_PATH):
        print(f'  not found — need to download from DepMap', file=sys.stderr); sys.exit(1)

    with open(IN_PATH, newline='') as f:
        rdr = csv.reader(f)
        header = next(rdr)
        # Find metadata column indices and the first gene column.
        meta_idx = {}
        for i, h in enumerate(header):
            if h in META_COLS or (i == 0 and not h):
                meta_idx[h or '_idx'] = i
        first_gene_col = max(meta_idx.values()) + 1
        gene_cols = header[first_gene_col:]

        # Parse gene symbols from "SYMBOL (entrez)" or "ENSG..._PAR_Y".
        # We drop PAR_Y and unnamed columns to keep the matrix to real
        # symbol-named genes.
        symbol_pat = re.compile(r'^(.+?)\s*\((\d+)\)$')
        gene_syms = []
        gene_col_idx = []  # csv column index for each kept gene
        seen = set()
        for j, col in enumerate(gene_cols):
            m = symbol_pat.match(col)
            if not m:
                continue  # skip PAR_Y and other unparseable entries
            sym = m.group(1).strip().upper()
            if sym in seen:
                continue
            seen.add(sym)
            gene_syms.append(sym)
            gene_col_idx.append(first_gene_col + j)
        print(f'  total cols: {len(header)}, kept genes: {len(gene_syms)}')

        # Stream rows, keeping only IsDefaultEntryForModel == 'Yes'.
        is_default_col = meta_idx.get('IsDefaultEntryForModel')
        model_col = meta_idx.get('ModelID')
        cl_ids = []
        rows_kept = []
        n_seen = 0
        for row in rdr:
            n_seen += 1
            if is_default_col is not None and row[is_default_col] != 'Yes':
                continue
            cl = row[model_col]
            if not cl: continue
            cl_ids.append(cl)
            rows_kept.append(row)
        print(f'  scanned {n_seen} rows, kept {len(rows_kept)} default-entry cell lines')

    n_genes = len(gene_syms)
    n_cl = len(cl_ids)
    print(f'  building matrix {n_genes} genes × {n_cl} cell lines = {n_genes * n_cl:,} entries')

    # gene-major: matrix[gi, ci]
    data = np.full((n_genes, n_cl), np.nan, dtype=np.float32)
    for ci, row in enumerate(rows_kept):
        for gi, col_i in enumerate(gene_col_idx):
            v = row[col_i].strip()
            if not v or v == 'NA': continue
            try:
                data[gi, ci] = float(v)
            except ValueError:
                pass

    valid = ~np.isnan(data)
    print(f'  valid values: {valid.sum():,} / {data.size:,} ({100*valid.mean():.1f}%)')
    print(f'  min={np.nanmin(data):.3f} max={np.nanmax(data):.3f} median={np.nanmedian(data):.3f}')
    above_clip = np.sum(data[valid] > 32767/SCALE_FACTOR)
    print(f'  would clip {above_clip} extreme amps (≥ {32767/SCALE_FACTOR:.2f})')

    print(f'  quantising to int16, scaleFactor={SCALE_FACTOR}')
    q = np.full(data.shape, NA_VALUE, dtype=np.int16)
    q[valid] = np.clip(np.round(data[valid] * SCALE_FACTOR).astype(np.int32), -32767, 32767).astype(np.int16)

    # Gene-major flatten for compact storage.
    flat = q.flatten()
    print(f'  writing {OUT_BIN} ...')
    with gzip.open(OUT_BIN, 'wb') as f:
        f.write(flat.tobytes())
    size_mb = os.path.getsize(OUT_BIN) / (1024 * 1024)
    print(f'  binary size: {size_mb:.1f} MB')

    meta = {
        '_doc': 'DepMap OmicsCNGeneWGS — relative copy-number values (1.0 = diploid). Int16, gene-major (matrix[gi * nCellLines + ci]). NaN encoded as -32768. Scale factor 3000 → divide raw value by 3000 to get CN.',
        'source': {'file': 'OmicsCNGeneWGS.csv', 'release': 'DepMap'},
        'genes': gene_syms,
        'cellLines': cl_ids,
        'nGenes': n_genes,
        'nCellLines': n_cl,
        'scaleFactor': SCALE_FACTOR,
        'naValue': NA_VALUE,
        'interpretation': '1.0 = diploid (2 copies). 3.0 = amplification (≈ 6 copies). 5.0 = strong amplification. 0.5 = single-copy loss. 0.3 = deep deletion. Range theoretically 0 to ~20 but typically clipped at 6.'
    }
    with open(OUT_META, 'w') as f:
        json.dump(meta, f, indent=1)
    print(f'  metadata: {OUT_META} ({os.path.getsize(OUT_META)} bytes)')
    print('done')

if __name__ == '__main__':
    main()
