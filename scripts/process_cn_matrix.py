#!/usr/bin/env python3
"""CN matrix from DepMap OmicsCNGene.csv (24Q4 release).

Source file:
  OmicsCNGene.csv (24Q4, ~1.4 GB)
    — DepMap's merged gene-level CN file. "Inferred from WGS or WES
      depending on the availability of the data type" — already hybrid
      by construction. ~1929 cell lines including lines never WGS'd
      (Jurkat, K562, ...). One row per cell line, ACH-... ID in
      column 0, gene columns to the right.

We use the 24Q4 OmicsCNGene values uniformly for every line and don't
tag per-line WGS-vs-WES provenance separately — DepMap's own per-line
choice is baked into the values, and the OmicsCNGeneWGS.csv 25Q3 file
that previously supplied the tag was a step ahead of the value vintage
(some lines tagged WGS by 25Q3 still had WES-derived 24Q4 values in our
binary), which created a small false-precision mismatch. Treating
everything as "DepMap OmicsCNGene 24Q4" is honest and uniform.

The gene list is intersected with the greenlisted CRISPR libraries
(Brunello, Gattinara, GeCKO v2, Jacquere, VBC, MinLibCas9, TKOv3,
Yusa v1) expanded through the human synonym table, so genes not
targetable by any major library are dropped — they have no actionable
CRISPR use case here.

Output:
  web_data/cn.bin.gz        — gzipped int16 matrix, nGenes × nCellLines
  web_data/cn_metadata.json — gene + cell-line lists, scaleFactor, naValue.

Scale factor: 3000 → range ±10.92 in CN units (1.0 = diploid).
"""
import csv, gzip, json, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HYBRID_PATH = os.path.join(ROOT, 'OmicsCNGene.csv')
OUT_BIN = os.path.join(ROOT, 'web_data', 'cn.bin.gz')
OUT_META = os.path.join(ROOT, 'web_data', 'cn_metadata.json')
GREENLISTED_LIB_DIR = '/Users/fredrikwermeling/Documents/greenlistedv2/libraries'
SCALE_FACTOR = 3000
NA_VALUE = -32768

META_COLS = {'SequencingID', 'ModelID', 'IsDefaultEntryForModel',
             'ModelConditionID', 'IsDefaultEntryForMC'}
SYMBOL_PAT = re.compile(r'^(.+?)\s*\((\d+)\)$')

def parse_symbol(col):
    m = SYMBOL_PAT.match(col)
    return m.group(1).strip().upper() if m else None

def load_crispr_library_symbols(lib_dir):
    """Return the set of uppercase gene symbols reachable through any human
    CRISPR library in greenlisted/libraries (Brunello, Gattinara, GeCKO v2,
    Jacquere, VBC, MinLibCas9, TKOv3, Yusa v1), expanded through the human
    synonym table so that an alias hit counts as a hit.

    Rationale: the only reason to include a gene in this CN matrix is for
    CRISPR target prioritisation. ~11k entries in OmicsCNGene (most
    pseudogenes, RNU/snoRNAs, IG / TR gene segments, and obscure lncRNAs)
    have no sgRNA in any major library, so they can't be screened. Dropping
    them shrinks the binary without removing any actionable target."""
    if not os.path.isdir(lib_dir):
        print(f'  greenlisted lib dir not found at {lib_dir} — skipping CRISPR-library cross-reference', file=sys.stderr)
        return None
    syms = set()
    # TSV libraries: (filename, 0-based column index of gene symbol)
    tsv_libs = [
        ('Brunello (human).txt', 1),
        ('Gattinara (human).txt', 1),
        ('GeCKO v2 (human) A+B.txt', 0),
        ('Jacquere (human).txt', 1),
        ('VBC (human).txt', 0),
    ]
    for fn, col in tsv_libs:
        p = os.path.join(lib_dir, fn)
        if not os.path.exists(p): continue
        with open(p, newline='') as f:
            rdr = csv.reader(f, delimiter='\t')
            next(rdr, None)
            for row in rdr:
                if len(row) > col and row[col].strip():
                    syms.add(row[col].strip().upper())
    # xlsx libraries
    try:
        import openpyxl
        xlsx_libs = [
            ('minlibcas9_raw.xlsx', 'Sheet1', 3),       # Approved_Symbol
            ('tkov3_raw.xlsx', 'TKOv3-Human-Library', 0),
            ('yusa_human_v1_raw.xlsx', 'Human v1', 1),
        ]
        for name, sheet, col in xlsx_libs:
            p = os.path.join(lib_dir, name)
            if not os.path.exists(p): continue
            wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
            ws = wb[sheet]
            for r in ws.iter_rows(min_row=2, values_only=True):
                if col < len(r) and r[col] is not None and str(r[col]).strip():
                    syms.add(str(r[col]).strip().upper())
    except ImportError:
        print('  openpyxl not installed — TKOv3 / Yusa / MinLibCas9 not crossreferenced', file=sys.stderr)
    n_raw = len(syms)
    # Expand through the bidirectional alias table so e.g. CCBL1 ↔ KYAT1.
    syn_path = os.path.join(lib_dir, 'human synonym.txt')
    if os.path.exists(syn_path):
        pairs = []
        with open(syn_path, newline='') as f:
            rdr = csv.reader(f, delimiter='\t')
            next(rdr, None)
            for row in rdr:
                if len(row) >= 2:
                    a, b = row[0].strip().upper(), row[1].strip().upper()
                    if a and b: pairs.append((a, b))
        for _ in range(5):
            grew = False
            for a, b in pairs:
                if a in syms and b not in syms: syms.add(b); grew = True
                if b in syms and a not in syms: syms.add(a); grew = True
            if not grew: break
    print(f'  CRISPR-library symbols: {n_raw} raw → {len(syms)} after synonym expansion')
    return syms

def main():
    if not os.path.exists(HYBRID_PATH):
        print(f'OmicsCNGene.csv not found at {HYBRID_PATH}', file=sys.stderr); sys.exit(1)

    # Drop obvious noise categories. IMPORTANT: pseudogene patterns
    # must be anchored to specific known prefixes (RPL, RPS, HBA, HBB,
    # HNRNP, EEF1, TUBB, etc.). A broad `P\d+$` regex catches TP53,
    # TP63, USP1..USP54, BMP1..BMP15, AKAP1..AKAP14, and many other
    # real gene symbols — never use it. Kept: protein-coding, named
    # lncRNAs, miRNA host genes, olfactory receptors (real genes),
    # immunoglobulin / TCR gene segments.
    NOISE_PATTERNS = [
        re.compile(r'^LOC\d+$'),                # LOC102345 — uncharacterized loci
        re.compile(r'^RNU\d'),                  # RNU6, RNU7 — small nuclear RNAs
        re.compile(r'^SNORA\d|^SNORD\d'),       # small nucleolar RNAs
        re.compile(r'^SCARNA\d'),               # small Cajal-body RNAs
        re.compile(r'-AS\d+$'),                 # antisense lncRNAs (PDE4B-AS1)
        re.compile(r'-OT\d+$'),                 # overlapping transcripts
        re.compile(r'-DT$'),                    # divergent transcripts
        # Narrow pseudogene patterns — match a specific known prefix
        # so we don't kill TP53 / USP1 / BMP1 / etc.
        re.compile(r'^RPL\d+[A-Z]?P\d+$'),      # ribosomal large pseudogenes (RPL32P26)
        re.compile(r'^RPS[\dA-Z]+P\d+$'),       # ribosomal small pseudogenes (RPS12P2)
        re.compile(r'^HNRNP[A-Z]?\d?P\d+$'),    # hnRNP pseudogenes
        re.compile(r'^(HBA|HBB|HBE|HBG)\d?P\d+$'),  # haemoglobin pseudogenes
        re.compile(r'^EEF1[A-Z]?\d?P\d+$'),     # EEF1 pseudogenes
        re.compile(r'^TUBB[A-Z]?\d?P\d+$'),     # tubulin-beta pseudogenes
        re.compile(r'^GAPDHP\d+$'),             # GAPDH pseudogenes
        re.compile(r'^Y_RNA'),                  # Y RNAs
        re.compile(r'^7SK$'),                   # 7SK ncRNA
    ]
    def is_noise(sym):
        for p in NOISE_PATTERNS:
            if p.search(sym): return True
        return False

    print('cross-referencing against greenlisted CRISPR libraries ...')
    crispr_syms = load_crispr_library_symbols(GREENLISTED_LIB_DIR)

    print(f'reading {HYBRID_PATH} ...')
    # Parse header to build gene_col_idx, deduplicating symbols.
    with open(HYBRID_PATH, newline='') as f:
        rdr = csv.reader(f)
        header = next(rdr)
        gene_syms = []
        gene_col_idx = []
        seen = set()
        n_filtered = 0
        n_no_sgrna = 0
        for j, col in enumerate(header):
            if j == 0: continue  # row-index column
            sym = parse_symbol(col)
            if not sym or sym in seen: continue
            if is_noise(sym):
                n_filtered += 1
                continue
            # Drop anything no CRISPR library can target — there's no
            # actionable point in showing CN for a gene with no sgRNAs.
            if crispr_syms is not None and sym not in crispr_syms:
                n_no_sgrna += 1
                continue
            seen.add(sym)
            gene_syms.append(sym)
            gene_col_idx.append(j)
        print(f'  total cols: {len(header)}, kept genes: {len(gene_syms)} (dropped {n_filtered} noise, {n_no_sgrna} no-sgRNA)')

        # Second pass: stream rows. One row per cell line (ACH-...).
        # Skip duplicates if any (use first occurrence).
        cl_ids = []
        rows_kept = []
        seen_cl = set()
        n_seen = 0
        for row in rdr:
            n_seen += 1
            if not row: continue
            cl = row[0].strip()
            if not cl or cl in seen_cl: continue
            seen_cl.add(cl)
            cl_ids.append(cl)
            rows_kept.append(row)
        print(f'  scanned {n_seen} rows, kept {len(cl_ids)} unique cell lines')

    n_genes = len(gene_syms)
    n_cl = len(cl_ids)
    print(f'  building matrix {n_genes} genes × {n_cl} cell lines = {n_genes * n_cl:,} entries')

    data = np.full((n_genes, n_cl), np.nan, dtype=np.float32)
    for ci, row in enumerate(rows_kept):
        for gi, col_i in enumerate(gene_col_idx):
            if col_i >= len(row): continue
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

    flat = q.flatten()
    print(f'  writing {OUT_BIN} ...')
    with gzip.open(OUT_BIN, 'wb') as f:
        f.write(flat.tobytes())
    size_mb = os.path.getsize(OUT_BIN) / (1024 * 1024)
    print(f'  binary size: {size_mb:.1f} MB')

    meta = {
        '_doc': 'CN matrix from DepMap OmicsCNGene.csv (24Q4 release) — DepMap\'s gene-level CN file that already chooses WGS or WES per line based on availability. Gene list is intersected with the greenlisted CRISPR libraries (Brunello, Gattinara, GeCKO v2, Jacquere, VBC, MinLibCas9, TKOv3, Yusa v1) expanded through the human synonym table; genes not reachable by any major library are dropped. Values are relative copy number (1.0 = the line\'s own modal baseline). Int16, gene-major (matrix[gi * nCellLines + ci]). NaN encoded as -32768. Scale factor 3000.',
        'source': 'OmicsCNGene.csv (DepMap 24Q4)',
        'genes': gene_syms,
        'cellLines': cl_ids,
        'nGenes': n_genes,
        'nCellLines': n_cl,
        'scaleFactor': SCALE_FACTOR,
        'naValue': NA_VALUE,
        'interpretation': '1.0 = the cell line\'s own modal baseline (≈ 2 copies for a diploid line, ≈ 4 for WGD). 3.0 = amplification. 5.0 = strong amplification. 0.5 = single-copy loss. 0.3 = deep deletion.'
    }
    with open(OUT_META, 'w') as f:
        json.dump(meta, f, indent=1)
    print(f'  metadata: {OUT_META} ({os.path.getsize(OUT_META)} bytes)')
    print('done')

if __name__ == '__main__':
    main()
