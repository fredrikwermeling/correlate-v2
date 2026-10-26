#!/usr/bin/env python3
"""Count viral reads per cell line in the UNMAPPED reads of each CCLE RNA BAM.

The CCLE hg19 BAMs are coordinate-sorted, so reads that matched nothing in
the human genome sit together at the end of the file. A binary search with
HTTP range reads finds where that tail begins, and only the tail is
downloaded. No index, no pysam, nothing written to disk.

Each unmapped read is cut into 25-mers (every 4th position) and looked up in
a table of every 25-mer of every reference genome (both strands). A read
counts for a virus when 2 or more of its 25-mers belong to that virus alone.
25-mers shared by two references (MLV / XMRV), or of low complexity, are
dropped. Breadth is the share of 50 equal bins of the genome that any read
touched: a real infection covers the genome, a vector fragment does not.

Not scanned: unmapped reads whose mate mapped (those are filed next to the
mate, inside the mapped part of the file).

Usage:
  python3 virus_scan.py lines.json out.json [--workers N] [--sample-mb N] [--names A,B]
  --sample-mb N   read only the last N MB of each file (a quick screen)
"""
import argparse
import glob
import json
import os
import struct
import sys
import time
import urllib.request
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K, STEP, MIN_HITS, NBINS, POSBIN = 25, 4, 2, 50, 50
CHUNK = 16 * 1024 * 1024
MAGIC = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC'

# --------------------------------------------------------------- references
_ASCII = np.full(256, 255, np.uint8)
for i, c in enumerate(b'ACGT'):
    _ASCII[c] = i
    _ASCII[c + 32] = i
# BAM 4-bit bases: A=1 C=2 G=4 T=8, anything else is not a base we use.
_NIB = np.full(16, 255, np.uint8)
_NIB[1], _NIB[2], _NIB[4], _NIB[8] = 0, 1, 2, 3


def _kmers(codes):
    """25-mer codes of a 1-D base array, with a validity mask."""
    n = len(codes) - K + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    w = np.lib.stride_tricks.sliding_window_view(codes, K)
    ok = (w != 255).all(1)
    out = np.zeros(n, np.uint64)
    for t in range(K):
        out = (out << np.uint64(2)) | (w[:, t] & 3).astype(np.uint64)
    # Low complexity: one base 80% or more of the k-mer, or a 2- or 3-base repeat.
    counts = np.stack([(w == b).sum(1) for b in range(4)], 1)
    ok &= counts.max(1) < 20
    ok &= ~(w[:, :-2] == w[:, 2:]).all(1)
    ok &= ~(w[:, :-3] == w[:, 3:]).all(1)
    return out, ok


def load_index():
    refs = json.load(open(os.path.join(HERE, 'references.json')))['references']
    keys, lens, allc, allv, allp = [], [], [], [], []
    for vi, r in enumerate(refs):
        lines = open(os.path.join(HERE, 'refs', f"{r['key']}.fa")).read().split('\n')
        seq = ''.join(lines[1:]).encode()
        fw = _ASCII[np.frombuffer(seq, np.uint8)]
        rc = np.where(fw == 255, 255, 3 - fw)[::-1].astype(np.uint8)
        L = len(fw)
        for strand, arr in ((0, fw), (1, rc)):
            c, ok = _kmers(arr)
            pos = np.arange(len(c))
            if strand:
                pos = L - K - pos
            allc.append(c[ok]); allv.append(np.full(ok.sum(), vi, np.int16)); allp.append(pos[ok].astype(np.int32))
        keys.append(r['key']); lens.append(L)
    codes = np.concatenate(allc); vid = np.concatenate(allv); pos = np.concatenate(allp)
    order = np.argsort(codes, kind='stable')
    codes, vid, pos = codes[order], vid[order], pos[order]
    # One entry per 25-mer; one found in two references belongs to neither.
    _, first, cnt = np.unique(codes, return_index=True, return_counts=True)
    multi = np.zeros(len(codes), bool)
    for f, n in zip(first[cnt > 1], cnt[cnt > 1]):
        if len(set(vid[f:f + n].tolist())) > 1:
            multi[f] = True
    codes, vid, pos, multi = codes[first], vid[first], pos[first], multi[first]
    vid = np.where(multi, -1, vid).astype(np.int16)
    return {'keys': keys, 'lens': np.array(lens), 'codes': codes, 'vid': vid, 'pos': pos}


IDX = None


def _init():
    global IDX
    IDX = load_index()


# ------------------------------------------------------------- HTTP / BGZF
def _get(url, a, b, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'Range': f'bytes={a}-{b}'})
            return urllib.request.urlopen(req, timeout=300).read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** k)


def _size(url):
    req = urllib.request.Request(url, method='HEAD')
    return int(urllib.request.urlopen(req, timeout=60).headers['Content-Length'])


def _inflate(raw, i):
    """Decompress whole BGZF blocks from raw[i:]; return (data, next offset)."""
    out = bytearray()
    while i + 18 <= len(raw) and raw[i:i + 4] == b'\x1f\x8b\x08\x04':
        bsize = struct.unpack_from('<H', raw, i + 16)[0] + 1
        if i + bsize > len(raw):
            break
        out += zlib.decompress(raw[i:i + bsize], 31)
        i += bsize
    return bytes(out), i


def _aligned(data, off, n=30):
    for _ in range(n):
        if off + 36 > len(data):
            return False
        bs, ref = struct.unpack_from('<ii', data, off)
        if not (32 < bs < 20000 and -1 <= ref < 100000):
            return False
        off += 4 + bs
    return True


def _first_refs(url, off):
    """refIDs of the first records after byte offset off (None if unreadable)."""
    raw = _get(url, off, off + 600_000)
    i = raw.find(MAGIC)
    if i < 0:
        return None
    data, _ = _inflate(raw, i)
    for o in range(min(len(data), 200_000)):
        if _aligned(data, o):
            refs = []
            for _ in range(30):
                bs, ref = struct.unpack_from('<ii', data, o)
                refs.append(ref)
                o += 4 + bs
            return refs
    return None


def find_unmapped_start(url, size):
    """Byte offset at or just before the first block of the unmapped tail."""
    lo, hi = 0, size
    while hi - lo > 8 * 1024 * 1024:
        mid = (lo + hi) // 2
        refs = _first_refs(url, mid)
        if refs is None:
            lo = mid  # unreadable spot: keep searching later in the file
        elif all(r == -1 for r in refs):
            hi = mid
        else:
            lo = mid
    return lo


# ------------------------------------------------------------- scanning
def _score(seqs, res):
    """seqs: list of uint8 base arrays (0-3, 255 = other). Updates res in place."""
    codes, vid, pos = IDX['codes'], IDX['vid'], IDX['pos']
    by_len = {}
    for s in seqs:
        by_len.setdefault(len(s), []).append(s)
    nv = len(IDX['keys'])
    for L, group in by_len.items():
        if L < K + STEP:
            continue
        M = np.stack(group)
        P = np.arange(0, L - K + 1, STEP)
        km = np.zeros((len(M), len(P)), np.uint64)
        bad = np.zeros((len(M), len(P)), bool)
        for t in range(K):
            col = M[:, P + t]
            bad |= col == 255
            km = (km << np.uint64(2)) | (col & 3).astype(np.uint64)
        flat = km.ravel()
        j = np.searchsorted(codes, flat)
        j = np.minimum(j, len(codes) - 1)
        hit = (codes[j] == flat) & ~bad.ravel()
        v = np.where(hit, vid[j], -2).reshape(km.shape)
        p = np.where(hit, pos[j], -1).reshape(km.shape)
        for vi in range(nv):
            m = v == vi
            reads = m.sum(1) >= MIN_HITS
            n = int(reads.sum())
            if not n:
                continue
            res['hits'][vi] += n
            hp = p[reads][m[reads]]
            res['bins'][vi] |= set((hp * NBINS // IDX['lens'][vi]).tolist())
            # Where each counted read sits: the mean of its matching positions,
            # in POSBIN-base bins, so reads can be assigned to the virus's genes.
            centre = np.nanmean(np.where(m[reads], p[reads], np.nan), axis=1)
            np.add.at(res['pos'][vi], (centre // POSBIN).astype(np.int64), 1)


def scan_one(rec, sample_mb):
    t0 = time.time()
    out = {'model_id': rec['model_id'], 'name': rec.get('name'), 'lineage': rec.get('lineage')}
    nv = len(IDX['keys'])
    res = {'hits': [0] * nv, 'bins': [set() for _ in range(nv)],
           'pos': [np.zeros(int(L) // POSBIN + 1, np.int64) for L in IDX['lens']]}
    try:
        url = rec['bam']
        size = _size(url)
        start = find_unmapped_start(url, size)
        out['unmapped_start'] = start
        if sample_mb:
            start = max(start, size - sample_mb * 1024 * 1024)
        out.update({'file_bytes': size, 'tail_start': start, 'mode': 'sample' if sample_mb else 'full'})
        raw, data, aligned = b'', b'', False
        unmapped = mapped_seen = 0
        seqs = []
        off = start
        while off < size:
            raw += _get(url, off, min(size, off + CHUNK) - 1)
            off = min(size, off + CHUNK)
            if not aligned:
                i = raw.find(MAGIC)
                if i < 0:
                    continue
                raw = raw[i:]
            dec, used = _inflate(raw, 0)
            raw = raw[used:]
            data += dec
            o = 0
            if not aligned:
                o = next((x for x in range(min(len(data), 300_000)) if _aligned(data, x)), None)
                if o is None:
                    continue
                aligned = True
            while o + 36 <= len(data):
                bs = struct.unpack_from('<i', data, o)[0]
                if o + 4 + bs > len(data):
                    break
                ref, _, l_name, _, _, n_cig, flag, l_seq = struct.unpack_from('<iiBBHHHi', data, o + 4)
                if flag & 4 and ref == -1:
                    unmapped += 1
                    p = o + 36 + l_name + 4 * n_cig
                    b = np.frombuffer(data, np.uint8, (l_seq + 1) // 2, p)
                    nib = np.empty(len(b) * 2, np.uint8)
                    nib[0::2] = b >> 4
                    nib[1::2] = b & 15
                    seqs.append(_NIB[nib[:l_seq]])
                    if len(seqs) >= 50_000:
                        _score(seqs, res)
                        seqs = []
                else:
                    mapped_seen += 1
                o += 4 + bs
            data = data[o:]
        if seqs:
            _score(seqs, res)
        out.update({'unmapped_scanned': unmapped, 'mapped_skipped': mapped_seen,
                    'bytes_read': size - start})
        out['hits'] = {IDX['keys'][i]: n for i, n in enumerate(res['hits']) if n}
        out['breadth'] = {IDX['keys'][i]: round(len(res['bins'][i]) / NBINS, 2) for i, n in enumerate(res['hits']) if n}
        out['posbin'] = POSBIN
        out['pos'] = {IDX['keys'][i]: [[int(b), int(c)] for b, c in enumerate(res['pos'][i]) if c] for i, n in enumerate(res['hits']) if n >= 3}
    except Exception as e:
        out['error'] = f'{type(e).__name__}: {e}'
    out['seconds'] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('lines'); ap.add_argument('out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 4)
    ap.add_argument('--sample-mb', type=int, default=0)
    ap.add_argument('--names', default='')
    a = ap.parse_args()
    lines = json.load(open(a.lines))
    if a.names:
        want = set(a.names.split(','))
        lines = [l for l in lines if l.get('name') in want or l['model_id'] in want]
    # Resume: keep finished lines; a line that errored is retried.
    done = {}
    try:
        for r in json.load(open(a.out)):
            if not r.get('error'):
                done[r['model_id']] = r
    except Exception:
        pass
    todo = [l for l in lines if l['model_id'] not in done]
    print(f"todo={len(todo)} (done={len(done)}) workers={a.workers} mode={'sample ' + str(a.sample_mb) + ' MB' if a.sample_mb else 'full tail'}", flush=True)
    res = list(done.values())
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers, initializer=_init) as ex:
        futs = [ex.submit(scan_one, l, a.sample_mb) for l in todo]
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            res.append(r)
            tag = r.get('error') or f"{r['seconds']}s unmapped={r['unmapped_scanned']:,} {r['bytes_read'] / 1e6:.0f} MB " + ' '.join(f"{k}={v}({r['breadth'][k]})" for k, v in r['hits'].items())
            print(f"[{i}/{len(todo)}] {str(r.get('name'))[:14]:14s} {tag}", flush=True)
            json.dump(res, open(a.out, 'w'))
    print(f"DONE {len(res)} lines in {(time.time() - t0) / 60:.1f} min -> {a.out}", flush=True)


if __name__ == '__main__':
    main()
