#!/usr/bin/env python3
"""Count retroelement transcription per cell line by STREAMING each CCLE BAM once.

Why streaming rather than the laptop's targeted fetches: DepMap published an
index for only 58 of 1,032 RNA BAMs, and without an index there is no random
access. Inside AWS us-east-1 that stops mattering - the bucket is in the same
region, so reading a whole BAM is fast and costs nothing in transfer. One pass
per file, no index, nothing written to disk.

The BAM is coordinate-sorted, and so is the panel, so interval matching is a
single advancing pointer per contig rather than a lookup per read.

Usage:  python3 te_quant_stream.py lines.json te_panel.json out.json [workers]
"""
import json, sys, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam

MIN_MAPQ = 20          # same threshold as the validated pilot
ACTIVE_CPM = 0.5       # an element counts as "on" above this

def load_panel(path):
    panel = json.load(open(path))
    by_chrom = defaultdict(list)
    for e in panel:
        by_chrom[str(e['chrom'])].append((int(e['start']), int(e['end']), e['id']))
    for c in by_chrom:
        by_chrom[c].sort()
    return panel, by_chrom

def count_one(rec, by_chrom):
    t0 = time.time()
    out = {'model_id': rec['model_id'], 'name': rec.get('name'),
           'lineage': rec.get('lineage')}
    try:
        bam = pysam.AlignmentFile(rec['bam'], "rb")     # no index needed
    except Exception as e:
        out['error'] = f'open: {e}'; return out
    counts_all = defaultdict(int); counts_uq = defaultdict(int)
    total = 0
    cur_chrom, iv, ptr, active, n_iv = None, [], 0, [], 0
    r_len = 300      # generous read-length pad for the fast reject
    try:
        for r in bam.fetch(until_eof=True):            # sequential stream
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            total += 1
            c = r.reference_name
            if c != cur_chrom:
                cur_chrom = c; iv = by_chrom.get(str(c), []); ptr = 0; active = []
                n_iv = len(iv)
                if not iv: continue
            if not iv: continue
            s = r.reference_start
            # Fast reject: reads are coordinate-ordered, so if this read starts
            # before the next unentered interval and nothing is currently open,
            # it cannot overlap anything. One comparison for the common case.
            if not active and ptr < n_iv and s < iv[ptr][0] - r_len:
                continue
            e = r.reference_end or (s + 1)
            while ptr < n_iv and iv[ptr][0] <= e:
                active.append(iv[ptr]); ptr += 1
            if active:
                if active[0][1] < s:                     # prune only when stale
                    active = [x for x in active if x[1] >= s]
                mq = r.mapping_quality >= MIN_MAPQ
                for st, en, eid in active:
                    if st <= e and en >= s:
                        counts_all[eid] += 1
                        if mq: counts_uq[eid] += 1
    except Exception as e:
        out['error'] = f'stream: {e}'
    bam.close()
    out.update({'total_mapped': total,
                'counts_all': dict(counts_all), 'counts_unique': dict(counts_uq),
                'seconds': round(time.time() - t0, 1)})
    return out

def main():
    lines = json.load(open(sys.argv[1]))
    panel, by_chrom = load_panel(sys.argv[2])
    outp = sys.argv[3]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    # resume: skip lines already done, so a dropped connection costs one file
    done = {}
    try:
        for r in json.load(open(outp)): done[r['model_id']] = r
    except Exception: pass
    todo = [l for l in lines if l['model_id'] not in done]
    print(f"panel={len(panel)} loci | todo={len(todo)} (done={len(done)}) | workers={workers}", flush=True)
    res = list(done.values())
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(count_one, l, by_chrom): l for l in todo}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result(); res.append(r)
            tag = r.get('error') or f"{r['seconds']}s mapped={r['total_mapped']:,}"
            print(f"[{i}/{len(todo)}] {str(r.get('name'))[:18]:18s} {tag}", flush=True)
            json.dump(res, open(outp, 'w'))
    print(f"DONE {len(res)} lines in {(time.time()-t0)/60:.1f} min -> {outp}", flush=True)

if __name__ == '__main__':
    main()
