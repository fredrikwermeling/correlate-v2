#!/usr/bin/env python3
"""Step one: does a CCLE BAM keep its unmapped reads, and how big is that part?

Reads the header (which aligner, which pipeline) and finds the unmapped tail
with a few HTTP range reads. Nothing large is downloaded.

Usage:  python3 check_one.py [cell line name]      (default RT4)
"""
import json
import os
import struct
import sys

import virus_scan as vs

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'RT4'
    rec = next(l for l in json.load(open(os.path.join(HERE, 'lines_all.json'))) if l['name'] == name)
    url = rec['bam']
    raw = vs._get(url, 0, 4_000_000)
    data, _ = vs._inflate(raw, 0)
    l_text = struct.unpack_from('<i', data, 4)[0]
    progs = sorted({f.split(':', 1)[1] for line in data[8:8 + l_text].decode(errors='replace').splitlines()
                    if line.startswith('@PG') for f in line.split('\t') if f.startswith('PN:')})
    size = vs._size(url)
    start = vs.find_unmapped_start(url, size)
    refs = vs._first_refs(url, max(start, size - 6_000_000)) or []
    has_tail = bool(refs) and all(r == -1 for r in refs)
    print(f'{name}: {size / 1e9:.1f} GB, programs: {", ".join(progs)}')
    print(f'unmapped tail: {"YES" if has_tail else "NO"}, {(size - start) / 1e9:.2f} GB ({100 * (size - start) / size:.1f}% of the file)')


if __name__ == '__main__':
    main()
