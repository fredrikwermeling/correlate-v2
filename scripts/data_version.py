#!/usr/bin/env python3
"""Keep the data-version token in index.html in step with web_data/.

The app adds this token to every data download, so a browser fetches the
data again only when the data itself has changed, not on every release.

  python3 scripts/data_version.py          check, exit 1 when out of date
  python3 scripts/data_version.py --write  update index.html

scripts/precheck.sh (the pre-push hook) runs the check, so a data change
cannot be pushed with a stale token.
"""
import hashlib
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INDEX = os.path.join(ROOT, "index.html")
META = re.compile(r'(<meta name="data-version" content=")([^"]*)(">)')


def current():
    files = subprocess.check_output(["git", "ls-files", "-z", "web_data"], cwd=ROOT).decode().split("\0")
    # Images carry their own ?v= in the markup and are not data downloads.
    files = sorted(f for f in files if f and not f.lower().endswith((".png", ".jpg", ".svg")))
    h = hashlib.sha1()
    for f in files:
        h.update(f.encode() + b"\0")
        with open(os.path.join(ROOT, f), "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        h.update(b"\0")
    return h.hexdigest()[:10]


def main():
    html = open(INDEX, encoding="utf-8").read()
    m = META.search(html)
    if not m:
        sys.exit('data_version: <meta name="data-version"> is missing from index.html')
    want = current()
    if m.group(2) == want:
        print("data_version: ok", want)
        return
    if "--write" in sys.argv:
        open(INDEX, "w", encoding="utf-8").write(META.sub(lambda x: x.group(1) + want + x.group(3), html, count=1))
        print("data_version: updated to", want)
        return
    sys.exit("data_version: web_data changed, run  python3 scripts/data_version.py --write")


if __name__ == "__main__":
    main()
