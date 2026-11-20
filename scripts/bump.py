#!/usr/bin/env python3
"""Bump the version in index.html: badge, app.js cache-bust and changelog.

  python3 scripts/bump.py 92.5 92.6 <<'END'
  <li>First change, in plain user language</li>
  <li>Second change</li>
  END
"""
import datetime
import os
import sys

old, new = sys.argv[1], sys.argv[2]
date = sys.argv[3] if len(sys.argv) > 3 else datetime.date.today().isoformat()
items = [l.strip() for l in sys.stdin.read().splitlines() if l.strip()]
path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "index.html")
s = open(path, encoding="utf-8").read()

badge_old, badge_new = f">v.{old}<", f">v.{new}<"
bust_old, bust_new = f'app.js?v=v{old.replace(".", "_")}"', f'app.js?v=v{new.replace(".", "_")}"'
head_old = f"                <p><strong>v.{old} ("
for needle in (badge_old, bust_old, head_old):
    if s.count(needle) != 1:
        sys.exit(f"bump: expected exactly one {needle!r}, found {s.count(needle)}")
block = (f"                <p><strong>v.{new} ({date})</strong></p><ul>\n"
         + "".join(f"                    {i}\n" for i in items) + "                </ul>\n")
s = s.replace(badge_old, badge_new).replace(bust_old, bust_new).replace(head_old, block + head_old)
open(path, "w", encoding="utf-8").write(s)
print(f"bump: v.{old} -> v.{new}, {len(items)} changelog lines")
