#!/bin/sh
# Commit and push V2, then rebuild, commit and push V1.
# usage: scripts/release.sh VERSION "commit subject" ["trailer lines"]
set -e
V="$1"; MSG="$2"; TRAILER="$3"
cd "$(dirname "$0")/.."
python3 scripts/data_version.py --write
scripts/precheck.sh
git add -A
git commit -q -m "v.$V: $MSG

$TRAILER"
git push -q
python3 scripts/build_v1.py | tail -2
cd "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app"
git add -A
git commit -q -m "Sync with Correlate V2 v.$V

$TRAILER"
git push -q
git log --oneline | head -1
cd - > /dev/null; git log --oneline | head -1; git status --short | head -3
