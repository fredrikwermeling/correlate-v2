#!/usr/bin/env bash
# How far along is the run? Read-only, safe at any time.
cd "$(dirname "$0")"
echo "=========== virus scan status  $(date -u '+%Y-%m-%d %H:%M UTC') ==========="
pgrep -f virus_scan.py >/dev/null && echo "PROCESS : running" || echo "PROCESS : not running"
TOT=$(python3 -c "import json;print(len(json.load(open('lines_all.json'))))")
for f in counts_all.json counts_sample.json counts_verify.json; do
  [ -f "$f" ] || continue
  python3 - "$f" "$TOT" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); tot = int(sys.argv[2])
if 'verify' in sys.argv[1]: tot = 5
ok = [r for r in d if not r.get('error')]
gb = sum(r.get('bytes_read', 0) for r in ok) / 1e9
print(f"RESULTS : {sys.argv[1]}  {len(ok)} of {tot} done, {len(d) - len(ok)} errors, {gb:.0f} GB read")
PY
done
echo "======================================================================="
