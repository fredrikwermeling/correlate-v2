#!/usr/bin/env bash
# Is the run alive, and how far along? Safe to run at any time, including
# while the job is running - it only reads.
cd "$(dirname "$0")"
echo "=========== TE run status  $(date -u '+%Y-%m-%d %H:%M UTC') ==========="

if pgrep -f te_quant_stream.py > /dev/null; then
  PID=$(pgrep -f te_quant_stream.py | head -1)
  echo "PROCESS : running (pid $PID)"
  CPU=$(ps -p "$PID" -o %cpu= 2>/dev/null | tr -d ' ')
  ET=$(ps -p "$PID" -o etime= 2>/dev/null | tr -d ' ')
  echo "          cpu ${CPU}%   running for ${ET}"
else
  echo "PROCESS : not running (finished, or never started)"
fi

# Network counter proves data is actually moving, not just a live process.
RX1=$(cat /sys/class/net/ens5/statistics/rx_bytes 2>/dev/null || echo 0)
sleep 5
RX2=$(cat /sys/class/net/ens5/statistics/rx_bytes 2>/dev/null || echo 0)
MBPS=$(( (RX2 - RX1) / 5 / 1048576 ))
echo "NETWORK : ${MBPS} MB/s inbound  $( [ "$MBPS" -gt 2 ] && echo '(streaming BAMs)' || echo '(idle - suspicious if the process is running)')"

for f in counts_all.json counts_verify.json; do
  [ -f "$f" ] || continue
  N=$(python3 -c "import json;print(len(json.load(open('$f'))))" 2>/dev/null || echo '?')
  TOT=$(python3 -c "import json;print(len(json.load(open('lines_all.json'))))" 2>/dev/null || echo 669)
  [ "$f" = counts_verify.json ] && TOT=5
  echo "RESULTS : $f  -> $N of $TOT cell lines done"
  if [ "$f" = counts_all.json ] && [ "$N" != '?' ] && [ "$N" -gt 0 ] 2>/dev/null; then
    W=$(nproc 2>/dev/null || echo 8)
    python3 - "$N" "$TOT" "$W" <<'PY'
import json, sys
n, tot, w = int(sys.argv[1]), int(sys.argv[2]), max(1, int(sys.argv[3]))
d = json.load(open('counts_all.json'))
secs = [r.get('seconds', 0) for r in d if r.get('seconds')]
if secs:
    per = sum(secs)/len(secs)
    # Per cell line PER WORKER. Wall clock is that divided by however many
    # run at once, which is one per core. This divisor was hardcoded to 12
    # from an earlier run and quietly under-reported the time left by a third.
    left = (tot - n) * per / w / 3600
    print(f"          mean {per:.0f}s per cell line per worker, {w} workers")
    print(f"          -> about {left:.1f} h remaining ({per/w:.0f}s per line of wall clock)")
errs = [r for r in d if r.get('error')]
print(f"          errors: {len(errs)}" + (f"  e.g. {errs[0]['error'][:60]}" if errs else ""))
PY
  fi
done
echo "=================================================================="
