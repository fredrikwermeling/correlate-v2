#!/usr/bin/env bash
# Run the virus scan. Works on a laptop or a fresh Ubuntu EC2 instance.
# From inside the virus_aws directory:
#   bash setup_and_run.sh verify          the five control cell lines
#   bash setup_and_run.sh full            every cell line, whole unmapped tail
#   bash setup_and_run.sh auto            verify, then full only if it passed
#   SAMPLE_MB=25 bash setup_and_run.sh full    quick screen: last 25 MB per file
set -euo pipefail
MODE="${1:-verify}"
WORKERS="${WORKERS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}"
SAMPLE="${SAMPLE_MB:-0}"

if ! python3 -c "import numpy" 2>/dev/null; then
  echo "== installing numpy =="
  sudo apt-get update -qq && sudo apt-get install -y -qq python3-pip tmux >/dev/null
  pip3 install --quiet numpy 2>/dev/null \
    || sudo apt-get install -y -qq python3-numpy >/dev/null 2>&1 \
    || pip3 install --quiet --break-system-packages numpy
fi
python3 fetch_refs.py

if [ "$MODE" = "auto" ]; then
  bash "$0" verify 2>&1 | tee verify.log
  if grep -q "VERIFIED" verify.log; then
    echo; echo "== verify passed, starting the full run =="
    exec bash "$0" full
  fi
  echo; echo "== VERIFY DID NOT PASS - full run NOT started; read verify.log =="
  exit 1
fi

case "$MODE" in
  verify) LINES=lines_verify.json; OUT=counts_verify.json; SAMPLE=25 ;;
  full)   LINES=lines_all.json;    OUT=$([ "$SAMPLE" -gt 0 ] && echo counts_sample.json || echo counts_all.json) ;;
  *) echo "usage: bash setup_and_run.sh [verify|full|auto]"; exit 1 ;;
esac

python3 virus_scan.py "$LINES" "$OUT" --workers "$WORKERS" --sample-mb "$SAMPLE"
if [ "$MODE" = "verify" ]; then python3 verify_controls.py "$OUT"; else python3 analyze.py "$OUT"; fi
