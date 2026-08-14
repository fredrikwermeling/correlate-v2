#!/usr/bin/env bash
# Bootstrap a fresh Ubuntu EC2 instance and run the retroelement quantification.
# Run from inside the te_aws directory:   bash setup_and_run.sh verify
#                                        bash setup_and_run.sh full
set -euo pipefail
MODE="${1:-verify}"
# One process per core: the counting loop is Python, so processes (not
# threads) are what actually use the cores, and more than one per core just
# adds contention.
WORKERS="${WORKERS:-$(nproc)}"

if ! python3 -c "import pysam" 2>/dev/null; then
  echo "== installing dependencies =="
  sudo apt-get update -qq
  sudo apt-get install -y -qq python3-pip tmux >/dev/null
  # Ubuntu 24.04 marks the system Python "externally managed" (PEP 668) and
  # refuses a plain pip install, so fall back through the options that work:
  # the distro package first, then pip with the override.
  pip3 install --quiet pysam 2>/dev/null \
    || sudo apt-get install -y -qq python3-pysam >/dev/null 2>&1 \
    || pip3 install --quiet --break-system-packages pysam
fi
python3 -c "import pysam" 2>/dev/null || { echo "ERROR: pysam did not install"; exit 1; }

# auto: verify first, and go on to the full run ONLY if it passed. Lets the
# whole thing be started in one go without watching the gate.
if [ "$MODE" = "auto" ]; then
  echo "== auto: verify, then full only if verify passes =="
  bash "$0" verify 2>&1 | tee verify.log
  if grep -q "VERIFIED" verify.log; then
    echo; echo "== verify passed, starting the full run =="
    exec bash "$0" full
  else
    echo; echo "== VERIFY DID NOT PASS - full run NOT started =="
    echo "== read verify.log before running anything else =="
    exit 1
  fi
fi

case "$MODE" in
  verify) LINES=lines_verify.json; OUT=counts_verify.json ;;
  full)   LINES=lines_all.json;    OUT=counts_all.json ;;
  *) echo "usage: bash setup_and_run.sh [verify|full|auto]"; exit 1 ;;
esac

echo "== running $MODE ($(python3 -c "import json;print(len(json.load(open('$LINES'))))") lines, $WORKERS workers) =="
python3 te_quant_stream.py "$LINES" te_panel.json "$OUT" "$WORKERS"

if [ "$MODE" = "verify" ]; then
  echo; echo "== comparing against the laptop pilot =="
  python3 verify_against_pilot.py "$OUT" pilot_reference_counts.json
else
  echo; echo "== validation =="
  python3 te_analyze.py "$OUT"
fi
