#!/usr/bin/env bash
# Bootstrap a fresh Ubuntu EC2 instance and run the retroelement quantification.
# Run from inside the te_aws directory:   bash setup_and_run.sh verify
#                                        bash setup_and_run.sh full
set -euo pipefail
MODE="${1:-verify}"
WORKERS="${WORKERS:-12}"

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

case "$MODE" in
  verify) LINES=lines_verify.json; OUT=counts_verify.json ;;
  full)   LINES=lines_all.json;    OUT=counts_all.json ;;
  *) echo "usage: bash setup_and_run.sh [verify|full]"; exit 1 ;;
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
