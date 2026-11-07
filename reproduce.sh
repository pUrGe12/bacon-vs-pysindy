#!/usr/bin/env bash
# Rerun experiments and check them against MANIFEST.json. Run ./setup.sh first.
#   ./reproduce.sh            everything (about 4.5 h on 20 cores; Test 10 is ~4 h)
#   ./reproduce.sh SET        one results folder, e.g. repro, rerun6b, test10
#   JOBS=8 ./reproduce.sh     number of parallel jobs (default: number of cores - 2)
set -euo pipefail
cd "$(dirname "$0")"
SET=${1:-}
JOBS=${JOBS:-$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))}

.venv/bin/python check_original.py
./experiments.sh run "$JOBS" $SET
.venv/bin/python manifest.py check $SET
if [ -z "$SET" ] || [ "$SET" = "repro" ]; then
  .venv/bin/python summarize.py repro
fi
