#!/usr/bin/env bash
# Build the exact environment every result was produced with (versions in pins.json).
#   ./setup.sh            .venv from requirements.lock + BACON at the pinned commit
# Needs: python3.10 (the results were made with 3.10.12), git, internet access once.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3.10}

pin() { python3 -c "import json; d=json.load(open('pins.json')); print($1)"; }
BACON_URL=$(pin 'd["bacon"]["url"]')
BACON_COMMIT=$(pin 'd["bacon"]["commit"]')

# 1. Python environment, every package pinned
if [ ! -x .venv/bin/python ]; then
  "$PY" -m venv .venv
fi
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -r requirements.lock
echo "python: $(.venv/bin/python --version) (results made with $(pin 'd["python"]'))"

# 2. Original BACON, at the pinned commit, in external/BACON
if [ ! -d external/BACON/.git ]; then
  git clone -q "$BACON_URL" external/BACON
fi
git -C external/BACON fetch -q origin
git -C external/BACON checkout -q "$BACON_COMMIT"
echo "BACON: $BACON_URL @ $(git -C external/BACON rev-parse --short HEAD)"

# 3. Check the files we actually use from BACON and ODEBench against their pinned hashes
.venv/bin/python - <<'EOF'
import hashlib, json, sys
pins = json.load(open("pins.json"))
bad = 0
def check(path, want):
    global bad
    got = hashlib.sha256(open(path, "rb").read()).hexdigest()
    ok = got == want
    bad += not ok
    print(f"  {'ok ' if ok else 'BAD'} {path}")
for f, want in pins["bacon"]["files_sha256"].items():
    check(f"external/BACON/{f}", want)
check(pins["odebench"]["converted_file"], pins["odebench"]["converted_sha256"])
sys.exit(1 if bad else 0)
EOF
echo "setup done. Next: ./reproduce.sh (see README, Reproducing from scratch)"
