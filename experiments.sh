#!/usr/bin/env bash
# Every experiment in the README, one command per line. Every setting that has ever
# changed is spelled out, so changing a default never changes an old experiment.
#   ./experiments.sh                 print all commands
#   ./experiments.sh run [N] [SET]   run them, N in parallel (default 8); SET filters by
#                                    results folder, e.g. "repro", "rerun2", "ablation2"
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-}; JOBS=${2:-8}; ONLY=${3:-}
PY=.venv/bin/python
B13="--bacon-fixes sign_aware signed_constancy bounded_retry"                      # Tests 1-3
B14="--bacon-fixes sign_aware signed_constancy bounded_retry intercept_vs_std"     # Test 4, rerun1
BOTH="--bacon-variants strict robust"

cmds=(
# --- Reproductions: each test exactly as it was originally run --------------------------
"$PY run_test.py --system decay      --out results/repro --tag test1       $B13 --bacon-variants strict --sindy-threshold 0.1 --sim-solver LSODA"
"$PY run_test.py --system logistic   --out results/repro --tag test2       $B13 $BOTH --sim-solver LSODA"
"$PY run_test.py --system oscillator --out results/repro --tag test3       $B13 $BOTH"
"$PY run_test.py --system oscillator --out results/repro --tag test4_1traj $B14 $BOTH"
"$PY run_test.py --system oscillator --out results/repro --tag test4_3traj $B14 $BOTH --n-traj 3"
)

SYSTEMS=("decay" "logistic" "oscillator" "oscillator --n-traj 3")
sfx() { [[ "$1" == *n-traj* ]] && echo "_3traj" || echo ""; }
without() { echo "$1" | tr ' ' '\n' | grep -vx "$2" | tr '\n' ' '; }

# --- Rerun 1 (after Test 4): fixes B1-B4, both variants -------------------------------
for s in "${SYSTEMS[@]}"; do
  t=$(sfx "$s"); cmds+=("$PY run_test.py --system $s --out results/rerun1 ${t:+--tag ${t#_}} $B14 $BOTH")
done
# --- Ablation 1: original BACON.1, and each of B1-B4 off from rerun 1 -------------------
for s in "${SYSTEMS[@]}"; do
  t=$(sfx "$s")
  cmds+=("$PY run_test.py --system $s --out results/ablation1 --tag original$t --bacon-fixes none $BOTH")
  for off in sign_aware signed_constancy bounded_retry intercept_vs_std; do
    cmds+=("$PY run_test.py --system $s --out results/ablation1 --tag no_$off$t --bacon-fixes $(without "sign_aware signed_constancy bounded_retry intercept_vs_std" $off) $BOTH")
  done
done

# --- Rerun 2 (current defaults): fixes B1, B3, B4; robust variant only ------------------
B134="sign_aware bounded_retry intercept_vs_std"
for s in "${SYSTEMS[@]}"; do
  t=$(sfx "$s"); cmds+=("$PY run_test.py --system $s --out results/rerun2 ${t:+--tag ${t#_}} --bacon-fixes $B134 --bacon-variants robust")
done
# --- Ablation 2: each of B1, B3, B4 off from rerun 2 -----------------------------------
# ("all off" = ablation1/original; "plus B2" = rerun1)
for s in "${SYSTEMS[@]}"; do
  t=$(sfx "$s")
  for off in sign_aware bounded_retry intercept_vs_std; do
    cmds+=("$PY run_test.py --system $s --out results/ablation2 --tag no_$off$t --bacon-fixes $(without "$B134" $off) --bacon-variants robust")
  done
done

# --- Rerun 3 (current defaults): fixes B1, B3, B4, B5; robust only; + power decay -------
B1345="sign_aware bounded_retry intercept_vs_std rate_seed"
SYSTEMS3=("${SYSTEMS[@]}" "powerdecay")
for s in "${SYSTEMS3[@]}"; do
  t=$(sfx "$s"); cmds+=("$PY run_test.py --system $s --out results/rerun3 ${t:+--tag ${t#_}} --bacon-fixes $B1345 --bacon-variants robust")
done
# Informed PySINDy baseline: polynomial library plus |x|^0.5, |x|^1.5, |x|^2.5
cmds+=("$PY run_test.py --system powerdecay --out results/rerun3 --tag halfpowers --bacon-fixes $B1345 --bacon-variants robust --sindy-library poly+halfpowers")
# --- Ablation 3: each of B1, B3, B4 off from rerun 3, every system ----------------------
# ("without B5" on the old systems = rerun 2; power decay gets it here, plus "all off")
for s in "${SYSTEMS3[@]}"; do
  t=$(sfx "$s")
  for off in sign_aware bounded_retry intercept_vs_std; do
    cmds+=("$PY run_test.py --system $s --out results/ablation3 --tag no_$off$t --bacon-fixes $(without "$B1345" $off) --bacon-variants robust")
  done
done
cmds+=("$PY run_test.py --system powerdecay --out results/ablation3 --tag no_rate_seed --bacon-fixes sign_aware bounded_retry intercept_vs_std --bacon-variants robust")
cmds+=("$PY run_test.py --system powerdecay --out results/ablation3 --tag original --bacon-fixes none --bacon-variants strict robust")

# --- Rerun 4 = Test 7 (current defaults): fixes B1, B3-B6; robust only -----------------
B13456="sign_aware bounded_retry intercept_vs_std rate_seed robust_line"
for s in "${SYSTEMS3[@]}"; do
  t=$(sfx "$s"); cmds+=("$PY run_test.py --system $s --out results/rerun4 ${t:+--tag ${t#_}} --bacon-fixes $B13456 --bacon-variants robust")
done
cmds+=("$PY run_test.py --system powerdecay --out results/rerun4 --tag halfpowers --bacon-fixes $B13456 --bacon-variants robust --sindy-library poly+halfpowers")
# --- Ablation 4: each of B1, B3, B4, B5 off from rerun 4 ("without B6" = rerun 3) -------
for s in "${SYSTEMS3[@]}"; do
  t=$(sfx "$s")
  for off in sign_aware bounded_retry intercept_vs_std rate_seed; do
    cmds+=("$PY run_test.py --system $s --out results/ablation4 --tag no_$off$t --bacon-fixes $(without "$B13456" $off) --bacon-variants robust")
  done
done

if [ "$MODE" = "run" ]; then
  mkdir -p results/logs
  printf '%s\n' "${cmds[@]}" | { [ -n "$ONLY" ] && grep -- "--out results/$ONLY " || cat; } |
    xargs -P "$JOBS" -I{} bash -c '
      cmd="$1"; log="results/logs/$(echo "$cmd" | md5sum | cut -c1-10).log"
      eval "$cmd" > "$log" 2>&1 && echo "done: $cmd" || echo "FAILED ($log): $cmd"' _ {}
else
  printf '%s\n' "${cmds[@]}"
fi
