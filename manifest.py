"""MANIFEST.json: every experiment, the exact command that makes it, its output files,
and a fingerprint of its results, so anyone can rerun and check.

    .venv/bin/python manifest.py make            # rebuild MANIFEST.json from results/
    .venv/bin/python manifest.py check [SET]     # do the result files on disk match it?
                                                 # SET = results folder, e.g. test10

The fingerprint is a SHA-256 of each results CSV's outcome columns (learned equations,
scores, errors), sorted by run, with run times left out and numbers written to 6
significant digits. Rerunning on the same setup gives the same fingerprint; see README,
"Reproducing from scratch", for what to expect on other machines.
"""
import datetime
import hashlib
import json
import os
import shlex
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
KEYS = ["n_traj", "noise", "seed", "deriv", "method"]
OUTCOME = ["equations", "structure_ok", "recovered", "coef_err", "sim_err", "bacon_invariant",
           "pred_err", "status"]  # last two: static laws (Test 12) only
# Results written by the original runs of Tests 1-4; experiments.sh reproduces each into
# results/repro/ and `summarize.py repro` compares them row by row.
ORIGINALS = {
    "results/decay.csv": "results/repro/decay_test1.csv",
    "results/logistic.csv": "results/repro/logistic_test2.csv",
    "results/oscillator.csv": "results/repro/oscillator_test3.csv",
    "results/oscillator_test4_1traj.csv": "results/repro/oscillator_test4_1traj.csv",
    "results/oscillator_test4_3traj.csv": "results/repro/oscillator_test4_3traj.csv",
}


def fingerprint(path):
    d = pd.read_csv(path)
    keys = [k for k in KEYS if k in d]
    cols = keys + [c for c in OUTCOME if c in d]
    d = d[cols].sort_values(keys).reset_index(drop=True)

    def fmt(v):
        if isinstance(v, float):
            return "nan" if v != v else f"{v:.6g}"
        return str(v)

    text = "\n".join(",".join(fmt(v) for v in row) for row in d.itertuples(index=False))
    return {"rows": len(d), "sha256": hashlib.sha256(text.encode()).hexdigest()}


def experiments():
    """Parse the commands printed by experiments.sh into (set, command, output csv)."""
    out = subprocess.run(["./experiments.sh"], cwd=HERE, capture_output=True, text=True,
                         check=True).stdout.splitlines()
    exps = []
    for cmd in out:
        a = shlex.split(cmd)
        get = lambda flag, default="": a[a.index(flag) + 1] if flag in a else default
        folder, system, tag = get("--out", "results"), get("--system"), get("--tag")
        name = system + (f"_{tag}" if tag else "")
        exps.append({"set": os.path.basename(folder), "command": cmd,
                     "results": f"{folder}/{name}.csv",
                     "summary": f"{folder}/{name}_summary.csv",
                     "run_manifest": f"{folder}/{name}_manifest.json"})
    return exps


def _git_head():
    r = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True)
    return r.stdout.strip() or None


def make():
    exps = experiments()
    for e in exps:
        full = os.path.join(HERE, e["results"])
        # A set defined in experiments.sh but never run (Test 11, stopped) has no results.
        e["fingerprint"] = fingerprint(full) if os.path.exists(full) else "not run"
    originals = {f: {"reproduced_by": r, "fingerprint": fingerprint(os.path.join(HERE, f))}
                 for f, r in ORIGINALS.items()}
    lock = open(os.path.join(HERE, "requirements.lock"), "rb").read()
    man = {
        "about": "Every experiment in this repo: command, outputs, result fingerprints. "
                 "Rebuild with `manifest.py make`, check with `manifest.py check`.",
        "made": datetime.datetime.now().isoformat(timespec="seconds"),
        "made_at_commit": _git_head(),
        "pins": json.load(open(os.path.join(HERE, "pins.json"))),
        "requirements_lock_sha256": hashlib.sha256(lock).hexdigest(),
        "how_to_reproduce": ["./setup.sh", "./reproduce.sh [SET]"],
        "fingerprint": "SHA-256 of the outcome columns " + ", ".join(OUTCOME) +
                       " sorted by " + ", ".join(KEYS) + "; numbers to 6 significant "
                       "digits; run times left out",
        "sets": sorted({e["set"] for e in exps}),
        "experiments": exps,
        "original_runs_tests_1_to_4": originals,
    }
    with open(os.path.join(HERE, "MANIFEST.json"), "w") as fh:
        json.dump(man, fh, indent=1)
    print(f"MANIFEST.json: {len(exps)} experiments in {len(man['sets'])} sets, "
          f"{len(originals)} original runs")


def check(only=None):
    man = json.load(open(os.path.join(HERE, "MANIFEST.json")))
    items = [(e["results"], e["fingerprint"]) for e in man["experiments"]
             if (only is None or e["set"] == only) and e["fingerprint"] != "not run"]
    if only is None:
        items += [(f, o["fingerprint"]) for f, o in man["original_runs_tests_1_to_4"].items()]
    same, differ, missing = 0, [], []
    for path, want in items:
        full = os.path.join(HERE, path)
        if not os.path.exists(full):
            missing.append(path)
            continue
        got = fingerprint(full)
        if got == want:
            same += 1
        else:
            differ.append(path)
    print(f"{same} of {len(items)} result files match MANIFEST.json"
          + (f"; {len(differ)} differ" if differ else "")
          + (f"; {len(missing)} missing" if missing else ""))
    for p in differ:
        print(f"  DIFFERS: {p}")
    for p in missing[:20]:
        print(f"  MISSING: {p}")
    return not differ and not missing


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make()
    elif sys.argv[1] == "check":
        sys.exit(0 if check(sys.argv[2] if len(sys.argv) > 2 else None) else 1)
    else:
        raise SystemExit(__doc__)
