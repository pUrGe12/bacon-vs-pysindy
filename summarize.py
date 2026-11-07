"""Build the README tables from the result files.

    .venv/bin/python summarize.py repro               # do reproductions match the original runs?
    .venv/bin/python summarize.py current rerun2      # one results folder, every system
    .venv/bin/python summarize.py ablation1           # each fix off, from rerun 1 (B1-B4)
    .venv/bin/python summarize.py ablation2           # each fix off, from rerun 2 (B1, B3, B4)
    .venv/bin/python summarize.py ablation3           # each fix off, from rerun 3 (B1, B3-B5)
    .venv/bin/python summarize.py ablation4           # each fix off, from rerun 4 (B1, B3-B6)
    .venv/bin/python summarize.py ablation5           # each fix off, from rerun 5 (B1, B3, B5-B7)
    .venv/bin/python summarize.py verify              # regenerated results == committed ones?
"""
import io
import os
import subprocess
import sys

import pandas as pd

KEY = ["noise", "seed", "deriv", "method"]
CHECK = ["equations", "recovered", "structure_ok"]
NUM = ["coef_err", "sim_err"]

# reproduction file -> file written by the original run (committed in git)
REPRO = {
    "Test 1 (decay)": ("results/repro/decay_test1.csv", "results/decay.csv"),
    "Test 2 (logistic)": ("results/repro/logistic_test2.csv", "results/logistic.csv"),
    "Test 3 (oscillator)": ("results/repro/oscillator_test3.csv", "results/oscillator.csv"),
    "Test 4 (oscillator, 1 traj)": ("results/repro/oscillator_test4_1traj.csv",
                                    "results/oscillator_test4_1traj.csv"),
    "Test 4 (oscillator, 3 traj)": ("results/repro/oscillator_test4_3traj.csv",
                                    "results/oscillator_test4_3traj.csv"),
}
SYSTEMS = [("decay", ""), ("logistic", ""), ("oscillator", ""), ("oscillator", "_3traj")]
# (label, file pattern) rows for each ablation table
ABLATION_ROWS = {
    "ablation1": [
        ("all fixes B1-B4 (rerun 1)", "results/rerun1/{s}{x}.csv"),
        ("no fixes (original BACON.1)", "results/ablation1/{s}_original{x}.csv"),
        ("without B1 sign_aware", "results/ablation1/{s}_no_sign_aware{x}.csv"),
        ("without B2 signed_constancy", "results/ablation1/{s}_no_signed_constancy{x}.csv"),
        ("without B3 bounded_retry", "results/ablation1/{s}_no_bounded_retry{x}.csv"),
        ("without B4 intercept_vs_std", "results/ablation1/{s}_no_intercept_vs_std{x}.csv"),
    ],
    "ablation5": [
        ("defaults B1, B3, B5, B6, B7 (rerun 5)", "results/rerun5/{s}{x}.csv"),
        ("no fixes (original BACON.1)", "results/{orig}/{s}_original{x}.csv"),
        ("without B1 sign_aware", "results/ablation5/{s}_no_sign_aware{x}.csv"),
        ("without B3 bounded_retry", "results/ablation5/{s}_no_bounded_retry{x}.csv"),
        ("without B5 rate_seed", "results/ablation5/{s}_no_rate_seed{x}.csv"),
        ("without B6 robust_line", "results/ablation5/{s}_no_robust_line{x}.csv"),
        ("without B7 own_rate_seed", "results/ablation5/{s}_no_own_rate_seed{x}.csv"),
        ("B1, B3-B6 with B4 (rerun 4)", "results/rerun4/{s}{x}.csv"),
    ],
    "ablation4": [
        ("defaults B1, B3-B6 (rerun 4)", "results/rerun4/{s}{x}.csv"),
        ("no fixes (original BACON.1)", "results/{orig}/{s}_original{x}.csv"),
        ("without B1 sign_aware", "results/ablation4/{s}_no_sign_aware{x}.csv"),
        ("without B3 bounded_retry", "results/ablation4/{s}_no_bounded_retry{x}.csv"),
        ("without B4 intercept_vs_std", "results/ablation4/{s}_no_intercept_vs_std{x}.csv"),
        ("without B5 rate_seed", "results/ablation4/{s}_no_rate_seed{x}.csv"),
        ("without B6 robust_line (rerun 3)", "results/rerun3/{s}{x}.csv"),
    ],
    "ablation3": [
        ("defaults B1, B3, B4, B5 (rerun 3)", "results/rerun3/{s}{x}.csv"),
        ("no fixes (original BACON.1)", "results/{orig}/{s}_original{x}.csv"),
        ("without B1 sign_aware", "results/ablation3/{s}_no_sign_aware{x}.csv"),
        ("without B3 bounded_retry", "results/ablation3/{s}_no_bounded_retry{x}.csv"),
        ("without B4 intercept_vs_std", "results/ablation3/{s}_no_intercept_vs_std{x}.csv"),
        ("without B5 rate_seed", "results/{norate}/{s}{x}.csv"),
    ],
    "ablation2": [
        ("defaults B1, B3, B4 (rerun 2)", "results/rerun2/{s}{x}.csv"),
        ("no fixes (original BACON.1)", "results/ablation1/{s}_original{x}.csv"),
        ("without B1 sign_aware", "results/ablation2/{s}_no_sign_aware{x}.csv"),
        ("without B3 bounded_retry", "results/ablation2/{s}_no_bounded_retry{x}.csv"),
        ("without B4 intercept_vs_std", "results/ablation2/{s}_no_intercept_vs_std{x}.csv"),
        ("plus B2 signed_constancy (rerun 1)", "results/rerun1/{s}{x}.csv"),
    ],
}


def repro():
    print("| Test | rows compared | identical | differing |")
    print("|---|---|---|---|")
    for name, (new, old) in REPRO.items():
        a, b = pd.read_csv(new), pd.read_csv(old)
        m = a.merge(b, on=KEY, suffixes=("_new", "_old"))
        same = _same(m)
        diff = [" ".join(map(str, r)) for r in m.loc[~same, KEY].itertuples(index=False)]
        print(f"| {name} | {len(m)} of {len(b)} | {int(same.sum())} | {', '.join(diff) or '—'} |")


def _same(m):
    same = pd.Series(True, index=m.index)
    for c in CHECK:
        same &= m[f"{c}_new"].astype(str) == m[f"{c}_old"].astype(str)
    for c in NUM:
        x, y = m[f"{c}_new"], m[f"{c}_old"]
        same &= ((x - y).abs() <= 1e-9 * (1 + y.abs())) | (x.isna() & y.isna()) | (x == y)
    return same


def verify():
    """Compare every committed results CSV with the file now on disk (ignores timing)."""
    files = [f for f in subprocess.run(["git", "ls-files", "results/*.csv"], capture_output=True,
                                       text=True, check=True).stdout.split()
             if not f.endswith("_summary.csv")]
    bad = 0
    for f in files:
        old = pd.read_csv(io.StringIO(subprocess.run(["git", "show", f"HEAD:{f}"],
                          capture_output=True, text=True, check=True).stdout))
        new = pd.read_csv(f)
        keys = [k for k in KEY + ["n_traj"] if k in old and k in new]
        m = new.merge(old, on=keys, suffixes=("_new", "_old"))
        ok = len(m) == len(old) == len(new) and bool(_same(m).all())
        bad += not ok
        if not ok:
            print(f"DIFFERS: {f} ({int((~_same(m)).sum())} rows)")
    print(f"{len(files) - bad} files identical, {bad} differ"
          if files else "no committed results")


def _counts(df, method):
    d = df[df.method == method]
    cell = []
    for noise in sorted(d.noise.unique()):
        dn = d[d.noise == noise]
        cell.append(" · ".join(f"{int(dn[dn.deriv == k].recovered.sum())}/{len(dn[dn.deriv == k])}"
                               for k in ("fd", "sfd")))
    return cell


def current(folder):
    systems = SYSTEMS + [("powerdecay", ""), ("powerdecay", "_halfpowers"), ("lotka_volterra", "")]
    for system, sfx in systems:
        if not os.path.exists(f"results/{folder}/{system}{sfx}.csv"):
            continue
        df = pd.read_csv(f"results/{folder}/{system}{sfx}.csv")
        noises = sorted(df.noise.unique())
        label = {"_3traj": " (3 trajectories)", "_halfpowers": " (PySINDy with half powers)"}
        print(f"\n**{system}{label.get(sfx, '')}** (fd · sfd)\n")
        print("| method | " + " | ".join(f"{n:.0%}" for n in noises) + " |")
        print("|---" * (len(noises) + 1) + "|")
        for method in ("BACON-dt", "BACON-dt-robust", "PySINDy"):
            if (df.method == method).any():
                print(f"| {method} | " + " | ".join(_counts(df, method)) + " |")


def _total(df, method):
    d = df[df.method == method]
    wrong = (~d.equations.str.contains("—")) & (~d.recovered)
    return f"{int(d.recovered.sum())} right · {int(wrong.sum())} wrong"


def ablation(which):
    print("Each cell: runs with the right law · runs that returned a wrong law, "
          "out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.\n")
    systems = SYSTEMS + ([("powerdecay", "")] if which in ("ablation3", "ablation4", "ablation5") else [])
    systems += [("lotka_volterra", "")] if which == "ablation5" else []
    cols = [f"{s}{' 3traj' if x else ''}" for s, x in systems]
    methods = ("BACON-dt", "BACON-dt-robust") if which == "ablation1" else ("BACON-dt-robust",)
    for method in methods:
        print(f"\n**{method}**\n")
        print("| config | " + " | ".join(cols) + " |")
        print("|---" * (len(cols) + 1) + "|")
        for label, pattern in ABLATION_ROWS[which]:
            cells = []
            for s, x in systems:
                new = s == "powerdecay"   # first run in rerun 3: its ablations live in ablation3
                orig = {"powerdecay": "ablation3", "lotka_volterra": "ablation5"}.get(s, "ablation1")
                path = pattern.format(s=s, x=x, orig=orig,
                                      norate="ablation3" if new else "rerun2")
                if new and "{norate}" in pattern:
                    path = "results/ablation3/powerdecay_no_rate_seed.csv"
                cells.append(_total(pd.read_csv(path), method) if os.path.exists(path) else "missing")
            print(f"| {label} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "repro":
        repro()
    elif cmd == "verify":
        verify()
    elif cmd == "current":
        current(sys.argv[2])
    else:
        ablation(cmd)
