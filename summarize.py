"""Build the README tables from the result files.

    .venv/bin/python summarize.py repro               # do reproductions match the original runs?
    .venv/bin/python summarize.py current rerun2      # one results folder, every system (rerun6b = Test 9b)
    .venv/bin/python summarize.py ablation1           # each fix off, from rerun 1 (B1-B4)
    .venv/bin/python summarize.py ablation2           # each fix off, from rerun 2 (B1, B3, B4)
    .venv/bin/python summarize.py ablation3           # each fix off, from rerun 3 (B1, B3-B5)
    .venv/bin/python summarize.py ablation4           # each fix off, from rerun 4 (B1, B3-B6)
    .venv/bin/python summarize.py ablation5           # each fix off, from rerun 5 (B1, B3, B5-B7)
    .venv/bin/python summarize.py verify              # regenerated results == committed ones?
    .venv/bin/python summarize.py tables              # rewrite every results table in README.md
                                                      # and CHANGES_FROM_BACON.md from the CSVs

Tables never pool runs: one column per noise level (with its run count), smoothed and
raw derivatives in separate tables.
"""
import io
import os
import re
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


# ---------------------------------------------------------------------------------------
# Results tables: one column per noise level, smoothed and raw derivatives separate.
NOISES = [(0.0, "no noise"), (0.01, "1% noise"), (0.05, "5% noise"), (0.1, "10% noise")]
DERIVS = [("sfd", "Smoothed derivatives"), ("fd", "Raw derivatives")]
CAPTION = {
    "exact": "runs that found the exact law",
    "sim": "runs whose learned equations simulate within 10% of the true trajectory",
    "rightwrong": "runs that found the exact law (right) or another law (wrong)",
    "outcome": "exact law / right terms but a coefficient more than 5% off / "
               "wrong terms / no law",
}
TITLE = {
    "decay": "Exponential decay, `x' = -0.5x`",
    "logistic": "Logistic growth, `x' = x - 0.1x²`",
    "oscillator": "Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`",
    "oscillator_3traj": "Harmonic oscillator, three trajectories pooled",
    "powerdecay": "Power decay, `x' = -0.5x^1.5` (PySINDy: plain polynomial terms)",
    "powerdecay_halfpowers": "Power decay, `x' = -0.5x^1.5` (PySINDy also given `x^0.5`, "
                             "`x^1.5`, `x^2.5`)",
    "lotka_volterra": "Lotka–Volterra, `x' = x - 0.1xy`, `y' = -1.5y + 0.075xy`",
}
LABEL = {"BACON-dt": "BACON-dt (strict)", "BACON-dt-robust": "BACON-dt-robust",
         "PySINDy": "PySINDy"}


def _cell(d, kind):
    if kind == "exact":
        return f"{int(d.recovered.sum())} of {len(d)}"
    if kind == "sim":
        return f"{int((d.sim_err < 0.1).sum())} of {len(d)}"
    none = d.equations.str.contains("—")
    if kind == "rightwrong":
        return f"{int(d.recovered.sum())} right, {int((~none & ~d.recovered).sum())} wrong"
    if kind == "outcome":
        close = d.structure_ok & ~d.recovered
        wrong = ~none & ~d.structure_ok
        return (f"{int(d.recovered.sum())} / {int(close.sum())} / "
                f"{int(wrong.sum())} / {int(none.sum())}")
    raise ValueError(kind)


def _load(path, method):
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path)
    return d[d.method == method]


def table(rows, deriv, kind, first="method"):
    """rows: [(label, DataFrame or None)] -> one markdown table for one derivative type."""
    have = [d for _, d in rows if d is not None]
    cols = [(n, name) for n, name in NOISES
            if any(((d.noise == n) & (d.deriv == deriv)).any() for d in have)]
    head = []
    for n, name in cols:
        k = len(have[0][(have[0].noise == n) & (have[0].deriv == deriv)])
        head.append(f"{name} ({k} run{'s' if k != 1 else ''})")
    dname = dict(DERIVS)[deriv]
    out = [f"{dname}: {CAPTION[kind]}", "",
           f"| {first} | " + " | ".join(head) + " |", "|---" * (len(cols) + 1) + "|"]
    for label, d in rows:
        if d is None:
            out.append(f"| {label} | " + " | ".join("not run" for _ in cols) + " |")
            continue
        cells = [_cell(d[(d.noise == n) & (d.deriv == deriv)], kind) for n, _ in cols]
        out.append(f"| {label} | " + " | ".join(cells) + " |")
    return "\n".join(out) + "\n"


def system_tables(title, rows, kinds=("exact",)):
    out = [f"**{title}**\n"] if title else []
    for deriv, _ in DERIVS:
        for kind in kinds:
            out.append(table(rows, deriv, kind))
    return "\n".join(out)


def folder_tables(folder, systems, methods=("BACON-dt", "BACON-dt-robust", "PySINDy")):
    out = []
    for f in systems:
        rows = [(LABEL[m], d) for m in methods
                for d in [_load(f"results/{folder}/{f}.csv", m)] if d is not None and len(d)]
        out.append(system_tables(TITLE[f], rows))
    return "\n".join(out)


OLD4 = ["decay", "logistic", "oscillator", "oscillator_3traj"]
ALL7 = OLD4 + ["powerdecay", "powerdecay_halfpowers", "lotka_volterra"]
T9 = [("BACON-dt-robust (same in every test)", "rerun5", "BACON-dt-robust"),
      ("PySINDy, Test 8 setup (raw columns, penalty 0.05)", "rerun5", "PySINDy"),
      ("PySINDy, 9a (normalized, penalty 0.05)", "rerun6", "PySINDy"),
      ("PySINDy, 9b (normalized, penalty 0)", "rerun6b", "PySINDy")]


def _versus(folder, systems):
    """BACON vs PySINDy, one table per derivative type, two rows per system."""
    out = []
    for deriv, dname in DERIVS:
        rows = []
        for f in systems:
            for m in ("BACON-dt-robust", "PySINDy"):
                rows.append((f"{TITLE[f].split(',')[0].split(' (')[0]}"
                             f"{' (½ powers)' if 'half' in f else ''}"
                             f"{' (3 traj.)' if '3traj' in f else ''} — {LABEL[m]}",
                             _load(f"results/{folder}/{f}.csv", m)))
        out.append(table(rows, deriv, "exact", first="system — method"))
    return "\n".join(out)


def _ablation_tables(which):
    systems = SYSTEMS + ([("powerdecay", "")] if which in ("ablation3", "ablation4", "ablation5") else [])
    systems += [("lotka_volterra", "")] if which == "ablation5" else []
    methods = ("BACON-dt", "BACON-dt-robust") if which == "ablation1" else ("BACON-dt-robust",)
    out = []
    for s, x in systems:
        new = s == "powerdecay"
        orig = {"powerdecay": "ablation3", "lotka_volterra": "ablation5"}.get(s, "ablation1")
        for method in methods:
            rows = []
            for label, pattern in ABLATION_ROWS[which]:
                path = pattern.format(s=s, x=x, orig=orig, norate="ablation3" if new else "rerun2")
                if new and "{norate}" in pattern:
                    path = "results/ablation3/powerdecay_no_rate_seed.csv"
                rows.append((label, _load(path, method)))
            title = TITLE[s + x] + (f" — {LABEL[method]}" if len(methods) > 1 else "")
            out.append(f"**{title}**\n")
            for deriv, _ in DERIVS:
                out.append(table(rows, deriv, "rightwrong", first="config"))
    return "\n".join(out)


def _b6_outcomes():
    out = []
    for f in ["decay", "logistic", "oscillator", "oscillator_3traj", "powerdecay"]:
        rows = [("B6 off (rerun 3)", _load(f"results/rerun3/{f}.csv", "BACON-dt-robust")),
                ("B6 on (rerun 4)", _load(f"results/rerun4/{f}.csv", "BACON-dt-robust"))]
        out.append(f"**{TITLE[f].split(' (')[0]}**\n")
        for deriv, _ in DERIVS:
            out.append(table(rows, deriv, "outcome", first="config"))
    return "\n".join(out)



# ---------------------------------------------------------------------------------------
# Test 10: ODEBench
def _odebench_setups():
    import json
    import sympy as sp
    from systems import ODEBENCH, ODEBENCH_FILE
    eqs = {e["id"]: e for e in json.load(open(ODEBENCH_FILE))["equations"]}
    out = []
    for name, make in ODEBENCH.items():
        sysm = make()
        eid, start = int(name[3:5]), int(name[-1])
        in_lib = True
        for r in sysm.true_rhs:
            try:
                in_lib &= sp.Poly(sp.expand(r), *sysm.symbols).total_degree() <= 3
            except sp.PolynomialError:
                in_lib = False
        out.append(dict(name=name, id=eid, start=start, dim=len(sysm.names), in_lib=in_lib,
                        label=f"#{eid} {eqs[eid]['eq_description']}, start {start}"))
    return sorted(out, key=lambda r: (r["id"], r["start"]))


def _odebench_counts(setups):
    """{(name, method, deriv, noise): (exact, runs)} from results/test10."""
    out = {}
    for st in setups:
        path = f"results/test10/{st['name']}.csv"
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path)
        for (m, k, n), g in d.groupby(["method", "deriv", "noise"]):
            out[(st["name"], m, k, n)] = (int(g.recovered.sum()), len(g))
    return out


def _odebench_table(rows, c, deriv):
    head = ["setup", "true law in PySINDy's terms?"]
    for n, nname in NOISES:
        runs = 1 if n == 0 else 20
        head += [f"{nname}: BACON ({runs} run{'s' if runs > 1 else ''})",
                 f"{nname}: PySINDy ({runs} run{'s' if runs > 1 else ''})"]
    out = ["| " + " | ".join(head) + " |", "|---" * len(head) + "|"]
    for st in rows:
        cells = [st["label"], "yes" if st["in_lib"] else "no"]
        for n, _ in NOISES:
            for m in ("BACON-dt-robust", "PySINDy"):
                v = c.get((st["name"], m, deriv, n))
                cells.append("not run" if v is None else f"{v[0]} of {v[1]}")
        out.append("| " + " | ".join(cells) + " |")
    out.append("")
    return out


def _odebench_full():
    setups = _odebench_setups()
    c = _odebench_counts(setups)
    out = []
    for deriv, dname in DERIVS:
        out.append(f"## {dname}\n")
        for dim in (1, 2, 3, 4):
            rows = [st for st in setups if st["dim"] == dim]
            out.append(f"### {dname}: systems with {dim} variable{'s' if dim > 1 else ''}"
                       f" ({len(rows)} setups)\n")
            out.append(f"{dname}: runs that found the exact law. Columns per noise level: "
                       "BACON, then PySINDy.\n")
            out += _odebench_table(rows, c, deriv)
    return "\n".join(out)


def _odebench_found():
    """Only the setups where at least one method found the exact law at least once."""
    setups = _odebench_setups()
    c = _odebench_counts(setups)
    out = []
    for deriv, dname in DERIVS:
        rows = [st for st in setups
                if any(c.get((st["name"], m, deriv, n), (0, 0))[0]
                       for m in ("BACON-dt-robust", "PySINDy") for n, _ in NOISES)]
        out.append(f"{dname}: runs that found the exact law, for the {len(rows)} setups where "
                   "either method found it at least once. Columns per noise level: BACON, "
                   "then PySINDy.\n")
        out += _odebench_table(rows, c, deriv)
    return "\n".join(out)


def _odebench_summary():
    setups = _odebench_setups()
    c = _odebench_counts(setups)
    groups = [("All setups", setups),
              ("Setups whose true law is inside PySINDy's terms (polynomial, degree ≤ 3)",
               [s for s in setups if s["in_lib"]]),
              ("Setups whose true law is outside PySINDy's terms",
               [s for s in setups if not s["in_lib"]])]
    out = []
    for deriv, dname in DERIVS:
        for gname, rows in groups:
            done = [s for s in rows if (s["name"], "PySINDy", deriv, 0.0) in c]
            out.append(f"{dname} — {gname} ({len(done)} of {len(rows)} setups run)\n")
            out.append("| | " + " | ".join(n for _, n in NOISES) + " |")
            out.append("|---" * (len(NOISES) + 1) + "|")
            lines = {k: [] for k in ("b_half", "p_half", "b_better", "p_better", "tie", "neither")}
            for n, _ in NOISES:
                cnt = dict.fromkeys(lines, 0)
                for s in done:
                    b, nb = c[(s["name"], "BACON-dt-robust", deriv, n)]
                    p_, np_ = c[(s["name"], "PySINDy", deriv, n)]
                    cnt["b_half"] += b > nb / 2
                    cnt["p_half"] += p_ > np_ / 2
                    if b == p_ == 0:
                        cnt["neither"] += 1
                    elif b > p_:
                        cnt["b_better"] += 1
                    elif p_ > b:
                        cnt["p_better"] += 1
                    else:
                        cnt["tie"] += 1
                for k in lines:
                    lines[k].append(f"{cnt[k]} of {len(done)}")
            labels = {"b_half": "BACON finds the exact law in more than half the runs",
                      "p_half": "PySINDy finds the exact law in more than half the runs",
                      "b_better": "BACON finds it in more runs than PySINDy",
                      "p_better": "PySINDy finds it in more runs than BACON",
                      "tie": "same number of runs, more than zero",
                      "neither": "neither finds it in any run"}
            for k, lab in labels.items():
                out.append(f"| {lab} | " + " | ".join(lines[k]) + " |")
            out.append("")
    return "\n".join(out)


BLOCKS = {
    "test1": lambda: system_tables(None, [
        ("BACON-dt (strict)", _load("results/decay.csv", "BACON-dt")),
        ("PySINDy (fixed threshold 0.1)", _load("results/decay.csv", "PySINDy"))]),
    "test2": lambda: system_tables(None, [(LABEL[m], _load("results/logistic.csv", m))
                                          for m in ("BACON-dt", "BACON-dt-robust", "PySINDy")]),
    "test3": lambda: system_tables(None, [(LABEL[m], _load("results/oscillator.csv", m))
                                          for m in ("BACON-dt", "BACON-dt-robust", "PySINDy")],
                                   kinds=("exact", "sim")),
    "test4": lambda: system_tables("One trajectory (same data as Test 3)", [
        ("BACON-dt-robust, Test 3 (before B4)", _load("results/oscillator.csv", "BACON-dt-robust")),
        ("BACON-dt (strict), Test 4", _load("results/oscillator_test4_1traj.csv", "BACON-dt")),
        ("BACON-dt-robust, Test 4", _load("results/oscillator_test4_1traj.csv", "BACON-dt-robust")),
        ("PySINDy", _load("results/oscillator_test4_1traj.csv", "PySINDy"))],
        kinds=("exact", "sim")) + "\n" + system_tables("Three trajectories pooled", [
        (LABEL[m], _load("results/oscillator_test4_3traj.csv", m))
        for m in ("BACON-dt", "BACON-dt-robust", "PySINDy")], kinds=("exact", "sim")),
    "rerun1": lambda: folder_tables("rerun1", OLD4),
    "rerun2": lambda: folder_tables("rerun2", OLD4),
    "test5": lambda: folder_tables("rerun3", OLD4),
    "test6": lambda: folder_tables("rerun3", ["powerdecay", "powerdecay_halfpowers"]),
    "test7": lambda: folder_tables("rerun4", OLD4 + ["powerdecay"]),
    "test8": lambda: folder_tables("rerun5", ["lotka_volterra"]),
    "test9": lambda: "\n".join(system_tables(TITLE[f], [
        (lab, _load(f"results/{folder}/{f}.csv", m)) for lab, folder, m in T9]) for f in ALL7),
    "test9_versus": lambda: _versus("rerun6b", ALL7),
    "ablation1": lambda: _ablation_tables("ablation1"),
    "ablation2": lambda: _ablation_tables("ablation2"),
    "ablation3": lambda: _ablation_tables("ablation3"),
    "ablation4": lambda: _ablation_tables("ablation4"),
    "ablation4_b6": _b6_outcomes,
    "ablation5": lambda: _ablation_tables("ablation5"),
    "odebench_full": _odebench_full,
    "test10_summary": _odebench_summary,
    "test10_found": _odebench_found,
}


def tables():
    pat = re.compile(r"(<!-- tables:(\w+) -->\n).*?(<!-- /tables -->)", re.S)
    for doc in ("README.md", "CHANGES_FROM_BACON.md", "ODEBENCH_RESULTS.md"):
        if not os.path.exists(doc):
            continue
        text = open(doc).read()
        names = pat.findall(text)
        text = pat.sub(lambda m: m.group(1) + "\n" + BLOCKS[m.group(2)]() + "\n" + m.group(3), text)
        open(doc, "w").write(text)
        print(f"{doc}: rewrote {len(names)} table blocks")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "repro":
        repro()
    elif cmd == "verify":
        verify()
    elif cmd == "current":
        current(sys.argv[2])
    elif cmd == "tables":
        tables()
    else:
        ablation(cmd)
