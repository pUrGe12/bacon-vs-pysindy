"""BACON vs PySINDy on a static law (BACON's home turf), across noise levels.

    .venv/bin/python run_static.py --system ideal_given --out results/test12 --seeds 20

Writes <system>.csv (one row per run), <system>_summary.csv and <system>_manifest.json,
like run_test.py. Laws and designs: static_laws.py.

Methods
  BACON            the original BACON.3 + BACON.1 pipeline, unchanged (bacon_static.py),
                   one frozen setting for every law and noise level: the code's defaults.
  BACON-authors    side run: the authors' own tuned args file for that dataset and noise
                   level, where one exists (given design, 0% and 1% noise only).
  PySINDy          PySINDy's STLSQ (as in Test 10: normalize_columns on, ridge alpha 0,
                   same threshold sweep and selection rule) fitting the target directly,
                   target = sum of c_k * term_k, with terms = every product of integer
                   powers of the inputs with total degree |e_1| + ... + |e_n| <= 3,
                   negative powers included (1, x, 1/x, x^2, x/y, x*y/z, 1/x^3, ...).

Scoring ("found the exact law"): both answers are solved for the target and put over a
common denominator, N(x)/D(x). Right terms = the same monomials in N and D as the true
law; then, after scaling so the true law's first denominator term matches, every
coefficient within 5%. `pred_err` = relative RMS error against the true law at 200 random
points inside the design's ranges (not on the grid).
"""
import argparse
import datetime
import json
import os
import sys
import time
import warnings
from importlib.metadata import version
from itertools import product

import numpy as np
import pandas as pd
import pysindy as ps
import sympy as sp

import bacon_static
import sindy_runner
from run_test import _git, _source_hashes
from static_laws import LAWS, SETUPS, add_noise, setup

SUCCESS_COEF_ERR = 0.05
LIB_DEGREE = 3
# Authors' tuned files (BACON repo, args/): dataset -> {noise: file}. Kepler's only file
# is BACON.6 with the answer's form written in ("j*nu**(2/3)"), so it is not used.
AUTHORS = {
    "ideal": {0.0: "args/ideal/ideal00.json", 0.01: "args/ideal/ideal10.json"},
    "ohm": {0.0: "args/ohm/ohm00.json", 0.01: "args/ohm/ohm10.json"},
    "black": {0.0: "args/black/black00.json", 0.01: "args/black/black10.json"},
    "synthetic": {0.0: "args/synthetic/synthetic.json", 0.01: "args/synthetic/synthetic.json"},
}


# --- PySINDy on a static table ------------------------------------------------------
def exponents(n, degree=LIB_DEGREE):
    """Integer exponent vectors with total degree <= degree, constant first."""
    out = [e for e in product(range(-degree, degree + 1), repeat=n)
           if sum(abs(v) for v in e) <= degree]
    return sorted(out, key=lambda e: (sum(abs(v) for v in e), [-v for v in e]))


def sindy_static(X, y, syms, normalize=True, alpha=0.0):
    exps = exponents(X.shape[1])
    theta = np.column_stack([np.prod(X ** np.array(e, dtype=float), axis=1) for e in exps])
    total_var = np.sum((y - y.mean()) ** 2)
    grid = sindy_runner.THRESHOLDS_NORM if normalize else sindy_runner.THRESHOLDS
    fits = []
    for thr in grid:
        opt = ps.STLSQ(threshold=thr, normalize_columns=normalize, alpha=alpha)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opt.fit(theta, y.reshape(-1, 1))
        coef = np.asarray(opt.coef_).ravel()
        k = int(np.count_nonzero(coef))
        if k > 0 and np.all(np.isfinite(coef)):
            fits.append((k, float(np.sum((y - theta @ coef) ** 2)), coef))
    if not fits:
        return None
    best = min(f[1] for f in fits)
    ok = [f for f in fits if f[1] <= best + sindy_runner.R2_TOL * total_var]
    coef = min(ok, key=lambda f: (f[0], f[1]))[2]
    terms = [sp.Mul(*[s ** int(p) for s, p in zip(syms, e)]) for e in exps]
    return sum((float(c) * t for c, t in zip(coef, terms) if c != 0), sp.Integer(0)), len(exps)


# --- scoring --------------------------------------------------------------------------
def _positive(expr, names):
    """BACON's symbols are only real=True; make every variable positive so sympy can
    simplify (P**2)**(1/3) to P**(2/3). All variables here are positive."""
    return expr.subs({s: sp.Symbol(s.name, positive=True) for s in expr.free_symbols
                      if s.name in names})


def _monomial_sum(expr):
    """{exponent dict: coef} if expr is a plain sum of c * x^a * y^b ..., else None."""
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, m = term.as_coeff_Mul()
        powers = m.as_powers_dict() if m != 1 else {}
        if not all(b.is_Symbol and x.is_Number for b, x in powers.items()):
            return None
        key = tuple(sorted((str(b), sp.Rational(x)) for b, x in powers.items()))
        out[key] = out.get(key, 0.0) + float(c)
    return out


def canonical(expr, names):
    """N/D over a common denominator: ({monomial: coef} for N, same for D)."""
    expr = _positive(expr, names)
    terms = _monomial_sum(expr)
    if terms is not None:
        # Fast path (every PySINDy answer): a sum of power products. Its common
        # denominator is one power product, so no symbolic cancelling is needed. Same
        # result as the general path below, which is very slow on 100-term sums.
        low = {}
        for key in terms:
            for b, x in key:
                low[b] = min(low.get(b, 0), x)
        shift = {b: -x for b, x in low.items() if x < 0}

        def mono(key):
            p = dict(key)
            for b, x in shift.items():
                p[b] = p.get(b, 0) + x
            return tuple(sorted((b, round(float(x), 6)) for b, x in p.items() if x != 0))
        num = {}
        for key, c in terms.items():
            if c != 0:
                num[mono(key)] = num.get(mono(key), 0.0) + c
        den = {tuple(sorted((b, round(float(x), 6)) for b, x in shift.items())): 1.0}
        return num, den
    e = sp.nsimplify(expr, rational=True)
    e = sp.cancel(sp.together(sp.powsimp(sp.powdenest(e, force=True), force=True)))
    num, den = sp.fraction(e)

    def terms(p):
        out = {}
        for term in sp.Add.make_args(sp.expand(p)):
            c, m = term.as_coeff_Mul()
            key = tuple(sorted((str(b), round(float(x), 6))
                               for b, x in (m.as_powers_dict().items() if m != 1 else [])))
            out[key] = out.get(key, 0.0) + float(c)
        return {k: v for k, v in out.items() if v != 0}
    return terms(num), terms(den)


def score(expr, law):
    names = set(law.inputs)
    tn, td = canonical(law.true_expr(), names)
    if expr is None:
        return {"structure_ok": False, "coef_err": np.nan}
    try:
        fn, fd = canonical(expr, names)
    except Exception:
        return {"structure_ok": False, "coef_err": np.nan}
    if set(fn) != set(tn) or set(fd) != set(td):
        return {"structure_ok": False, "coef_err": np.nan}
    ref = sorted(td)[0]
    scale = fd[ref] / td[ref]
    errs = [abs(f[k] / scale - t[k]) / abs(t[k]) for f, t in ((fn, tn), (fd, td)) for k in t]
    return {"structure_ok": True, "coef_err": float(max(errs))}


def pred_err(expr, law):
    if expr is None:
        return np.nan
    pts = law.held_out()
    true = law.evaluate(pts)
    syms = [sp.Symbol(v, positive=True) for v in law.inputs]
    try:
        f = sp.lambdify(syms, _positive(expr, set(law.inputs)), "numpy")
        with np.errstate(all="ignore"):
            got = np.asarray(f(*[pts[v] for v in law.inputs]), dtype=complex)
            got = np.broadcast_to(got, true.shape)
        if np.any(np.abs(got.imag) > 1e-9 * np.abs(got.real) + 1e-12):
            return np.inf
        rel = (got.real - true) / true
        return float(np.sqrt(np.mean(rel ** 2))) if np.all(np.isfinite(rel)) else np.inf
    except Exception:
        return np.inf


def fmt(expr):
    return "—" if expr is None else str(sp.N(expr, 4))


# --- run -------------------------------------------------------------------------------
def check_given_matches_bacon(law):
    """The given design must be exactly BACON's datasets.py data (clean)."""
    data, symbols = bacon_static.bacon_datasets.allowed_data()[law.name](0)
    names = [str(s) for s in symbols]
    ours = law.data("given")
    assert names == law.inputs + [law.target], (names, law.inputs)
    for v, col in zip(names, data):
        assert np.allclose(ours[v], np.asarray(col, dtype=float), rtol=1e-12, atol=0), v


def run(system, noises, seeds, out_dir, authors=True):
    law, design = setup(system)
    if design == "given":
        check_given_matches_bacon(law)
    clean = law.data(design)
    X = np.column_stack([clean[v] for v in law.inputs])
    psyms = [sp.Symbol(v, positive=True) for v in law.inputs]
    author_files = AUTHORS.get(law.name, {}) if (authors and design == "given") else {}

    meta = {
        "command": " ".join([os.path.basename(sys.executable)] + sys.argv),
        "started": datetime.datetime.now().isoformat(timespec="seconds"),
        "this_repo": _git(os.path.dirname(os.path.abspath(__file__))),
        "source_sha256": _source_hashes(),
        "bacon_repo": _git(bacon_static.BACON_PATH),
        "packages": {p: version(p) for p in ("numpy", "scipy", "sympy", "pandas", "pysindy")},
        "python": sys.version.split()[0],
        "law": {"name": law.name, "description": law.description, "inputs": law.inputs,
                "target": law.target, "true": law.true, "design": design,
                "rows": len(X), "levels": {v: sorted(set(map(float, clean[v])))
                                           for v in law.inputs}},
        "bacon": {"pipeline": "data_space + layer + BACON.1 (original code, unchanged)",
                  "frozen_config": bacon_static.FROZEN,
                  "code_defaults": bacon_static.CODE_DEFAULTS,
                  "time_limit_s": bacon_static.TIME_LIMIT,
                  "authors_files": {str(k): v for k, v in author_files.items()}},
        "pysindy": {"optimizer": "STLSQ", "normalize_columns": True, "ridge_alpha": 0.0,
                    "library": f"products of integer powers, total degree <= {LIB_DEGREE}",
                    "n_terms": len(exponents(len(law.inputs))),
                    "sweep_grid": list(map(float, sindy_runner.THRESHOLDS_NORM)),
                    "r2_tol": sindy_runner.R2_TOL},
        "harness": {"noise_levels": noises, "seeds": seeds,
                    "noise_model": "Gaussian on the target, sd = noise * |value|",
                    "success_coef_err": SUCCESS_COEF_ERR, "held_out_points": 200},
    }

    rows = []
    for noise in noises:
        for seed in (range(seeds) if noise > 0 else [0]):
            rng = np.random.default_rng(seed)
            y = add_noise(clean[law.target], noise, rng)
            cols = {**{v: clean[v] for v in law.inputs}, law.target: y}
            methods = ["BACON", "PySINDy"] + (["BACON-authors"] if noise in author_files else [])
            for method in methods:
                t0 = time.perf_counter()
                invariants, status = "", "ok"
                if method.startswith("BACON"):
                    cfg = (bacon_static.FROZEN if method == "BACON" else
                           bacon_static.authors_config(
                               os.path.join(bacon_static.BACON_PATH, author_files[noise])))
                    expr, inv, status = bacon_static.discover(cols, law.inputs, law.target, cfg)
                    invariants = "; ".join(f"{e.lhs} = {sp.N(e.rhs, 6)}" for e in inv)
                else:
                    res = sindy_static(X, y, psyms)
                    expr = res[0] if res else None
                    status = "ok" if res else "no law"
                secs = time.perf_counter() - t0
                s = score(expr, law)
                rows.append({
                    "system": system, "law": law.name, "design": design,
                    "noise": noise, "seed": seed, "method": method,
                    "equations": f"{law.target} = {fmt(expr)}",
                    "structure_ok": s["structure_ok"], "coef_err": s["coef_err"],
                    "recovered": bool(s["structure_ok"] and s["coef_err"] < SUCCESS_COEF_ERR),
                    "pred_err": pred_err(expr, law),
                    "status": status, "bacon_invariant": invariants, "seconds": secs,
                })
                r = rows[-1]
                print(f"noise={noise:<5} seed={seed:<2} {method:<13} recovered={str(r['recovered']):<5} "
                      f"pred_err={r['pred_err']:.3g}  {r['equations']}  [{status}] {invariants}",
                      flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, f"{system}.csv"), index=False)
    summary = (df.groupby(["noise", "method"])
                 .agg(recovered=("recovered", "mean"), structure_ok=("structure_ok", "mean"),
                      pred_err_median=("pred_err", "median"), runs=("recovered", "size"))
                 .reset_index())
    print(f"\n=== {law.description}, {design} design ({len(X)} rows) ===")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    summary.to_csv(os.path.join(out_dir, f"{system}_summary.csv"), index=False)
    meta["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    with open(os.path.join(out_dir, f"{system}_manifest.json"), "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    return df


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--system", required=True, choices=SETUPS)
    p.add_argument("--noise", type=float, nargs="+", default=[0, 0.01, 0.05, 0.1])
    p.add_argument("--seeds", type=int, default=20)
    p.add_argument("--out", default="results/test12")
    p.add_argument("--authors", default="on", choices=["on", "off"],
                   help="also run BACON with the authors' tuned args files where they exist")
    a = p.parse_args()
    run(a.system, a.noise, a.seeds, a.out, authors=a.authors == "on")
