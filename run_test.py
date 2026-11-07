"""Run BACON-dt and PySINDy on one system across noise levels and derivative methods.

    .venv/bin/python run_test.py --system decay
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
import sympy as sp

import bacon_dt
import sindy_runner
from derivatives import differentiate, METHODS
from scoring import score_equation, simulation_error
from systems import SYSTEMS, add_noise

SUCCESS_COEF_ERR = 0.05   # "recovered" = right terms AND coefficients within 5%
EDGE = 2                  # drop this many samples at each end (finite-difference edges)


def fmt(expr):
    return "—" if expr is None else str(sp.N(expr, 4))


def run(system_name, noises, seeds, out_dir):
    system = SYSTEMS[system_name]()
    X_clean = system.simulate()
    t, syms = system.t, system.symbols
    rows = []

    for noise in noises:
        for seed in (range(seeds) if noise > 0 else [0]):
            X = add_noise(X_clean, noise, np.random.default_rng(seed))
            for deriv in METHODS:
                Xdot = differentiate(X, t, deriv)
                sl = slice(EDGE, -EDGE)
                Xf, Xdf, tf = X[sl], Xdot[sl], t[sl]

                for method in ("BACON-dt", "BACON-dt-robust", "PySINDy"):
                    t0 = time.perf_counter()
                    if method.startswith("BACON"):
                        found = bacon_dt.discover(Xf, Xdf, syms, robust=method.endswith("robust"))
                        rhs = [r["rhs"] for r in found]
                        extra = "; ".join(
                            f"{r['invariant']} = {r['constant']:.4g} (delta {r['delta']:.3f})"
                            if r["rhs"] is not None else "no invariant" for r in found)
                    else:
                        rhs = sindy_runner.discover(Xf, Xdf, tf, syms)
                        extra = ""
                    secs = time.perf_counter() - t0

                    scores = [score_equation(r, tr) for r, tr in zip(rhs, system.true_rhs)]
                    structure_ok = all(s["structure_ok"] for s in scores)
                    coef_err = max(s["coef_err"] for s in scores)
                    rows.append({
                        "system": system_name, "noise": noise, "seed": seed,
                        "deriv": deriv, "method": method,
                        "equations": " | ".join(f"D{s} = {fmt(r)}" for s, r in zip(syms, rhs)),
                        "structure_ok": structure_ok,
                        "coef_err": coef_err,
                        "sim_err": simulation_error(rhs, syms, system.x0, t, X_clean),
                        "recovered": bool(structure_ok and coef_err < SUCCESS_COEF_ERR),
                        "seconds": secs,
                        "bacon_invariant": extra,
                    })
                    r = rows[-1]
                    print(f"noise={noise:<5} seed={seed} {deriv:<3} {method:<15} "
                          f"recovered={str(r['recovered']):<5} coef_err={r['coef_err']:.3f} "
                          f"sim_err={r['sim_err']:.3f}  {r['equations']}  {extra}")

    df = pd.DataFrame(rows)
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, f"{system_name}.csv"), index=False)

    summary = (df.groupby(["noise", "deriv", "method"])
                 .agg(recovered=("recovered", "mean"),
                      coef_err_median=("coef_err", "median"),
                      sim_err_median=("sim_err", "median"),
                      runs=("recovered", "size"))
                 .reset_index())
    print(f"\n=== {system.description} ===")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    summary.to_csv(os.path.join(out_dir, f"{system_name}_summary.csv"), index=False)
    return df, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--system", default="decay", choices=list(SYSTEMS))
    p.add_argument("--noise", type=float, nargs="+", default=[0, 0.01, 0.05, 0.1])
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--out", default="results")
    a = p.parse_args()
    run(a.system, a.noise, a.seeds, a.out)
