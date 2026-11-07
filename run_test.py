"""Run BACON-dt and PySINDy on one system across noise levels and derivative methods.

    .venv/bin/python run_test.py --system decay

Every run writes <name>.csv (one row per run), <name>_summary.csv, and
<name>_manifest.json (full config, BACON fixes, git commits, package versions,
command line), so any result can be traced to the exact code that made it.
README.md lists the command that reproduces each test.
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time
from importlib.metadata import version

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
VARIANTS = {"strict": "BACON-dt", "robust": "BACON-dt-robust"}


def fmt(expr):
    return "—" if expr is None else str(sp.N(expr, 4))


def _git(path):
    try:
        head = subprocess.run(["git", "-C", path, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", path, "status", "--porcelain"],
                                    capture_output=True, text=True).stdout.strip())
        return {"path": path, "commit": head, "uncommitted_changes": dirty}
    except Exception as e:
        return {"path": path, "commit": None, "error": str(e)}


def _source_hashes():
    # Identifies the exact code even when the run was made from uncommitted changes.
    here = os.path.dirname(os.path.abspath(__file__))
    files = sorted(f for f in os.listdir(here) if f.endswith(".py"))
    files += [os.path.join(bacon_dt.BACON_PATH, "space_of_laws/laws_methods/bacon1.py")]
    out = {}
    for f in files:
        with open(os.path.join(here, f), "rb") as fh:
            out[os.path.relpath(os.path.join(here, f), here)] = hashlib.sha256(fh.read()).hexdigest()[:16]
    return out


def manifest(system, cfg):
    return {
        "command": " ".join([os.path.basename(sys.executable)] + sys.argv),
        "started": datetime.datetime.now().isoformat(timespec="seconds"),
        "this_repo": _git(os.path.dirname(os.path.abspath(__file__))),
        "source_sha256": _source_hashes(),
        "bacon_repo": _git(bacon_dt.BACON_PATH),
        "packages": {p: version(p) for p in ("numpy", "scipy", "sympy", "pandas", "pysindy")},
        "python": sys.version.split()[0],
        "system": {"name": system.name, "description": system.description,
                   "true_rhs": [str(r) for r in system.true_rhs],
                   "initial_conditions": system.initial_conditions[:cfg["n_traj"]],
                   "t_start": float(system.t[0]), "t_end": float(system.t[-1]),
                   "samples_per_trajectory": len(system.t)},
        "bacon": {"variants": [VARIANTS[v] for v in cfg["bacon_variants"]],
                  "fixes_on": list(cfg["bacon_fixes"]),
                  "fixes_off": [f for f in bacon_dt.ALL_FIXES if f not in cfg["bacon_fixes"]],
                  "params": bacon_dt.DEFAULT_PARAMS},
        "pysindy": {"poly_degree": 3, "optimizer": "STLSQ",
                    "threshold": cfg["sindy_threshold"],
                    "library_kind": cfg["sindy_library"],
                    "sweep_grid": list(map(float, sindy_runner.THRESHOLDS)),
                    "r2_tol": sindy_runner.R2_TOL},
        "harness": {"noise_levels": cfg["noises"], "seeds": cfg["seeds"],
                    "n_traj": cfg["n_traj"], "derivatives": list(METHODS),
                    "edge_trim": EDGE, "success_coef_err": SUCCESS_COEF_ERR,
                    "sim_solver": cfg["sim_solver"]},
    }


def run(system_name, noises, seeds, out_dir, n_traj=1, tag="",
        bacon_fixes=bacon_dt.DEFAULT_FIXES, bacon_variants=("robust",),
        sindy_threshold="sweep", sim_solver="RK45", sindy_library="poly"):
    system = SYSTEMS[system_name]()
    cfg = dict(noises=noises, seeds=seeds, n_traj=n_traj, bacon_fixes=tuple(bacon_fixes),
               bacon_variants=tuple(bacon_variants), sindy_threshold=sindy_threshold,
               sindy_library=sindy_library,
               sim_solver=sim_solver)
    meta = manifest(system, cfg)

    x0s = system.initial_conditions[:n_traj]
    if len(x0s) < n_traj:
        raise ValueError(f"{system_name} only defines {len(x0s)} initial conditions")
    clean = [system.simulate(x0) for x0 in x0s]
    t, syms = system.t, system.symbols
    sl = slice(EDGE, -EDGE)
    methods = [VARIANTS[v] for v in bacon_variants] + ["PySINDy"]
    fixes_str = ",".join(bacon_fixes) or "none"
    rows = []

    for noise in noises:
        for seed in (range(seeds) if noise > 0 else [0]):
            rng = np.random.default_rng(seed)
            noisy = [add_noise(Xc, noise, rng) for Xc in clean]
            for deriv in METHODS:
                # Derivatives per trajectory, then pool the rows. Both methods get the same pool.
                Xf = np.vstack([X[sl] for X in noisy])
                Xdf = np.vstack([differentiate(X, t, deriv)[sl] for X in noisy])
                dt = t[1] - t[0]

                for method in methods:
                    t0 = time.perf_counter()
                    if method.startswith("BACON"):
                        found = bacon_dt.discover(Xf, Xdf, syms, fixes=bacon_fixes,
                                                  robust=method.endswith("robust"))
                        rhs = [r["rhs"] for r in found]
                        extra = "; ".join(
                            f"{r['invariant']} = {r['constant']:.4g} (delta {r['delta']:.3f})"
                            if r["rhs"] is not None else "no invariant" for r in found)
                    else:
                        rhs = sindy_runner.discover(Xf, Xdf, dt, syms, threshold=sindy_threshold,
                                                    library=sindy_library)
                        extra = ""
                    secs = time.perf_counter() - t0

                    scores = [score_equation(r, tr) for r, tr in zip(rhs, system.true_rhs)]
                    structure_ok = all(s["structure_ok"] for s in scores)
                    coef_err = max(s["coef_err"] for s in scores)
                    rows.append({
                        "system": system_name, "n_traj": n_traj, "noise": noise, "seed": seed,
                        "deriv": deriv, "method": method,
                        "equations": " | ".join(f"D{s} = {fmt(r)}" for s, r in zip(syms, rhs)),
                        "structure_ok": structure_ok,
                        "coef_err": coef_err,
                        "sim_err": float(np.mean([
                            simulation_error(rhs, syms, x0, t, Xc, solver=sim_solver)
                            for x0, Xc in zip(x0s, clean)])),
                        "recovered": bool(structure_ok and coef_err < SUCCESS_COEF_ERR),
                        "seconds": secs,
                        "bacon_invariant": extra,
                        "bacon_fixes": fixes_str if method.startswith("BACON") else "",
                        "sindy_threshold": sindy_threshold if method == "PySINDy" else "",
                        "sindy_library": sindy_library if method == "PySINDy" else "",
                    })
                    r = rows[-1]
                    print(f"noise={noise:<5} seed={seed} {deriv:<3} {method:<15} "
                          f"recovered={str(r['recovered']):<5} coef_err={r['coef_err']:.3f} "
                          f"sim_err={r['sim_err']:.3f}  {r['equations']}  {extra}")

    df = pd.DataFrame(rows)
    os.makedirs(out_dir, exist_ok=True)
    name = system_name + (f"_{tag}" if tag else "")
    df.to_csv(os.path.join(out_dir, f"{name}.csv"), index=False)

    summary = (df.groupby(["noise", "deriv", "method"])
                 .agg(recovered=("recovered", "mean"),
                      coef_err_median=("coef_err", "median"),
                      sim_err_median=("sim_err", "median"),
                      runs=("recovered", "size"))
                 .reset_index())
    print(f"\n=== {system.description}, {n_traj} trajectory(ies), BACON fixes: {fixes_str} ===")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    summary.to_csv(os.path.join(out_dir, f"{name}_summary.csv"), index=False)

    meta["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    with open(os.path.join(out_dir, f"{name}_manifest.json"), "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    return df, summary


def _threshold(v):
    return v if v == "sweep" else float(v)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--system", default="decay", choices=list(SYSTEMS))
    p.add_argument("--noise", type=float, nargs="+", default=[0, 0.01, 0.05, 0.1])
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--out", default="results")
    p.add_argument("--tag", default="", help="suffix for the results file names")
    p.add_argument("--n-traj", type=int, default=1,
                   help="number of trajectories (different initial conditions) to pool")
    p.add_argument("--bacon-fixes", nargs="+", default=list(bacon_dt.DEFAULT_FIXES),
                   choices=list(bacon_dt.ALL_FIXES) + ["none"],
                   help="BACON fixes to switch on (default: sign_aware bounded_retry "
                        "intercept_vs_std). 'none' = original BACON.1")
    p.add_argument("--bacon-variants", nargs="+", default=["robust"], choices=list(VARIANTS),
                   help="robust (default) and/or strict (Tests 1-4 ran both)")
    p.add_argument("--sindy-threshold", type=_threshold, default="sweep",
                   help="'sweep' (default) or a fixed STLSQ threshold, e.g. 0.1")
    p.add_argument("--sim-solver", default="RK45", choices=["RK45", "LSODA"])
    p.add_argument("--sindy-library", default="poly", choices=["poly", "poly+halfpowers"],
                   help="PySINDy term library (default: polynomial degree 3)")
    a = p.parse_args()
    fixes = [] if a.bacon_fixes == ["none"] else a.bacon_fixes
    run(a.system, a.noise, a.seeds, a.out, a.n_traj, a.tag,
        fixes, a.bacon_variants, a.sindy_threshold, a.sim_solver, a.sindy_library)
