"""Check that SignedBACON1 with every fix off makes the same decisions as the
original BACON_1, using BACON's original default tolerances.

    .venv/bin/python check_original.py
"""
import sys
import warnings

import numpy as np
import pandas as pd
import sympy as sp

from bacon_dt import BACON_1, SignedBACON1
from derivatives import differentiate
from systems import decay, logistic, oscillator

warnings.filterwarnings("ignore")
sys.setrecursionlimit(5000)  # the original retries by recursion
ORIGINAL_DEFAULTS = dict(epsilon=0.001, delta=0.1, epsilon_scale=1.1, delta_scale=1.05, c_val=0.02)


def cases():
    n = np.arange(1, 7, dtype=float)
    P, D = sp.Symbol("P"), sp.Symbol("D")
    yield "kepler (BACON's own dataset)", pd.DataFrame({P: n**3, D: n**2})
    for make in (decay, logistic, oscillator):
        s = make()
        X = s.simulate()
        Xd = differentiate(X, s.t, "fd")
        target = sp.Symbol(f"D{s.symbols[0]}")
        for j, xj in enumerate(s.symbols):
            yield f"{s.name}: {xj} vs {target}", pd.DataFrame({xj: X[2:-2, j], target: Xd[2:-2, 0]})


def describe(out):
    return None if out is None else (str(out[1]), round(float(np.mean(out[0])), 6))


ok = True
for name, df in cases():
    syms = list(df.columns)
    orig = describe(BACON_1(df, syms, **ORIGINAL_DEFAULTS).bacon_iterations())
    ours = describe(SignedBACON1(df, syms, fixes=(), robust=False, max_delta=1e9,
                                 **ORIGINAL_DEFAULTS).bacon_iterations())
    same = orig == ours
    ok &= same
    print(f"{'SAME' if same else 'DIFF'}  {name:<28} original: {orig}   fixes off: {ours}")
print("\nAll match." if ok else "\nMISMATCH: fixes-off mode is not the original.")
sys.exit(0 if ok else 1)
