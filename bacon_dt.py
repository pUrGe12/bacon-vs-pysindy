"""BACON with a "differential" notion, for time-series data.

Reuses BACON.1 from the original repo (~/Desktop/BACON, override with BACON_PATH)
and adds:
  * derivative columns D(x) = dx/dt, built from differences between adjacent samples,
    which BACON then pairs with state columns to look for invariants;
  * fixes for three things that break BACON.1 on ODE data:
      1. sign-blind heuristic: it picks product vs ratio from the sign of the
         correlation, so for Dx = -k*x it multiplies forever. We decide using
         the trend of the magnitudes |a|, |b| instead.
      2. constancy test `M(1-d) < v < M(1+d)` is never true when M < 0.
         We use |v - M| < d*|M|, which is the same test for M > 0.
      3. on failure it loosens its tolerances and retries without limit, so it
         always "finds" something. We cap the loosening and can return nothing.

Optional `robust=True` swaps BACON's "every point within delta of the mean" for
"the typical point within delta of the median" (median absolute deviation), and
reports the constant as the median. One bad noisy point no longer vetoes a law.
"""
import os
import sys
from statistics import fmean, median

import numpy as np
import pandas as pd
import sympy as sp
from scipy.stats import linregress

sys.path.append(os.environ.get("BACON_PATH", os.path.expanduser("~/Desktop/BACON")))
from space_of_laws.laws_methods.bacon1 import BACON_1, nu, eta  # noqa: E402


class SignedBACON1(BACON_1):
    """BACON.1 with sign-aware heuristics and bounded tolerance loosening."""

    def __init__(self, df, all_found_symbols, epsilon=0.001, delta=0.01,
                 epsilon_scale=1.5, delta_scale=1.5, c_val=0.02,
                 max_delta=0.3, robust=False, verbose=False):
        super().__init__(df, all_found_symbols, epsilon=epsilon, delta=delta,
                         epsilon_scale=epsilon_scale, delta_scale=delta_scale,
                         c_val=c_val, verbose=verbose)
        self.max_delta = max_delta
        self.robust = robust

    def center(self, data):
        return median(data) if self.robust else fmean(data)

    def _reset(self):
        self.data = [self.initial_df[c] for c in self.init_symbols]
        self.symbols = [nu, eta]

    def _one_pass(self):
        # Same loop as BACON_1.bacon_iterations, without its unbounded retry.
        self.lin_data = None
        self.update = ""
        j = 0
        while self.update != "constant" and j < 5:
            sy_start = len(self.symbols)
            self.bacon_instance(0, -1)
            if self.update in ("product", "division"):
                self.bacon_instance(1, -2)
            elif self.update not in ("linear", "constant"):
                self.bacon_instance(1, -1)
            if sy_start == len(self.symbols) or self.update == "linear":
                break
            j += 1
        return self.update in ("constant", "linear")

    def bacon_iterations(self):
        """Returns (data, expr, lin_data, delta_used), or None if nothing found."""
        while self.delta <= self.max_delta:
            self._reset()
            if self._one_pass():
                return self.data[-1], self.subs_expr(self.symbols[-1]), self.lin_data, self.delta
            if self.verbose:
                print(f"BACON-dt: nothing at delta={self.delta:.3f}, loosening")
            self.epsilon *= self.epsilon_scale
            self.delta *= self.delta_scale
        return None

    def bacon_instance(self, start, finish):
        a, b = self.data[start], self.data[finish]
        a_, b_ = self.symbols[start], self.symbols[finish]

        m, c, r, _, _ = linregress(a, b)
        self.check_constant(b_, b)
        if self.update == "constant":
            return

        if 1 - abs(r) < self.epsilon and abs(c / fmean(b)) > self.c_val:
            self.linear(a_, b_, a, b, m)
            return

        # Fix 1: compare magnitudes, so negative derivatives don't flip the choice.
        r_mag = np.corrcoef(np.abs(a), np.abs(b))[0, 1]
        if r_mag > 0:
            sy = sp.simplify(a_ / b_)
            if self.new_term(sy):
                self.division(a_, b_, a, b)
        elif r_mag < 0:
            sy = sp.simplify(a_ * b_)
            if self.new_term(sy):
                self.product(a_, b_, a, b)

    def check_constant(self, symbol, data):
        # Same early exit as BACON_1: a raw input column can't be "the law".
        if len(self.symbols) == 2:
            real_symb = self.subs_expr(symbol)
            if len(real_symb.free_symbols) == 1:
                return
        M = self.center(data)
        if self.robust:
            ok = median(abs(v - M) for v in data) < self.delta * abs(M)
        else:
            # Fix 2: sign-safe version of BACON's "every point within delta of the mean".
            ok = all(abs(v - M) < self.delta * abs(M) for v in data)
        if ok:
            self.update = "constant"
            if self.verbose:
                print(f"BACON-dt: {self.subs_expr(symbol)} is constant (delta={self.delta:.3f})")


def _fit_r2(rhs, target_sym, cols, y):
    syms = sorted(rhs.free_symbols, key=str)
    if not syms:
        pred = np.full_like(y, float(rhs))
    else:
        pred = sp.lambdify(syms, rhs, "numpy")(*[cols[s] for s in syms])
    pred = np.asarray(pred, dtype=float)
    if not np.all(np.isfinite(pred)):
        return -np.inf
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def discover(X, Xdot, symbols, robust=False, verbose=False):
    """For each state x_i, pair D(x_i) with every state x_j and let BACON look for
    an invariant. Solve the invariant for D(x_i) to get the ODE right-hand side.

    Returns a list (one per state) of dicts: rhs (sympy or None), invariant,
    constant, delta, partner, r2.
    """
    cols = {s: X[:, j] for j, s in enumerate(symbols)}
    results = []
    for i, xi in enumerate(symbols):
        D = sp.Symbol(f"D{xi}", real=True)
        y = Xdot[:, i]
        best = {"rhs": None, "invariant": None, "constant": None,
                "delta": None, "partner": None, "r2": -np.inf}
        for xj in symbols:
            df = pd.DataFrame({xj: cols[xj], D: y})
            bacon = SignedBACON1(df, [xj, D], robust=robust, verbose=verbose)
            out = bacon.bacon_iterations()
            if out is None:
                continue
            data, expr, _, delta = out
            C = bacon.center(data)
            if D not in expr.free_symbols:
                continue
            for root in sp.solve(sp.Eq(expr, C), D):
                if not root.is_real and root.has(sp.I):
                    continue
                r2 = _fit_r2(root, D, cols, y)
                # Prefer the law found at the tightest tolerance, then best fit.
                key = (-delta, r2)
                if best["rhs"] is None or key > (-best["delta"], best["r2"]):
                    best = {"rhs": root, "invariant": expr, "constant": C,
                            "delta": delta, "partner": xj, "r2": r2}
        results.append(best)
    return results
