"""BACON with a "differential" notion, for time-series data.

Reuses BACON.1 from the original repo (~/Desktop/BACON, override with BACON_PATH).
Every deviation from the original BACON.1 code is a named, switchable fix, so any
test can be rerun exactly as it was. CHANGES_FROM_BACON.md documents each one
(what the original does, what we do, why, and which test introduced it).

    B1 sign_aware        product-vs-ratio decided from the trend of |a|, |b|
    B2 signed_constancy  constancy tested on signed values (original: on |b|)
    B3 bounded_retry     stop loosening tolerances at max_delta and report nothing
    B4 intercept_vs_std  "significant intercept" judged against std(b), not mean(b)

Defaults (DEFAULT_FIXES): B1, B3, B4. B2 was on in Tests 1-4 and the first rerun, then
switched off by default: the ablation showed no benefit, so it is one less change from
the original. It stays available so those runs can be reproduced.

With every fix switched off, SignedBACON1 makes the same decisions as BACON_1
(apart from a hard stop after HARD_STOP loosenings, so a run can't hang).

`robust=True` (the BACON-dt-robust variant) replaces the constancy test with a
median-based one; see CHANGES_FROM_BACON.md, R1.
"""
import os
import sys
from statistics import fmean, median

import numpy as np
import pandas as pd
import sympy as sp
from scipy.stats import linregress

BACON_PATH = os.environ.get("BACON_PATH", os.path.expanduser("~/Desktop/BACON"))
sys.path.append(BACON_PATH)
from space_of_laws.laws_methods.bacon1 import BACON_1, nu, eta, dummies  # noqa: E402

ALL_FIXES = ("sign_aware", "signed_constancy", "bounded_retry", "intercept_vs_std")
DEFAULT_FIXES = ("sign_aware", "bounded_retry", "intercept_vs_std")
HARD_STOP = 200

# Tolerance schedule (CHANGES_FROM_BACON.md, P1). Original BACON_1 defaults:
# epsilon=0.001, delta=0.1, epsilon_scale=1.1, delta_scale=1.05, c_val=0.02.
DEFAULT_PARAMS = dict(epsilon=0.001, delta=0.01, epsilon_scale=1.5, delta_scale=1.5,
                      c_val=0.02, max_delta=0.3)


class SignedBACON1(BACON_1):
    """BACON.1 with switchable fixes for ODE data."""

    def __init__(self, df, all_found_symbols, epsilon=0.001, delta=0.01,
                 epsilon_scale=1.5, delta_scale=1.5, c_val=0.02,
                 max_delta=0.3, robust=False, fixes=DEFAULT_FIXES, verbose=False):
        super().__init__(df, all_found_symbols, epsilon=epsilon, delta=delta,
                         epsilon_scale=epsilon_scale, delta_scale=delta_scale,
                         c_val=c_val, verbose=verbose)
        unknown = set(fixes) - set(ALL_FIXES)
        if unknown:
            raise ValueError(f"unknown fixes: {unknown}")
        self.fixes = set(fixes)
        self.max_delta = max_delta
        self.robust = robust

    def center(self, data):
        return median(data) if self.robust else fmean(data)

    def _reset(self):
        self.data = [self.initial_df[c] for c in self.init_symbols]
        self.symbols = [nu, eta]

    def _one_pass(self):
        # Same loop as BACON_1.bacon_iterations, minus its retry (handled below).
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
        for _ in range(HARD_STOP):
            # B3 off: original behaviour, keep loosening until something passes.
            if "bounded_retry" in self.fixes and self.delta > self.max_delta:
                return None
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
        # B2: original passes abs(b) here.
        self.check_constant(b_, b if "signed_constancy" in self.fixes else abs(b))
        if self.update == "constant":
            return

        # B4: original compares the intercept with mean(b), ~0 for oscillating signals.
        scale = np.std(b) if "intercept_vs_std" in self.fixes else abs(fmean(b))
        if 1 - abs(r) < self.epsilon and abs(c) > self.c_val * scale:
            self.linear(a_, b_, a, b, m)
            return

        # B1: original uses the sign of r, which flips for negative derivatives.
        trend = np.corrcoef(np.abs(a), np.abs(b))[0, 1] if "sign_aware" in self.fixes else r
        if trend > 0:
            sy = sp.simplify(a_ / b_)
            if self.new_term(sy):
                self.division(a_, b_, a, b)
        elif trend < 0:
            sy = sp.simplify(a_ * b_)
            if self.new_term(sy):
                self.product(a_, b_, a, b)

    def check_constant(self, symbol, data):
        # Same early exit as BACON_1: a raw input column can't be "the law".
        if len(self.symbols) == 2:
            real_symb = self.subs_expr(symbol)
            if len(real_symb.free_symbols) == 1 and str(real_symb) not in dummies:
                return
        M = self.center(data)
        if self.robust:
            # R1: the typical point (median absolute deviation) within delta of the median.
            ok = median(abs(v - M) for v in data) < self.delta * abs(M)
        elif "signed_constancy" in self.fixes:
            ok = all(abs(v - M) < self.delta * abs(M) for v in data)
        else:
            ok = all(M * (1 - self.delta) < v < M * (1 + self.delta) for v in data)  # original
        if ok:
            self.update = "constant"
            if self.verbose:
                print(f"BACON-dt: {self.subs_expr(symbol)} is constant (delta={self.delta:.3f})")


def _fit_r2(rhs, cols, y):
    syms = sorted(rhs.free_symbols, key=str)
    if not syms:
        pred = np.full_like(y, float(rhs))
    else:
        pred = sp.lambdify(syms, rhs, "numpy")(*[cols[s] for s in syms])
    pred = np.asarray(pred, dtype=float)
    if not np.all(np.isfinite(pred)):
        return -np.inf
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def discover(X, Xdot, symbols, robust=True, fixes=DEFAULT_FIXES, verbose=False):
    """Time-series adaptation (CHANGES_FROM_BACON.md, A2/A3): for each state x_i, pair
    D(x_i) with every state x_j and let BACON.1 look for an invariant. Solve the
    invariant for D(x_i) to get the ODE right-hand side.

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
            bacon = SignedBACON1(df, [xj, D], robust=robust, fixes=fixes,
                                 verbose=verbose, **DEFAULT_PARAMS)
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
                r2 = _fit_r2(root, cols, y)
                # Prefer the law found at the tightest tolerance, then best fit.
                key = (-delta, r2)
                if best["rhs"] is None or key > (-best["delta"], best["r2"]):
                    best = {"rhs": root, "invariant": expr, "constant": C,
                            "delta": delta, "partner": xj, "r2": r2}
        results.append(best)
    return results
