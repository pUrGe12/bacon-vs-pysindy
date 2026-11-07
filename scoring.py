"""Compare a learned ODE right-hand side with the ground truth."""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def _clean_powers(expr):
    """Turn float exponents like x**1.0 / x**1.5 into exact rationals."""
    return expr.replace(
        lambda e: e.is_Pow and e.exp.is_Float,
        lambda e: e.base ** sp.nsimplify(e.exp, tolerance=1e-6, rational=True))


def terms(expr):
    """{monomial: coefficient} of an expanded expression.

    Nested powers are merged as if variables were positive, so BACON's sqrt(x**3)
    and the true x**(3/2) count as the same term (CHANGES_FROM_BACON.md, H6).
    """
    out = {}
    expr = sp.powdenest(_clean_powers(sp.sympify(expr)), force=True)
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, mono = term.as_coeff_Mul()
        if abs(float(coeff)) > 1e-10:
            out[mono] = out.get(mono, 0.0) + float(coeff)
    return out


def score_equation(learned, true):
    """structure_ok: same set of terms. coef_err: relative L2 error over all terms."""
    if learned is None:
        return {"structure_ok": False, "coef_err": np.nan}
    lt, tt = terms(learned), terms(true)
    keys = set(lt) | set(tt)
    diff = np.sqrt(sum((lt.get(k, 0) - tt.get(k, 0)) ** 2 for k in keys))
    norm = np.sqrt(sum(v ** 2 for v in tt.values()))
    return {"structure_ok": set(lt) == set(tt), "coef_err": diff / norm}


def simulation_error(rhs_list, symbols, x0, t, X_clean, solver="RK45"):
    """Integrate the learned ODE from the true x0; relative L2 error vs clean data.

    solver: Tests 1-2 used LSODA, Tests 3+ use RK45 (LSODA hung when the
    blow-up guard raised inside its Fortran callback).
    """
    if any(r is None for r in rhs_list):
        return np.nan
    f = sp.lambdify(symbols, rhs_list, "numpy")
    limit = 1e3 * np.abs(X_clean).max()
    calls = [0]

    def rhs(_t, s):
        # Bad learned equations can be stiff, complex-valued or explode;
        # give up on those instead of letting the solver crawl.
        calls[0] += 1
        if calls[0] > 20000 or np.any(np.abs(s) > limit):
            raise RuntimeError("diverged")
        out = np.asarray(f(*s)).ravel()
        if np.iscomplexobj(out) or not np.all(np.isfinite(out)):
            raise RuntimeError("non-finite")
        return out.astype(float)

    try:
        sol = solve_ivp(rhs, (t[0], t[-1]), x0, t_eval=t, method=solver,
                        rtol=1e-8, atol=1e-10)
        if not sol.success or sol.y.shape[1] != len(t):
            return np.inf
        return float(np.linalg.norm(sol.y.T - X_clean) / np.linalg.norm(X_clean))
    except Exception:
        return np.inf
