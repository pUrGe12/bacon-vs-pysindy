"""PySINDy baseline: polynomial library + STLSQ, fed the same derivatives as BACON.

The STLSQ threshold is not fixed. A fixed threshold silently drops any true term
whose coefficient is below it (e.g. the -0.1 x^2 of the logistic equation), and
tuning it per system would use knowledge of the answer. Instead we sweep a grid
and keep the sparsest model whose R^2 is within R2_TOL of the best one.
(AIC/BIC were tried too: they keep tiny junk terms on noise-free data.)
Pass a number as `threshold` to use one fixed threshold instead (Test 1 used 0.1).

`library="poly+halfpowers"` adds |x|^0.5, |x|^1.5, |x|^2.5 for each state to the
polynomial library. This is an *informed* baseline (we know the answer has a
half-integer power); the default is the plain polynomial library.
`normalize=True` turns on STLSQ's `normalize_columns` (each library column divided by
its L2 norm before thresholding, so the threshold no longer depends on each term's
scale). The threshold then applies to coefficient x column norm, which is much larger
than the raw coefficient, so the sweep uses THRESHOLDS_NORM: the same density (4 per
decade) over 0.001-1000. Same selection rule.

`alpha` is STLSQ's ridge penalty (PySINDy default 0.05). STLSQ thresholds the *ridge*
coefficients and only refits without the penalty at the end. On raw columns (lengths
~10-6000 here) 0.05 is negligible; on normalized columns (length 1) it is a 5% penalty
that spreads weight over correlated terms (x, x^2, x^3), so junk terms are never
thresholded out. Test 9 found this; alpha=0 is plain least squares, as in the original
SINDy algorithm (Brunton et al. 2016).
See CHANGES_FROM_BACON.md, S1-S6.
"""
import warnings

import numpy as np
import pysindy as ps
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

THRESHOLDS = np.logspace(-3, 0, 13)
THRESHOLDS_NORM = np.logspace(-3, 3, 25)
R2_TOL = 0.01


HALF_POWERS = (0.5, 1.5, 2.5)
# One plain one-argument function per power: CustomLibrary counts each function's
# arguments to decide how many inputs it takes, so `lambda x, p=p: ...` gets dropped.
_HALF_FUNCS = [lambda x: np.abs(x) ** 0.5, lambda x: np.abs(x) ** 1.5,
               lambda x: np.abs(x) ** 2.5]
_HALF_NAMES = [lambda x: f"{x}^0.5", lambda x: f"{x}^1.5", lambda x: f"{x}^2.5"]


def _library(kind, degree):
    poly = ps.PolynomialLibrary(degree=degree)
    if kind == "poly":
        return poly
    if kind == "poly+halfpowers":
        half = ps.CustomLibrary(library_functions=_HALF_FUNCS, function_names=_HALF_NAMES)
        return ps.GeneralizedLibrary([poly, half])
    raise ValueError(kind)


def discover(X, Xdot, dt, symbols, degree=3, threshold="sweep", library="poly",
             normalize=False, alpha=0.05):
    names = [str(s) for s in symbols]
    grid = THRESHOLDS_NORM if normalize else THRESHOLDS
    thresholds = grid if threshold == "sweep" else [float(threshold)]
    total_var = np.sum((Xdot - Xdot.mean(axis=0)) ** 2)
    fits = []
    for thr in thresholds:
        model = ps.SINDy(optimizer=ps.STLSQ(threshold=thr, normalize_columns=normalize, alpha=alpha),
                         feature_library=_library(library, degree))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # "threshold eliminated all coefficients"
            model.fit(X, t=dt, x_dot=Xdot, feature_names=names)
        rss = np.sum((Xdot - model.predict(X)) ** 2)
        k = int(np.count_nonzero(model.coefficients()))
        if k > 0:
            fits.append((k, rss, model))
    best_rss = min(f[1] for f in fits)
    ok = [f for f in fits if f[1] <= best_rss + R2_TOL * total_var]
    model = min(ok, key=lambda f: (f[0], f[1]))[2]

    local = {str(s): s for s in symbols}
    feats = [parse_expr(f.replace("^", "**").replace(" ", "*"), local_dict=local)
             for f in model.get_feature_names()]
    rhs = []
    for row in model.coefficients():
        rhs.append(sum((float(c) * f for c, f in zip(row, feats) if c != 0), sp.Integer(0)))
    return rhs
