"""Derivative estimates shared by both methods, so the comparison isolates the search."""
import pysindy as ps

METHODS = {
    # Plain central differences: PySINDy's own default.
    "fd": lambda: ps.FiniteDifference(order=2),
    # Savitzky-Golay smoothing, then differences.
    "sfd": lambda: ps.SmoothedFiniteDifference(
        smoother_kws={"window_length": 11, "polyorder": 3}),
}


def differentiate(X, t, method):
    return METHODS[method]()(X, t)
