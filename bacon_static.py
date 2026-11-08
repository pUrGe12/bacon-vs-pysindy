"""The original BACON pipeline (BACON.3 layers over BACON.1) on a static table, unchanged.

This is Miller & Banerjee's own code path, the one `main.py --dataset ideal` runs:
data_space -> layer -> BACON.1, then utils.losses combines the per-layer invariants and
solves for the dependent variable. None of our time-series fixes (B1-B7) is involved.

Two things are done around it, neither changes a decision BACON makes:
  * Their pipeline prints the final law (df_helper.score); we capture it instead, using
    the same two lines score() runs (simplify_eqns(...).iterate_through_dummys()).
  * laws_manager imports PySR and layer_select imports GPy at module level. Neither is
    used on this path (no PySR laws, no GP denoising or GP ranking), and PySR needs Julia,
    so if they are not installed an empty stand-in module is registered for them.
A run that raises (e.g. "No relationships found compatible with BACON.3") or exceeds the
time limit counts as "no law found".
"""
import contextlib
import io
import signal
import sys
import types

import pandas as pd
import sympy as sp

from bacon_dt import BACON_PATH

for _mod, _attrs in {"pysr": ["PySRRegressor"], "GPy": []}.items():
    try:
        __import__(_mod)
    except ImportError:
        _stub = types.ModuleType(_mod)
        for _a in _attrs:
            setattr(_stub, _a, None)
        sys.modules[_mod] = _stub

if BACON_PATH not in sys.path:
    sys.path.insert(0, BACON_PATH)

import datasets as bacon_datasets                                  # noqa: E402
from space_of_laws.laws_manager import laws_main                  # noqa: E402
from space_of_data.layer_manager import layer_main                # noqa: E402
from space_of_data.space_methods import data_space_manager as dsm  # noqa: E402
from utils import losses                                           # noqa: E402

# Frozen setting (Test 12 on): BACON.3 layer selection, every BACON.1 and data-space
# parameter at the code's own default (laws_manager.bacon_1, data_space.__init__), i.e.
# what main.py runs with no --args file.
FROZEN = {"layer_method": "bacon.3", "layer_args": {"verbose": False},
          "laws_method": "bacon.1", "laws_args": {}, "data_space_args": {}}
CODE_DEFAULTS = {"epsilon": 0.001, "delta": 0.1, "epsilon_scale": 1.1,
                 "delta_scale": 1.05, "c_val": 0.02, "Delta": 0.1}
TIME_LIMIT = 300  # seconds per run


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def discover(cols, inputs, target, config=FROZEN, time_limit=TIME_LIMIT):
    """Returns (law or None, invariants found, status). law is a sympy expression for the
    target in terms of the inputs (BACON's own symbols)."""
    syms = [sp.Symbol(v, real=True) for v in inputs + [target]]   # as in datasets.py
    df = pd.DataFrame({s: cols[v] for s, v in zip(syms, inputs + [target])})
    captured = {}

    def score(init_df, eqns):  # replaces df_helper.score: same computation, no printing
        captured["invariants"] = list(eqns)
        captured["law"] = losses.simplify_eqns(init_df, eqns, init_df.columns[-1]) \
                                .iterate_through_dummys()

    old_score, dsm.df_helper.score = dsm.df_helper.score, score
    old_handler = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(time_limit)
    status = "ok"
    try:
        layer = config.get("layer_method") or None
        laws = laws_main(config["laws_method"], config.get("laws_args") or {})
        with contextlib.redirect_stdout(io.StringIO()):
            if len(inputs) == 1 and not layer:
                layer = "bacon.3"   # never used: with two columns there are no layers
            ds = dsm.data_space(df, layer_main(layer, config.get("layer_args") or {}), laws,
                                **(config.get("data_space_args") or {}))
            ds.run_iterations()
    except _Timeout:
        status = "time limit"
    except Exception as e:  # BACON gave up (no relationship, invalid layer, ...)
        status = f"{type(e).__name__}: {str(e)[:80]}"
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
        dsm.df_helper.score = old_score

    law = captured.get("law")
    if law is None:
        return None, captured.get("invariants", []), status if status != "ok" else "no law"
    # BACON's final form is Eq(expression, target)
    expr = law.lhs if law.rhs == syms[-1] else sp.solve(law, syms[-1])[-1]
    return expr, captured["invariants"], status


def authors_config(path):
    import json
    with open(path) as fh:
        a = json.load(fh)
    return {k: a.get(k) for k in ("layer_method", "layer_args", "laws_method",
                                  "laws_args", "data_space_args")}
