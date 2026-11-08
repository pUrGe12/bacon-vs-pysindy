"""BACON's home turf: static laws measured in controlled experiments (Test 12 on).

No time, no derivatives: a table where each independent variable is set to a few levels
in every combination (a full factorial design) and one dependent variable is measured.

The laws are BACON's own datasets (datasets.py in the BACON repo), with the same
variables, levels and constants. Each law comes in two designs:
  given  exactly the authors' levels (3 per variable; Kepler 6 points, synthetic 5).
         run_static.py checks these match BACON's datasets.py value for value.
  dense  more levels over the same ranges: 6 per variable, or 12 points when there is
         a single input, so regression methods get more rows than they have terms.
Two of BACON's datasets are left out: `birthday` (no closed-form law to recover) and
`lv` (a time series; Tests 1-11 cover those).

Noise, as in BACON's datasets.py: Gaussian on the dependent variable only, with standard
deviation `noise` times its value.
"""
from dataclasses import dataclass, field
from itertools import product

import numpy as np
import sympy as sp

DESIGNS = ("given", "dense")


@dataclass
class Law:
    name: str
    description: str
    inputs: list            # independent variable names, in BACON's column order
    target: str             # dependent variable name (BACON's last column)
    levels: dict            # given design: name -> list of levels
    true: str               # target = this expression (sympy syntax)
    dense_levels: dict = field(default_factory=dict)
    # Kepler: P and D are both made from a planet index n; this maps n to the inputs.
    from_index: object = None

    def symbols(self):
        return {v: sp.Symbol(v, positive=True) for v in self.inputs + [self.target]}

    def true_expr(self):
        return sp.sympify(self.true, locals=self.symbols())

    def design(self, kind):
        """Rows of inputs, in BACON's order (last input varies fastest)."""
        if self.from_index is not None:
            n = (np.arange(1, 7, dtype="float64") if kind == "given"
                 else np.linspace(1, 6, 12))
            return {v: col for v, col in zip(self.inputs, self.from_index(n))}
        if kind == "given":
            levels = self.levels
        else:
            k = 12 if len(self.inputs) == 1 else 6
            levels = {v: np.linspace(min(l), max(l), k) for v, l in self.levels.items()}
        rows = np.array(list(product(*[levels[v] for v in self.inputs])), dtype="float64")
        return {v: rows[:, i] for i, v in enumerate(self.inputs)}

    def evaluate(self, cols):
        syms = self.symbols()
        f = sp.lambdify([syms[v] for v in self.inputs], self.true_expr(), "numpy")
        return np.asarray(f(*[cols[v] for v in self.inputs]), dtype="float64")

    def data(self, kind):
        cols = self.design(kind)
        cols[self.target] = self.evaluate(cols)
        return cols

    def held_out(self, n=200, seed=12345):
        """Random points inside the design's box, for checking a law away from the grid."""
        rng = np.random.default_rng(seed)
        given = self.design("given")
        return {v: rng.uniform(given[v].min(), given[v].max(), n) for v in self.inputs}


def add_noise(y, noise, rng):
    return y + rng.normal(0.0, 1.0, len(y)) * noise * np.abs(y) if noise > 0 else y.copy()


LAWS = {
    "ideal": Law("ideal", "Ideal gas law, V = M(T + 273)/P",
                 ["M", "T", "P"], "V",
                 {"M": [1, 2, 3], "T": [10, 20, 30], "P": [1000, 2000, 3000]},
                 "M*(T + 273)/P"),
    "ohm": Law("ohm", "Ohm's law for a wire, I = T D^2 / (2 (L + 3))",
               ["T", "D", "L"], "I",
               {"T": [100, 120, 140], "D": [0.01, 0.02, 0.03], "L": [0.5, 1, 1.5]},
               "T*D**2/(2*(L + 3))"),
    "black": Law("black", "Black's law of mixtures, T_f = (M_1 T_1 + M_2 T_2)/(M_1 + M_2)",
                 ["M_1", "M_2", "T_1", "T_2"], "T_f",
                 {"M_1": [1, 2, 3], "M_2": [1, 2, 3], "T_1": [50, 60, 70],
                  "T_2": [50, 60, 70]},
                 "(M_1*T_1 + M_2*T_2)/(M_1 + M_2)"),
    "kepler": Law("kepler", "Kepler's third law, D^3 = P^2 (so D = P^(2/3))",
                  ["P"], "D", {}, "P**(2/3)",
                  from_index=lambda n: [n ** 3]),
    "synthetic": Law("synthetic", "BACON's synthetic test, Y = 3/(X + 2)",
                     ["X"], "Y", {"X": [1, 2, 3, 4, 5]}, "3/(X + 2)"),
}

SETUPS = [f"{law}_{d}" for law in LAWS for d in DESIGNS]


def setup(name):
    law, design = name.rsplit("_", 1)
    return LAWS[law], design
