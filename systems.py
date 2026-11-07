"""Time-series test systems: ground-truth ODEs, clean trajectories, and noise."""
from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


@dataclass
class System:
    name: str
    names: list            # state variable names, e.g. ["x"]
    f: callable            # f(t, x) -> dx/dt, used to simulate the clean data
    true_rhs: list         # sympy expressions for each dx_i/dt, in terms of `names`
    x0: list
    t: np.ndarray
    description: str = ""
    extra_x0s: list = field(default_factory=list)  # more starting points, for --n-traj > 1
    symbols: list = field(init=False)

    def __post_init__(self):
        self.symbols = [sp.Symbol(n, real=True) for n in self.names]

    @property
    def initial_conditions(self):
        return [self.x0] + self.extra_x0s

    def simulate(self, x0=None):
        x0 = self.x0 if x0 is None else x0
        sol = solve_ivp(self.f, (self.t[0], self.t[-1]), x0, t_eval=self.t,
                        rtol=1e-10, atol=1e-12, method="LSODA")
        return sol.y.T  # shape (n_samples, n_states)


def add_noise(X, level, rng):
    """Gaussian noise with std = level * std of each state column."""
    if level == 0:
        return X.copy()
    return X + rng.normal(0, level * X.std(axis=0), size=X.shape)


def decay():
    k = 0.5
    x = sp.Symbol("x", real=True)
    return System(
        name="decay",
        names=["x"],
        f=lambda t, s: [-k * s[0]],
        true_rhs=[-k * x],
        x0=[3.0],
        t=np.linspace(0, 5, 101),
        description="Exponential decay  x' = -0.5 x",
    )


def logistic():
    r, K = 1.0, 10.0
    x = sp.Symbol("x", real=True)
    return System(
        name="logistic",
        names=["x"],
        f=lambda t, s: [r * s[0] * (1 - s[0] / K)],
        true_rhs=[r * x - (r / K) * x**2],
        x0=[0.5],
        t=np.linspace(0, 8, 101),
        description="Logistic growth  x' = x(1 - x/10)  =  x - 0.1 x^2",
    )


def oscillator():
    w2 = 4.0  # omega = 2, period ~3.14
    x, v = sp.symbols("x v", real=True)
    return System(
        name="oscillator",
        names=["x", "v"],
        f=lambda t, s: [s[1], -w2 * s[0]],
        true_rhs=[v, -w2 * x],
        x0=[1.0, 0.0],
        extra_x0s=[[0.5, 1.0], [-1.5, 0.5]],  # energies x^2 + v^2/4 = 1, 0.5, 2.31
        t=np.linspace(0, 10, 201),
        description="Harmonic oscillator  x' = v,  v' = -4x  (position and velocity observed)",
    )


def powerdecay():
    k = 0.5
    x = sp.Symbol("x", real=True)
    return System(
        name="powerdecay",
        names=["x"],
        f=lambda t, s: [-k * abs(s[0]) ** 1.5],
        true_rhs=[-k * x ** sp.Rational(3, 2)],
        x0=[4.0],
        t=np.linspace(0, 5, 101),  # x falls from 4 to 0.33, staying clear of 0 under noise
        description="Fractional power decay  x' = -0.5 x^1.5",
    )


def lotka_volterra():
    a, b, c, d = 1.0, 0.1, 1.5, 0.075
    x, y = sp.symbols("x y", real=True)
    return System(
        name="lotka_volterra",
        names=["x", "y"],
        f=lambda t, s: [a * s[0] - b * s[0] * s[1], -c * s[1] + d * s[0] * s[1]],
        true_rhs=[a * x - b * x * y, -c * y + d * x * y],
        x0=[10.0, 5.0],
        extra_x0s=[[5.0, 5.0], [15.0, 8.0]],
        t=np.linspace(0, 15, 301),  # period ~5, so ~3 cycles; dt = 0.05
        description="Lotka-Volterra  x' = x - 0.1xy,  y' = -1.5y + 0.075xy",
    )


SYSTEMS = {
    "decay": decay,
    "logistic": logistic,
    "oscillator": oscillator,
    "powerdecay": powerdecay,
    "lotka_volterra": lotka_volterra,
}


# --- ODEBench (Test 10): all 63 systems, each starting point as its own setup ----------
# Equations, constants and starting points exactly as in ODEBench (see import_odebench.py);
# time grid as ODEBench's own trajectories: t = 0..10, 150 samples.
import json as _json
import os as _os

ODEBENCH_FILE = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                              "data", "odebench_equations.json")


def _odebench_system(eq, start, time):
    n = eq["dim"]
    names = [f"x{i}" for i in range(n)]
    syms = sp.symbols(" ".join(names), real=True)
    syms = [syms] if n == 1 else list(syms)
    local = {f"x_{i}": s for i, s in enumerate(syms)}
    local.update({f"c_{i}": sp.Float(c) for i, c in enumerate(eq["consts"][0])})
    local.update({"sin": sp.sin, "cos": sp.cos, "exp": sp.exp, "log": sp.log,
                  "cot": sp.cot, "abs": sp.Abs})
    rhs = [sp.sympify(part.strip().replace("^", "**"), locals=local)
           for part in eq["eq"].split("|")]
    f_np = sp.lambdify(syms, rhs, "numpy")
    return System(
        name=f"ode{eq['id']:02d}_start{start + 1}",
        names=names,
        f=lambda t, s, f_np=f_np: np.asarray(f_np(*s), dtype=float).ravel(),
        true_rhs=rhs,
        x0=[float(v) for v in eq["init"][start]],
        t=np.linspace(time["t_start"], time["t_end"], time["n_samples"]),
        description=f"ODEBench #{eq['id']} ({eq['eq_description']}), start {start + 1}: "
                    + " | ".join(f"D{n_} = {r}" for n_, r in zip(names, rhs)),
    )


def odebench_systems():
    with open(ODEBENCH_FILE) as fh:
        data = _json.load(fh)
    out = {}
    for eq in data["equations"]:
        for start in range(len(eq["init"])):
            name = f"ode{eq['id']:02d}_start{start + 1}"
            out[name] = (lambda eq=eq, start=start, time=data["time"]:
                         _odebench_system(eq, start, time))
    return out


ODEBENCH = odebench_systems() if _os.path.exists(ODEBENCH_FILE) else {}
SYSTEMS.update(ODEBENCH)
