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


SYSTEMS = {
    "decay": decay,
    "logistic": logistic,
    "oscillator": oscillator,
}
