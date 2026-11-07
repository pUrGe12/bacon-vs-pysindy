# BACON vs PySINDy on time-series data

Can BACON (Langley's heuristic law finder, as rebuilt in
[Miller & Banerjee 2025](https://www.researchgate.net/publication/389562521)) discover
ODEs if we give it a notion of derivatives? We compare it against PySINDy.

## Setup

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_test.py --system decay
```

BACON is imported from `~/Desktop/BACON` (set `BACON_PATH` to change). That repo is not modified.

## How it works

| File | What it does |
|---|---|
| `systems.py` | Ground-truth ODEs, clean trajectories, noise (std = level × std of each state) |
| `derivatives.py` | Derivative estimates given to **both** methods: `fd` (central differences, PySINDy's default) and `sfd` (Savitzky–Golay smoothed, window 11) |
| `bacon_dt.py` | BACON with derivative columns `Dx`, plus fixes (below) |
| `sindy_runner.py` | PySINDy: polynomial library degree 3, STLSQ. Threshold picked automatically: sweep 0.001–1, keep the sparsest model whose R² is within 0.01 of the best (since Test 2; Test 1 used a fixed 0.1) |
| `scoring.py` | Term structure match, coefficient error, error when simulating the learned ODE |
| `run_test.py` | Runs a system across noise levels × seeds × derivative methods |

**BACON-dt** pairs each derivative `D(x_i)` with each state `x_j` and runs BACON.1 to look for
an invariant (e.g. `x/Dx = const`), then solves the invariant for `Dx` to get the ODE.

Fixes to BACON.1 needed for ODE data:
1. **Sign-blind heuristic.** BACON.1 picks product vs ratio from the sign of the
   correlation. For `Dx = -kx` the correlation is negative, so it multiplies (`Dx·x`,
   `Dx·x²`, …) and never tries the ratio. BACON-dt decides from the trend of `|a|`, `|b|`.
2. **Constancy test broken for negative values.** `M(1-δ) < v < M(1+δ)` is never true when
   `M < 0`. Replaced with `|v - M| < δ|M|` (identical for `M > 0`).
3. **Unbounded retries.** On failure BACON.1 loosens its tolerances and retries without end,
   so it always returns *something*. BACON-dt loosens δ from 1% by ×1.5 up to 30%, then
   reports "no invariant".

**BACON-dt-robust** (since Test 2) also replaces "every point within δ of the mean" with
"the typical point (median absolute deviation) within δ of the median", and reports the
constant as the median.

Framing: BACON was never built for time series, so losing here is expected. Later we
will also adapt PySINDy to BACON's home turf (static laws like Kepler, ideal gas) to make
the comparison fair in both directions.

"Recovered" = exactly the right terms **and** coefficients within 5%.

## Results log

### Test 1: Exponential decay, `x' = -0.5x`
x0 = 3, t ∈ [0, 5], 101 samples, 5 noise seeds (1 at zero noise).

Recovery rate (out of 5 runs; 1 run at zero noise):

| noise | BACON-dt (fd) | BACON-dt (sfd) | PySINDy (fd) | PySINDy (sfd) |
|---|---|---|---|---|
| 0%  | 1/1 | 1/1 | 1/1 | 1/1 |
| 1%  | 0/5 | 0/5 | 5/5 | 5/5 |
| 5%  | 0/5 | 0/5 | 4/5 | 5/5 |
| 10% | 0/5 | 0/5 | 2/5 | 3/5 |

**Takeaways**
- On clean data both get `x' = -0.5001x`. BACON-dt finds it as the invariant `x/Dx = -2`.
- With any noise BACON-dt finds nothing. The invariant is there: the median of `x/Dx` is
  -2.005 (1% noise, sfd) and -1.94 (5%, sfd). But BACON's test requires **every** point within δ
  of the mean, and one bad point rules the law out:

  | noise, deriv | typical point off by | worst point off by |
  |---|---|---|
  | 1%, sfd | 4% | 43% |
  | 1%, fd  | 46% | 3288% |
  | 5%, sfd | 23% | 3172% |

- One near-miss (1% noise, sfd, seed 4): at δ = 5% BACON accepted a convoluted invariant,
  `x(1.97Dx + x)/Dx = -0.051`, i.e. `Dx = -x²/(1.97x + 0.051)`. It simulates almost
  perfectly (0.5% error) but has the wrong structure.
- PySINDy is robust because least squares averages noise over all points. Its failures at
  high noise are extra terms (`x²`, `x³`, constant) that stay above the 0.1 threshold.

**Next step to try:** make BACON's constancy test robust, e.g. use a median and a
spread measure (MAD), or "most points within δ" instead of "all points within δ", and
estimate the constant with the median, not the mean.

### Test 2: Logistic growth, `x' = x(1 - x/10) = x - 0.1x²`
x0 = 0.5, t ∈ [0, 8], 101 samples, 5 noise seeds (1 at zero noise).

| noise | BACON-dt (fd / sfd) | BACON-dt-robust (fd / sfd) | PySINDy (fd / sfd) |
|---|---|---|---|
| 0%  | 0/1 · 0/1 | 0/1 · 0/1 | 1/1 · 1/1 |
| 1%  | 0/5 · 0/5 | 0/5 · 0/5 | 4/5 · 5/5 |
| 5%  | 0/5 · 0/5 | 0/5 · 0/5 | 1/5 · 3/5 |
| 10% | 0/5 · 0/5 | 0/5 · 0/5 | 1/5 · 3/5 |

**Takeaways**
- BACON-dt finds nothing, even on clean data, and the robust variant doesn't help. The
  cause is BACON's search heuristic, not its law test: `Dx` rises then falls as `x` grows,
  but BACON assumes steady (monotonic) relationships. It sees a mostly downward trend,
  tries `Dx·x`, and runs out of moves. It never tries `Dx/x`.
- Given `Dx/x` as a column, BACON finds `Dx/x + 0.1x = 1` (exactly the true law) on clean
  data and at 1% noise (δ = 7.6%). It fails at 5% noise.
- PySINDy gets the `x` and `x²` terms right in every run. All its failures are one extra
  constant term (+0.06 to +0.68) that makes the structure "wrong", although the simulated
  trajectory stays within 1–8%.
- This test exposed a bug in my PySINDy setup: a fixed threshold of 0.1 deleted the
  true `-0.1x²` term, so PySINDy failed even on clean data. The threshold is now chosen
  automatically, the same way for every system. Test 1 still used the fixed 0.1, but its
  PySINDy coefficient (0.5) is far above it. A quick check with the new rule (smoothed
  derivatives) gives the same story for decay (5/5, 5/5, 3/5 at 1/5/10%).

**Next step to try:** a derivative-aware search for BACON-dt: always try the per-unit rate
`Dx_i/x_j` (and `Dx_i/x_j²`, …) as a starting term, instead of relying only on the
monotonic up/down heuristic.

### Test 3: Harmonic oscillator, `x' = v`, `v' = -4x`
Position and velocity both observed (noisy). x0 = 1, v0 = 0, t ∈ [0, 10] (~3 periods),
201 samples, 5 noise seeds (1 at zero noise). First test with two state variables.

Cells: exactly recovered / runs (in brackets: runs whose learned ODE simulates within 10%
of the true trajectory).

| noise | BACON-dt (fd · sfd) | BACON-dt-robust (fd · sfd) | PySINDy (fd · sfd) |
|---|---|---|---|
| 0%  | 1/1 (1) · 1/1 (1) | 1/1 (1) · 1/1 (1) | **0/1 (1) · 0/1 (0)** |
| 1%  | 0/5 (0) · 0/5 (0) | 0/5 (3) · 0/5 (2) | 5/5 (5) · 5/5 (5) |
| 5%  | 0/5 (0) · 0/5 (1) | 0/5 (0) · 0/5 (4) | 3/5 (5) · 5/5 (5) |
| 10% | 0/5 (0) · 0/5 (0) | 0/5 (0) · 1/5 (2) | 2/5 (2) · 2/5 (2) |

**Takeaways**
- **First BACON win: clean data.** BACON-dt finds `v/Dx = 1.002` and `x/Dv = -0.2504`, i.e.
  `x' = 0.998v`, `v' = -3.993x`. (The 0.2% shortfall is the usual finite-difference error on
  oscillations, sin(ω·dt)/(ω·dt); PySINDy shows it too.)
- **PySINDy fails on clean data because its term library is degenerate.** A single
  noise-free trajectory stays on one energy level, `x² + v²/4 = 1`, so the library columns
  `x`, `x³` and `x·v²` are linearly dependent: `x = x³ + x·v²/4`. Any multiple of
  `x(1 - x² - v²/4)` can be added for free. With smoothed derivatives it returned
  `v' = 4.975e5·x - 4.975e5·x³ - 1.244e5·x·v²` (= `-4x` plus a huge multiple of zero), which
  blows up as soon as it leaves that energy level. Noise breaks the exact dependence,
  which is why PySINDy does better with noise. This is a known SINDy pitfall; the usual fix
  is data from several initial conditions. BACON is immune because it only ever compares
  two columns at a time.
- **Under noise, BACON's misses are mostly near-misses.** The robust variant often gets the
  dynamics right in simulation (4/5 at 5% with sfd) but fails the exact check because of
  a tiny constant, e.g. `Dx = 1.0v - 0.0038`. Cause: BACON.1 takes the "linear" branch
  (`Dx - m·v = const`) when the intercept is large compared with the **mean** of `Dx`.
  An oscillating signal has a mean near zero, so any noise-level intercept looks large, and
  BACON writes a line instead of the ratio `v/Dx = 1`. Same root cause as Test 1's sign
  problem: BACON.1 assumes positive quantities.
- **Zero crossings hurt strict BACON.** `x/Dv` blows up wherever `x` and `Dv` cross zero, so
  the strict "every point" test fails (e.g. 1% noise, seed 3: no law found for `v'`). The
  robust median test survives these and finds `x/Dv = -0.2512`.
- Solving some odd invariants for `Dx` gives monsters (cube roots, coefficients of 1e16),
  e.g. `v(v - 1.009·Dx)/Dx = -0.003`. They are close to `Dx = v` but have the wrong structure.

Also fixed: `scoring.py` now simulates with RK45. LSODA (Fortran) hung when our "this
equation is blowing up" guard raised an error inside it.

**Next steps to try**
1. BACON's linear branch: compare the intercept with the spread (std) of the data, not its
   mean, so zero-mean signals don't trigger it.
2. PySINDy: give it several trajectories from different initial conditions. That is the
   standard practice, and it removes the library degeneracy.
