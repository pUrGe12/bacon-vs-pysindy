# BACON vs PySINDy: analysis of Tests 1–10 (ODE time series) and Test 12 (static laws)

**Answer: for discovering ODEs from time series, PySINDy is the better tool.** On an
outside benchmark that neither method was tuned for (ODEBench, Test 10), PySINDy finds the
exact law in more setups at every noise level, and its lead grows with noise. BACON,
adapted with derivative columns and five fixes, is competitive only on simple laws, where
a variable's rate of change depends on a single variable. It wins only where PySINDy hits
known weak spots.

This was PySINDy's home turf: time series and derivatives are what SINDy was built for,
and BACON (Langley's 1970s–80s law finder, in the 2025 rebuild by Miller & Banerjee) never
was. The reverse comparison, both methods on BACON's own static laws (Test 12), is at the
end: there BACON finds laws PySINDy cannot express, but only its two-quantity search
survives noise.

Tables in this file are generated from the result files by
`.venv/bin/python summarize.py tables`. Details of every test: [README.md](README.md).
Every change to either method: [CHANGES_FROM_BACON.md](CHANGES_FROM_BACON.md).

## What was compared

- **Task:** given noisy samples of a system's state over time, find the equations
  `x' = f(x)` that generate them.
- **Same inputs for both:** each method gets the same noisy data and the same derivative
  estimates, either **smoothed** (slope of a small curve fitted over 11 points, which
  averages noise out) or **raw** (slope between neighbouring points, which magnifies noise).
  So the comparison is between their search methods, not their differentiation.
- **BACON-dt:** BACON.1 looks for a combination of two quantities that stays constant (like
  Kepler's `D³/P²`). We give it a derivative column `Dx`; an invariant such as `x/Dx = -2`
  is solved for `Dx`, giving `x' = -0.5x`.
- **PySINDy:** fits each derivative as a sparse sum of candidate terms (here all
  polynomials up to degree 3), keeping only terms above a threshold.
- **"Found the exact law":** exactly the right terms, and every coefficient within 5%.
- **Runs:** every run uses the same system, starting point and time grid; runs differ in
  noise level (0, 1, 5, 10%), noise draw (5 per level in Tests 1–9, 20 in Test 10) and
  derivative type.

## Finding 1: on systems neither method was tuned for, PySINDy wins

**Test 10:** all 63 ODEBench systems (mostly from Strogatz's textbook), each from both of
its starting points, so 126 setups. 20 noise draws per level, ODEBench's own constants and
time grid. Both methods frozen before the run, with the rules written down first.

56 of the 126 setups have laws that PySINDy's terms cannot express (`sin`, `exp`, `log`,
fractions, powers like `x^1.2`); **neither method ever found any of them**. The table below
is for the other 70, whose laws are polynomials of degree ≤ 3. Counts are of **setups**: a
method counts for a setup when it finds the exact law in more than half its runs (the 1
clean run, or more than 10 of the 20 noisy ones).

<!-- tables:analysis_odebench -->

Smoothed derivatives — Setups whose true law is inside PySINDy's terms (polynomial, degree ≤ 3) (70 of 70 setups run)

| | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|
| BACON finds the exact law in more than half the runs | 17 of 70 | 16 of 70 | 6 of 70 | 4 of 70 |
| PySINDy finds the exact law in more than half the runs | 23 of 70 | 22 of 70 | 15 of 70 | 10 of 70 |
| BACON finds it in more runs than PySINDy | 5 of 70 | 9 of 70 | 4 of 70 | 4 of 70 |
| PySINDy finds it in more runs than BACON | 11 of 70 | 13 of 70 | 18 of 70 | 13 of 70 |
| same number of runs, more than zero | 12 of 70 | 6 of 70 | 3 of 70 | 1 of 70 |
| neither finds it in any run | 42 of 70 | 42 of 70 | 45 of 70 | 52 of 70 |

Raw derivatives — Setups whose true law is inside PySINDy's terms (polynomial, degree ≤ 3) (70 of 70 setups run)

| | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|
| BACON finds the exact law in more than half the runs | 18 of 70 | 14 of 70 | 2 of 70 | 0 of 70 |
| PySINDy finds the exact law in more than half the runs | 29 of 70 | 24 of 70 | 15 of 70 | 7 of 70 |
| BACON finds it in more runs than PySINDy | 2 of 70 | 6 of 70 | 3 of 70 | 3 of 70 |
| PySINDy finds it in more runs than BACON | 13 of 70 | 16 of 70 | 23 of 70 | 18 of 70 |
| same number of runs, more than zero | 16 of 70 | 6 of 70 | 0 of 70 | 1 of 70 |
| neither finds it in any run | 39 of 70 | 42 of 70 | 44 of 70 | 48 of 70 |

<!-- /tables -->

- PySINDy leads at every noise level, with both derivative types. The gap is smallest on
  clean data and grows with noise.
- In 42 of the 70 setups neither method finds the law even on clean data (smoothed
  derivatives; 39 with raw). This includes
  every 3- and 4-variable system (e.g. Lorenz, Rössler): 150 samples over 10 time units is
  too coarse for their derivatives. That is a limit of this benchmark's sampling as much as
  of either method.
- Full per-setup tables: [ODEBENCH_RESULTS.md](ODEBENCH_RESULTS.md).

## Finding 2: BACON can only find laws where a derivative depends on one variable

BACON.1 relates two quantities at a time, and our adaptation pairs each derivative with
one state variable (or a per-unit rate with one other variable). So it can express
`x' = -0.5x`, `x' = x - 0.1x²` or `x' = x(a - by)`, but never a law that needs two other
variables in separate terms.

The evidence is consistent across every test:

- **Every law of that kind failed in every run.** On ODEBench: the damped oscillator (#25:
  BACON finds `v' = -4.29x` and drops the damping term `-0.43v`), Van der Pol (#37, #38),
  Duffing (#40), the Brusselator (#49), the falling object (#5). PySINDy finds each of these
  in at least one setup, mostly at 0–1% noise.
- **Every ODEBench law BACON did find has this form:** the RC circuit (#1), growth models
  (#2, #3, #6, #9, #12), `x' = -x³` (#11), the undamped oscillator (#24), Lotka–Volterra
  (#27).

This is a property of the search, not of noise, so more data would not change it. A
damped oscillator was on the list of things to try after Test 8 for exactly this reason.

## Finding 3: BACON wins only where PySINDy has known weak spots

Where BACON beats PySINDy, the cause is always on PySINDy's side:

1. **One clean trajectory.** Along a single noise-free orbit, PySINDy's candidate terms can
   be exact combinations of each other (on the oscillator, `x = x³ + x·v²/4`), so any
   multiple of zero can be added and it returns huge cancelling coefficients. Undamped
   oscillator on clean data: BACON finds it, PySINDy returns `1.3e6·x0³ …` (Test 3, and
   ODEBench #24). BACON is immune because it only ever compares two columns. The standard
   remedy, data from several starting points, makes PySINDy recover it (Test 4). Test 11
   would have applied it to ODEBench; it was stopped before any results.
2. **Junk terms on Lotka–Volterra.** PySINDy keeps the right terms plus small extra ones.
   ODEBench #27, start 2, smoothed: BACON finds the exact law in 20, 19 and 10 of 20 runs at
   1/5/10% noise, PySINDy in 0 at every level (Test 10). On our own Lotka–Volterra
   system BACON also wins at 10% noise (Test 9).
3. **Data that barely moves.** ODEBench #9, start 2, starts at 0.55 next to its resting
   point 0.533; PySINDy's columns are nearly identical and it finds nothing, BACON 17 of 20
   at 1% noise. (Our noise is scaled to how much a variable moves, so it is tiny here too.)

On our own systems (Tests 1–9, final settings of both methods, 5 noise draws per level):

<!-- tables:test9_versus -->

Smoothed derivatives: runs that found the exact law

| system — method | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| Exponential decay — BACON-dt-robust | 1 of 1 | 5 of 5 | 4 of 5 | 0 of 5 |
| Exponential decay — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 3 of 5 |
| Logistic growth — BACON-dt-robust | 1 of 1 | 5 of 5 | 4 of 5 | 1 of 5 |
| Logistic growth — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 4 of 5 |
| Harmonic oscillator — BACON-dt-robust | 1 of 1 | 5 of 5 | 5 of 5 | 4 of 5 |
| Harmonic oscillator — PySINDy | 0 of 1 | 5 of 5 | 5 of 5 | 5 of 5 |
| Harmonic oscillator (3 traj.) — BACON-dt-robust | 1 of 1 | 5 of 5 | 5 of 5 | 5 of 5 |
| Harmonic oscillator (3 traj.) — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 5 of 5 |
| Power decay — BACON-dt-robust | 1 of 1 | 4 of 5 | 0 of 5 | 0 of 5 |
| Power decay — PySINDy | 0 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Power decay (½ powers) — BACON-dt-robust | 1 of 1 | 4 of 5 | 0 of 5 | 0 of 5 |
| Power decay (½ powers) — PySINDy | 1 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Lotka–Volterra — BACON-dt-robust | 1 of 1 | 5 of 5 | 5 of 5 | 5 of 5 |
| Lotka–Volterra — PySINDy | 1 of 1 | 5 of 5 | 4 of 5 | 0 of 5 |

Raw derivatives: runs that found the exact law

| system — method | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| Exponential decay — BACON-dt-robust | 1 of 1 | 3 of 5 | 0 of 5 | 0 of 5 |
| Exponential decay — PySINDy | 1 of 1 | 5 of 5 | 4 of 5 | 2 of 5 |
| Logistic growth — BACON-dt-robust | 1 of 1 | 5 of 5 | 1 of 5 | 0 of 5 |
| Logistic growth — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 4 of 5 |
| Harmonic oscillator — BACON-dt-robust | 1 of 1 | 5 of 5 | 1 of 5 | 0 of 5 |
| Harmonic oscillator — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 2 of 5 |
| Harmonic oscillator (3 traj.) — BACON-dt-robust | 1 of 1 | 5 of 5 | 0 of 5 | 0 of 5 |
| Harmonic oscillator (3 traj.) — PySINDy | 1 of 1 | 5 of 5 | 5 of 5 | 5 of 5 |
| Power decay — BACON-dt-robust | 1 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Power decay — PySINDy | 0 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Power decay (½ powers) — BACON-dt-robust | 1 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Power decay (½ powers) — PySINDy | 1 of 1 | 0 of 5 | 0 of 5 | 0 of 5 |
| Lotka–Volterra — BACON-dt-robust | 1 of 1 | 5 of 5 | 1 of 5 | 0 of 5 |
| Lotka–Volterra — PySINDy | 1 of 1 | 5 of 5 | 2 of 5 | 1 of 5 |

<!-- /tables -->

One BACON strength showed up only on our own data: on power decay `x' = -0.5x^1.5`
(Test 6), BACON found the fractional power by itself (4 of 5 runs at 1% noise, smoothed),
which PySINDy's polynomial terms cannot express at all. Given half-integer powers as extra
terms, PySINDy also finds it, but then we told it the answer's form. On ODEBench, none of
the 56 laws outside PySINDy's terms was found by BACON either, so this advantage is not
confirmed on outside data.

## Finding 4: the original BACON cannot do time series without fixes

With every fix off (the original BACON.1 code), BACON returned a **wrong law in every run
on every system**, clean data included. Each fix is a named switch, measured by switching
it off on its own (ablations 1–5 in CHANGES_FROM_BACON.md):

| Fix | What it changes | Without it |
|---|---|---|
| B1 `sign_aware` | Handles negative derivatives | Nothing found on any system (before B5), clean data included |
| B3 `bounded_retry` | Caps BACON's retry loop, so it can say "nothing found" | Most noisy runs return wrong laws instead of nothing |
| B5 `rate_seed` | Also starts from per-unit rates `Dx/x` | Logistic growth never found (smoothed, 1% noise: 5 of 5 → 0 of 5) |
| B6 `robust_line` | Robust straight-line test (Theil–Sen slope, same tolerance as other laws) | Logistic at 5% noise, smoothed: 4 of 5 → 0 of 5; Lotka–Volterra at 5–10%: 5 of 5 → 0 of 5 |
| B7 `own_rate_seed` | Pairs `Dx/x` with each *other* variable | Lotka–Volterra never found (0 at every noise level) |
| Robust variant | A law counts when the *typical* point is within tolerance, not every point | Decay at 1% noise: 5 of 5 → 0 of 5; the original "every point" test never finds the exact law under any noise |

B2 and B4 were tried and dropped: B2 changed one cell in all of ablation 1; B4 had no
effect once B6 was on. `check_original.py` confirms that with every fix off our code makes
exactly the same decisions as the original BACON.1.

Caveat: these fixes were built while looking at Tests 1–8, which is why Test 10 used an
outside benchmark with BACON frozen.

## Finding 5: how PySINDy is configured matters as much as the method

- **Threshold:** a fixed threshold of 0.1 deleted the true `-0.1x²` logistic term, so
  PySINDy failed even on clean data (Test 2). We pick the threshold automatically, the same
  way for every system.
- **Normalized columns (Test 9):** turning on PySINDy's `normalize_columns` with its
  default ridge penalty (0.05) broke it: it found nothing on logistic and Lotka–Volterra
  at any noise level, clean data included, because on normalized columns the penalty keeps correlated junk terms alive.
  With the penalty off, PySINDy improved on logistic, the oscillator and Lotka–Volterra and
  got worse nowhere on its default terms. That setting (chosen after seeing 9a, which is
  noted) is the one frozen for Test 10.

Without these, Tests 1–8 understated PySINDy. Any comparison like this one depends on
giving the baseline a fair setup.

## What these results do not establish

- **That PySINDy is better overall.** Findings 1–5 are time series, PySINDy's home turf.
  On BACON's static laws (Test 12) the picture reverses, though only on BACON's own five
  datasets.
- **Anything about 3- and 4-variable systems.** Nobody found them, partly because of
  ODEBench's coarse sampling.
- **Results under other conditions:** one noise model (Gaussian, scaled to each variable's
  spread), one sampling grid per system, and PySINDy with polynomial terms only.
- **Small-sample precision in Tests 1–9.** Those used 5 noise draws per level, so gaps like
  4 of 5 vs 5 of 5 are within chance. Test 10 used 20.
- **Reproducibility across machines.** A fresh clone reproduced every checked result
  exactly, but on the same machine; see README, "Reproducing from scratch".

## BACON's home turf (Test 12)

BACON was built for static laws from controlled experiments: tables where each variable
is set to a few values in every combination, with no time and no derivatives. Test 12 ran
both methods on BACON's own five datasets (ideal gas, Ohm, Black, Kepler, synthetic), as
the authors give them and with more levels, 20 noise draws per level, both frozen
beforehand. Full tables: README, Test 12.

- **BACON found laws that PySINDy could not express at all.** PySINDy found none of the
  five laws in any run. For four of them that was certain: a sum in the denominator, or
  the power 2/3, is outside its candidate terms.
- **BACON is robust on two-quantity laws:** Kepler exact in 18, 13 and 5 of 20 runs at
  1%, 5% and 10% noise.
- **BACON breaks on multi-variable laws under even 1% noise:** Ohm and Black exact on clean
  data, then no law in any of the 120 noisy runs per design. Its BACON.3 layers require every
  slice of the table to agree on one expression.
- **The ideal gas law defeated both** with frozen settings, for the same reason: `T + 273`
  varies too little over the authors' range, and both methods' tolerances drop it.
- The authors' tuned settings recover some laws at 1% noise (ideal gas 6 of 20, Ohm 7 of
  20), not Black.

So neither study crowns a winner in general. On time series, PySINDy is better. On static
tables, BACON's search can express laws that a fixed list of terms cannot, but only its
two-quantity search survives noise.

## Possible next steps

- **An outside set of static laws** (e.g. the Feynman equations with up to 4 variables),
  chosen by a rule, so the laws are not BACON's own. Including laws inside PySINDy's terms
  would make it a contest rather than a test of expressibility.
- **A PySINDy that can express ratios:** PySINDy's implicit variant (SINDy-PI) fits
  `T_f·(M₁ + M₂) = M₁T₁ + M₂T₂`-type laws. That is the fair counterpart to BACON's search.
