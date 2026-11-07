# BACON vs PySINDy on time-series data

Can BACON (Langley's heuristic law finder, as rebuilt in
[Miller & Banerjee 2025](https://www.researchgate.net/publication/389562521)) discover
ODEs if we give it a notion of derivatives? We compare it against PySINDy.

## Setup

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_test.py --system decay          # one system, current code
./experiments.sh run 14                              # every experiment, 14 in parallel
./experiments.sh run 8 rerun2                        # only one set (repro, rerun2, ablation2, ...)
```

BACON is imported from `~/Desktop/BACON` (set `BACON_PATH` to change). That repo is not modified.

## How BACON was adapted to time series

BACON looks for a combination of two variables that stays constant (e.g. Kepler's
`D³/P²`). We give it a **derivative column** `Dx = dx/dt`, computed from differences between
adjacent samples, and pair it with each state variable. If BACON finds an invariant such as
`x/Dx = -2`, we solve it for `Dx` and get the ODE `x' = -0.5x`.

- Only BACON's core law finder (BACON.1) is used. Its multi-variable layers (BACON.3/5,
  MCTS) assume controlled experiments where all variables but one are held fixed, which
  time series never are.
- Both methods get the **same derivative estimates** (`fd` raw, `sfd` smoothed), so the
  comparison is between their search methods.
- BACON.1 needed fixes for ODE data (negative derivatives, oscillating signals, a retry
  loop that always "finds" something). Each one is a named switch (`--bacon-fixes`), and
  `--bacon-fixes none` runs the original algorithm. Current defaults: B1, B3, B5, B6, B7.
- B5 also starts BACON's search from per-unit rates (`Dx_i/x_j`), because BACON.1's
  up/down heuristic never proposes them when a derivative rises and then falls.
- B6 makes BACON's straight-line test robust: a Theil–Sen slope, and the line must pass the
  same constancy test as every other law, with no leftover trend.
- B7 also pairs each variable's own per-capita rate (`Dx/x`) with every *other* variable,
  the form of interaction laws like Lotka–Volterra (`Dx/x = a - b·y`).
- The default BACON variant is **BACON-dt-robust**: a law counts as constant when the
  *typical* point is within tolerance (median-based), not *every* point. The strict
  variant (`--bacon-variants strict`) is BACON's original test; it only ever worked on
  clean data.

**[CHANGES_FROM_BACON.md](CHANGES_FROM_BACON.md)** lists every difference from the
original BACON and from PySINDy's defaults, with the test that introduced it and its
measured effect (an ablation over all systems).

Framing: BACON was never built for time series, so losing here is expected. Later we
will also adapt PySINDy to BACON's home turf (static laws like Kepler, ideal gas) to make
the comparison fair in both directions.

## Files

| File | What it does |
|---|---|
| `systems.py` | Ground-truth ODEs (decay, logistic, oscillator, power decay, Lotka–Volterra), clean trajectories, noise (std = level × std of each state) |
| `derivatives.py` | Derivative estimates given to **both** methods: `fd` (central differences, PySINDy's default) and `sfd` (Savitzky–Golay smoothed, window 11) |
| `bacon_dt.py` | BACON.1 subclass with derivative columns, switchable fixes B1–B7 (default B1, B3, B5–B7) and the robust variant (default) |
| `sindy_runner.py` | PySINDy: polynomial library degree 3, STLSQ, threshold chosen by a sweep |
| `scoring.py` | Term structure match, coefficient error, error when simulating the learned ODE |
| `run_test.py` | Runs one system across noise levels × seeds × derivative methods; writes results + manifest |
| `experiments.sh` | The exact command for every experiment in this README |
| `summarize.py` | Builds the tables below from the result files |
| `check_original.py` | Checks that "all fixes off" gives the same output as the original BACON.1 |

## Reproducibility

- Every run writes `results/.../<name>.csv` (one row per run, including which BACON fixes
  were on), `<name>_summary.csv` and `<name>_manifest.json`. The manifest records the
  command, all settings, git commits of this repo and of BACON (with an "uncommitted
  changes" flag), a hash of each source file, and package versions.
- `results/<system>.csv` and `results/oscillator_test4_*.csv` are the files written by the
  original runs of Tests 1–4. `experiments.sh` reproduces each test with its original
  settings into `results/repro/`. All 448 of 448 rows are identical:
  `.venv/bin/python summarize.py repro`. Every command in `experiments.sh` spells out its
  settings, so changing a default never changes an old experiment.

| Test | Config (details in CHANGES_FROM_BACON.md) | Reproduce |
|---|---|---|
| 1 | fixes B1–B3, strict only, PySINDy threshold 0.1, LSODA | `run_test.py --system decay --bacon-fixes sign_aware signed_constancy bounded_retry --bacon-variants strict --sindy-threshold 0.1 --sim-solver LSODA` |
| 2 | fixes B1–B3, threshold sweep, LSODA | `run_test.py --system logistic --bacon-fixes sign_aware signed_constancy bounded_retry --sim-solver LSODA` |
| 3 | fixes B1–B3, threshold sweep, RK45 | `run_test.py --system oscillator --bacon-fixes sign_aware signed_constancy bounded_retry` |
| 4 | fixes B1–B4, both variants | `run_test.py --system oscillator --bacon-fixes sign_aware signed_constancy bounded_retry intercept_vs_std --bacon-variants strict robust` (and `--n-traj 3`) |
| rerun 1 | as Test 4, every system | see `experiments.sh` (`results/rerun1/`) |
| rerun 2 | fixes B1, B3, B4, robust only | see `experiments.sh` (`results/rerun2/`) |
| 5, 6 (rerun 3) | fixes B1, B3, B4, B5, robust only | see `experiments.sh` (`results/rerun3/`) |
| 7 (rerun 4) | fixes B1, B3–B6, robust only | see `experiments.sh` (`results/rerun4/`) |
| 8 (rerun 5) | current defaults: fixes B1, B3, B5, B6, B7, robust only | `run_test.py --system <name>` (`results/rerun5/`); PySINDy side run: `--system powerdecay --sindy-library poly+halfpowers` |

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

> **Correction (after Test 4):** the fixes list written at the time said the original
> constancy test "is never true when `M < 0`". That was wrong: the original applies the test
> to `abs(b)`, so negative constants already worked. Fix B2 instead makes the test use
> signed values. See CHANGES_FROM_BACON.md, B2. The results below are unaffected.

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

### Test 4: Oscillator again, with both fixes from Test 3
Same system and noise as Test 3. Two changes: BACON fix 4 (intercept vs std, both BACON
variants), and an option to pool several trajectories. Run twice:
- **1 trajectory:** identical data to Test 3, so any change comes from the BACON fix
  (`results/oscillator_test4_1traj.csv`).
- **3 trajectories:** starting points (1, 0), (0.5, 1), (-1.5, 0.5), i.e. energies
  `x² + v²/4` = 1, 0.5, 2.31. Both methods get the same pooled 3× data
  (`results/oscillator_test4_3traj.csv`).

Exactly recovered / runs (in brackets: simulates within 10%). Both derivative types shown
for the robust variant; strict BACON-dt is 0/5 under noise in every configuration.

| noise | Test 3 robust (fd · sfd) | 1 traj robust (fd · sfd) | 1 traj PySINDy (fd · sfd) | 3 traj robust (fd · sfd) | 3 traj PySINDy (fd · sfd) |
|---|---|---|---|---|---|
| 0%  | 1/1 · 1/1 | 1/1 · 1/1 | 0/1 · 0/1 | 1/1 · 1/1 | **1/1 · 1/1** |
| 1%  | 0/5 · 0/5 | **5/5 · 5/5** | 5/5 · 5/5 | 5/5 · 5/5 | 5/5 · 5/5 |
| 5%  | 0/5 · 0/5 | 0/5 · **5/5** (3) | 3/5 · 5/5 | 0/5 · 5/5 (3) | 5/5 · 5/5 |
| 10% | 0/5 · 1/5 | 0/5 · **3/5** (0) | 2/5 · 2/5 | 0/5 · **5/5** (0) | 5/5 · 5/5 (0) |

**Takeaways**
- **The intercept fix was the main thing holding BACON back.** On the same data as Test 3,
  BACON-dt-robust with smoothed derivatives goes from 0/5 to 5/5 at 1% and 5% noise, and
  to 3/5 at 10%, ahead of PySINDy's 2/5 on one trajectory.
- **Multiple trajectories fix PySINDy's clean-data failure.** With 3 energy levels the
  library columns are no longer dependent, and PySINDy gets `x' = 0.998v`,
  `v' = -3.992x` on clean data, and 5/5 at every noise level.
- **With 3 trajectories and smoothed derivatives it's a tie:** both methods 5/5 at every
  noise level.
- **PySINDy still wins with raw derivatives** (fd): it gets 5/5 at 5% and 10% noise, where
  BACON-dt-robust finds nothing. Raw differences amplify noise, so the typical `v/Dx` point
  is too scattered for even the median test. Least squares averages it out.
- **Strict BACON-dt is 0/5 under noise in every configuration.** `x/Dv` blows up at every
  zero crossing, so the "every point within δ" test cannot pass. For signals that cross
  zero, the robust test is required.
- **The "simulates within 10%" column is a weak signal for oscillators.** A 2–3% frequency
  error becomes a large phase error after 3 periods, so at 10% noise no method passes it,
  not even the exactly recovered runs.
- Speed: PySINDy ~0.07 s per run (13-threshold sweep), BACON-dt-robust ~0.4 s.

**Next steps to try**
1. Rerun Tests 1–2 with the current code (fix 4, threshold sweep). Decay should now work for
   BACON-dt-robust; logistic still needs a search that tries `Dx/x`.
2. Derivative-aware search for BACON (always try `Dx_i/x_j` as a starting term), the
   fix that logistic growth needs.

### Rerun 1: every system with the Test 4 code
Config: BACON fixes B1–B4, both variants, PySINDy threshold sweep, RK45. Commands in
`experiments.sh`, results in `results/rerun1/`, tables from `summarize.py current rerun1`.
Cells: exactly recovered / runs, raw (fd) · smoothed (sfd) derivatives.

**decay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt | 1/1 · 1/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 4/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 4/5 · 5/5 | 2/5 · 3/5 |

**logistic** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt | 0/1 · 0/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| BACON-dt-robust | 0/1 · 0/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 4/5 · 5/5 | 1/5 · 3/5 | 1/5 · 3/5 |

**oscillator** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt | 1/1 · 1/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 3/5 |
| PySINDy | 0/1 · 0/1 | 5/5 · 5/5 | 3/5 · 5/5 | 2/5 · 2/5 |

**oscillator (3 trajectories)** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt | 1/1 · 1/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 5/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 5/5 · 5/5 | 5/5 · 5/5 |

**Takeaways**
- **Decay changed most.** BACON-dt-robust (which didn't exist in Test 1) now recovers
  `x' = -0.5x` in 5/5 runs at 1% noise with either derivative, and 4/5 at 5% (smoothed).
  The ablation shows this comes from the robust variant, not fix B4. PySINDy is still
  ahead at 5–10% noise.
- **Logistic is unchanged:** BACON finds nothing, PySINDy as before. The missing piece is
  still the search order (trying `Dx/x`), not the law test.
- **Oscillator results match Test 4 exactly,** as they should (same code and settings).
- **Strict BACON-dt only works on clean data.** Every noisy recovery comes from the robust
  variant, so the robust variant should be the default BACON-dt going forward.
- **Ablation** (CHANGES_FROM_BACON.md): B1 is essential everywhere; B4 is decisive for
  oscillations; B3 (retry cap) prevents most wrong laws (the original BACON returns a
  wrong law in every run on every system); B2 has almost no effect.

**Decision after rerun 1:** drop B2 (switched off by default) and make the robust variant
the only default BACON variant.

### Rerun 2: B2 off, robust variant only
Config: BACON fixes B1, B3, B4, BACON-dt-robust only, PySINDy threshold sweep, RK45.
`./experiments.sh run 8 rerun2`, results in `results/rerun2/`, tables from
`summarize.py current rerun2`. Cells: exactly recovered / runs, fd · sfd.

**decay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 4/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 4/5 · 5/5 | 2/5 · 3/5 |

**logistic** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 0/1 · 0/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 4/5 · 5/5 | 1/5 · 3/5 | 1/5 · 3/5 |

**oscillator** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 3/5 |
| PySINDy | 0/1 · 0/1 | 5/5 · 5/5 | 3/5 · 5/5 | 2/5 · 2/5 |

**oscillator (3 trajectories)** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 5/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 5/5 · 5/5 | 5/5 · 5/5 |

**Takeaways**
- **Same recoveries as rerun 1** on every system. Dropping B2 changed 4 of 128 runs, all
  at 5% noise with raw derivatives on the oscillator: where BACON used to find nothing, it
  now returns a loose law, one of them the right form with coefficients 6% off. Details in
  CHANGES_FROM_BACON.md, ablation 2.
- **Where things stand:**
  - **Decay:** BACON ties PySINDy up to 1% noise, then falls behind (4/5 vs 5/5 at 5% smoothed).
  - **Logistic:** BACON finds nothing.
  - **Oscillator, one trajectory:** BACON is ahead on clean data and at 10% noise (smoothed).
  - **Oscillator, three trajectories:** a tie with smoothed derivatives.
  - **Raw derivatives:** PySINDy is ahead at 5% and 10% noise on every system.

**Next:** logistic growth, which needs a derivative-aware search order (try `Dx/x`).

### Test 5: logistic growth fixed (B5 rate seed), every system rerun
New fix **B5 `rate_seed`**: for each derivative, BACON also starts its search from the
per-unit rate `Dx_i/x_j`, not only from `Dx_i`. Config: fixes B1, B3, B4, B5, robust variant,
PySINDy threshold sweep, RK45. `./experiments.sh run 8 rerun3`, results in `results/rerun3/`.
Cells: exactly recovered / runs, fd · sfd.

**decay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 4/5 · 5/5 | 2/5 · 3/5 |

**logistic** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 0/5 · 5/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 4/5 · 5/5 | 1/5 · 3/5 | 1/5 · 3/5 |

**oscillator** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 4/5 |
| PySINDy | 0/1 · 0/1 | 5/5 · 5/5 | 3/5 · 5/5 | 2/5 · 2/5 |

**oscillator (3 trajectories)** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 5/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 5/5 · 5/5 | 5/5 · 5/5 |

**Takeaways**
- **Logistic now works for BACON:** clean data and 5/5 at 1% noise (smoothed), as the
  invariant `Dx/x + 0.1x = 1`, i.e. exactly `x' = x - 0.1x²`. Before B5: nothing, ever.
- **It stops at 1% noise.** At 5% the `Dx/x` versus `x` data is too scattered for
  BACON.1's straight-line test, which needs `1 - |r|` below about 0.026 at the loosest
  tolerance. Measured `1 - |r|` at 5% noise: 0.07–0.28. Early on `x` is small (0.5), so
  dividing by it amplifies the noise. A robust line test would be the next fix.
- **B5 helped elsewhere too:** decay at 5% noise (smoothed) 4/5 → 5/5, oscillator at 10%
  (smoothed) 3/5 → 4/5. It made nothing worse (ablation 3 in CHANGES_FROM_BACON.md).
- PySINDy is unchanged; it still leads on logistic at 5–10% noise.

### Test 6: fractional power decay, `x' = -0.5x^1.5`
x0 = 4, t ∈ [0, 5] (x falls to 0.33), 101 samples, 5 noise seeds. Expected to be BACON's
home turf: its ratio search finds rational powers (like Kepler's `D³/P²`), while PySINDy's
polynomial term list cannot express `x^1.5`. As an informed side run, PySINDy also gets a
term list with half-integer powers (`--sindy-library poly+halfpowers`).

**powerdecay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 0/5 · 5/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 0/1 · 0/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |

**powerdecay (PySINDy with half powers)** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 0/5 · 5/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 3/5 · 4/5 | 0/5 · 1/5 |

(BACON's rows are identical in both tables; only PySINDy's term list differs.)

**Takeaways**
- **BACON finds the exact law** on clean data and at 1% noise (smoothed, 5/5), as the
  invariant `x³/Dx² = 4.0`, i.e. `x' = -0.5·√(x³)`. Same kind of law as Kepler's.
- **PySINDy with its default polynomial list can't, at any noise level.** It fits
  `x' = -0.32x - 0.18x²`, which simulates within 1–2% of the data but is not the law.
- **Given half powers, PySINDy wins with noise:** 5/5 at 1%, 3–4/5 at 5%, 0–1/5 at 10%,
  where BACON finds nothing above 1%. But it only works because we put `x^1.5` in its list,
  i.e. we told it the answer's form. BACON found the form on its own.
- **Why BACON stops at 1%:** the invariant `Dx²/x³` squares the derivative and divides by
  `x³`, which amplifies noise. Its typical spread is 10% of its value at 1% noise, 45% at
  5% and 71% at 10%, above BACON's 30% tolerance cap.
- **With raw derivatives at 1% noise BACON settles for a straight line,** `x' = 0.36 - 0.94x`,
  found at a loose tolerance (δ = 0.26). Counted wrong, but it simulates within 4%.

Also fixed for this test:
- **Scoring (H6):** `√(x³)` and `x^(3/2)` now count as the same term.
- **PySINDy's custom term list (S5):** written so PySINDy keeps all three half-power
  terms. My first version was silently dropped by PySINDy, because it counts each
  function's arguments.

**Check:** all 50 committed results files regenerate identically with the new code
(`summarize.py verify`), so neither B5 (off in old commands) nor the scoring change
altered any earlier result.

**Next steps to try**
1. A robust straight-line test for BACON (e.g. fit on the median trend, not every point):
   the current bottleneck on logistic and power decay above 1% noise.
2. Next system: Lotka–Volterra (two interacting species; `Dx/x` is a straight line in `y`).

### Test 7: robust straight-line test (B6), every system rerun
BACON.1 accepts a straight-line law `b - m·a = c` only if the least-squares correlation is
almost perfect (`1 - |r|` < 0.026 at its loosest), which noise breaks even when the line is
clearly there. New fix **B6 `robust_line`**:
1. Slope by **Theil–Sen** (median of all pairwise slopes; up to ~29% of points can be way
   off without moving it).
2. The leftover `b - m·a` must pass **the same constancy test as every other law**, at the
   current tolerance δ, so lines, ratios and products are held to one standard.
3. **No leftover trend:** split the data into thirds by `a`. If one third's median leftover
   is both significantly off (more than 3 standard errors) and more than δ/2 off, it's a
   curve, not a line. Without this, the first version accepted power decay's curved
   `x/Dx` as a line; without the δ/2 floor, it rejected true lines on clean data.

Config: fixes B1, B3, B4, B5, B6, robust variant, PySINDy threshold sweep, RK45.
`./experiments.sh run 8 rerun4`, results in `results/rerun4/`. Cells: exactly recovered /
runs, fd · sfd. (The PySINDy half-power side run on power decay is unchanged from Test 6.)

**decay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 3/5 · 5/5 | 0/5 · 4/5 | 0/5 · 0/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 4/5 · 5/5 | 2/5 · 3/5 |

**logistic** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 1/5 · 4/5 | 0/5 · 1/5 |
| PySINDy | 1/1 · 1/1 | 4/5 · 5/5 | 1/5 · 3/5 | 1/5 · 3/5 |

**oscillator** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 1/5 · 5/5 | 0/5 · 4/5 |
| PySINDy | 0/1 · 0/1 | 5/5 · 5/5 | 3/5 · 5/5 | 2/5 · 2/5 |

**oscillator (3 trajectories)** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 0/5 · 5/5 | 0/5 · 5/5 |
| PySINDy | 1/1 · 1/1 | 5/5 · 5/5 | 5/5 · 5/5 | 5/5 · 5/5 |

**powerdecay** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 0/5 · 4/5 | 0/5 · 0/5 | 0/5 · 0/5 |
| PySINDy | 0/1 · 0/1 | 0/5 · 0/5 | 0/5 · 0/5 | 0/5 · 0/5 |

**Takeaways**
- **Logistic jumped:** 5% noise 0/5 → 4/5 (smoothed) and 1/5 (raw), 10% noise 0/5 → 1/5
  (smoothed), and raw derivatives at 1% noise 0/5 → 5/5. BACON now matches or beats PySINDy
  on logistic at 0–5% noise. Every other logistic law it returns has the right form, with
  coefficients 5–8% off.
- **Oscillator improved slightly** (5% raw 0/5 → 1/5) and returns no wrong-form laws.
- **Decay and power decay lost a little.** Decay at 1% raw derivatives went 5/5 → 3/5 and at 5%
  smoothed 5/5 → 4/5; power decay at 1% smoothed went 5/5 → 4/5. In each case a false
  straight line (e.g. `0.55·Dx + x/Dx = -2.49`) passed at a tighter tolerance than the true
  law. At 5% noise power decay now returns wrong laws instead of nothing.
- **The trade:** on BACON-dt-robust over all 5 setups, exact recoveries 74 → 82 of 160 and
  wrong-form laws 8 → 11. Details in CHANGES_FROM_BACON.md, ablation 4.
- **With B6 on, B4 has no effect any more** (B6's constancy test already rejects lines with
  `c ≈ 0`). A candidate to drop, like B2.
- **Check:** all 73 committed results files regenerate identically (`summarize.py verify`).

**Next steps to try**
1. Stop false lines winning on tolerance alone: when a line and a ratio law both pass,
   prefer the one with fewer fitted constants (a ratio has one, a line two), or compare
   them on fit, not just on the δ at which they passed.
2. Drop B4 (no effect with B6).
3. Next system: Lotka–Volterra.

### Test 8: Lotka–Volterra (and B4 dropped)
Predator–prey: `x' = x - 0.1xy` (prey), `y' = -1.5y + 0.075xy` (predators). x0 = 10, y0 = 5,
t ∈ [0, 15] (~3 cycles), 301 samples, one trajectory, 5 noise seeds.

Two changes:
- **B4 dropped** from the defaults (no effect with B6; confirmed row for row on every
  older system).
- **New fix B7 `own_rate_seed`:** pair each variable's own per-capita rate `Dx_i/x_i` with
  every *other* variable. The prey law is `Dx/x = 1 - 0.1y`, a straight line in the
  *predator*; B5 only ever divides by the partner variable, so it never built that pair.

Config: fixes B1, B3, B5, B6, B7, robust variant, PySINDy threshold sweep, RK45.
`./experiments.sh run 8 rerun5`, results in `results/rerun5/`. The other systems are
unchanged from Test 7 (identical rows), so only Lotka–Volterra is shown. Cells: exactly
recovered / runs, fd · sfd.

**lotka_volterra** (fd · sfd)

| method | 0% | 1% | 5% | 10% |
|---|---|---|---|---|
| BACON-dt-robust | 1/1 · 1/1 | 5/5 · 5/5 | 1/5 · 5/5 | 0/5 · 5/5 |
| PySINDy | 1/1 · 0/1 | 4/5 · 3/5 | 1/5 · 0/5 | 0/5 · 0/5 |

**Takeaways**
- **BACON wins Lotka–Volterra clearly.** With smoothed derivatives it finds both laws
  exactly at every noise level, 10% included: `Dx/x + 0.100y = 1.0` and
  `Dy/y - 0.075x = -1.5`. With raw derivatives it holds up to 1% noise.
- **PySINDy gets the right terms every time but keeps junk.** Every miss is the true
  `x`, `xy`, `y` terms plus a small constant (0.03–0.25) or a small `y` term in the prey
  equation. Simulated, they stay within 0.2–21% of the data (median 4%), but fail the exact check. Likely cause: its
  threshold applies to raw coefficients, and with `x` reaching 40 the smallest true
  coefficient (0.075) is no bigger than the junk. PySINDy's `normalize_columns` option
  exists for exactly this; we haven't used it so far.
- **Ablation 5:** B7 is what makes it work (23 → 0 right without it), B6 is essential
  (23 → 10), and B7 changes nothing on any other system.
- **Check:** all 99 committed results files regenerate identically.

**Next steps to try**
1. **Fairness for PySINDy:** add `normalize_columns` as a PySINDy switch and rerun every
   system. Lotka–Volterra suggests our PySINDy setup, not PySINDy itself, may be what loses.
2. Stop false lines winning on tolerance alone (decay and power decay regressions from Test 7).
3. Next system: damped oscillator (`v' = -ω²x - γv`, three terms; BACON compares two at a time).
