# Changes from the original BACON

Every way this project differs from the original BACON code, and from PySINDy's
defaults, with the test that introduced each change and its measured effect.

**Original:** `~/Desktop/BACON` at commit `4dbdc957`, the rebuild of Langley's BACON by
Miller & Banerjee (2025). That repo is **not modified**. We import its `BACON_1` class
from `space_of_laws/laws_methods/bacon1.py` and subclass it in `bacon_dt.py`.

**Check:** `.venv/bin/python check_original.py` runs our class with every fix off and
BACON's original default tolerances, next to the original `BACON_1`, on 5 inputs (Kepler
plus the decay, logistic and oscillator data). All 5 give identical output, so "fixes off"
really is the original algorithm.

---

## What we use from the original, and what we don't

| Used, unchanged | Not used |
|---|---|
| `BACON_1`: the product / ratio / linear search, `new_term`, `product`, `division`, `linear`, `subs_expr` | `main.py`: it imports PySR (needs Julia) and GPy at startup, and `requirements.txt` lists neither |
| | The "space of data" layers (BACON.3 ranking, BACON.5, MCTS) and `data_space`: they assume factorial experiments (all variables but one held fixed), which time series never are |
| | BACON.6 (fits a user-supplied functional form), BACON.7, the PySR law method |

---

## A. Adapting BACON to time series (new code, not in BACON)

| ID | What | Since |
|---|---|---|
| A1 | **Derivative columns.** `D(x_i) = dx_i/dt` from differences between adjacent samples: `fd` = central differences (`pysindy.FiniteDifference(order=2)`), `sfd` = Savitzky–Golay smoothing (window 11, order 3), then differences. Two samples are dropped at each end. PySINDy gets the same derivatives. | Test 1 |
| A2 | **Pairing.** For each state `x_i`, run BACON.1 on a two-column table `(x_j, D(x_i))` for every state `x_j` (`x_j` is BACON's `nu`, `D(x_i)` its `eta`). BACON.1 then searches products, ratios and one linear step for a constant. | Test 1 |
| A3 | **Invariant → ODE.** BACON returns `expression = C`. We solve it for `D(x_i)` with sympy, drop complex roots, and per state keep the law found at the tightest tolerance δ (ties: best R² against `D(x_i)`). `C` is the mean (strict) or median (robust). | Test 1 |
| A4 | **Several trajectories** (`--n-traj`). Derivatives are computed per trajectory, then rows are pooled. Both methods get the same pool. | Test 4 |

Known limit of A2: each derivative is explained by one partner variable plus at most one
linear step, so laws needing sums of terms in different variables (e.g. a damped
oscillator, `v' = -ω²x - γv`) are out of reach unless the linear step happens to chain.

---

## B. Fixes to BACON.1 (switchable: `--bacon-fixes`)

Default (since rerun 5 / Test 8): **B1, B3, B5, B6, B7 on; B2 and B4 off.** `--bacon-fixes none`
gives the original BACON.1 decisions. Tests 1–3 ran with B1–B3, Test 4 and rerun 1 with
B1–B4, rerun 2 with B1, B3, B4, rerun 3 (Tests 5–6) with B1, B3, B4, B5, rerun 4 (Test 7)
with B1, B3–B6.

| ID | Flag | Original BACON.1 | BACON-dt | Why | Since |
|---|---|---|---|---|---|
| B1 | `sign_aware` | Ratio if `r > 0`, product if `r < 0`, where `r` = correlation of `a`, `b` | Same rule on the correlation of `|a|`, `|b|` | For `Dx = -kx` the correlation is negative, so the original multiplies (`Dx·x`, `Dx·x²`, …) and never tries the ratio `x/Dx` | Test 1 |
| B2 | `signed_constancy` | `check_constant(b_, abs(b))`, then every `v` in `M(1-δ) < v < M(1+δ)`, i.e. on **absolute values** | Every `v` within `|v - M| < δ|M|`, on **signed values** | Stops a quantity that flips sign (e.g. ±2) from counting as constant. **Off by default since rerun 2:** the ablation showed no benefit, so we keep the original behaviour. *Correction:* Tests 1–4 described this as "the original fails when `M < 0`". That was wrong: the original passes `abs(b)`, so negative constants already worked. | Test 1 |
| B3 | `bounded_retry` | If nothing is found, build a new `BACON_1` with `ε × epsilon_scale`, `δ × delta_scale` and recurse, **without limit** | Stop once `δ > 0.3` and report "no invariant" | Without a cap, BACON always returns *something* once the tolerance is loose enough. With the cap off, a hard stop at 200 loosenings keeps runs from hanging. | Test 1 |
| B4 | `intercept_vs_std` | Take the linear branch if `1 - |r| < ε` and `|c / mean(b)| > c_val` | `... and |c| > c_val · std(b)` | Oscillating signals have mean ≈ 0, so any noise-level intercept passes and BACON writes `Dx - m·v = c` instead of the ratio `v/Dx`. **Off by default since rerun 5:** with B6 on it changes nothing (ablation 4; rerun 5 confirms every result row is identical without it). | Test 4 |
| B5 | `rate_seed` | The search starts from the two input columns `(x_j, D(x_i))` only | Also run BACON.1 from `(x_j, D(x_i)/x_j)`, the per-unit rate (skipped where `x_j` has zeros). Implemented in `discover()` (the pairing layer, A2), not inside BACON.1; all candidates compete under the same rule (tightest δ, then best R²). | BACON.1 picks its next term from whether one variable goes up or down with the other. In logistic growth `Dx` rises, then falls as `x` grows, so BACON never proposes `Dx/x`, although `Dx/x = r - (r/K)x` is a straight line in `x`. Per-unit rates are a standard building block of ODE models, not a logistic-specific trick. | Test 5 |
| B6 | `robust_line` | Straight-line law `b - m·a = c`: least-squares slope `m`, accepted if `1 - |r| < ε` (ε ≤ 0.026 at the loosest tolerance) and the intercept is significant (B4). The leftover `b - m·a` is never tested for constancy. | Slope by **Theil–Sen** (median of pairwise slopes). Accepted if the leftover `b - m·a` passes **the same constancy test as every other law** at the current δ, **and** shows no leftover trend: split the data into thirds by `a`; a trend counts if some third's median leftover is more than 3 standard errors (z, from the leftover's own spread) from the overall median **and** more than δ/2 of it. Intercept check (B4) still applies. ε is no longer used. | With noise, `1 - |r|` collapses (logistic `Dx/x` vs `x` at 5% noise: 0.07–0.28, limit 0.026), mostly because noise is larger where `x` is small, although the line is clearly there. The trend check stops the robust test from accepting a curve as a line (first version without it turned power decay's `x/Dx` curve into a false line at 1% noise); the δ/2 floor stops it from rejecting true lines on noise-free data, where tiny finite-difference errors are statistically "significant". | Test 7 |
| B7 | `own_rate_seed` | (as B5) | Also run BACON.1 from `(x_j, D(x_i)/x_i)`: each variable's own per-capita rate, paired with every *other* variable `x_j` (skipped where `x_i` has zeros). In `discover()`, like B5. | In Lotka–Volterra the law is `Dx/x = a - b·y`: the prey's per-capita rate is a straight line in the *predator*. B5 divides by the partner (`Dx/y` vs `y`), so it never builds that pair. Per-capita rates against other species are the textbook form of interaction models. Has no effect on one-variable systems; on the oscillator, every result row is identical with and without it. | Test 8 |

### R. Robust variant (BACON-dt-robust, since Test 2; the only default variant since rerun 2)

| ID | Original | BACON-dt-robust | Why |
|---|---|---|---|
| R1 | Constant if **every** point is within δ of the **mean**; the constant is reported as the mean | Constant if the **median** absolute deviation from the **median** is below δ·|median|; the constant is reported as the median | With noisy derivatives, one bad point (often at a zero crossing) vetoes an otherwise clear law |

### P. Parameter choices

| ID | Original | Ours | Why | Since |
|---|---|---|---|---|
| P1 | `BACON_1` defaults: ε = 0.001, δ = 0.1, ε × 1.1 and δ × 1.05 per retry, c_val = 0.02 (the repo's args files use δ = 0.03–0.04) | ε = 0.001, δ starts at 0.01, both × 1.5 per retry, cap δ ≤ 0.3, c_val = 0.02 | Start strict and record the δ a law needed (a confidence signal); ×1.5 keeps runs fast | Test 1 |
| P2 | `check_constant` skips a raw input column unless its name is one of BACON's dummy letters a–n | Same. Tests 1–4 had dropped the "unless" part by mistake; restored after Test 4. No effect on any result: no variable was named a–n. | Faithfulness | restored after Test 4 |

---

## Measured effect of each fix (ablation)

### Ablation 1: from the rerun-1 settings (B1–B4, both variants)

Each fix switched off on its own, plus all off (= original BACON.1), on every system.
Made by `./experiments.sh run 8 ablation1`; tables from `.venv/bin/python summarize.py ablation1`.

Each cell: runs with the right law · runs that returned a wrong law, out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.


**BACON-dt**

| config | decay | logistic | oscillator | oscillator 3traj |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 2 right · 4 wrong | 0 right · 0 wrong | 2 right · 0 wrong | 2 right · 0 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without sign_aware | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong |
| without signed_constancy | 2 right · 4 wrong | 0 right · 0 wrong | 2 right · 0 wrong | 2 right · 0 wrong |
| without bounded_retry | 3 right · 29 wrong | 0 right · 32 wrong | 10 right · 22 wrong | 4 right · 28 wrong |
| without intercept_vs_std | 2 right · 3 wrong | 0 right · 0 wrong | 2 right · 13 wrong | 2 right · 14 wrong |

**BACON-dt-robust**

| config | decay | logistic | oscillator | oscillator 3traj |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 16 right · 1 wrong | 0 right · 0 wrong | 20 right · 3 wrong | 22 right · 0 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without sign_aware | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong |
| without signed_constancy | 16 right · 1 wrong | 0 right · 0 wrong | 20 right · 4 wrong | 22 right · 0 wrong |
| without bounded_retry | 16 right · 16 wrong | 0 right · 32 wrong | 22 right · 10 wrong | 22 right · 10 wrong |
| without intercept_vs_std | 16 right · 1 wrong | 0 right · 0 wrong | 3 right · 20 wrong | 2 right · 20 wrong |

What each fix is worth (ablation 1):
- **B1 is essential.** Without it nothing is found on any system, clean data included.
- **B2 has almost no effect:** without it, one extra wrong law in 128 robust runs
  (oscillator), nothing else changes. A candidate for removal, to stay closer to the original.
- **B3 trades a little recall for honesty.** Without the cap BACON always returns a law:
  slightly more right ones (oscillator strict 2 → 10), but far more wrong ones (0 → 22).
  The original BACON returns a wrong law in every single run.
- **B4 matters only for oscillating signals,** where it is decisive: robust on the
  oscillator goes from 3 right · 20 wrong without it to 20 right · 3 wrong with it.

### Ablation 2: from the rerun-2 defaults (B1, B3, B4, robust variant)

Made by `./experiments.sh run 12 ablation2`; tables from `.venv/bin/python summarize.py ablation2`.
"No fixes" reuses ablation 1; "plus B2" is rerun 1.

Each cell: runs with the right law · runs that returned a wrong law, out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.


**BACON-dt-robust**

| config | decay | logistic | oscillator | oscillator 3traj |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 16 right · 1 wrong | 0 right · 0 wrong | 20 right · 4 wrong | 22 right · 0 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without B1 sign_aware | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong | 0 right · 0 wrong |
| without B3 bounded_retry | 16 right · 16 wrong | 0 right · 32 wrong | 22 right · 10 wrong | 22 right · 10 wrong |
| without B4 intercept_vs_std | 16 right · 1 wrong | 0 right · 0 wrong | 4 right · 20 wrong | 2 right · 20 wrong |
| plus B2 signed_constancy (rerun 1) | 16 right · 1 wrong | 0 right · 0 wrong | 20 right · 3 wrong | 22 right · 0 wrong |

- **Dropping B2 changes 4 of 128 runs,** all at 5% noise with raw derivatives on the
  oscillator, all accepted only at a loose tolerance (δ = 0.26). Before, they returned
  nothing; now 3 return a law for `v'` only (`v' ≈ -4.1x`, right form) and 1 returns both
  laws with the right terms but coefficients 6% off (`x' = 1.059v`, `v' = -4.266x`), which
  counts as wrong under the 5% rule. No junk laws.
- B1, B3 and B4 matter exactly as in ablation 1.

### Ablation 3: from the rerun-3 defaults (B1, B3, B4, B5, robust variant)

Made by `./experiments.sh run 17 ablation3`; tables from `.venv/bin/python summarize.py ablation3`.
"Without B5" on the old systems is rerun 2; power decay gets its own "all off" and
"without B5" runs in `results/ablation3/`.

Each cell: runs with the right law · runs that returned a wrong law, out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.


**BACON-dt-robust**

| config | decay | logistic | oscillator | oscillator 3traj | powerdecay |
|---|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 17 right · 0 wrong | 7 right · 0 wrong | 21 right · 4 wrong | 22 right · 1 wrong | 7 right · 5 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without B1 sign_aware | 17 right · 0 wrong | 7 right · 0 wrong | 22 right · 1 wrong | 18 right · 4 wrong | 0 right · 12 wrong |
| without B3 bounded_retry | 20 right · 12 wrong | 10 right · 22 wrong | 27 right · 5 wrong | 30 right · 2 wrong | 7 right · 25 wrong |
| without B4 intercept_vs_std | 17 right · 0 wrong | 7 right · 0 wrong | 6 right · 19 wrong | 4 right · 20 wrong | 7 right · 5 wrong |
| without B5 rate_seed | 16 right · 1 wrong | 0 right · 0 wrong | 20 right · 4 wrong | 22 right · 0 wrong | 2 right · 10 wrong |

- **B5 helps everywhere and hurts nowhere.** Logistic 0 → 7 right; power decay 2 → 7 right
  and 10 → 5 wrong; decay 16 → 17 right and 1 → 0 wrong; oscillator 20 → 21 right.
- **With B5, B1 is only essential for power decay** (7 → 0 right without it). On decay,
  logistic and the oscillator, the rate seed hands BACON the right ratio directly, so the
  up/down heuristic matters less; on the oscillator with 3 trajectories B1 still helps
  (18 → 22 right).
- **B3 and B4 behave as in ablations 1–2:** without the cap, BACON returns many more wrong
  laws (logistic 0 → 22); without B4, the oscillator collapses (21 → 6 right).

### Ablation 4: from the rerun-4 defaults (B1, B3, B4, B5, B6, robust variant)

Made by `./experiments.sh run 20 ablation4`; tables from `.venv/bin/python summarize.py ablation4`.
"Without B6" is rerun 3; "no fixes" reuses ablations 1 and 3.

Each cell: runs with the right law · runs that returned a wrong law, out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.


**BACON-dt-robust**

| config | decay | logistic | oscillator | oscillator 3traj | powerdecay |
|---|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 14 right · 3 wrong | 18 right · 9 wrong | 22 right · 2 wrong | 22 right · 0 wrong | 6 right · 8 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without B1 sign_aware | 17 right · 1 wrong | 18 right · 9 wrong | 22 right · 1 wrong | 22 right · 0 wrong | 0 right · 7 wrong |
| without B3 bounded_retry | 16 right · 16 wrong | 18 right · 14 wrong | 29 right · 3 wrong | 31 right · 1 wrong | 6 right · 26 wrong |
| without B4 intercept_vs_std | 14 right · 3 wrong | 18 right · 9 wrong | 22 right · 2 wrong | 22 right · 0 wrong | 6 right · 8 wrong |
| without B5 rate_seed | 12 right · 5 wrong | 0 right · 0 wrong | 20 right · 4 wrong | 22 right · 0 wrong | 5 right · 9 wrong |
| without B6 robust_line (rerun 3) | 17 right · 0 wrong | 7 right · 0 wrong | 21 right · 4 wrong | 22 right · 1 wrong | 7 right · 5 wrong |

Splitting "wrong" into *right form, coefficients more than 5% off* and *wrong form*
(B6 off = rerun 3 → B6 on = rerun 4, out of 32 runs each):

| system | exact | right form, coef > 5% off | wrong form | nothing |
|---|---|---|---|---|
| decay | 17 → 14 | 0 → 0 | 0 → 3 | 15 → 15 |
| logistic | 7 → **18** | 0 → 9 | 0 → 0 | 25 → 5 |
| oscillator | 21 → 22 | 2 → 2 | 2 → 0 | 7 → 8 |
| oscillator, 3 traj | 22 → 22 | 0 → 0 | 1 → 0 | 9 → 10 |
| power decay | 7 → 6 | 0 → 0 | 5 → 8 | 20 → 18 |

- **B6 does what it was built for:** logistic goes from 7 to 18 exact, and every other
  logistic law it returns has the right form. The oscillator also loses its wrong-form laws.
- **Its cost is false lines at loose tolerance:** decay loses 3 exact runs to lines like
  `0.55·Dx + x/Dx = -2.49`, accepted at δ = 0.17 before the true `x/Dx` passes (BACON keeps
  the law found at the tightest δ). Power decay loses 1 exact run and returns wrong
  rational laws at 5% noise instead of nothing (they simulate within 3–5%).
- **With B6, B4 no longer changes anything:** an oscillating line `Dx - m·v = c` has `c ≈ 0`,
  so its leftover can never pass the constancy test relative to `|c|`. B4 still matters
  whenever B6 is off. A candidate to drop, as B2 was.
- **B1 matters only for power decay** (6 → 0 right without it), as in ablation 3.

### Ablation 5: from the current defaults (B1, B3, B5, B6, B7, robust variant)

Made by `./experiments.sh run 20 ablation5`; tables from `.venv/bin/python summarize.py ablation5`.
"Without B7" is the rerun-4 setting minus B4, so comparing it with rerun 4 (last row)
measures dropping B4: identical, row for row, on every system.

Each cell: runs with the right law · runs that returned a wrong law, out of 32 (clean + noisy, fd and sfd pooled). The rest returned nothing.


**BACON-dt-robust**

| config | decay | logistic | oscillator | oscillator 3traj | powerdecay | lotka_volterra |
|---|---|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 14 right · 3 wrong | 18 right · 9 wrong | 22 right · 2 wrong | 22 right · 0 wrong | 6 right · 8 wrong | 23 right · 1 wrong |
| no fixes (original BACON.1) | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong | 0 right · 32 wrong |
| without B1 sign_aware | 17 right · 1 wrong | 18 right · 9 wrong | 22 right · 1 wrong | 22 right · 0 wrong | 0 right · 7 wrong | 23 right · 1 wrong |
| without B3 bounded_retry | 16 right · 16 wrong | 18 right · 14 wrong | 29 right · 3 wrong | 31 right · 1 wrong | 6 right · 26 wrong | 26 right · 6 wrong |
| without B5 rate_seed | 12 right · 5 wrong | 0 right · 0 wrong | 20 right · 4 wrong | 22 right · 0 wrong | 5 right · 9 wrong | 23 right · 0 wrong |
| without B6 robust_line | 17 right · 0 wrong | 7 right · 0 wrong | 6 right · 19 wrong | 4 right · 20 wrong | 7 right · 5 wrong | 10 right · 0 wrong |
| without B7 own_rate_seed | 14 right · 3 wrong | 18 right · 9 wrong | 22 right · 2 wrong | 22 right · 0 wrong | 6 right · 8 wrong | 0 right · 4 wrong |
| B1, B3-B6 with B4 (rerun 4) | 14 right · 3 wrong | 18 right · 9 wrong | 22 right · 2 wrong | 22 right · 0 wrong | 6 right · 8 wrong | missing |

- **B7 is what makes Lotka–Volterra work:** 23 right without → 0 right · 4 wrong with it off.
  On every other system it changes nothing (identical rows).
- **Dropping B4 changes nothing** while B6 is on (identical rows on all 5 older setups).
- **B6 now also carries the oscillator:** without B6 *and* without B4, the oscillator
  collapses (22 → 6 right). B6 subsumed what B4 did.
- **B6 is essential for Lotka–Volterra** (23 → 10 right): both laws are straight lines.
- **B3** behaves as before: without the cap, more right laws but many more wrong ones.

---

## S. PySINDy configuration (differences from PySINDy defaults)

| ID | PySINDy default | Ours | Since |
|---|---|---|---|
| S1 | Polynomial library degree 2 | Degree 3, STLSQ optimizer | Test 1 |
| S2 | STLSQ threshold 0.1 | Test 1: fixed 0.1. Test 2+: sweep 13 log-spaced thresholds 0.001–1 and keep the sparsest model whose R² is within 0.01 of the best (`--sindy-threshold`). Reason: a fixed 0.1 deleted the true `-0.1x²` logistic term. AIC/BIC were tried and kept tiny junk terms on clean data. | Test 2 |
| S3 | Differentiates the data itself (default `FiniteDifference`) | Given the same `x_dot` as BACON. `fd` equals PySINDy's default. | Test 1 |
| S4 | Multiple trajectories passed as a list | Rows concatenated and fitted with the given `x_dot` (same result, since no differentiation happens inside) | Test 4 |
| S5 | — | Optional `--sindy-library poly+halfpowers`: polynomial degree 3 plus `|x|^0.5`, `|x|^1.5`, `|x|^2.5` per state. An **informed** baseline (it assumes we already suspect half-integer powers), used only as a labelled side run on power decay. Default stays polynomial. | Test 6 |

## H. Harness and scoring history

| ID | Change | Since |
|---|---|---|
| H1 | Simulation solver: LSODA in Tests 1–2, RK45 from Test 3 (LSODA hung when our blow-up guard raised inside its Fortran callback). `--sim-solver` reproduces either. Only affects `sim_err`. | Test 3 |
| H2 | Simulation guard: give up after 20,000 evaluations or once the state exceeds 1000× the data range | Test 1 |
| H3 | "Recovered" = exactly the right terms and coefficient error < 5%. Unchanged throughout. | Test 1 |
| H4 | Noise: Gaussian, std = level × std of each state column. Unchanged throughout. | Test 1 |
| H5 | Every run writes `<name>_manifest.json`: command, settings, BACON fixes on/off, git commits of both repos (with an "uncommitted changes" flag), a hash of each source file, package versions | after Test 4 |
| H6 | Term matching merges nested powers as if variables were positive, so `sqrt(x**3)` (how BACON writes it) and `x**(3/2)` count as the same term. Checked: no committed result changed. | Test 6 |

## Configuration of each test

| Test | System | BACON fixes | BACON variants | PySINDy threshold | Solver | Trajectories |
|---|---|---|---|---|---|---|
| 1 | decay | B1–B3 | strict | 0.1 | LSODA | 1 |
| 2 | logistic | B1–B3 | strict, robust | sweep | LSODA | 1 |
| 3 | oscillator | B1–B3 | strict, robust | sweep | RK45 | 1 |
| 4 | oscillator | B1–B4 | strict, robust | sweep | RK45 | 1 and 3 |
| rerun 1 | all | B1–B4 | strict, robust | sweep | RK45 | 1 (oscillator also 3) |
| rerun 2 | all | B1, B3, B4 | robust | sweep | RK45 | 1 (oscillator also 3) |
| rerun 3 = Tests 5 & 6 | all + power decay | B1, B3, B4, B5 | robust | sweep (+ half-power side run) | RK45 | 1 (oscillator also 3) |
| rerun 4 = Test 7 | all + power decay | B1, B3, B4, B5, B6 | robust | sweep (+ half-power side run) | RK45 | 1 (oscillator also 3) |
| rerun 5 = Test 8 (current defaults) | all + power decay + Lotka–Volterra | B1, B3, B5, B6, B7 | robust | sweep (+ half-power side run) | RK45 | 1 (oscillator also 3) |

Exact commands for each are in `experiments.sh`. Reproducing Tests 1–4 with those commands
gives results identical to the original runs (448 of 448 rows; `summarize.py repro`), and
still does after the default changes of reruns 2 and 3, because every old command spells
out its settings. `summarize.py verify` regenerates-and-compares every committed results
file: after Test 6, all 50 were identical (timing columns aside); after Test 7, all 73;
after Test 8, all 99.
