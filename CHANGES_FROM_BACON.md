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
| **Test 12 on (static laws):** the whole original pipeline, unchanged: `data_space`, the BACON.3 layer ranking, `BACON_1`, and `utils.losses` to combine the invariants and solve for the measured quantity; the datasets in `datasets.py`; the authors' args files (side run). None of the fixes B1–B7 is used there. | Test 12: `main.py` itself (we call the same functions it calls; see H8) |

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

Each cell counts runs that found the exact law ("right") and runs that returned some other law ("wrong"); the remaining runs at that noise level found no law for at least one variable.

<!-- tables:ablation1 -->

**Exponential decay, `x' = -0.5x` — BACON-dt (strict)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 1 right, 4 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 3 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Exponential decay, `x' = -0.5x` — BACON-dt-robust**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²` — BACON-dt (strict)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²` — BACON-dt-robust**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x` — BACON-dt (strict)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 3 wrong | 0 right, 3 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 2 right, 3 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 3 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x` — BACON-dt-robust**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 1 right, 4 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 2 right, 3 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 1 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled — BACON-dt (strict)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 2 right, 3 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 3 wrong | 0 right, 3 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled — BACON-dt-robust**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| all fixes B1-B4 (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B2 signed_constancy | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

<!-- /tables -->

What each fix is worth (ablation 1):
- **B1 is essential.** Without it nothing is found on any system, clean data included.
- **B2 has almost no effect:** without it, exactly one cell changes (robust, oscillator,
  raw derivatives, 5% noise: 1 wrong → 2 wrong). A candidate for removal, to stay closer
  to the original.
- **B3 trades a little recall for honesty.** Without the cap BACON returns a law in almost
  every run. A few more are right (strict, oscillator, smoothed, 1% noise: 0 → 5 of 5
  right), but most extra laws are wrong: at 5% and 10% noise, cells go from "0 right,
  0 wrong" to "0 right, 5 wrong" on every system. The original BACON returns a wrong law
  in every single run.
- **B4 matters only for oscillating signals,** where it is decisive. Robust, oscillator,
  smoothed derivatives without → with B4: 1% noise 0 → 5 right, 5% noise 0 → 5 right,
  10% noise 1 → 3 right; with three trajectories, 0 → 5 right at every noise level.

### Ablation 2: from the rerun-2 defaults (B1, B3, B4, robust variant)

Made by `./experiments.sh run 12 ablation2`; tables from `.venv/bin/python summarize.py ablation2`.
"No fixes" reuses ablation 1; "plus B2" is rerun 1.

Each cell counts runs that found the exact law ("right") and runs that returned some other law ("wrong"); the remaining runs at that noise level found no law for at least one variable.

<!-- tables:ablation2 -->

**Exponential decay, `x' = -0.5x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 2 right, 3 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 2 right, 3 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4 (rerun 2) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| plus B2 signed_constancy (rerun 1) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

<!-- /tables -->

- **Dropping B2 changes 4 runs, all in one cell** (oscillator, raw derivatives, 5% noise), all accepted only at a loose tolerance (δ = 0.26). Before, they returned
  nothing; now 3 return a law for `v'` only (`v' ≈ -4.1x`, right form) and 1 returns both
  laws with the right terms but coefficients 6% off (`x' = 1.059v`, `v' = -4.266x`), which
  counts as wrong under the 5% rule. No junk laws.
- B1, B3 and B4 matter exactly as in ablation 1.

### Ablation 3: from the rerun-3 defaults (B1, B3, B4, B5, robust variant)

Made by `./experiments.sh run 17 ablation3`; tables from `.venv/bin/python summarize.py ablation3`.
"Without B5" on the old systems is rerun 2; power decay gets its own "all off" and
"without B5" runs in `results/ablation3/`.

Each cell counts runs that found the exact law ("right") and runs that returned some other law ("wrong"); the remaining runs at that noise level found no law for at least one variable.

<!-- tables:ablation3 -->

**Exponential decay, `x' = -0.5x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 2 right, 3 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 4 right, 1 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 2 right, 3 wrong | 4 right, 1 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 2 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 4 right, 1 wrong | 4 right, 1 wrong | 4 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 2 right, 3 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Power decay, `x' = -0.5x^1.5` (PySINDy: plain polynomial terms)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B4, B5 (rerun 3) | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

<!-- /tables -->

- **B5 helps and hurts nowhere** (without → with B5, smoothed unless stated). Logistic:
  clean data 0 → 1 of 1 (both derivative types) and 1% noise 0 → 5 of 5. Power decay, 1%
  noise: 0 right, 5 wrong → 5 right, 0 wrong. Decay, 5% noise: 4 → 5 right. Oscillator,
  10% noise: 3 → 4 right. The only other change is one extra wrong law (oscillator, raw, 5%).
- **With B5, B1 is only essential for power decay:** without it, clean data 1 → 0 of 1
  and 1% noise (smoothed) 5 → 0 of 5, all returned as wrong laws. On decay and logistic
  it changes nothing, because the rate seed hands BACON the right ratio directly. With 3
  oscillator trajectories it still helps (smoothed, 1/5/10% noise: 4 → 5 right each).
- **B3 and B4 behave as in ablations 1–2:** without the cap, high-noise cells turn into
  "5 wrong" (e.g. logistic, raw, 5% and 10% noise); without B4 the oscillator collapses
  (smoothed, 1% and 5% noise: 5 → 0 right).

### Ablation 4: from the rerun-4 defaults (B1, B3, B4, B5, B6, robust variant)

Made by `./experiments.sh run 20 ablation4`; tables from `.venv/bin/python summarize.py ablation4`.
"Without B6" is rerun 3; "no fixes" reuses ablations 1 and 3.

Each cell counts runs that found the exact law ("right") and runs that returned some other law ("wrong"); the remaining runs at that noise level found no law for at least one variable.

<!-- tables:ablation4 -->

**Exponential decay, `x' = -0.5x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 2 right, 3 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 2 right, 3 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 4 right, 1 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 3 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 5 right, 0 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong | 0 right, 0 wrong |

**Power decay, `x' = -0.5x^1.5` (PySINDy: plain polynomial terms)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3-B6 (rerun 4) | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B4 intercept_vs_std | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line (rerun 3) | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

<!-- /tables -->

What B6 changed, run by run (rerun 3 = B6 off, rerun 4 = B6 on, BACON-dt-robust). Each cell: runs with the exact law / the right terms but a coefficient more than 5% off / wrong terms / no law.

<!-- tables:ablation4_b6 -->

**Exponential decay, `x' = -0.5x`**

Smoothed derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 4 / 0 / 1 / 0 | 0 / 0 / 0 / 5 |

Raw derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 3 / 0 / 2 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |

**Logistic growth, `x' = x - 0.1x²`**

Smoothed derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 4 / 1 / 0 / 0 | 1 / 4 / 0 / 0 |

Raw derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 1 / 4 / 0 / 0 | 0 / 0 / 0 / 5 |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`**

Smoothed derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 4 / 1 / 0 / 0 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 4 / 1 / 0 / 0 |

Raw derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 1 / 2 / 2 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 1 / 1 / 0 / 3 | 0 / 0 / 0 / 5 |

**Harmonic oscillator, three trajectories pooled**

Smoothed derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 |

Raw derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 1 / 4 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |

**Power decay, `x' = -0.5x^1.5`**

Smoothed derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 5 / 0 / 0 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 4 / 0 / 1 / 0 | 0 / 0 / 3 / 2 | 0 / 0 / 0 / 5 |

Raw derivatives: exact law / right terms but a coefficient more than 5% off / wrong terms / no law

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| B6 off (rerun 3) | 1 / 0 / 0 / 0 | 0 / 0 / 5 / 0 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |
| B6 on (rerun 4) | 1 / 0 / 0 / 0 | 0 / 0 / 4 / 1 | 0 / 0 / 0 / 5 | 0 / 0 / 0 / 5 |

<!-- /tables -->

- **B6 does what it was built for:** logistic with smoothed derivatives goes 0 → 4 of 5
  exact at 5% noise and 0 → 1 of 5 at 10%; with raw derivatives 0 → 5 of 5 at 1% noise.
  Every other logistic law it returns has the right terms. The oscillator also loses its wrong-form laws.
- **Its cost is false lines at loose tolerance:** decay loses 3 exact runs (raw 1% noise
  5 → 3, smoothed 5% noise 5 → 4) to lines like
  `0.55·Dx + x/Dx = -2.49`, accepted at δ = 0.17 before the true `x/Dx` passes (BACON keeps
  the law found at the tightest δ). Power decay loses 1 exact run (smoothed, 1% noise, 5 → 4) and returns wrong
  rational laws at 5% noise instead of nothing (they simulate within 3–5%).
- **With B6, B4 no longer changes anything:** an oscillating line `Dx - m·v = c` has `c ≈ 0`,
  so its leftover can never pass the constancy test relative to `|c|`. B4 still matters
  whenever B6 is off. A candidate to drop, as B2 was.
- **B1 matters only for power decay:** without it, clean data 1 → 0 of 1 and 1% noise
  (smoothed) 4 → 0 of 5, as in ablation 3. On decay, switching B1 *off* now slightly helps
  (raw 1% noise 3 → 5 right, smoothed 5% noise 4 → 5 right), because it avoids the same
  false lines that cost B6 those runs.

### Ablation 5: from the current defaults (B1, B3, B5, B6, B7, robust variant)

Made by `./experiments.sh run 20 ablation5`; tables from `.venv/bin/python summarize.py ablation5`.
"Without B7" is the rerun-4 setting minus B4, so comparing it with rerun 4 (last row)
measures dropping B4: identical, row for row, on every system.

Each cell counts runs that found the exact law ("right") and runs that returned some other law ("wrong"); the remaining runs at that noise level found no law for at least one variable.

<!-- tables:ablation5 -->

**Exponential decay, `x' = -0.5x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 2 right, 3 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 2 right, 3 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Logistic growth, `x' = x - 0.1x²`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 1 right, 4 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 5 wrong |
| without B5 rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 4 wrong | 0 right, 0 wrong |

**Harmonic oscillator, one trajectory, `x' = v`, `v' = -4x`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 3 right, 2 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 4 right, 1 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 4 right, 1 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |

**Harmonic oscillator, three trajectories pooled**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 2 right, 3 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 5 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 2 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Power decay, `x' = -0.5x^1.5` (PySINDy: plain polynomial terms)**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 3 right, 2 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 3 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 4 right, 1 wrong | 0 right, 3 wrong | 0 right, 0 wrong |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 0 right, 1 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 0 right, 5 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | 1 right, 0 wrong | 0 right, 4 wrong | 0 right, 0 wrong | 0 right, 0 wrong |

**Lotka–Volterra, `x' = x - 0.1xy`, `y' = -1.5y + 0.075xy`**

Smoothed derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong | 5 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 5 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 2 wrong | 0 right, 2 wrong |
| B1, B3-B6 with B4 (rerun 4) | not run | not run | not run | not run |

Raw derivatives: runs that found the exact law (right) or another law (wrong)

| config | no noise (1 run) | 1% noise (5 runs) | 5% noise (5 runs) | 10% noise (5 runs) |
|---|---|---|---|---|
| defaults B1, B3, B5, B6, B7 (rerun 5) | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| no fixes (original BACON.1) | 0 right, 1 wrong | 0 right, 5 wrong | 0 right, 5 wrong | 0 right, 5 wrong |
| without B1 sign_aware | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 1 wrong | 0 right, 0 wrong |
| without B3 bounded_retry | 1 right, 0 wrong | 5 right, 0 wrong | 4 right, 1 wrong | 0 right, 5 wrong |
| without B5 rate_seed | 1 right, 0 wrong | 5 right, 0 wrong | 1 right, 0 wrong | 0 right, 0 wrong |
| without B6 robust_line | 1 right, 0 wrong | 3 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| without B7 own_rate_seed | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong | 0 right, 0 wrong |
| B1, B3-B6 with B4 (rerun 4) | not run | not run | not run | not run |

<!-- /tables -->

- **B7 is what makes Lotka–Volterra work.** With B7, smoothed derivatives: 5 of 5 right at
  1%, 5% and 10% noise. Without it: 0 right at every noise level and derivative type
  (the rest find nothing, or 2 wrong laws at 5% and 10%). On every other system B7
  changes nothing (identical rows).
- **Dropping B4 changes nothing** while B6 is on (identical rows on all 5 older setups).
- **B6 now also carries the oscillator:** without B6 *and* without B4, the oscillator
  collapses (smoothed, 1% and 5% noise: 5 → 0 right, all wrong laws). B6 subsumed what B4 did.
- **B6 is essential for Lotka–Volterra** (smoothed, 5% and 10% noise: 5 → 0 right):
  both laws are straight lines.
- **B3** behaves as before: without the cap, more right laws but many more wrong ones.

---

## S. PySINDy configuration (differences from PySINDy defaults)

| ID | PySINDy default | Ours | Since |
|---|---|---|---|
| S1 | Polynomial library degree 2 | Degree 3, STLSQ optimizer | Test 1 |
| S2 | STLSQ threshold 0.1 | Test 1: fixed 0.1. Test 2+: sweep 13 log-spaced thresholds 0.001–1 and keep the sparsest model whose R² is within 0.01 of the best (`--sindy-threshold`). Reason: a fixed 0.1 deleted the true `-0.1x²` logistic term. AIC/BIC were tried and kept tiny junk terms on clean data. | Test 2 |
| S3 | Differentiates the data itself (default `FiniteDifference`) | Given the same `x_dot` as BACON. `fd` equals PySINDy's default. | Test 1 |
| S4 | Multiple trajectories passed as a list | Rows concatenated and fitted with the given `x_dot` (same result, since no differentiation happens inside) | Test 4 |
| S6 | `normalize_columns=False` | Optional `--sindy-normalize on`: each library column divided by its L2 length before thresholding. The threshold then applies to coefficient × column length, so the sweep grid becomes 25 log-spaced thresholds 0.001–1000 (same density as S2). Default **on** since Test 9 (off before). Effect: Test 9a (with the default ridge penalty, PySINDy collapses: 0/32 on logistic and Lotka–Volterra) and 9b. | Test 9 |
| S7 | STLSQ ridge penalty `alpha=0.05` | `--sindy-alpha`, default **0** since Test 9 (0.05 before). Test 9b uses 0 (plain least squares inside STLSQ, as in the original SINDy algorithm). Reason: STLSQ thresholds the ridge coefficients; on normalized columns 0.05 is a 5% penalty that keeps correlated junk terms alive. **Chosen after seeing 9a.** Effect with S6, smoothed derivatives (runs with the exact law, out of 5 per noise level, Test 8 setup → 9b): logistic 5% noise 3 → 5, 10% 3 → 4; one-trajectory oscillator 10% 2 → 5; Lotka–Volterra 1% 3 → 5, 5% 0 → 4, 10% 0 → 0; half-power side run 1% 5 → 0, 5% 4 → 0 (cancelling coefficients); decay, three-trajectory oscillator and plain power decay unchanged. Full tables: README, Test 9. | Test 9 |
| S8 | Fits `x' = f(x)` from time series | Static laws (Test 12, `run_static.py`): STLSQ (as S2, S6, S7) fitting the measured quantity directly from the inputs, no derivatives. Candidate terms: every product of integer powers of the inputs with total degree ≤ 3, negative powers included (7 terms for one input, 63 for three, 129 for four), built as a matrix and passed to `ps.STLSQ`. Only the ideal gas law among BACON's datasets can be written in these terms; stated before running. | Test 12 |
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
| H8 | Static laws (Test 12): BACON's own pipeline called from `bacon_static.py` with the same calls `main.py` makes. Its final `df_helper.score` (which prints the law) is swapped for a function doing the same two computations and keeping the result. PySR and GPy, imported at module level but unused on this path, get empty stand-in modules if not installed. A raised error or more than 300 s = "no law". Noise as in BACON's `datasets.py` (Gaussian on the measured quantity, sd = level × value), **not** H4. "Exact law" compares the numerator and denominator terms of both laws over a common denominator, coefficients within 5% (H3); answers that are plain sums of power products (all of PySINDy's) are put over their common denominator directly instead of by sympy's `cancel`, which gives the same verdict and is much faster. Data: `static_laws.py`; the given design is checked against `datasets.py`. | Test 12 |
| H7 | ODEBench setups (Test 10): equations, constants and starting points from ODEBench; time grid t = 0..10 with 150 samples (ODEBench's own); each starting point run as its own one-trajectory setup; 20 noise draws per level. Clean data simulated with our usual LSODA at rtol 1e-10 (ODEBench uses rtol 1e-5). Manifests also hash `data/*.json`. | Test 10 |

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
| rerun 5 = Test 8 | all + power decay + Lotka–Volterra | B1, B3, B5, B6, B7 | robust | sweep (+ half-power side run) | RK45 | 1 (oscillator also 3) |
| rerun 6 = Test 9a | as rerun 5 | B1, B3, B5, B6, B7 | robust | sweep 0.001–1000, normalized columns (S6) | RK45 | 1 (oscillator also 3) |
| rerun 6b = Test 9b (current defaults) | as rerun 5 | B1, B3, B5, B6, B7 | robust | as 9a, ridge penalty 0 (S7) | RK45 | 1 (oscillator also 3) |
| Test 10 | ODEBench: 63 systems × 2 starting points, 20 noise draws | B1, B3, B5, B6, B7 | robust | as 9b | RK45 | 1 per setup |
| Test 12 | BACON's static datasets (ideal, Ohm, Black, Kepler, synthetic) × given and dense designs, 20 noise draws | none: original BACON.3 + BACON.1 pipeline at code defaults (+ authors' args files as side run) | — | as 9b, static terms (S8) | — | — |

Exact commands for each are in `experiments.sh`. Reproducing Tests 1–4 with those commands
gives results identical to the original runs (448 of 448 rows; `summarize.py repro`), and
still does after the default changes of reruns 2 and 3, because every old command spells
out its settings. `summarize.py verify` regenerates-and-compares every committed results
file: after Test 6, all 50 were identical (timing columns aside); after Test 7, all 73;
after Test 8, all 99.
