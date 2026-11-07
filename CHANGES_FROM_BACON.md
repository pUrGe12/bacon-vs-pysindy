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

Default (since rerun 2): **B1, B3, B4 on, B2 off.** `--bacon-fixes none` gives the original
BACON.1 decisions. Tests 1–3 ran with B1–B3, Test 4 and rerun 1 with B1–B4.

| ID | Flag | Original BACON.1 | BACON-dt | Why | Since |
|---|---|---|---|---|---|
| B1 | `sign_aware` | Ratio if `r > 0`, product if `r < 0`, where `r` = correlation of `a`, `b` | Same rule on the correlation of `|a|`, `|b|` | For `Dx = -kx` the correlation is negative, so the original multiplies (`Dx·x`, `Dx·x²`, …) and never tries the ratio `x/Dx` | Test 1 |
| B2 | `signed_constancy` | `check_constant(b_, abs(b))`, then every `v` in `M(1-δ) < v < M(1+δ)`, i.e. on **absolute values** | Every `v` within `|v - M| < δ|M|`, on **signed values** | Stops a quantity that flips sign (e.g. ±2) from counting as constant. **Off by default since rerun 2:** the ablation showed no benefit, so we keep the original behaviour. *Correction:* Tests 1–4 described this as "the original fails when `M < 0`". That was wrong: the original passes `abs(b)`, so negative constants already worked. | Test 1 |
| B3 | `bounded_retry` | If nothing is found, build a new `BACON_1` with `ε × epsilon_scale`, `δ × delta_scale` and recurse, **without limit** | Stop once `δ > 0.3` and report "no invariant" | Without a cap, BACON always returns *something* once the tolerance is loose enough. With the cap off, a hard stop at 200 loosenings keeps runs from hanging. | Test 1 |
| B4 | `intercept_vs_std` | Take the linear branch if `1 - |r| < ε` and `|c / mean(b)| > c_val` | `... and |c| > c_val · std(b)` | Oscillating signals have mean ≈ 0, so any noise-level intercept passes and BACON writes `Dx - m·v = c` instead of the ratio `v/Dx` | Test 4 |

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

### Ablation 2: from the current defaults (B1, B3, B4, robust variant)

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

---

## S. PySINDy configuration (differences from PySINDy defaults)

| ID | PySINDy default | Ours | Since |
|---|---|---|---|
| S1 | Polynomial library degree 2 | Degree 3, STLSQ optimizer | Test 1 |
| S2 | STLSQ threshold 0.1 | Test 1: fixed 0.1. Test 2+: sweep 13 log-spaced thresholds 0.001–1 and keep the sparsest model whose R² is within 0.01 of the best (`--sindy-threshold`). Reason: a fixed 0.1 deleted the true `-0.1x²` logistic term. AIC/BIC were tried and kept tiny junk terms on clean data. | Test 2 |
| S3 | Differentiates the data itself (default `FiniteDifference`) | Given the same `x_dot` as BACON. `fd` equals PySINDy's default. | Test 1 |
| S4 | Multiple trajectories passed as a list | Rows concatenated and fitted with the given `x_dot` (same result, since no differentiation happens inside) | Test 4 |

## H. Harness and scoring history

| ID | Change | Since |
|---|---|---|
| H1 | Simulation solver: LSODA in Tests 1–2, RK45 from Test 3 (LSODA hung when our blow-up guard raised inside its Fortran callback). `--sim-solver` reproduces either. Only affects `sim_err`. | Test 3 |
| H2 | Simulation guard: give up after 20,000 evaluations or once the state exceeds 1000× the data range | Test 1 |
| H3 | "Recovered" = exactly the right terms and coefficient error < 5%. Unchanged throughout. | Test 1 |
| H4 | Noise: Gaussian, std = level × std of each state column. Unchanged throughout. | Test 1 |
| H5 | Every run writes `<name>_manifest.json`: command, settings, BACON fixes on/off, git commits of both repos (with an "uncommitted changes" flag), a hash of each source file, package versions | after Test 4 |

## Configuration of each test

| Test | System | BACON fixes | BACON variants | PySINDy threshold | Solver | Trajectories |
|---|---|---|---|---|---|---|
| 1 | decay | B1–B3 | strict | 0.1 | LSODA | 1 |
| 2 | logistic | B1–B3 | strict, robust | sweep | LSODA | 1 |
| 3 | oscillator | B1–B3 | strict, robust | sweep | RK45 | 1 |
| 4 | oscillator | B1–B4 | strict, robust | sweep | RK45 | 1 and 3 |
| rerun 1 | all | B1–B4 | strict, robust | sweep | RK45 | 1 (oscillator also 3) |
| rerun 2 (current defaults) | all | B1, B3, B4 | robust | sweep | RK45 | 1 (oscillator also 3) |

Exact commands for each are in `experiments.sh`. Reproducing Tests 1–4 with those commands
gives results identical to the original runs (448 of 448 rows; `summarize.py repro`), and
still does after the default changes of rerun 2, because every old command spells out its
settings.
