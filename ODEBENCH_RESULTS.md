# Test 10: ODEBench results, setup by setup

All 63 ODEBench systems, each from its two starting points (126 setups). Rules and
summary: README, Test 10. Generated from `results/test10/` by
`.venv/bin/python summarize.py tables`; do not edit the tables by hand.

How to read: each cell counts runs that found the exact law (exactly the right terms,
every coefficient within 5%). With no noise there is 1 run per setup; at 1%, 5% and 10%
noise there are 20 runs, one per noise draw. "True law in PySINDy's terms?" says whether
the true law is a polynomial of degree ≤ 3, i.e. something PySINDy's term list can
express at all. Smoothed and raw derivatives are in separate sections.

<!-- tables:odebench_full -->

## Smoothed derivatives

### Smoothed derivatives: systems with 1 variable (46 setups)

Smoothed derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #1 RC-circuit (charging capacitor), start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 18 of 20 | 0 of 20 | 5 of 20 | 0 of 20 | 0 of 20 |
| #1 RC-circuit (charging capacitor), start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 19 of 20 | 8 of 20 | 8 of 20 | 0 of 20 | 1 of 20 |
| #2 Population growth (naive), start 1 | yes | 1 of 1 | 1 of 1 | 16 of 20 | 20 of 20 | 6 of 20 | 18 of 20 | 0 of 20 | 11 of 20 |
| #2 Population growth (naive), start 2 | yes | 1 of 1 | 1 of 1 | 16 of 20 | 20 of 20 | 6 of 20 | 18 of 20 | 0 of 20 | 11 of 20 |
| #3 Population growth with carrying capacity, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 12 of 20 | 12 of 20 |
| #3 Population growth with carrying capacity, start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 12 of 20 | 14 of 20 |
| #4 RC-circuit with non-linear resistor (charging capacitor), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #4 RC-circuit with non-linear resistor (charging capacitor), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #5 Velocity of a falling object with air resistance, start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 16 of 20 | 0 of 20 | 6 of 20 | 0 of 20 | 4 of 20 |
| #5 Velocity of a falling object with air resistance, start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #6 Autocatalysis with one fixed abundant chemical, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 1 of 20 | 20 of 20 | 0 of 20 | 0 of 20 |
| #6 Autocatalysis with one fixed abundant chemical, start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 19 of 20 | 0 of 20 | 10 of 20 | 0 of 20 | 0 of 20 |
| #7 Gompertz law for tumor growth, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #7 Gompertz law for tumor growth, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #8 Logistic equation with Allee effect, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #8 Logistic equation with Allee effect, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #9 Language death model for two languages, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 19 of 20 | 12 of 20 | 15 of 20 | 8 of 20 | 6 of 20 |
| #9 Language death model for two languages, start 2 | yes | 1 of 1 | 1 of 1 | 17 of 20 | 0 of 20 | 9 of 20 | 0 of 20 | 6 of 20 | 0 of 20 |
| #10 Refined language death model for two languages, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #10 Refined language death model for two languages, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #11 Naive critical slowing down (statistical mechanics), start 1 | yes | 1 of 1 | 0 of 1 | 5 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #11 Naive critical slowing down (statistical mechanics), start 2 | yes | 1 of 1 | 0 of 1 | 13 of 20 | 11 of 20 | 0 of 20 | 11 of 20 | 0 of 20 | 11 of 20 |
| #12 Photons in a laser (simple), start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 16 of 20 | 0 of 20 | 4 of 20 | 0 of 20 | 0 of 20 |
| #12 Photons in a laser (simple), start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 7 of 20 | 19 of 20 | 0 of 20 | 2 of 20 |
| #13 Overdamped bead on a rotating hoop, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #13 Overdamped bead on a rotating hoop, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #14 Budworm outbreak model with predation, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #14 Budworm outbreak model with predation, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #15 Budworm outbreak with predation (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #15 Budworm outbreak with predation (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #16 Landau equation (typical time scale tau = 1), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #16 Landau equation (typical time scale tau = 1), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #17 Logistic equation with harvesting/fishing, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #17 Logistic equation with harvesting/fishing, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #18 Improved logistic equation with harvesting/fishing, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #18 Improved logistic equation with harvesting/fishing, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #19 Improved logistic equation with harvesting/fishing (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #19 Improved logistic equation with harvesting/fishing (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #20 Autocatalytic gene switching (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #20 Autocatalytic gene switching (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #21 Dimensionally reduced SIR infection model for dead people (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #21 Dimensionally reduced SIR infection model for dead people (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #22 Hysteretic activation of a protein expression (positive feedback, basal promoter expression), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #22 Hysteretic activation of a protein expression (positive feedback, basal promoter expression), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #23 Overdamped pendulum with constant driving torque/fireflies/Josephson junction (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #23 Overdamped pendulum with constant driving torque/fireflies/Josephson junction (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Smoothed derivatives: systems with 2 variables (56 setups)

Smoothed derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #24 Harmonic oscillator without damping, start 1 | yes | 1 of 1 | 0 of 1 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 18 of 20 |
| #24 Harmonic oscillator without damping, start 2 | yes | 1 of 1 | 0 of 1 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 17 of 20 |
| #25 Harmonic oscillator with damping, start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 |
| #25 Harmonic oscillator with damping, start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 |
| #26 Lotka-Volterra competition model (Strogatz version with sheeps and rabbits), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #26 Lotka-Volterra competition model (Strogatz version with sheeps and rabbits), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #27 Lotka-Volterra simple (as on Wikipedia), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #27 Lotka-Volterra simple (as on Wikipedia), start 2 | yes | 1 of 1 | 0 of 1 | 20 of 20 | 0 of 20 | 19 of 20 | 0 of 20 | 10 of 20 | 0 of 20 |
| #28 Pendulum without friction, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #28 Pendulum without friction, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #29 Dipole fixed point, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #29 Dipole fixed point, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #30 RNA molecules catalyzing each others replication, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #30 RNA molecules catalyzing each others replication, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #31 SIR infection model only for healthy and sick, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 4 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #31 SIR infection model only for healthy and sick, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #32 Damped double well oscillator, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 3 of 20 | 0 of 20 | 0 of 20 |
| #32 Damped double well oscillator, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #33 Glider (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #33 Glider (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #34 Frictionless bead on a rotating hoop (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #34 Frictionless bead on a rotating hoop (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #35 Rotational dynamics of an object in a shear flow, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #35 Rotational dynamics of an object in a shear flow, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #36 Pendulum with non-linear damping, no driving (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #36 Pendulum with non-linear damping, no driving (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #37 Van der Pol oscillator (standard form), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 9 of 20 | 0 of 20 | 0 of 20 |
| #37 Van der Pol oscillator (standard form), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 19 of 20 | 0 of 20 | 11 of 20 |
| #38 Van der Pol oscillator (simplified form from Strogatz), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 10 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #38 Van der Pol oscillator (simplified form from Strogatz), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 1 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #39 Glycolytic oscillator, e.g., ADP and F6P in yeast (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #39 Glycolytic oscillator, e.g., ADP and F6P in yeast (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #40 Duffing equation (weakly non-linear oscillation), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 5 of 20 |
| #40 Duffing equation (weakly non-linear oscillation), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 18 of 20 | 0 of 20 | 9 of 20 |
| #41 Cell cycle model by Tyson for interaction between protein cdc2 and cyclin (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #41 Cell cycle model by Tyson for interaction between protein cdc2 and cyclin (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #42 Reduced model for chlorine dioxide-iodine-malonic acid rection (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #42 Reduced model for chlorine dioxide-iodine-malonic acid rection (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #43 Driven pendulum with linear damping / Josephson junction (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #43 Driven pendulum with linear damping / Josephson junction (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #44 Driven pendulum with quadratic damping (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #44 Driven pendulum with quadratic damping (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #45 Isothermal autocatalytic reaction model by Gray and Scott 1985 (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #45 Isothermal autocatalytic reaction model by Gray and Scott 1985 (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #46 Interacting bar magnets, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #46 Interacting bar magnets, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #47 Binocular rivalry model (no oscillations), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #47 Binocular rivalry model (no oscillations), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #48 Bacterial respiration model for nutrients and oxygen levels, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #48 Bacterial respiration model for nutrients and oxygen levels, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #49 Brusselator: hypothetical chemical oscillation model (dimensionless), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 1 of 20 | 0 of 20 | 0 of 20 |
| #49 Brusselator: hypothetical chemical oscillation model (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #50 Chemical oscillator model by Schnackenberg 1979 (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #50 Chemical oscillator model by Schnackenberg 1979 (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #51 Oscillator death model by Ermentrout and Kopell 1990, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #51 Oscillator death model by Ermentrout and Kopell 1990, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Smoothed derivatives: systems with 3 variables (20 setups)

Smoothed derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #52 Maxwell-Bloch equations (laser dynamics), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #52 Maxwell-Bloch equations (laser dynamics), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #53 Model for apoptosis (cell death), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #53 Model for apoptosis (cell death), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #54 Lorenz equations in well-behaved periodic regime, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #54 Lorenz equations in well-behaved periodic regime, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #55 Lorenz equations in complex periodic regime, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #55 Lorenz equations in complex periodic regime, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #56 Lorenz equations standard parameters (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #56 Lorenz equations standard parameters (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #57 Rössler attractor (stable fixed point), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #57 Rössler attractor (stable fixed point), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #58 Rössler attractor (periodic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #58 Rössler attractor (periodic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #59 Rössler attractor (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #59 Rössler attractor (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #60 Aizawa attractor (chaotic), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #60 Aizawa attractor (chaotic), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #61 Chen-Lee attractor; system for gyro motion with feedback control of rigid body (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #61 Chen-Lee attractor; system for gyro motion with feedback control of rigid body (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Smoothed derivatives: systems with 4 variables (4 setups)

Smoothed derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #62 Binocular rivalry model with adaptation (oscillations), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #62 Binocular rivalry model with adaptation (oscillations), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #63 SEIR infection model (proportions), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #63 SEIR infection model (proportions), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

## Raw derivatives

### Raw derivatives: systems with 1 variable (46 setups)

Raw derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #1 RC-circuit (charging capacitor), start 1 | yes | 1 of 1 | 1 of 1 | 0 of 20 | 16 of 20 | 0 of 20 | 2 of 20 | 0 of 20 | 0 of 20 |
| #1 RC-circuit (charging capacitor), start 2 | yes | 1 of 1 | 1 of 1 | 17 of 20 | 17 of 20 | 0 of 20 | 4 of 20 | 0 of 20 | 0 of 20 |
| #2 Population growth (naive), start 1 | yes | 1 of 1 | 1 of 1 | 19 of 20 | 20 of 20 | 0 of 20 | 16 of 20 | 0 of 20 | 13 of 20 |
| #2 Population growth (naive), start 2 | yes | 1 of 1 | 1 of 1 | 19 of 20 | 20 of 20 | 0 of 20 | 16 of 20 | 0 of 20 | 14 of 20 |
| #3 Population growth with carrying capacity, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 9 of 20 | 18 of 20 | 0 of 20 | 17 of 20 |
| #3 Population growth with carrying capacity, start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 17 of 20 | 13 of 20 | 9 of 20 | 9 of 20 |
| #4 RC-circuit with non-linear resistor (charging capacitor), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #4 RC-circuit with non-linear resistor (charging capacitor), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #5 Velocity of a falling object with air resistance, start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 14 of 20 | 0 of 20 | 1 of 20 | 0 of 20 | 1 of 20 |
| #5 Velocity of a falling object with air resistance, start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #6 Autocatalysis with one fixed abundant chemical, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 10 of 20 | 20 of 20 | 5 of 20 | 1 of 20 |
| #6 Autocatalysis with one fixed abundant chemical, start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 19 of 20 | 4 of 20 | 7 of 20 | 0 of 20 | 2 of 20 |
| #7 Gompertz law for tumor growth, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #7 Gompertz law for tumor growth, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #8 Logistic equation with Allee effect, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #8 Logistic equation with Allee effect, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #9 Language death model for two languages, start 1 | yes | 1 of 1 | 1 of 1 | 17 of 20 | 19 of 20 | 1 of 20 | 10 of 20 | 0 of 20 | 4 of 20 |
| #9 Language death model for two languages, start 2 | yes | 1 of 1 | 1 of 1 | 14 of 20 | 0 of 20 | 8 of 20 | 0 of 20 | 6 of 20 | 0 of 20 |
| #10 Refined language death model for two languages, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #10 Refined language death model for two languages, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #11 Naive critical slowing down (statistical mechanics), start 1 | yes | 1 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 6 of 20 | 0 of 20 | 9 of 20 |
| #11 Naive critical slowing down (statistical mechanics), start 2 | yes | 1 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 18 of 20 | 0 of 20 | 7 of 20 |
| #12 Photons in a laser (simple), start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 9 of 20 | 6 of 20 | 4 of 20 | 0 of 20 | 1 of 20 |
| #12 Photons in a laser (simple), start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 20 of 20 | 11 of 20 | 20 of 20 | 8 of 20 | 4 of 20 |
| #13 Overdamped bead on a rotating hoop, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #13 Overdamped bead on a rotating hoop, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #14 Budworm outbreak model with predation, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #14 Budworm outbreak model with predation, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #15 Budworm outbreak with predation (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #15 Budworm outbreak with predation (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #16 Landau equation (typical time scale tau = 1), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #16 Landau equation (typical time scale tau = 1), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #17 Logistic equation with harvesting/fishing, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #17 Logistic equation with harvesting/fishing, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #18 Improved logistic equation with harvesting/fishing, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #18 Improved logistic equation with harvesting/fishing, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #19 Improved logistic equation with harvesting/fishing (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #19 Improved logistic equation with harvesting/fishing (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #20 Autocatalytic gene switching (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #20 Autocatalytic gene switching (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #21 Dimensionally reduced SIR infection model for dead people (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #21 Dimensionally reduced SIR infection model for dead people (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #22 Hysteretic activation of a protein expression (positive feedback, basal promoter expression), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #22 Hysteretic activation of a protein expression (positive feedback, basal promoter expression), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #23 Overdamped pendulum with constant driving torque/fireflies/Josephson junction (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #23 Overdamped pendulum with constant driving torque/fireflies/Josephson junction (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Raw derivatives: systems with 2 variables (56 setups)

Raw derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #24 Harmonic oscillator without damping, start 1 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 19 of 20 | 1 of 20 | 18 of 20 | 0 of 20 | 9 of 20 |
| #24 Harmonic oscillator without damping, start 2 | yes | 1 of 1 | 1 of 1 | 19 of 20 | 17 of 20 | 2 of 20 | 16 of 20 | 0 of 20 | 11 of 20 |
| #25 Harmonic oscillator with damping, start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 |
| #25 Harmonic oscillator with damping, start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 19 of 20 |
| #26 Lotka-Volterra competition model (Strogatz version with sheeps and rabbits), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #26 Lotka-Volterra competition model (Strogatz version with sheeps and rabbits), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #27 Lotka-Volterra simple (as on Wikipedia), start 1 | yes | 1 of 1 | 0 of 1 | 3 of 20 | 3 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #27 Lotka-Volterra simple (as on Wikipedia), start 2 | yes | 1 of 1 | 1 of 1 | 20 of 20 | 18 of 20 | 7 of 20 | 12 of 20 | 0 of 20 | 5 of 20 |
| #28 Pendulum without friction, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #28 Pendulum without friction, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #29 Dipole fixed point, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #29 Dipole fixed point, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #30 RNA molecules catalyzing each others replication, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #30 RNA molecules catalyzing each others replication, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #31 SIR infection model only for healthy and sick, start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 3 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #31 SIR infection model only for healthy and sick, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #32 Damped double well oscillator, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 3 of 20 | 0 of 20 | 0 of 20 |
| #32 Damped double well oscillator, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #33 Glider (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #33 Glider (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #34 Frictionless bead on a rotating hoop (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #34 Frictionless bead on a rotating hoop (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #35 Rotational dynamics of an object in a shear flow, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #35 Rotational dynamics of an object in a shear flow, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #36 Pendulum with non-linear damping, no driving (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #36 Pendulum with non-linear damping, no driving (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #37 Van der Pol oscillator (standard form), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 10 of 20 | 0 of 20 | 1 of 20 |
| #37 Van der Pol oscillator (standard form), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 18 of 20 | 0 of 20 | 7 of 20 |
| #38 Van der Pol oscillator (simplified form from Strogatz), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #38 Van der Pol oscillator (simplified form from Strogatz), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #39 Glycolytic oscillator, e.g., ADP and F6P in yeast (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #39 Glycolytic oscillator, e.g., ADP and F6P in yeast (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #40 Duffing equation (weakly non-linear oscillation), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 20 of 20 | 0 of 20 | 9 of 20 |
| #40 Duffing equation (weakly non-linear oscillation), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 19 of 20 | 0 of 20 | 11 of 20 |
| #41 Cell cycle model by Tyson for interaction between protein cdc2 and cyclin (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #41 Cell cycle model by Tyson for interaction between protein cdc2 and cyclin (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #42 Reduced model for chlorine dioxide-iodine-malonic acid rection (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #42 Reduced model for chlorine dioxide-iodine-malonic acid rection (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #43 Driven pendulum with linear damping / Josephson junction (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #43 Driven pendulum with linear damping / Josephson junction (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #44 Driven pendulum with quadratic damping (dimensionless), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #44 Driven pendulum with quadratic damping (dimensionless), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #45 Isothermal autocatalytic reaction model by Gray and Scott 1985 (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #45 Isothermal autocatalytic reaction model by Gray and Scott 1985 (dimensionless), start 2 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #46 Interacting bar magnets, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #46 Interacting bar magnets, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #47 Binocular rivalry model (no oscillations), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #47 Binocular rivalry model (no oscillations), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #48 Bacterial respiration model for nutrients and oxygen levels, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #48 Bacterial respiration model for nutrients and oxygen levels, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #49 Brusselator: hypothetical chemical oscillation model (dimensionless), start 1 | yes | 0 of 1 | 1 of 1 | 0 of 20 | 20 of 20 | 0 of 20 | 3 of 20 | 0 of 20 | 0 of 20 |
| #49 Brusselator: hypothetical chemical oscillation model (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #50 Chemical oscillator model by Schnackenberg 1979 (dimensionless), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #50 Chemical oscillator model by Schnackenberg 1979 (dimensionless), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #51 Oscillator death model by Ermentrout and Kopell 1990, start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #51 Oscillator death model by Ermentrout and Kopell 1990, start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Raw derivatives: systems with 3 variables (20 setups)

Raw derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #52 Maxwell-Bloch equations (laser dynamics), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #52 Maxwell-Bloch equations (laser dynamics), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #53 Model for apoptosis (cell death), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #53 Model for apoptosis (cell death), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #54 Lorenz equations in well-behaved periodic regime, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #54 Lorenz equations in well-behaved periodic regime, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #55 Lorenz equations in complex periodic regime, start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #55 Lorenz equations in complex periodic regime, start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #56 Lorenz equations standard parameters (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #56 Lorenz equations standard parameters (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #57 Rössler attractor (stable fixed point), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #57 Rössler attractor (stable fixed point), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #58 Rössler attractor (periodic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #58 Rössler attractor (periodic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #59 Rössler attractor (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #59 Rössler attractor (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #60 Aizawa attractor (chaotic), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #60 Aizawa attractor (chaotic), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #61 Chen-Lee attractor; system for gyro motion with feedback control of rigid body (chaotic), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #61 Chen-Lee attractor; system for gyro motion with feedback control of rigid body (chaotic), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

### Raw derivatives: systems with 4 variables (4 setups)

Raw derivatives: runs that found the exact law. Columns per noise level: BACON, then PySINDy.

| setup | true law in PySINDy's terms? | no noise: BACON (1 run) | no noise: PySINDy (1 run) | 1% noise: BACON (20 runs) | 1% noise: PySINDy (20 runs) | 5% noise: BACON (20 runs) | 5% noise: PySINDy (20 runs) | 10% noise: BACON (20 runs) | 10% noise: PySINDy (20 runs) |
|---|---|---|---|---|---|---|---|---|---|
| #62 Binocular rivalry model with adaptation (oscillations), start 1 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #62 Binocular rivalry model with adaptation (oscillations), start 2 | no | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #63 SEIR infection model (proportions), start 1 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |
| #63 SEIR infection model (proportions), start 2 | yes | 0 of 1 | 0 of 1 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 | 0 of 20 |

<!-- /tables -->
