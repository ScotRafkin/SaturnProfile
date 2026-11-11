# REPORT_21 Step 1: the general error budget, and hypothesis 4

**Status: reported, not accepted. Check 5 did not finish.** SPEC_21 v0.5, after REVIEW_21_step1_preliminary_v04.
Nothing is committed. Latitudes are planetographic with planetocentric in parentheses.

Logs: `reports/step21_1/run_main.log` (the full run from the start: 1a to 4, then a crash in check 5) and
`reports/step21_1/run_resume.log` (from the resume point after the fix in section 6: the sections, checks 2 to 4
again, then a second crash in check 5). Figures: `reports/figures/step21_1_*.png`. There is no `summary.json`.

## 1. Where it stands

| Check | Result |
|---|---|
| 1a resolution, hypothesis 1 as truth | reported; 4 passes |
| 1b the balanced equatorial layer | reported; 5 passes |
| 1c statistics, 20 prior draws | FAIL on the bound, IRIS 290 mbar only (0.77); sampling, as accepted |
| 1d the prior's draws | reported |
| 2 hypothesis 4, absolute first | reported; first fit 4 passes, second 3 |
| 3 the checks | reported; overall ratio 0.14 |
| 4 the solve | PASS |
| 5 the reference does not matter | not finished: the forward model stopped at its wind grid's top (section 7) |

Every fit converged under v0.5. The stability constraint no longer holds up convergence.

## 2. The v0.5 stability constraint (section 1 item 1b), as built

- **Layers** are the intervals between the wind grid's pressure nodes (121 levels, 0.05 decade, 0.115 in `ln p`).
  Each target's delivered `ln T` is taken at the nodes, linear in `ln p` between its levels. A layer's slope is
  the difference over the interval.
- **Rows** come from `H_forward`, the delivered rows with the shear as the forward run sees it: the wind linear
  in `ln p` between nodes, so the shear is constant per cell (`Basis.x_matrix_forward`).
- **Relative to the source:**
  - `z = (slope - kappa) - max(slope_anchor - kappa_anchor, 0)`, with `kappa = R / c_p`, normal hydrogen, `c_p`
    from the composition.
  - Excess already in the anchor is reported, not penalized.
- **Penalty:** `s = width log(1 + exp(z / width))`, every layer always in, `y = 0`, `R = sigma^2`. The control
  file sets `sigma = 0.005` and `width = 0.005`; its comment gives the basis: about 1.5 percent of `kappa`, about
  0.07 K over one layer.
- **Exemptions:** the constraint is exempt from Huber's rule, and is left out of the consistency check, 1c and
  leave-one-out.
- **"Active"** means `z > 0`.
- **No finer wind grid.**

**Smoke check of the rows** (hypothesis 4 start; the 1000 to 1122 mbar layer; two perturbed forward runs;
`scratchpad/diag_rows2.py`, not in the repository). Actual over predicted change in the layer slopes:

| Perturbation | Target layer, 31.2 S | Median actual / predicted (range), per target |
|---|---|---|
| plain, along the row | -0.00714 against -0.00700 | 1.03 (0.99 to 1.05), 1.02 (0.99 to 1.05), 1.03 (1.00 to 1.06) |
| prior-weighted | -0.00716 against -0.00700 | 1.00 (0.86 to 1.07), 1.00 (0.89 to 1.17), 0.99 (0.92 to 1.13) |

Under v0.4 the rows overpredicted by about 8 times. They now match the forward run.

## 3. Where the constraint was active, and by how much

Excess over the allowed value, `z` in `d ln T / d ln p`, by pass. Every pass not listed had none. The Lindal anchor
has no excess in any layer at wind-grid resolution, so nothing was exempted in any Lindal-referenced fit.

| Fit | Pass | Active |
|---|---|---|
| 1a | 1 to 4 | none |
| 1b | 2 | 31.2 S 1000 to 1122 mbar +0.0068; 74.2 S 794 to 891 +0.0020; 3.0 S 1000 to 1122 +0.0022 |
| hypothesis 4, first fit | 1 | 31.2 S 1000 to 1122 +0.0471, 1122 to 1259 +0.0031; 74.2 S 1000 to 1122 +0.0346; 3.0 S 1000 to 1122 +0.0407 |
| | 2 | 3.0 S 1000 to 1122 +0.0001 |
| hypothesis 4, second fit | 1 to 3 | none |
| check 4, from the prior mean, first fit | 2 | 31.2 S 1000 to 1122 +0.0062 |
| check 5, Voyager 2 egress reference | 1 | 36.3 N 708 to 794 +0.0035; 74.2 S 794 to 891 +0.0006 |

**The constraint acts only on the opening passes, and nowhere at any final fit.** At the hypothesis 4 start the
1000 to 1122 mbar layer was 0.035 to 0.047 beyond the adiabat; one pass brought it back.

**Voyager 2 egress has excess of its own.** As the check 5 anchor it carries +0.0282 at 794 to 891 mbar. The
delivered profiles there carry +0.0234 (36.3 N), +0.0288 (74.2 S) and +0.0280 (3.0 S). That layer is reported
and not penalized, as v0.5 rules.

## 4. Synthetic checks (1a to 1d)

**1a, hypothesis 1 as the truth:**
- RMS shear error 15.5 against a true RMS of 8.5 m/s per unit `ln p`.
- 99.7 percent of the shear within 1 posterior sd; 100 percent within 2.
- Misfit over expected: 31.2 S 0.98, 74.2 S 0.91, 3.0 S 0.72, IRIS 110 0.26, 290 0.59, 730 0.21, all 0.55.

**1b, a balanced 20 K layer at the equator near 15 mbar:**
- RMS shear error 19.6 against a true RMS of 12.6.
- At 15 mbar the fit does not recover the layer's shear: the truth is -94 to -117 m/s per unit `ln p` within
  5 deg of the equator, and the fit gives +2 to +18. The posterior sd there is 230 to 430, so the fit states that
  it does not know.
- Misfit ratios 0.12 to 1.16, all 0.61.

**1c, 20 prior draws:**
- Pooled within 2 sd: 95.4 percent (bound 90).
- Mean ratio by group (spread): 31.2 S 0.94 (0.16), 74.2 S 1.01 (0.16), 3.0 S 0.96 (0.21), IRIS 110 1.15
  (1.00), IRIS 290 **0.77** (0.32), IRIS 730 0.96 (0.50), all 0.99 (0.11).
- Effective degrees of freedom: 61.2, 57.0, 56.5, 4.0, 4.5, 11.6. A 20-draw mean for IRIS 290 has an sd of 0.15
  under the model.
- Over 500 draws: 1.02, 1.00, 1.00, 0.97, **1.04** (its standard error 0.03), 1.02, all 1.00.
- Condition number of S: 3.78e5.

The failure is the sampling effect accepted in REVIEW_21_step1_preliminary_v04.

**1d, the prior's draws:**
- Largest `|T'|`, against `sigma_T`:
  - 1000 mbar: 1.7 to 4.0 K, against 2 K;
  - 110 mbar: 7.3 to 12.2 K, against 6 K;
  - 10 mbar: 24.6 to 31.0 K, against 15 K;
  - 1 mbar: 19.2 to 39.8 K, against 15 K.
- At the equator:
  - the largest shear is 1374 m/s per unit `ln p`;
  - the wind change from 1000 to 10 mbar ranges from -1861 to +1519 m/s over the 20 draws.
- **The deformation radius used at the equator,** as a latitude length:
  - 5.00 deg at 613 mbar, rising to 9.26 deg at 6.7 mbar;
  - 7.6 deg at 0.2 mbar;
  - the floor binds within 0 to 14.2 deg of the equator.
- **Stability:** not active at either hypothesis 4 fit, final state; no anchor excess.

## 5. Hypothesis 4 (check 2) and the checks (check 3)

**The matrix** (`step21_1_matrix_absolute.png`, hypothesis 4 solid, 3 dashed, both adiabats drawn).

At 110 mbar, hypothesis 4 brings the southern profiles together. Delivered minus observed (K):

| Pair | Hypothesis 3 | Hypothesis 4 |
|---|---|---|
| 31.2 S from 74.2 S | +6.74 | +2.71 |
| 74.2 S from 31.2 S | -6.71 | -2.67 |
| 31.2 S from Lindal | +3.27 | -0.81 |
| 74.2 S from Lindal | -3.58 | -3.62 |

At 10 mbar, hypothesis 3's ±10 to 12 K misses involving 3.0 S fall to 0.1 to 2.1 K. The full table is in
`run_main.log` under check 2.

Largest remaining misses:
- **74.2 S from Lindal:** +3.57 K at 730 mbar, -3.62 K at 110 mbar. Its stated sigma is 5.8 and 4.5 K.
- **3.0 S at 1 mbar:** +4.1 to +6.2 K. Its stated sigma is 10.7 K.

**Shear at 110 mbar near the equator** (m/s per unit `ln p`, planetocentric 5 N to 5 S):
- hypothesis 4: 84, 98, 65, 40, 31;
- hypothesis 3: 96, 131, 135, 110, 70;
- hypothesis 1: 66, 63, 58, 54, 52;
- hypothesis 4's sd: 58 to 223.

**IRIS residuals at the fit** (RMS against median stated, K/deg):
- 110 mbar: 0.045 against 0.074;
- 290 mbar: 0.016 against 0.076;
- 730 mbar: 0.028 against 0.111.

**Sections** (`step21_1_sections_temperature.png`): the zonal wind is shaded under labeled temperature contours,
with no hatching and no IRIS dots.

**Consistency (check 3).** Misfit over expected, before and after IRIS's systematic parts:

| Group | Before | After |
|---|---|---|
| 31.2 S | 0.12 | 0.12 |
| 74.2 S | 0.10 | 0.10 |
| 3.0 S | 0.16 | 0.16 |
| IRIS 110 mbar | 2.08 | 1.12 |
| IRIS 290 mbar | 0.54 | 0.62 |
| IRIS 730 mbar | 0.27 | 0.28 |
| all | 0.15 | **0.14** |

**The overall ratio stays well below 1, at 0.14.** The radio groups are 0.10 to 0.16. Only IRIS 110 mbar needed
a systematic part: 0.0527 per radian, 20 deg Gaussian. IRIS 290 and 730 took none. No value was tuned.

The low radio ratios say the stated radio errors, the representativeness and the prior together are larger than
the misfits, by about a factor of 2.5 to 3 in sd. That is the opposite of the concern in section 2, that they
may be too optimistic.

**Leave-one-out** (weighted error over rows, 1 if right):

| Left out | 31.2 S | 74.2 S | 3.0 S | IRIS 110 | IRIS 290 | IRIS 730 |
|---|---|---|---|---|---|---|
| Weighted error | 0.13 | 0.11 | 0.21 | 0.14 | 0.10 | 0.17 |

Huber's rule downweighted 0 of 523 rows.

## 6. The solve (check 4), and a departure

**Check 4: PASS.**
- Passes: first fit 4, second 3. The changes in delivered T fell from 6.7 to 11.0 K to under 0.03 K.
- From the prior mean: within 0.015 K (bound 0.1 K).
- Pressure identity: voyager2_egress 5.34e-4, voyager1_ingress 7.40e-4, voyager1_egress 4.35e-3.

The resume run repeated the sections and checks 2 to 4 with the same numbers.

**Departure: `radio.temperature_on_wind_grid`.** It read the isobars' refractivity as the hard-coded
`ln_refractivity_lindal`. A transfer run names it after its anchor, so the Voyager 2 egress-referenced fit
crashed (`KeyError`, first check 5 attempt). It now takes the run's single `ln_refractivity_` variable and raises
if there is not exactly one. For every Lindal-referenced run it picks the same variable, so no reported value
changes. The code predates SPEC_21; no earlier run reached check 5.

## 7. Check 5: stopped at the wind grid's top

> **Correction (11 November 2026).** This section is superseded by `reports/REPORT_21_step1_check5_diagnosis.md`.
> The fit did not run one pass and then stop. Its solve asked for about 15 K in every pass while the forward runs
> did not follow, and the forward run reaching the wind grid's top was a symptom of that divergence, not a
> geometry edge to be fixed for its own sake. The cause was that the fit was asked to use structure below the
> forward model's resolution. With every radio profile averaged into the forward wind grid's layers, the swap fit
> converges, and check 5 is accepted on the materiality rule (REVIEW_21_step1_check5). The text below stands as
> first written.

With the fix, the Voyager 2 egress-referenced fit ran its first pass (stability as in section 3). The solve
predicted 15.0 K; delivered T moved 9.9 K at 36.3 N, 13.4 K at 74.2 S and 13.7 K at 3.0 S.

The forward run for the next pass then raised: at 22.7 S planetocentric, geopotential 3.8e6 m2/s2, the isobar map
asked for wind at 0.9973 Pa. The wind file's grid ends at 1.0 Pa and is not extended (SPEC_04 Step 1
deliverable 2). The script exited with code 1.

This is a geometry edge of the forward run, not a failure of the fit to converge. The column's top rose about
0.3 percent in pressure past the grid. The Lindal-referenced fits never reached it.

**Possible remedies, none tried:**
1. Raise the wind file's top above 1 Pa for the fit's runs.
2. Damp the first step of a fit.

Either changes how check 5 runs; for decision.

## 8. What the earlier reviews accepted, carried

These are unchanged:
- the deformation-radius floor (v0.4);
- the diagnosis of the v0.4 non-convergence;
- the section mesh's duplicate-node fix;
- IRIS 290 mbar in 1c as sampling.

Lindal's registered files are untouched. No regression was run: the model does not change.
