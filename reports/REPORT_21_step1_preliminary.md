# REPORT_21 Step 1, preliminary: the general error budget and hypothesis 4

**Status: preliminary. The acceptance did not finish; this is not the Step 1 report.** SPEC_21 v0.3, REVIEW_21_step0.
STATE.md is unchanged and nothing is committed. Latitudes are planetographic with planetocentric in parentheses,
except where a table says planetocentric.

## 1. Where it stands

The acceptance (`tests/step21_1/accept_step21_1.py`, log `reports/step21_1/run.log`) ran the synthetic checks
1a to 1d, the hypothesis 4 fit, its three Lindal runs and the nine `_h4` matrix runs, and drew the matrix figure.
It stopped at the first hypothesis 4 section (section 2). Checks 2 (sections, shear and wind, IRIS residuals),
3 (consistency, leave-one-out, Huber), 4 (start from the prior mean, the pressure identity) and 5 (reference swap)
have not run. There is no `summary.json`.

| Check | Result |
|---|---|
| 1a resolution, hypothesis 1 as truth | reported |
| 1b balanced 20 K equatorial layer near 15 mbar | reported: not recovered at 15 mbar |
| 1c linear statistics, 20 prior draws | **FAIL**: IRIS 290 mbar mean ratio 0.73 |
| 1d physical draws | reported: equatorial wind in the draws unphysical |
| 2 to 5 | not run: the Lindal section's southern march did not converge |

Three items need a decision. Two are mine to diagnose with the author's approval (section 4.1 and 4.2); one
is a ruling on SPEC_21 section 2's prior (section 4.3).

## 2. What Step 1 built

As REVIEW_21_step0 ruled, with every amplitude, length, form and threshold in the control file with its basis:

* **Ingest.** An `[uncertainty]` table in each Voyager scalars file; kind T now carries each anchor's stated
  error. Voyager 2 egress uses its own curve; both Voyager 1 profiles use Lindal's Fig. 4 curve times 2, and
  Voyager 1 ingress adds 2 K in quadrature from 600 to 1200 mbar. At 110 mbar: 31.2 S (26.1 S) 2.27 K,
  74.2 S (71.1 S) and 3.0 S (2.4 S) 4.54 K. Lindal's Fig. 4 curve on his 66 levels is a new file beside his,
  `occul_data/lindal/lindal_temperature_uncertainty.csv`. Lindal's registered files are unchanged (md5).
* **Prior stated in temperature.** A Matérn 5/2 temperature covariance in latitude, its gradient covariance
  conditioned to zero gradient at the equator and both poles, divided by `w` at the base wind and Lindal's mean
  temperature, formed on a fine grid and projected onto the 2.5 deg knots (1533 unknowns). Built in about 82 s
  with `balance_weight` in parallel.
* **Radio.** Each profile's stated error with a 0.3 K floor, plus representativeness 0.5 K (1 bar/p)^0.5,
  capped at 15 K, Matérn 3/2 over half a scale height, for every profile including the reference anchor.
* **IRIS.** The equatorial omission removed and the gap rule kept. The 110 mbar placement as a rank-one term per
  level, `d dᵀ/12` with `d` the operator's change from 110 to 150 mbar, recomputed each pass (ruling 8). A
  systematic part, Gaussian over 20 deg, sized after the first fit.
* **Huber** at 2 on whitened residuals per group; the reference wind held at zero spread.
* **Controls.** `fits/voyager_fit_2/voyager_fit_2.toml` and `voyager_fit_2_reference_31s.toml`.

## 3. Results so far

### 3.1 The hypothesis 4 fit (`reports/step21_1/command_fit.log`)

Converged in 4 passes; the last plain step changed the delivered temperature by at most 0.016 K. IRIS needed no
systematic part, so the second fit is the first. Huber's rule downweighted 0 of 298 rows.

| Group | 31.2 S | 74.2 S | 3.0 S | IRIS 110 | IRIS 290 | IRIS 730 | all |
|---|---|---|---|---|---|---|---|
| misfit / expected | 0.13 | 0.10 | 0.17 | 0.79 | 0.36 | 0.27 | 0.14 |

IRIS 110's ratio was 2.32 under SPEC_19. The wind at the equator for 1000, 110, 10 and 1 mbar is 428, 381, 222
and 96 m/s, against hypothesis 3's 496, 270, 98 and 89. At 75 N planetocentric it is 97, 88, 81 and 80, against
104 at 110 mbar and 127 at 10 mbar under hypothesis 3. Over 1.3 bar to 0.2 mbar the wind spans -39 to 468 m/s.

**The budget is generous.** In 1a the data were drawn from the stated R, and the radio groups came out at 0.97,
0.90 and 0.69. The real data's radio ratios are 0.10 to 0.17, so the stated radio errors (the anchors' curves and
the representativeness term) are larger than the real scatter by a factor of 2.5 to 3 in sd. The overall 0.14 is
`innovation · S⁻¹ innovation / m`, so S as a whole is too large for these data, which section 4.3 supports.

### 3.2 The anchor matrix (`reports/figures/step21_1_matrix_absolute.png`)

All nine `_h4` runs finished. Absolute profiles, hypothesis 4 solid, hypothesis 3 dashed, both hydrogen adiabats.
Hypothesis 4 removes hypothesis 3's warm bump near 10 mbar in the transfers from 3.0 S (2.4 S), and the transfers
sit inside the stated bands over more of the column at every target. The numeric matrix was not printed before
the stop; it comes with the full run.

### 3.3 Checks 1a and 1b (`reports/figures/step21_1_synthetic.png`)

| | RMS shear error | true RMS | within 1 sd | within 2 sd | all ratio |
|---|---|---|---|---|---|
| 1a, hypothesis 1 | 24.75 | 8.47 | 99.7% | 100.0% | 0.55 |
| 1b, equatorial layer | 27.14 | 12.61 | 99.5% | 100.0% | 0.61 |

(Shear in m/s per unit ln p.) The posterior sd covers the error almost everywhere, but it is wide: at the equator
about 490 at 110 mbar and 1260 at 15 mbar. Projected on the 2.5 deg knots the 110 mbar sd is about 500 against
about 600 for the exact prior, a smoothing to note. In 1b the 15 mbar layer is not recovered: the true equatorial
shear there is -117, the fit +18. The data at the Voyager latitudes and IRIS's levels do not reach it, and the
prior lets the fit place a large shear anywhere near the equator, so the recovered field also shows a spurious
positive shear column from about 30 to 0.3 mbar in both 1a and 1b.

### 3.4 Check 1c: FAIL

Pooled, 95.5% of the shear lies within 2 posterior sd (bound 90%). Each group's mean ratio over 20 draws, with
its spread:

| 31.2 S | 74.2 S | 3.0 S | IRIS 110 | IRIS 290 | IRIS 730 | all |
|---|---|---|---|---|---|---|
| 0.94 (0.16) | 1.01 (0.16) | 0.97 (0.21) | 1.20 (1.56) | **0.73 (0.29)** | 0.98 (0.48) | 0.98 (0.10) |

IRIS 110 is at the bound, but its spread is 1.56, because the rank-one placement term dominates it, so 1.20 is
consistent with 1 over 20 draws. IRIS 290 is about four standard errors low. In 1c the truth comes from the prior
and the noise from the same R the check uses, so the expected value is exact and the ratio should average 1.
Section 4.2.

### 3.5 Check 1d

The largest |T′| of the 20 draws, its mean along the isobar removed, by the balance:

| | 36.3 N (30.8 N) | 31.2 S (26.1 S) | 74.2 S (71.1 S) | 3.0 S (2.4 S) | 87.5 N | 87.5 S | σ_T |
|---|---|---|---|---|---|---|---|
| 1000 mbar | 3.4 | 3.2 | 3.5 | 1.6 | 3.9 | 3.7 | 2.0 |
| 110 mbar | 11.4 | 7.5 | 9.9 | 9.0 | 12.9 | 11.6 | 6.0 |
| 10 mbar | 26.7 | 29.7 | 24.6 | 25.0 | 30.1 | 20.8 | 15.0 |
| 1 mbar | 38.7 | 31.7 | 33.3 | 38.3 | 24.7 | 20.1 | 15.0 |

(87.5 N and S are planetocentric.) The temperatures are physical: a maximum over 20 draws near 2σ_T, as expected.

**The equatorial wind is not.** The largest |shear| of the draws at the equator is 3861 m/s per unit ln p. The
wind change from 1000 to 10 mbar, by draw, in m/s:
-2299, +2976, +605, -1207, +203, +2431, +1689, +120, -447, -3708, +3538, -3102, -1594, -2275, -961, -3413,
-4493, -3404, -568, +593.

## 4. Open items

### 4.1 The southern march does not converge (stopped the acceptance)

`casspian-section` on `sections/lindal_voyager_fit_2` finished its northern march and failed in the southern one.
The transfer's outer loop reached about 10⁻⁶ in |d ln p| by pass 4, then wandered between 1×10⁻⁶ and 1×10⁻⁵ for
50 passes against a tolerance of 10⁻⁸, and jumped to 8×10⁻⁴ at passes 38 and 39. The jump looks like a discrete
switch, not slow convergence. The same tool converged under hypothesis 3, so hypothesis 4's wind is the trigger.
The error does not name the latitude.

Proposed diagnosis, for approval: find the stalled latitude and what changes at pass 38, then fix whichever of
the transfer or the hand-off of the 2.5 deg wind it is. No tolerance is to be loosened.

### 4.2 IRIS 290 mbar in 1c

Proposed diagnosis, for approval: at the final linearization, 500 draws for the ratios alone (no posterior sd,
so cheap), the effective degrees of freedom of each group's ratio, and the condition number of S. This tells
sampling and heavy tails apart from numerical loss in `cho_solve`.

### 4.3 The prior at the equator: a ruling for the reviewer

Dividing the temperature-gradient covariance by `w`, which goes to zero at the equator, gives the shear prior no
physical limit there. Section 3.5 shows the result: equatorial winds changing by thousands of m/s across the
column. The real fit is not harmed (its equatorial wind is 428 to 96 m/s), but the prior's breadth is why 1b
cannot recover an equatorial layer near 15 mbar, why the synthetic fits show a spurious equatorial shear column,
and much of why the overall misfit ratio is 0.14. Any limit, for example a cap on the shear sd or a wind-based
term near the equator, changes SPEC_21 section 2's prior, so it is not mine to choose.

## 5. Departures and notes

* `fits/voyager_fit_1`'s control file uses SPEC_19's format and no longer loads with the new `casspian-wind-fit`.
  Its products are unchanged and read as they are; only a refit of it would need the old tool.
* `tests/step20_1/build_matrix.py` was generalized (a suffix, the Lindal row's run names, the hypothesis label)
  so SPEC_21 makes the `_h4` row with it. SPEC_20's behavior with its defaults is unchanged.
* The acceptance runs from the start each time, about 2.5 hours to the point where it stopped. I propose a resume
  point after the matrix runs, so a rerun of sections 2 to 5 does not refit.
* No regression was run; no computed value of an accepted product changed.
