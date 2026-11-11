# REPORT_21 Step 1, preliminary under v0.4: the stability constraint stops the fits from converging

**Status: preliminary. The v0.4 acceptance did not finish; this is not the Step 1 report.** SPEC_21 v0.4,
REVIEW_21_step1_preliminary. STATE.md is unchanged and nothing is committed. Latitudes are planetographic with
planetocentric in parentheses.

## 1. Where it stands

The two diagnoses REVIEW_21_step1_preliminary approved (section 4.1, the southern march; section 4.2, IRIS 290 mbar
in 1c) are done and the mesh fix is in; their results are as reported in chat and are repeated in section 5 for the
record. The two v0.4 rules (section 1 items 1a and 1b) are built as ruled. The rerun of Step 1
(`tests/step21_1/accept_step21_1.py`, log `reports/step21_1/run.log`; the first run's log kept as `run_v03.log`)
failed:

| Fit | v0.3 run | v0.4 run |
|---|---|---|
| 1a synthetic, hypothesis 1 as truth | converged, 7 passes | not converged in 20 passes |
| 1b synthetic, the balanced equatorial layer | converged, 4 passes | not converged in 20 passes |
| hypothesis 4, first fit | converged, 4 passes | not converged in 20 passes; the script stopped (exit 1) |

1c, 1d, checks 2 to 5 and the figures did not run. There is no `summary.json`.

In each v0.4 fit the solve's predicted change in delivered temperature falls to about 0.13 to 0.2 K by pass 4 and
stays there; the iteration ends only on a plain step under 0.1 K, so it never ends. The Anderson weights swing
(as far as -26.5) and the delivered temperatures jump by 1 to 2 K on the passes where a pair enters or leaves the
stability constraint's active set. One pair, 1045 to 1095 mbar, is active at every target from pass 3 on, with an
excess over the adiabat that does not shrink (+0.005 to +0.007 in `d ln T / d ln p`, every pass).

## 2. Which v0.4 change: the stability constraint

The hypothesis 4 first fit, same start, 15 passes at most (`scratchpad/diag_v04.py`, not in the repository):

| Run | Deformation floor | Stability constraint | Result |
|---|---|---|---|
| A | on | off | converged in 3 passes, last change 0.044 K |
| B | off | on | not converged in 15 passes |

Run B repeats the failure: from pass 4 on the solve asks for 0.13 to 0.2 K every pass; on pass 10, as pairs at
500 to 575 mbar entered, the delivered temperature moved 2.1 K. The deformation-radius floor (item 1a) does not
affect convergence.

The active pairs are neighboring delivered levels about 0.047 apart in ln p. With sigma 0.005 in `d ln T / d ln p`,
a row holds the pair's temperature difference to 0.020 to 0.032 K, far tighter than any radio row.

## 3. Why: the row predicts a slope change the forward run does not make

At the hypothesis 4 start I stepped the state along the 31.2 S (26.1 S) row of the 1045 to 1095 mbar pair and
reran the three forward runs (`scratchpad/diag_rows.py`):

| Step along the row | Aimed slope change at 31.2 S (26.1 S) | Actual | Fraction delivered |
|---|---|---|---|
| prior-weighted, full | -0.0070 | -0.00090 | 0.13 |
| prior-weighted, half | -0.0035 | -0.00046 | 0.13 |
| least-squares (identity-weighted) | -0.0070 | -0.00016 | 0.02 |

The half step delivers half the change, so the response is linear; the mismatch is not a large-step effect. The
other two targets, 74.2 S (71.1 S) and 3.0 S (2.4 S), give the same fractions.

The rows are right elsewhere. Over each delivered profile below 100 mbar the predicted and actual temperature
changes agree to 0.05 K or better for the prior-weighted steps (largest change about 0.67 K) and to 0.09 K for the
least-squares step. The disagreement is only in the difference between two close levels, which the row divides
by 0.047.

The forward run moves the two levels together. In the full prior-weighted step at 31.2 S (26.1 S):

| Level (mbar) | Predicted dT (K) | Actual dT (K) |
|---|---|---|
| 997.9 | +0.376 | +0.429 |
| 1044.6 | +0.338 | +0.330 |
| 1094.7 | +0.299 | +0.329 |
| 1199.8 | +0.230 | +0.239 |

Every step and every target show the same: the actual changes at 1044.6 and 1094.7 mbar agree to a few
thousandths of a kelvin. The base run's wind grid
(`forward/lindal_iris_v_sl_ii_110/inputs/lindal_iris_v_sl_ii_110_wind.nc`, 121 levels, 0.05 decade, 0.115 in ln p)
has nodes at 1000 and 1122 mbar, so both levels sit inside one interval.

My reading, inferred from these numbers and not yet traced through the forward code: the forward run takes the
wind between grid nodes as linear in ln p, so the shear, and the temperature change it balances, is constant
across an interval; two delivered levels in one interval move together and their slope cannot change. The
row, built from the B-spline shear evaluated at each level's label (`radio.delivered`), sees a shear that varies
within the interval. For a single level the two agree closely; for a slope over 0.047 in ln p they do not.

The pinned pair is therefore one whose excess comes from the base run's own profile, which the fit cannot
remove. Each solve pushes hard on it, through a row with about eight times the leverage the forward model has;
the push moves the rest of the profile by about 0.15 K, the rerun leaves the pair's slope where it was, and the
iteration does not settle. Pairs that switch in and out (500 to 575 mbar, 630 to 830 mbar) add the larger jumps.

## 4. For decision

This is a ruling on item 1b's form; no change has been made. Options, not exclusive:

1. **Slopes at the forward grid's resolution.** Take the constraint between levels on either side of a wind-grid
   node, or over the wind grid's spacing, so a pair's slope is one the forward run can move.
2. **Rows that see the wind as the forward run does.** Form the stability rows from the shear interpolated on the
   wind grid, so their leverage matches the forward run's. A pair inside one interval then has almost no leverage
   and the constraint can only report it.
3. **A finer wind grid** in the base run. This changes the forward runs and would need its own regression.

Whichever form is chosen, a pair whose excess is set by the base run's profile, not the shear, cannot be removed
by this fit; the ruling may want such pairs reported rather than penalized.

## 5. The approved diagnoses, for the record

**Section 4.1, the southern march.** The section mesh had 29 duplicate latitude nodes, 2.5 deg multiples placed
twice about 3e-17 rad apart. Differences across the zero-width gaps turned rounding into jumps of exactly 8 or 16 in
`d ln p / d phi` every pass, holding the loop near 1e-6 with occasional kicks (pass 38). `lib.mesh.build_mesh` now
drops a uniform node within 1e-6 of a spacing of a required latitude. The southern march converges in 6 passes
(3.5e-10), the northern in 5. No forward run's mesh has such a pair, so no run's values change; no tolerance was
loosened. The SPEC_20 hypothesis 3 Lindal section remade with the fix moves by at most 0.034 K in temperature,
38 m in height and nothing in wind; those section products are not registered and the SPEC_20 figures are not
remade.

**Section 4.2, IRIS 290 mbar in 1c.** Over 500 draws at the final linearization of the first run each group's
mean ratio is 0.96 to 1.01 (IRIS 290 mbar 0.982 +/- 0.030; all 0.999). IRIS 290 mbar has 4.7 effective degrees of
freedom (IRIS 110 mbar 2.1, IRIS 730 mbar 12.2, the radio groups 57 to 61), so its 20-draw mean has sd 0.15 and the
first run's 0.73 was 1.8 sd low. S's condition number is 3.8e5 (1.9e4 equilibrated). The 20-draw bound stands;
the acceptance now reports each group's effective degrees of freedom, the 500-draw means and S's condition number
beside it.

## 6. What v0.4 built (unchanged by this failure)

* **Deformation floor (item 1a).** `prior.temperature_prior` warps latitude so the length is
  `max(l0, L_R / a)`, `L_R = c / sqrt(f^2 + beta c)`, `c` = N times the prior's vertical length in height, N from
  Lindal's profile with normal-hydrogen `c_p`. The covariance stays positive definite; with the floor off it
  matches v0.3 to 2e-10. At the equator the length used is 7 to 9.3 deg in the stratosphere, smaller near
  0.8 bar; no floor beyond 20 deg.
* **No superadiabatic delivered profile (item 1b).** `estimation/stability.py`: one-sided rows against the
  normal-hydrogen dry adiabat, `c_p` from the run's composition, sigma 0.005 in the control file with its basis;
  active only where violated and reported each pass; exempt from Huber's rule and the consistency check.
* **Acceptance.** 1c's added lines (section 5); 1d reports the equatorial length used and where the constraint
  was active; a resume point after the matrix runs (`--resume`). `read_control` takes 80 to 206 s with the floor.
