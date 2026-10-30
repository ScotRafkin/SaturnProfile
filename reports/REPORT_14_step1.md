# REPORT 14, Step 1. Before acceptance

Coding agent, 30 October 2026. Specification: `docs/specs/SPEC_14_Reversed_Shear.md` v0.2. Working
tree on `main` at `e2d683a`, uncommitted. No regression (SPEC_14 §0).

Acceptance: `tests/step14_1/accept_step14_1.py` **2 of 3**. Check 1 misses its "value for value"
by 4.5e-3 m/s, which the diagnosis of §3 traces to the construction's stopping rule. Output in
`reports/step14_1/output.txt`; figures and tables beside it.

**In short.**
- **Case (iii)** is in `tools/lindal/lindal_wind.py` as `above_top = "reversed"`, with `hold_Pa` and
  `cross_Pa` read for that case only, and the integral in closed form. The closed form agrees with
  a numerical integral of `s_top f(p)` to 1.3e-10 m/s (§1).
- **The run** `lindal_iris_v_sl_iii_110` was made by `casspian-new-run` from
  `lindal_iris_v_sl_ii_110` and differs in `[shear]` only. Both runs were built and run on the
  same code here.
- **Delivered minus egress at 31.2 S (26.06 S planetocentric):** +0.20, -0.77, -2.17 and +2.88 K at
  10, 4.5, 1 and 0.3 mbar, against -1.90, +0.10, +5.24 and +15.10 K for the base. Every value is
  inside the egress's sigma_T of 2.88, 3.47, 5.42 and 8.74 K. At 110, 290 and 730 mbar the two
  runs agree to 0.09 K (+5.09, +2.30 and +0.40 K).
- **Against SPEC_14 §2 ruling 2's expectation** of about -0.6, -1.0 and +4.1 K at 10, 1 and 0.3
  mbar: +0.20, -2.17 and +2.88 K. These differ by 0.8, 1.2 and 1.2 K. The ruling expected that a
  line in `ln p` would not fit both 1 and 0.3 mbar, and it does not.
- **The equatorial wind** at 0.2 mbar is 472 m/s against 269 in the base. The largest wind aloft
  is 1052 m/s at 0.01 mbar, 5.6 deg planetographic. This is ruling 1's accepted consequence of
  holding f above 0.2 mbar.
- **Check 1.** Below 110 mbar the wind differs from the base by at most 4.5e-3 m/s, and the
  delivered temperature at and below 110 mbar by at most 0.030 K. Both are within the bounds the
  check also states, but not "value for value". The construction stops when the wind changes
  nowhere by more than 0.1 m/s. The larger change aloft takes case (iii) one iteration more, so
  the shear at the three levels comes from a later iterate. With the tolerance lowered to 1e-6
  m/s in a diagnostic process, the difference falls to 1.8e-7 m/s (§3). Nothing was tuned. The
  author rules (§6).

## 1. Case (iii)

`f(p)` multiplies the top level's shear above `top_level_Pa`, with `x = ln(hold_Pa / p)` and
`x_c = ln(hold_Pa / cross_Pa)`:
- `f = 1` from the top level up to `hold_Pa`;
- `f = 1 - x / x_c` above `hold_Pa`, through zero at `cross_Pa`, up to `REVERSED_TOP_Pa = 20 Pa`
  (0.2 mbar), a constant of the case beside `ABOVE_SCALE_HEIGHTS`;
- above 0.2 mbar, `f` held at its 0.2 mbar value, -2.586 for the first guesses.

The integral from the top level, closed form (`integrated_shear`, the `reversed` branch):
`s_top (x_h - x_1)` down to `hold_Pa`; then `- s_top (X - X^2 / (2 x_c))` with `X` the depth above
`hold_Pa`, capped at 0.2 mbar; then `s_top f_t (ln p - ln 20 Pa)` above. Below the top level the
code is untouched and `integrated_shear` returns the relaxed case's values bit for bit (checked in
memory on a test shear).

**Check of the closed form.** On a shear of 2, 1 and 0.5 m/s at the three levels, the closed form
against a trapezoid of `s_top f(p)` on 400,001 points from 1 Pa to 110 mbar: largest difference
1.3e-10 m/s.

`reversed_factor(p, hold_Pa, cross_Pa)` is `f(p)`, used by the integral's check and the figures.
`construct` refuses `above_top = "reversed"` without both values, values out of the order
`top_level_Pa >= hold_Pa > cross_Pa > 20 Pa` (the closed form assumes it), and either value with
another case (a key no code reads is refused as a typo). The `lindal_iris` case of
`casspian-wind-shear` accepts the two keys as optional and passes them through.

`vertical_structure` names the case, its two values and that they are tuned to the Voyager 2
egress, so the egress then checks consistency, not independence.

## 2. The run

`casspian-new-run lindal_iris_v_sl_iii_110 --from forward/lindal_iris_v_sl_ii_110`, then in
`[shear]` only: `above_top = "reversed"`, `hold_Pa = 1500.0`, `cross_Pa = 450.0`, and its comment,
which says in capitals that the two values are tuned to the Voyager 2 egress. With the run names
equated, the namelists are identical and the build files differ in those lines only.

Outer loop: 7 passes (the base 6). Pressure identity: largest |T x residual| 0.168 K at 3.99 mbar,
RMS 0.064 K over 66 levels (the base 0.064 K); every level in `pressure_identity.csv`.

## 3. Check 1, measured

| | Shipped tolerance, 0.1 m/s | Diagnostic, 1e-6 m/s |
|---|---|---|
| iterations, base and case (iii) | 5 and 6 | 12 and 13 |
| last change, base and case (iii) | 0.044, 0.031 m/s | 7.4e-7, 5.9e-7 m/s |
| largest wind difference below 110 mbar | 4.5e-3 m/s | 1.8e-7 m/s |
| largest shear difference at the three levels | 3.5e-3 m/s | 1.3e-7 m/s |

(`tests/step14_1/diagnose_below_top.py`, `reports/step14_1/diagnose_below_top.txt`; the tolerance
was lowered in that process only.) The shear at the three levels depends only on the wind at the
levels, which case (iii) does not reach. But the stopping test takes the largest change anywhere
on the grid, and case (iii)'s change aloft is about three times the base's at every iteration (577
against 177 m/s at the first). The difference shrinks with the residual, so it is the iteration's
and not the case's.

The other parts of check 1: `u_reference_ms`, `u_total_uncertainty_ms` and the latitudes are equal
(the first run of the suite compared the uncertainty without `equal_nan` and reported it unequal;
the suite was corrected). Delivered temperature at and below 110 mbar, new minus base: largest
-0.030 K, at 114.8 mbar. That is within the check's 0.1 K and needs no explanation beyond the above.

## 4. For the author

| Level (mbar) | 730 | 290 | 110 | 10 | 4.5 | 1 | 0.3 |
|---|---|---|---|---|---|---|---|
| base minus egress (K) | +0.40 | +2.32 | +5.00 | -1.90 | +0.10 | +5.24 | +15.10 |
| case (iii) minus egress (K) | +0.40 | +2.30 | +5.09 | +0.20 | -0.77 | -2.17 | +2.88 |
| egress sigma_T (K) | 2.73 | 2.29 | 2.27 | 2.88 | 3.47 | 5.42 | 8.74 |
| base density vs egress | -0.3 % | -2.3 % | -5.5 % | +1.8 % | +0.9 % | -3.7 % | -10.0 % |
| case (iii) density vs egress | -0.3 % | -2.3 % | -5.6 % | -0.2 % | +1.5 % | +1.6 % | -1.8 % |

The IRIS fit (smooth window) at 31.2 S minus the egress: +2.48, +3.34 and -1.10 K at 110, 290 and
730 mbar. The 110 mbar difference of about +5 K is the radio and IRIS disagreement of REPORT_12_step1
and REPORT_13_step3; case (iii) leaves it as it was.

**Winds** (`winds.png`, `largest_wind_aloft.csv` for every level above 110 mbar):
- The equatorial wind at 10, 1 and 0.2 mbar is +244, +285 and +472 m/s, against +301, +276 and
  +269 m/s in the base.
- The largest |u| is +468 m/s at 10 mbar (12.3 deg planetographic), +454 at 1 mbar, +522 at 0.2
  mbar (5.6 deg) and +1052 at 0.01 mbar (5.6 deg). The base stays between 425 and 461 m/s.

**f(p) and the shear** (`shear.png`). At the anchor, 36.3 N (30.81 N planetocentric), the 110 mbar
shear is negative, about -5.8 m/s per e-fold, and the reversed shear above 0.2 mbar is +15. At the
target, 31.2 S (26.06 S planetocentric), it is about +7.9, and -20.4 above 0.2 mbar.

**Lapse rate** -dT/dz above 110 mbar at 31.2 S (`lapse_rate.png`): -0.55 to +0.44 K/km for case
(iii), against -0.62 to +0.35 K/km for the base. At 10, 4.5, 1 and 0.3 mbar it is -0.55, -0.35,
+0.14 and +0.25 K/km, against -0.62, -0.44, +0.05 and +0.17.

**Profiles** (`profiles.png`): temperature and density at 31.2 S, base and case (iii), with the
anchor, the egress and its sigma_T band, the IRIS fit's values at the two labels, and delivered
minus egress.

## 5. Files

- `src/casspian/tools/lindal/lindal_wind.py`: case (iii), `REVERSED_TOP_Pa`, `reversed_factor`,
  the reversed branch of `integrated_shear`, the two keys in `construct`, `vertical_structure`.
- `src/casspian/tools/wind/shear.py`: `hold_Pa` and `cross_Pa` among the parameter keys, optional
  for `lindal_iris`, passed through.
- `forward/lindal_iris_v_sl_iii_110/`: the two control files (inputs and product not committed).
- `tests/step14_1/accept_step14_1.py` and `tests/step14_1/diagnose_below_top.py`.

## 6. For the author

Check 1 asks for the wind below 110 mbar "value for value". It differs by 4.5e-3 m/s, a remainder
of the construction's stopping rule (§3), and the delivered temperature by 0.03 K. Two ways
forward, either one revertible:
1. **Accept as measured.** The check's wording becomes "equal to within the construction's
   tolerance", and the suite asks for the difference below 110 mbar to be smaller than
   `TOLERANCE_MS`. No code change.
2. **Stop on the levels.** The stopping test takes the change at the three levels only, since only
   those feed back into the balance. The wind below 110 mbar is then the base's bit for bit,
   because both cases stop on the same iterate. That changes when every Lindal-line wind stops,
   so the three existing runs would shift within 0.1 m/s and their suites would need rerunning.

Recommended: 1.
