# REPORT 13, Step 2. The three runs

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_13_SL_Wind.md` v0.3. Working tree on
`main` at `94a3090`, uncommitted, with Step 1's changes (REPORT_13_step1). Acceptance:
`tests/step13_2/accept_step13_2.py`, **2 of 2** (both reported); output in
`reports/step13_2/output.txt`, figures `profiles.png` and `winds.png`. No regression.

**In short.**
- **Two runs, each one change from the one before.** Both were made by `casspian-new-run` and
  converge, in 5 and 6 passes, with pressure identities within 0.078 K.
- **Part (A) of REPORT_11_step2 vanishes to its linearization.** At 110 mbar it falls from -1.21 K
  (slope) to -0.17 K (derivative), of which -0.16 K is `T_anchor ln(T_t/T_a)` against `T_t - T_a`.
  So the round trip at 110 mbar closes to -0.16 K, against -1.50 K with the slope.
- **Against the egress the derivative does worse at 110 mbar**: +5.26 K, against +3.91 K. The
  transfer now carries the IRIS fit's full warming of +5.2 K from 36.3 N to 31.2 S, which the radio
  occultations do not show (REPORT_12_step1 §3).
- **The Sanchez-Lavega cloud wind changes the delivered profile by at most 0.14 K at the three
  levels.** It changes the wind by up to 100 m/s, and the shear changes with it, mostly within
  15 deg of the equator, as ruling 3 expects.

## 1. The runs

- **`forward/lindal_iris_v_ii_110`**, from `lindal_iris_ii_110`. Its one change is
  `gradient = "value"` in `[shear]`. Its namelist comment and description say so.
- **`forward/lindal_iris_v_sl_ii_110`**, from `lindal_iris_v_ii_110`. Its one change is its
  `[wind]`:
  - `curve_format = "table"`, `uncertainty_source = "table"`, `gap_rule = "pchip_bridge"` and
    `ring_gap_deg = [-10.7, -2.1]`;
  - the Sanchez-Lavega table and `sanchezlavega2000.toml` in place of the Ingersoll and Pollard
    curve, points and properties;
  - the bins and the join window removed.

  Its wind is built into its own `inputs/`.
- **`lindal_iris_ii_110`** is unchanged. Its product, from SPEC_12's rerun at the same target, is
  read and not rerun. REPORT_13_step1 acceptance 0 shows that the edited tools reproduce its wind.
- **The two new runs.**

  | Run | Outer loop | Mesh | Pressure identity |
  |---|---|---|---|
  | `lindal_iris_v_ii_110` | 5 passes | 1142 by 82 | at most 0.078 K |
  | `lindal_iris_v_sl_ii_110` | 6 passes | 1142 by 82 | at most 0.075 K |

  Their inputs and products are not committed, as for every unregistered run.

## 2. Acceptance 1: the round trip

The round trip is SPEC_11 Step 2's: the delivered `T(target) - T(anchor)` against the IRIS fit's
change in value. The anchor is at 36.30 N (30.806 N planetocentric). The target, 26.061 S
planetocentric, is 31.11 S planetographic by the wind file's own conversion, the one SPEC_11 Step 2
used (it is 31.2 S by the anchor's fixed point). So the IRIS change here is +5.18 K at 110 mbar,
against +5.21 K at 31.2 S in REPORT_12_step1.

| Level | Run | Delivered | IRIS change | Difference | (A) | (A1) | (A2) | (B) | (C) | (D) |
|---|---|---|---|---|---|---|---|---|---|---|
| 110 mbar | slope, IP | +3.67 | +5.18 | **-1.50** | -1.21 | -1.152 | -0.061 | -0.11 | -0.09 | -0.09 |
| | value, IP | +5.02 | +5.18 | **-0.16** | -0.17 | -0.016 | -0.156 | +0.18 | -0.12 | -0.05 |
| | value, SL | +5.08 | +5.18 | **-0.10** | -0.17 | -0.016 | -0.156 | +0.26 | -0.12 | -0.07 |
| 290 mbar | slope, IP | +1.31 | +0.61 | **+0.70** | +0.31 | +0.308 | +0.000 | +0.44 | -0.21 | +0.16 |
| | value, IP | +1.22 | +0.61 | **+0.60** | -0.02 | -0.014 | -0.002 | +0.60 | -0.19 | +0.21 |
| | value, SL | +1.30 | +0.61 | **+0.69** | -0.02 | -0.014 | -0.002 | +0.72 | -0.19 | +0.18 |
| 730 mbar | slope, IP | -5.26 | -4.38 | **-0.88** | -0.07 | +0.013 | -0.082 | -0.88 | +0.18 | -0.10 |
| | value, IP | -5.08 | -4.38 | **-0.70** | +0.07 | +0.151 | -0.077 | -0.85 | +0.18 | -0.11 |
| | value, SL | -4.94 | -4.38 | **-0.56** | +0.07 | +0.151 | -0.077 | -0.64 | +0.17 | -0.14 |

All values are in K. IP is Ingersoll and Pollard, SL Sanchez-Lavega. The parts add to the
difference.

- **The parts** are REPORT_11_step2 §3's, with the gradient the wind was built from in place of the
  slope:
  - (A) is `T_anchor ∫ d ln T` of that gradient against the IRIS change. It is split into (A1), the
    plain `∫ dT` against the change, and (A2), the linearization.
  - (B) and (C) are the wind's implied change against the gradient's, inside and outside the band
    `|phi_c| < 5`.
  - (D) is the delivered change against the wind's implied one.
- **One correction to the SPEC_11 diagnosis.** Its integrals ran over the grid nodes inside the
  path and left out the partial cells at the two ends. Here each integral runs from the anchor's
  latitude to the target's exactly. With the old integral (A1) was -0.10 K at 110 mbar for the
  derivative, and that was the truncation, not the gradient.
- **(A) vanishes as the specification expects, to its linearization.**
  - At 110 and 290 mbar (A1) is -0.016 and -0.014 K.
  - (A2), the `-(dT)^2 / (2 T_anchor)` of reading a change in `ln T` as one in `T`, is -0.156 K
    at 110 mbar. It is the same with either gradient, and REPORT_11 counted it in (A) without
    separating it.
  - At 730 mbar (A1) is +0.15 K. The target, 31.11 S, lies in the IRIS gap (-27.7 to -38.5),
    where `fit` widens its window and the 730 mbar value bends sharply (the spike at -31.6 deg in
    REPORT_13_step1's figure). The trapezoid of the differenced value across that bend, read to
    an end point inside it, is good to about 0.15 K there.
- **What is left is the equatorial band, (B)**: +0.60 to +0.72 K at 290 mbar and -0.64 to -0.88 K
  at 730 mbar. It is the PCHIP bridge of SPEC_11 ruling 2. The derivative leaves (B) as it was; the
  Sanchez-Lavega wind moves it by up to 0.2 K.

## 3. Acceptance 2: the profiles and the winds

**The three profiles against the egress** (`profiles.png`), at 31.2 S (26.06 S planetocentric):

| Level | slope, IP | value, IP | value, SL | Egress sigma_T |
|---|---|---|---|---|
| 110 mbar | +3.91 K (-4.3 %) | +5.26 K (-5.7 %) | +5.32 K (-5.8 %) | 2.27 K |
| 290 mbar | +2.37 K (-2.4 %) | +2.27 K (-2.3 %) | +2.36 K (-2.4 %) | 2.29 K |
| 730 mbar | -0.95 K (+0.8 %) | -0.77 K (+0.7 %) | -0.63 K (+0.5 %) | 2.73 K |

Each entry is the run minus the egress, with the density difference relative to the egress in
parentheses.

- **The derivative moves the run toward the IRIS change and away from the radio one.** With the
  wind now standing for the fitted temperatures, the transfer carries the IRIS change from 36.3 N
  to 31.2 S almost in full: +5.0 K at 110 mbar, against the fit's +5.2. The two radio profiles
  differ there by +0.03 K (REPORT_12_step1 §3). So the run moves from 1.7 to 2.3 sigma_T above the
  egress at 110 mbar.
  - The difference between the run and the egress at 110 mbar is now almost all the IRIS and
    radio disagreement, and almost none of it is the model's.
  - At 290 and 730 mbar the derivative changes the run by 0.1 to 0.2 K.
- **The Sanchez-Lavega wind** moves the profile by at most 0.14 K at the three levels, and the two
  value-gradient curves lie almost on top of each other.

**The winds** (`winds.png`): `u_total` of the three runs, then each minus the one before.
- **Derivative minus slope.** Largest -98.1 m/s, at -13.6 deg near 0.01 mbar, with an RMS of
  8.8 m/s. The reference wind is unchanged. The differences are narrow columns at the latitudes of
  REPORT_13_step1's spikes, -48, -32 and -13 deg, together with the broad equatorial difference.
  They lie mostly above 398 mbar, carried upward by the relaxation, with one column below it at
  -32 deg, where the 730 mbar spike is.
- **Sanchez-Lavega minus Ingersoll and Pollard, both with the derivative.**
  - The total: largest -99.6 m/s, RMS 22.3 m/s.
  - The reference: RMS 21.8 m/s, nearly all of the change.
  - **The shear changes too (ruling 3)**: largest -45.7 m/s at -10.5 deg, with an RMS of 15.2 m/s
    within 15 deg of the equator against 1.0 m/s beyond. It changes through `(du/dphi)_p` and the
    isobar geometry, and the third panel of the bottom row shows it. The implementation does not
    hold the shear fixed.

## 4. For the author

- **The derivative closes the round trip and opens the comparison.** With the wind built from the
  fitted temperatures, the model reproduces the IRIS change at 110 mbar to 0.16 K, about the
  linearization. The comparison with the egress then shows plainly that the IRIS and the radio
  disagree on this change, by about 5 K at 110 mbar.
- **The spikes.** REPORT_13_step1 §2 shows that the derivative carries the fit's bends where its
  window widens. They become narrow wind columns of up to about 100 m/s at the top of the grid.
  The target's path crosses them, and their effect on the delivered profile is not separated here.
- **The band remains.** After (A), the largest part is (B), up to 0.88 K at 730 mbar, set by the
  bridge rather than the data.
- **`lindal_iris_ii_110`** was not rerun. Its product is from SPEC_12's acceptance, and acceptance 0
  of Step 1 shows that its wind is reproduced exactly. If the author wants all three products from
  one tree, the script reruns it when its product is absent.

## 5. Files

- **New:**
  - `forward/lindal_iris_v_ii_110/` (its two control files)
  - `forward/lindal_iris_v_sl_ii_110/` (its two control files)
  - `tests/step13_2/accept_step13_2.py`
- **Not committed:** the two runs' inputs and products.
