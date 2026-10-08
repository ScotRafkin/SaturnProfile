# REPORT 13, Step 1. The gradient from the fitted temperatures, and the Sanchez-Lavega wind

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_13_SL_Wind.md` v0.3. Working tree on
`main` at `94a3090`, uncommitted. Acceptance: `tests/step13_1/accept_step13_1.py`, **3 of 3** (two
reported items, and a check that the existing paths are unchanged); output in
`reports/step13_1/output.txt`, figures `shear_slope_value.png` and `cloud_winds.png`. No regression.

**In short.**
- **The wind's uncertainty in a run** (§3 ruling 2, asked before the runs): nothing in a run reads
  it. A NaN there does nothing to the delivered profile (§1).
- **The existing paths are unchanged.** The Ingersoll and Pollard cloud wind rebuilt by the edited
  tool is identical in every variable and attribute but the build stamps. The `lindal_iris` wind at
  the default gradient is identical to `lindal_iris_ii_110`'s.
- **The derivative differs from the slope mainly near the fit's widened windows.** Elsewhere the
  two agree closely. Near about -48 deg, the IRIS gap (-27.7 to -38.5) and -13 deg, the derivative
  has narrow spikes up to 1.4 K per degree, and the shear moves by up to 40 m/s (§2).
- **The Sanchez-Lavega cloud wind** passes through all 260 rows. It runs 6 m/s slower than
  Ingersoll and Pollard on average, with an RMS difference of 21 m/s. In the ring gap it is 379 to
  436 m/s, where the reflected curve gives 395 to 491 (§3).

## 1. Whether a run reads the wind's uncertainty

No step of a run reads `u_total_uncertainty_ms`; a NaN in it does nothing to the delivered
profile. In detail:
- **The model** reads `u_total_ms` only. The run's loader (`control._admit_wind`) checks
  `u_total_ms` and the poles, and the schema checks only that the companion exists (SPEC_00 §5
  allows NaN). No value of the companion is read.
- **The product's** `u_column_uncertainty_ms` is NaN by construction, "not propagated", whatever
  the wind's uncertainty is (`forward/production.py`, SPEC_03 decision 5).
- **The `lindal_iris` case** copies the cloud wind's uncertainty at its reference level to every
  pressure, and computes nothing from it. The NaNs of the ring gap pass into the run's wind file
  unchanged.
- **The run's figures** (`figures_product`) take the reference wind and discard its uncertainty.
- **NaN is not new to a run.** The current Ingersoll and Pollard wind already has NaN uncertainty
  in 82 of its 361 latitude rows, poleward of the data at both ends.
- **One cosmetic effect, outside a run.** `figures_inputs.panel_wind_latitude`, which plots an
  input file on request, fills the ±1 sigma band through the finite points only. Across the ring
  gap the band is drawn straight from edge to edge, not left open. It is a figure of the input, not
  of any run, and it is not changed here.

## 2. The gradient (§1 item 1)

- **`lindal_wind.construct`** takes `gradient = "slope"` (the default, SPEC_11) or `"value"`.
  `derivative_of_value` differences the fit's `value_K` on the planetographic grid: central inside
  the fit's span, one-sided at its ends, NaN outside. The span must be one contiguous run of the
  grid, or the function refuses. The derivative is then converted per planetocentric radian by
  the same `dphi_g/dphi_c` as the slope. `fit` is unchanged.
- **The case** `lindal_iris` of `casspian-wind-shear` accepts an optional `gradient`, `"slope"`
  when absent. For `"value"` the wind's `method` attribute says so; for `"slope"` it is as before.
- **Acceptance 1** (`shear_slope_value.png`): derivative minus slope.

  | Level | dT/dphi RMS | dT/dphi largest | Shear RMS | Shear largest |
  |---|---|---|---|---|
  | 110 mbar | 0.152 K/deg | 1.42 K/deg | 4.3 m/s | +39.6 m/s at -13.6 |
  | 290 mbar | 0.117 K/deg | 1.32 K/deg | 2.6 m/s | +21.0 m/s at -47.3 |
  | 730 mbar | 0.145 K/deg | 1.41 K/deg | 3.0 m/s | +31.2 m/s at -31.6 |

  K/deg is K per planetocentric degree. The construction settles in 6 iterations with the
  derivative, against 5 with the slope.
- **Where they differ, and why.** The fit's value is continuous but its derivative is not smooth.
  Where the data thin out, `fit` widens its window to keep 3 effective points, and the width bends
  in latitude. The local slope does not see the bend, but the derivative of the value does. So the
  two agree to a few hundredths of a K per degree over most of the globe. They differ in narrow
  spikes near -48 deg, across the IRIS gap from -27.7 to -38.5 deg, and near -13 deg. They also
  differ broadly between about -10 and +5 deg at 110 and 730 mbar, where the shear is
  interpolated across the equatorial band from its edge values.

  This is what the change asks for: a shear that stands for the temperatures the fit gives, bends
  included. Nothing is smoothed. Whether the spikes are wanted is for the author.

## 3. The Sanchez-Lavega cloud wind (§1 item 2, §3 rulings 1 and 2)

- **`tools/wind/curve.py`**: `read_table` and `TableCurve`.
  - The table is read from its second (planetographic) column, `u` and `u_rms`.
  - The two rows marked "fake data added here" are dropped.
  - `gap_rule = "pchip_bridge"` is one PCHIP through all 260 rows, and the caps follow the
    existing polar rule.
  - The declared ring gap's edges must be rows, and no row may lie inside it.
- **Provenance.**
  - Observed among the rows.
  - Interpolated inside any interval wider than `GAP_FACTOR` (2) times the median row spacing:
    the ring gap and the five small gaps (-29.6 to -28.1, -22.0 to -20.4, -20.4 to -18.9, -17.4 to
    -14.8, 7.2 to 8.7).
  - Extrapolated in the caps.

  The ruling names the ring gap; the small gaps take the same rule. The normal row spacing is 0.45
  to 0.62 deg planetographic, and every gap is 1.5 deg or wider.
- **Uncertainty.** `u_rms` linear in latitude. NaN in the ring gap (14 grid latitudes) and
  poleward of the rows, 83 NaN latitudes in all.
- **`tools/wind/build_wind.py`.**
  - Three new keys, `curve_format`, `uncertainty_source` and `ring_gap_deg`. The fourth, the ring
    gap, declares what ruling 2 names, so that no number is in the code.
  - `points_source`, `bin_width_deg`, `min_count_per_bin` and `join_window_deg` are now required
    only with `curve_format = "segments"`, the default. A key the format does not use is refused.
  - The texts of a tabulated source come from a `[source]` table in its data properties file.
    The Ingersoll and Pollard texts stay in the module, so that file is written as before
    (acceptance 0).
- **`data_static/winds/sanchezlavega2000.toml`** (new, for the author: data_static is normally
  written by hand).
  - The observation level, the epoch and the season are carried from `smith1982_fig4.toml`: 1 bar
    assigned, 1981-08-26, Ls 18.2.
  - The convention, the frame and `u_rms` are from the paper and SPEC_13 §0.
  - `vasavada_saturn_winds.note.md` gains one line citing SPEC_13 (ruling 4).
- **Acceptance 2** (`cloud_winds.png`).
  - The curve passes through every row exactly.
  - Within the rows and outside the ring gap, Sanchez-Lavega minus Ingersoll and Pollard is
    -6.1 m/s on average, with an RMS of 20.8 m/s.
  - The peak is 467.1 m/s at +5.6 deg. The current wind's 490.5 m/s at -7.4 deg is in its
    reflected gap fill.
  - In the ring gap the new wind is 379 to 436 m/s, against 395 to 491 for the reflection.
- **Step 1's control file** is `tests/step13_1/sl_wind_build.toml`. The third run's `[wind]`
  section is the same but for its paths.

## 4. Files

- **Changed:**
  - `src/casspian/tools/lindal/lindal_wind.py`
  - `src/casspian/tools/wind/shear.py`
  - `src/casspian/tools/wind/curve.py`
  - `src/casspian/tools/wind/build_wind.py`
  - `data_static/winds/vasavada_saturn_winds.note.md`
- **New:**
  - `data_static/winds/sanchezlavega2000.toml`
  - `tests/step13_1/accept_step13_1.py`
  - `tests/step13_1/sl_wind_build.toml`
