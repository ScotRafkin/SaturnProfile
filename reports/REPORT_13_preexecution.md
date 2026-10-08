# REPORT 13, pre-execution. The gradient from the fitted temperatures, and the Sanchez-Lavega wind

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_13_SL_Wind.md` v0.2. `main` at
`94a3090`.

**In short.** The spec's status line reads "draft for the author", so no work has started. Two
points need a ruling before Step 1 (items 1 and 2); the rest are proposals that can stand unless
ruled otherwise.

## 1. The curve path does not take this table as it stands (blocking)

`tools/wind/curve.py` and `build_wind.py` were built for Ingersoll and Pollard, and three of their
fixed parts do not fit a single tabulation that carries its own `u_rms`:

- **The input.** `read_curve` reads a CSV with two named segments, `solid_north` and
  `solid_south`. The Sanchez-Lavega table is one fixed-width text file in four columns.
- **The gap.** `gap_rule` offers reflection of the northern segment, with a linear join blended
  onto the southern segment, or reflection followed by the Smith bins. SPEC_13 §1 item 2 asks for a
  PCHIP bridge between the data edges at -2.1 and -10.7, with no reflection and no join.
- **The uncertainty.** It is the RMS of a second source, the Smith points, about the curve per
  2 degree bin. SPEC_13 asks for the table's own `u_rms`.

Proposal: give the curve path one more input form and one more value of each rule, chosen in the
control file:
- `curve_format = "table"`, with the two columns read as one segment and the ±90 "fake data" rows
  dropped;
- `gap_rule = "pchip_bridge"`, one PCHIP through every row, which covers the ring gap and the small
  gaps alike, with no join;
- `uncertainty_source = "table"`, `u_rms` interpolated linearly in latitude.

Ingersoll and Pollard keep their current path and values, so no existing file changes. This adds
options rather than replacing any, so it needs the author's leave under "replace, never add".

## 2. The uncertainty and the provenance inside the gaps (blocking)

The table has no `u_rms` in the ring gap (-2.1 to -10.7), nor in the small gaps (-14.8 to -17.4,
-18.9 to -22.0, -28.1 to -29.6, +7.2 to +8.7). Two choices:
- (a) interpolate `u_rms` across each gap. This is continuous, but in the ring gap it understates
  what is not known.
- (b) NaN in the ring gap, which SPEC_00 §5 permits, and interpolated across the small gaps.

I propose (b). Wind values in the ring gap would carry provenance INTERPOLATED, and so would any
latitude farther than one table spacing (0.5 deg planetocentric) from a row. Rows are OBSERVED,
and the caps beyond 80.7 and -70.9 are EXTRAPOLATED by the polar rule.

## 3. Proposals that can stand unless ruled otherwise

- **`gradient = "value"`.** Central differences (`np.gradient`) of `fit`'s `value_K` on the
  planetographic grid on which `construct` already evaluates `fit`, then converted to
  per-planetocentric-radian by the same `dg_dc` as the slope. At the data's ends, one-sided
  differences. Where `value_K` is NaN, so is the gradient, as now.
- **The reference wind moves with `[wind]`.** In `lindal_iris_v_sl_ii_110` the IRIS wind is
  referenced to the Sanchez-Lavega cloud wind at 398 mbar. So both the reference and the total
  change between the second and third runs, and the shear does not. This is what "identical except
  for that wind" gives.
- **The fit's grid** is the cloud wind's latitude grid, so it changes between the second and third
  runs as well. The fitted values are the same function evaluated at other points, and the
  difference is at the level of the differencing.
- **The epoch, the frame and `u_rms`.** `vasavada_saturn_winds.note.md` lists these as to be
  confirmed. SPEC_13 §0 settles them (Voyager 1 and 2, System III, Eq. 6), and I would add one line
  to the note citing it.
- **The output.** The Sanchez-Lavega wind is built into the third run's own inputs
  (`forward/lindal_iris_v_sl_ii_110/inputs/`), as `lindal_transfer_wind.nc` is built into its run.
  SPEC_13 §2 requires this: every kind file a run is driven by stays in that run's build. Step 1's
  figures use a copy built under `reports/step13_1/`.

## 4. Not blocking

- `casspian-new-run` exists (`tools.run.new_run`).
- The table's two latitude columns are consistent with the reduction's flattening
  (0.0983 ± 0.0012, against 0.09822; from the note).
