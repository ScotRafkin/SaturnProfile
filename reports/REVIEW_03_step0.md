# REVIEW 03, Step 0. The pressure grid and Eq. B1 along the local vertical

Review of `reports/REPORT_03_step0.md` against SPEC_03 v0.5 Step 0, SPEC_01 v0.23, SPEC_02
v0.10 and SPEC_00 v0.16. Reviewer: S. Rafkin, 14 September 2026. **Disposition: accepted.**
The two failing checks are specification errors, both confirmed; the numbers are restated in
SPEC_03 v0.6 and SPEC_01 v0.24 and the checks pass as restated. Step 1 proceeds after the
acceptance commit and the sweep.

## What was verified independently

Opened with netCDF4 and numpy, without the acceptance script.

- The raw bundle (`a2a276ba...`) carries both pressure columns, the index and the flag. The
  `k` sequence starts at −70 and runs 10 apart to 130, 8 apart to 194, 6 apart to 236, 4
  apart to 268, 2 apart to 310; 65 rows snapped, the last as printed. Every snapped value,
  rounded to the decimals its CSV row prints, equals the printed value (66 rows checked from
  the CSV). The top eight values are 19.9526, 25.1189, 31.6228, 39.8107, 50.1187, 63.0957,
  79.4328, 100.0 Pa.
- The candidate product (`e72011d0...`): `radius_m` at the anchor level equals `r0`
  (58,516,188.288884744 m) exactly and at the top level equals `r0 + (h_top − h_ref) cos ψ`
  to the last digit (58,801,571.048 m); the bottom is 58,412,566.575 m. `n` at the top is
  1.041934e22, at the second level 1.270497e22, at 1 bar 5.373123528e25, each equal to
  `p / (k_B T)` from the file's own thermo group. `phi_c`, `psi`, both polar radii and the
  asymmetry are the SPEC_02 Step 6 values to the last printed figure. The record carries the
  projection rule, a residual of 12.44 m and a drift of +0.02694 and −0.00975 degrees; the
  reviewer's own estimates before the step were about 12 m, +0.027 and −0.010.
- Finding 2: the change of the effective gravity at the top level was recomputed from the
  embedded G, R and W as 5.284e-5, the report's figure.
- The hashes of the raw bundle, kind T and kind C on disk match section 7 of the report.

## Rulings on the findings

1. **The `k` run endpoints.** The specification was wrong at both boundaries, as the report
   says: 199.53 mbar is `k` 230 and 229.09 mbar is `k` 236 (spacing 6 ends at 236); 478.63 mbar
   is `k` 268 and 501.19 is 270 (spacing 4 ends at 268). Verified from the bundle. Restated in
   SPEC_03 v0.6 and SPEC_01 v0.24. Check 1 passes as restated; rerun it so the report shows it.
2. **The gravity bound.** The specification quoted the Newtonian estimate. The report's
   decomposition is right and was reproduced: 4.488e-5 divided by `g_eff / g_N` gives
   5.02e-5, and the centrifugal part adds 2.6e-6. The bound is now 6e-5 with the
   decomposition stated (SPEC_03 v0.6). Check 10 passes as restated; rerun it.
3. **Keep-by-content.** Correct: `raw_bundle` is already compared as content, and adding
   `input_hashes` would defeat SPEC_02 Step 6 decision 4. The specification sentence is
   restated (SPEC_03 v0.6). No code change.
4. **The ammonia fill at the two filled levels.** Recorded. The change (below 4e-9 in `ℛ̄`) is
   the correct consequence of interpolating in the grid pressures and needs no action.
5. **`-dirty` refusal against acceptance before a commit.** The author's ruling (the candidate
   built in memory, the relaxation confined to the acceptance script and named in its output)
   is now the general rule in SPEC_03 v0.6 §0 for any step that changes an input file, so the
   next such step does not have to ask.

## On the decisions

1. The in-memory candidate: accepted, and made the rule (finding 5).
2. The snap where the rule left room (candidate search window, half-to-even rounding from the
   exact binary value, "two nearest fixed rows" read as the two adjacent rows, passes in
   increasing pressure, the exact-decade conversion, the refusal of a transcription without
   `[pressure_grid]`): accepted as reported. The tie rule is never exercised on Table I, which
   the report says; good.
3. `pressure_grid_index` as float64 with NaN for the excluded row: accepted.
4. No uncertainty companion on `pressure_printed_Pa`, schema version not bumped: accepted.
5. Gravity and rotation not rebuilt: accepted. They do not read the raw bundle, their content
   and hashes are unchanged, and a rebuild would change writer attributes only. SPEC_03 v0.6
   says so.
6. The projection residual and drift as defined (trapezoid on the tabulated levels, ψ at each
   level from the projected radius): accepted. The definitions are in the record's note, which
   is where they belong.
7. The DataTree guard fixed before the acceptance ran: noted, and the right way to report it.

## Order of work

1. Commit the author's documents (SPEC_03 v0.6, SPEC_01 v0.24, this review) in their own
   commit.
2. Rerun checks 1 and 10 against the restated numbers; fill in section 2; commit the step
   (the acceptance commit for Step 0). Push.
3. Sweep on the clean tree: `casspian-lindal-raw`, then the wind and composition tools, then
   `casspian-lindal-inputs`, then `casspian-refrac` with figures; gravity and rotation kept.
   Read back every commit attribute and refuse any `-dirty`. Record the SHA-256 of every
   product in section 7 in the row format the Step 02_1 suite reads, replacing the
   provisional table.
4. Run Steps 02_1, 02_4, 02_5 and 02_6 on the clean products; restore the committed Step 5
   figures afterward as the handoff says; record the results in section 8. Attach F2 and F4
   from the clean product to the report. Commit the record ("Record the Step 0 sweep").
5. Set `STATE.md` to accepted with the acceptance commit. Step 1 then proceeds.
