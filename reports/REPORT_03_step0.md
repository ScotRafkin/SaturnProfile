# REPORT 03, Step 0. The pressure grid and Eq. B1 along the local vertical

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 14 September 2026. Specification: `SPEC_03_Forward_Production.md` v0.5, Step 0, with
`SPEC_01_Lindal_Tool_Chain.md` v0.23 (Steps 2, 8, 9), `SPEC_02_Refrac_and_Diagnostics.md` v0.10
(Steps 3, 4) and `SPEC_00_Architecture_and_Data_Files.md` v0.16 (sections 3.3, 6.1, 6.7).
Status: **accepted (REVIEW_03_step0, 14 September 2026)**, with findings 1 and 2 restated in
SPEC_03 v0.6 and SPEC_01 v0.24.

**Thirteen of thirteen acceptance checks pass** as restated. As first run against v0.5, eleven of
thirteen passed: the two failures were specification numbers, not code (findings 1 and 2), and the
checks were left as written until the review restated them. Everything the anchoring decision
depends on is unchanged bit for bit (check 9).

The acceptance ran on a **candidate** kind N built in memory (author ruling of 14 September
2026, option 1; decision 1). The registered product under `occul_data/lindal/` is still the clean
`2cf4457` product and is rebuilt only at the sweep after the acceptance commit.

## 0. Before this step

| Commit | What |
|---|---|
| `97f07e7` | The author's documents: SPEC_03 v0.4, SPEC_00 v0.16, SPEC_01 v0.22, SPEC_02 v0.10, `STATE.md`, the handoffs of 13 and 14 September, and `Claude outputs/` added to `.gitignore`. |
| `5529459` | SPEC_03 v0.5 and SPEC_01 v0.23 (the snap rule resolved from the fixed side, the basis count 48), `STATE.md`. |

The status lines of SPEC_03 v0.5 and `STATE.md` were read before starting. The pre-execution
review of v0.3 and v0.4 (three ambiguous rows, the basis count, the twelfth Step 3 case, the
compensated composition edit, the gauge match) is recorded in SPEC_03 section 8.

---

## 1. What was built

1. **`tools/lindal/build_raw.py`, the pressure grid** (SPEC_01 v0.23 Step 2).
   - `grid_candidates(text, denominator)`: every integer `k` whose grid value `10^(k/100)` mbar,
     rounded to the decimals printed in the row, equals the printed value.
   - `resolve_grid(...)`: one candidate fixes a row; ambiguous rows are settled one at a time
     from the fixed side by continuing the spacing of the two nearest fixed rows, both sides
     agreeing when both are fixed; a resolved row counts as fixed; anything left ambiguous, or a
     row with no candidate that is not excluded, raises `PressureGridError` with the row named.
     No spacing number is in the code.
   - `table1` now carries `pressure_printed_Pa` (the CSV value, `value_source = "Table I"`),
     `pressure_Pa` on the grid, `pressure_grid_index` (`k`, NaN where kept as printed) and the
     flag `pressure_grid_applied` (int8, `as_printed snapped_to_grid`). The rule, the counts, each
     resolved row with its candidates, and the specification's note on 0.32 mbar are in
     `history`.
2. **`raw/lindal_scalars.toml`** gains `[pressure_grid]` with the v0.5 text, basis count 48.
   **`raw/notes.md`** gains a section stating the same in prose. `lindal_table1.csv` is
   unchanged.
3. **`lib.schema`**: kind T allows `pressure_printed_Pa`, and refuses it without the global
   `pressure_grid_rule` (SPEC_00 v0.16 section 6.1).
4. **`tools/lindal/build_inputs.py`**: kind T carries `pressure_printed_Pa` and
   `pressure_grid_rule` (the rule, the exclusion, the basis and the value source quoted), and
   refuses a raw bundle without them. The keep-by-content rule was not changed (finding 3).
5. **`lib.reduction.absolute_radius(h, h_ref, r0, psi)`** = `r0 + (h − h_ref) cos ψ`
   (SPEC_02 v0.10 Step 3). `psi` is required; there is no default.
6. **`refrac.reduce`**: B1 with the frozen anchor tilt; `field_line_projection` measures what that
   neglects, with ψ at every level from `lib.gravity.g_eff_vector` at the projected radius,
   `phi_c` and `u(phi_c)`. `Reduction` gains `projection`.
7. **`refrac.product`**: `reduction_record` gains `radius_projection_rule`,
   `radius_projection_residual_m`, `latitude_drift_neglected_deg` and a `radius_projection_note`
   (decision 6). The `long_name` of `radius_m` and of `height_above_anchor_isobar_m` say what
   each now is.
8. **Acceptance suites updated** for the new values, other checks unchanged:
   - `step2`: check 3 reads the top row on the grid (19.9526 Pa) and as printed (20.0 Pa); check
     5 locates the first and last ammonia rows by printed pressure.
   - `step9`: kind T's variable set includes `pressure_printed_Pa`.
   - `step02_3`: levels off the exact decades (794.33, 0.20, 0.25 mbar) are found by printed
     pressure; check 3 expects `n` 1.041934e22 at the top and 1.270497e22 at the second level;
     check 4 expects `r0 − h_ref cos ψ` at 1 bar and `r0 + (h_top − h_ref) cos ψ` at the top.
   - `step02_1`: reads the hash table of this report as well.

## 2. Acceptance results

Run by `reports/step03_0/accept_step03_0.py`; full output in `reports/step03_0/output.txt`.
"Was" values are read from copies of every product taken before the rebuild
(`reports/step03_0/before/`).

| Check | Measured |
|---|---|
| 1. Raw bundle: both pressures, the flag, the `k` sequence | **Pass, as restated at SPEC_03 v0.6** (finding 1; against v0.5 it failed on the runs (6, 230), (4, 270)). `k` starts at −70; runs (spacing, last `k`): **(10, 130), (8, 194), (6, 236), (4, 268), (2, 310)**, as expected. Everything else in the check holds: 65 of 66 rows snapped, the last row as printed (129,848.0 Pa, `k` NaN, flag 0); top eight to six figures 19.9526, 25.1189, 31.6228, 39.8107, 50.1187, 63.0957, 79.4328, 100.0; the 1 bar row 100,000.0 Pa; `pressure_printed_Pa` equals the previous bundle's `pressure_Pa` bit for bit; temperature, ammonia and height unchanged bit for bit; the basis carries 48. `history`: 62 rows with one candidate, 3 resolved, 1 excluded; 0.32 mbar took −50, then 0.25 took −60, then 0.20 took −70, each from the fixed rows below. |
| 2. A CSV copy with one printed value no grid value rounds to is refused with the row named | **Pass.** 1000.00 edited to 1001.00 (candidates: none): "row 60 of lindal_table1.csv (1001.00 mbar): no value of the declared grid rounds to the printed 1001.00 mbar, and the row is not in excluded_printed_values_mbar. Refused rather than kept as printed." |
| 3. The excluded row moved onto the grid is written with 66 snapped rows | **Pass.** 1298.48 edited to 1318.26 (candidates `[312]`), exclusion list unchanged: 66 snapped, last row `k` 312, `pressure_Pa` 131,825.674 Pa. |
| 4. Kind T reads back with `pressure_grid_rule` and `pressure_printed_Pa` | **Pass.** Both present; `pressure_Pa` and `pressure_printed_Pa` equal the bundle's bit for bit. |
| 4b. Beyond the specification: T and D rebuilt, not kept | **Pass.** `created_at` 2026-09-11 to 2026-09-14 for both; `raw_bundle` hash `3054cc6c...` to `a2a276ba...`. |
| 5. Beyond the specification: kind C on kind T's levels | **Pass.** Levels equal bit for bit. `x_NH3` changes at two levels only: 1047.13 mbar (the interpolated gap) by −3.8e-6 relative, 1298.48 mbar (the extrapolated end) by +4.29e-5 relative, because the linear fill now runs through grid pressures. |
| 6. `radius_m` at the anchor equals `r0` exactly and at the top `r0 + (h_top − h_ref) cos ψ` to 1e-6 m | **Pass.** Anchor level 29: 58,516,188.288884744 m, equal to `r0`. Top: 58,801,571.048 m, difference from the formula 0.0 m. Top 58,801,571.0 m (was 58,802,888.3, change −1,317.2); bottom 58,412,566.6 m (was 58,412,088.3, change +478.3); `cos ψ` 0.995406. `height_above_anchor_isobar_m` unchanged bit for bit. |
| 7. `n` at the top and second level; 1 bar unchanged; `N` moves with `p` | **Pass.** Top 1.041934e22 (was 1.044408e22); second 1.270497e22 (was 1.264485e22); 1 bar 5.373123528e25, bit for bit unchanged. `n_new / n_old` equals `p_grid / p_printed` to 2.2e-16 at every level; `N` likewise at the 64 levels where `R_bar` is unchanged. At 1047.13 and 1298.48 mbar `R_bar` moves by +6.1e-11 and −3.4e-9. Largest pressure change 1.18e-2, at 0.32 mbar. |
| 8. The record carries the projection rule, its residual and the drift | **Pass.** `radius_projection_rule = "cos psi at the anchor"`; residual **12.44 m** (about 12); drift **[+0.02694, −0.00975] deg** (about +0.027 and −0.010); ψ from 5.475835 deg at the bottom to 5.546374 at the top. The tolerances, 3 m and 0.0015 deg, are this script's, stated in the output. |
| 9. `r0`, `phi_c`, `psi`, polar radii, asymmetry and every `reduction_record` number unchanged | **Pass, bit for bit.** `r0` 58,516,188.288884744 m; `phi_c` 30.80556842739218 deg; ψ 5.494431541512493 deg; `r0` uncertainty 11,345.547 m; `phi_c` uncertainty 0.109191 deg; north 54,420,643.464 m, south 54,449,383.835 m, asymmetry 28,740.370 m. 57 record attributes compared, none differ; the commit and the input hashes not compared; four added, none removed. `N` fractional uncertainty 0.023318 at every level. |
| 10. F2's gravity along the profile changes by under 6e-5 (restated at SPEC_03 v0.6) | **Pass** (finding 2; against v0.5's 5e-5 it failed). Largest `\|g(r_new) / g(r_old) − 1\|` = **5.284e-5**, at the top level. |
| 11a. Beyond the specification: the wind file identical in content | **Pass.** Rebuilt because it records the raw bundle hash; identical after dropping the writer and provenance attributes. |
| 11b. Beyond the specification: commit and hash of every file now | **Pass.** Provisional; section 7. |

The candidate's F1 to F4 were rendered to `reports/step03_0/candidate/figures/` and looked at: F2
and F4 are as before apart from the moved radii. F2 and F4 from the clean product of the sweep
(footer `2149b64`, SHA-256 `64c5d01a...`) are attached as REVIEW_03_step0 asked:
`reports/figures/step03_0_lindal_diag_F2_gravity_profile.png` and
`reports/figures/step03_0_lindal_diag_F4_product.png`.

## 3. Findings

REVIEW_03_step0 confirmed findings 1 to 3 and restated the specifications (SPEC_03 v0.6, SPEC_01
v0.24): the `k` runs end at 236 and 268, the gravity bound is 6e-5 with the decomposition stated,
and the keep-by-content sentence names `raw_bundle`. Finding 4 is recorded with no action;
finding 5 is now the rule of SPEC_03 v0.6 section 0. The acceptance was rerun in full after the
review: thirteen of thirteen (`reports/step03_0/output.txt`).

1. **The `k` runs in the specification are off by one step at two boundaries.** SPEC_03 v0.5
   Step 0 (expected values) and SPEC_01 v0.23 Step 2 (acceptance) say 6 apart to 230 and 4 apart
   to 270. Table I says 6 apart to **236** and 4 apart to **268**:
   - 199.53 mbar is `k` 230 and 229.09 mbar is `k` 236 (`10^2.36` = 229.087), each with exactly
     one candidate; 251.19 mbar is `k` 240, so the step of 6 ends at 236.
   - 478.63 mbar is `k` 268 (`10^2.68` = 478.630) and 501.19 mbar is `k` 270, so the step of 4
     ends at 268 and 268 to 270 is already a step of 2.

   The rule, the three resolved rows, the top eight values and the basis string ("10, 8, 6, 4, 2
   hundredths of a decade in succession") are all right; only the two endpoints quoted for the
   runs are wrong. Proposed: the expected sequence reads 10 apart to 130, 8 to 194, 6 to 236, 4 to
   268, 2 to 310, in both specifications. The check was not changed.

2. **The gravity bound of 5e-5 is the Newtonian estimate; the effective gravity moves more.**
   The top level drops 1,317.24 m, `dr/r` = −2.240e-5. The Newtonian part changes by 4.488e-5,
   which is `−2 dr/r`. Measured against `g_eff`, which is smaller than `g_N` by the factor 1.1183,
   that becomes 5.019e-5, and the rigid and wind centrifugal part, which falls with `r`, adds
   2.65e-6: total **5.284e-5**, at the top level where the radius moves most. Measured
   independently of the acceptance script with `lib.gravity.g_newton` and `g_eff_radial`.
   Proposed: "under 6e-5", or the expected value 5.3e-5 stated with this decomposition. The check
   was not changed.

3. **The keep-by-content rule needed no change to rebuild T and D.** SPEC_03 Step 0 asks for
   "content compared including `input_hashes`". Kinds T and D already carry the raw bundle hash in
   the `raw_bundle` attribute, which is compared as content, so a new raw bundle rebuilds both
   (check 4b). Adding `input_hashes` to the comparison would also rebuild them on every edit of
   `lindal_build.toml`, because `input_hashes` lists the control file; that is what SPEC_02 Step 6
   decision 4 excluded it to prevent. Proposed: the specification sentence says the raw bundle
   hash is compared through `raw_bundle`, and `input_hashes` stays excluded.

4. **The ammonia fill moves at the two filled levels** (check 5). Linear interpolation in
   pressure through grid values instead of printed values moves the interior gap at 1047.13 mbar
   by −3.8e-6 relative and the downward extrapolation at 1298.48 mbar by +4.3e-5 relative (79.253
   to 79.256 ppm). The effect on `R_bar` is below 4e-9. Recorded because the specification lists
   the levels whose pressures change but not this.

5. **The refusal of `-dirty` inputs and the Step 0 acceptance cannot both hold before a commit.**
   Every input the step rebuilds changes content, so `refrac` refuses them until they are
   rebuilt on a clean tree, and the specification's acceptance needs a product built from them.
   Resolved for this step by the author's ruling (decision 1). A later step that changes an input
   will meet the same thing; the specification may want to say how such a step is accepted.

## 4. Decisions

1. **The candidate product is built in memory** (author ruling, option 1). `reports/step03_0/
   relaxed.py` wraps `lib.control._load_into_memory` in the acceptance process only, renaming the
   `-dirty` suffix of each input's commit to "-dirty (relaxed in the Step 03_0 acceptance)". The
   manifest and the six inputs are copied to `reports/step03_0/candidate/` and
   `refrac.product.build_product` writes the candidate there. No code, rule or file under
   `occul_data/` is affected. The four relaxed inputs are listed in the output: thermo,
   composition, geodesy, wind.
2. **The snap, where the rule leaves room.**
   - Candidates are searched over `k` from `floor(100 log10(p − half unit))` to
     `ceil(100 log10(p + half unit))`, widened by one each side so float logarithms cannot drop an
     end; the grid value is rounded from its exact binary value, half to even (no grid value lies
     on a half unit of Table I's precision, so the tie rule is never exercised).
   - "The two nearest fixed rows on its fixed side" is read as the two adjacent rows on that side,
     both fixed. Each pass settles the first resolvable row in increasing pressure, then starts
     again.
   - The grid value in Pa is the float `10^(k/100)` mbar times 100 in exact decimal, then rounded
     once, as the tool already converts printed values. The exact decades come out exact
     (10,000.0 and 100,000.0 Pa).
   - The rule string must read `10^(k/{denominator}) {unit}` with the unit Table I prints, or the
     declaration is refused; a transcription without `[pressure_grid]` is refused.
3. **`pressure_grid_index` is carried in the raw bundle** so the `k` sequence the acceptance
   checks is in the file, float64 with NaN for the excluded row rather than an integer with a
   fill value.
4. **`pressure_printed_Pa` in kind T has no uncertainty companion**: it is the printed number,
   and the companion stays on `pressure_Pa`. The kind T schema version is not bumped, since the
   variable is optional and every existing file still validates.
5. **Only products that read the raw bundle were rebuilt in the working tree**: the raw bundle,
   wind, composition, T and D (the manifest rewrite gave the same hash, `97e542d6...`). Gravity and
   rotation do not read the raw bundle and keep `ece58d2`. **Question for the review:** SPEC_03
   Step 0 says "every product under `occul_data/lindal/` is rebuilt" at the sweep. Rebuilding
   gravity and rotation would change only their writer attributes and hashes. Proposed: keep them.
6. **The projection residual and drift, defined.** The residual is the largest magnitude over the
   levels of the difference between radii integrated along the local vertical with ψ at each
   level (`dr = cos ψ dh`) and radii projected with the anchor ψ; the drift is `sin ψ dh / r`
   integrated from the anchor, reported at the top and the bottom. Both by the trapezoid on the
   tabulated levels. A `radius_projection_note` in the record says so.
7. **A guard I wrote was wrong and was fixed before the acceptance ran.** The first rebuild of T
   and D refused a good raw bundle because `"scalars/pressure_grid" in tree` does not see nested
   groups of an `xarray.DataTree`; the check now looks the group up and catches the `KeyError`.

## 5. Suites not run before the commit

Steps 02_1, 02_4, 02_5 and 02_6 are run after the acceptance commit and the sweep:
- 02_1 compares input and manifest hashes with the report tables, which this report fills only at
  the sweep, and includes a `-dirty` refusal case the relaxation would defeat;
- 02_4, 02_5 and 02_6 run `casspian-refrac` in a subprocess on the on-disk inputs, which refuses
  them until they are clean, and 02_4 to 02_6 rewrite `occul_data/lindal/lindal_refractivity.nc`.

They were run after the sweep; results in section 8.

## 6. Manuscript notes carried (SPEC_03 decision 3)

For the handoff's list: the notation list should stop saying z and r are interchangeable; B1
gains `cos ψ` and a sentence on the field-line altitude; B2 is restated along the local vertical
with the magnitude of the effective gravity, or kept radial with the levels placed on the field
line, but not the present mixture; the latitude drift along the field line (+0.027 deg at the top,
−0.010 at the bottom) is mentioned once; B3.1 or B2 states the pressure grid inference and that
Table I is printed to two figures in its top decade; the Draft 8 passage in B3.2 arguing for polar
anchoring is rewritten to SPEC_02 decision 9.

## 7. Hashes

The sweep, after the acceptance commit `2149b64`, on a clean tree (`reports/step03_0/sweep.sh`,
log in `reports/step03_0/sweep.txt`): `casspian-lindal-raw`, the wind and composition tools,
`casspian-lindal-inputs` (T and D rebuilt, the manifest rewritten with the same hash), then
`casspian-refrac` with its figures. Every rebuilt product carries commit `2149b64` and none
`-dirty`; gravity and rotation are kept at `ece58d2` (REVIEW_03_step0, decision 5). These rows
supersede REPORT_01_step9 section 6, REPORT_01_step8 and REPORT_02_step6 section 7 for every file
they name, and the Step 02_1 suite reads them. The provisional working tree hashes of the review
are superseded.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_reduction.toml` | manifest | `97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee` |
| `lindal_raw.nc` | raw | `07e0727c9ef07b468722ad9d62f565309dedfa28e7c6bcef84e395a9d428c713` |
| `lindal_thermo.nc` | thermo | `5705f1dc153b9c7b956e2414c21b3fe2e14e3e30ba2c768598402b1edbab73f4` |
| `lindal_geodesy.nc` | geodesy | `51c545a6c6e8843b81e82ba2bb3ebf13f93c44179e2b22eed4b98d952d2d100a` |
| `lindal_gravity.nc` | gravity | `6db0129cefd01bc9c151020e604cf4563a82381e61c97381d0a7ee383826911e` |
| `lindal_rotation.nc` | rotation | `a5c72017a423c50f2a952ddcbc4ccf777c1e94683fa6965ee96f7f43d949a424` |
| `lindal_wind.nc` | wind | `26fe83a3d90fa9eab97b82c94a89c923c13664a645621955fb0bd80434e0fd57` |
| `lindal_composition.nc` | composition | `6af9dbfbcc7599161dce5ad892716182ae96b47ce75f69c017f34121e9a23717` |
| `lindal_refractivity.nc` | refractivity | `64c5d01a98692fba49c264b890c88949503ee16e98c45abaa74095ac92f0a1ab` |

## 8. Regression

Run by `reports/step03_0/run_regression.sh` after the library, tool and `refrac` changes, on the
rebuilt working tree inputs; results in `reports/step03_0/regression.txt`, each suite's output
beside it. **All pass:** 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 7 and 9 for SPEC_02
Steps 2 and 3.
- The SPEC_01 suites ran as they are, `step2` and `step9` with the updates of section 1. `step8`
  check 3 (the gap at 1047.13 mbar, 15.9 +- 0.1 ppm) passes on the rebuilt composition.
- Steps 02_2 and 02_3 ran through `relaxed.py`, since `lib.control` refuses the `-dirty` inputs
  they load; both pin `mean_polar_radius` in memory as before. The Step 02_3 check 4 measures
  `r0 + (h_top − h_ref) cos ψ` under that rule.
- Steps 02_1, 02_4, 02_5 and 02_6: after the sweep (section 5).

No suite wrote into `occul_data/`; the SPEC_01 Step 3 suite writes its gravity file into its own
directory.

**After the sweep** (`reports/step03_0/post_sweep_suites.sh`, results in
`reports/step03_0/post_sweep.txt`), on the clean products of section 7: **Step 02_1 6 of 6,
Step 02_4 9 of 9, Step 02_5 7 of 7, Step 02_6 7 of 7.**
- Step 02_1 matched every input and the manifest against the rows of section 7.
- Step 02_6 first read **6 of 7**: its check 1 compared the hashes of kinds T and D with
  REPORT_01_step9 only, which the Step 0 rebuild supersedes, although the tool kept both files by
  content as it should. Its `recorded_hashes` now reads the rows of REPORT_01_step9,
  REPORT_01_step8, REPORT_02_step6 and REPORT_03_step0, a later row overriding an earlier one, as
  the Step 02_1 suite does; rerun, 7 of 7. The check's "was" figure for the manifest now shows the
  latest recorded hash, the same file.
- Steps 02_4 to 02_6 rewrite `lindal_refractivity.nc` and its figures from the working tree, and
  Step 02_5 recopies the committed Step 5 report figures. The swept product and figures were
  copied aside before the suites and restored after them; the product's SHA-256 after the restore
  is `64c5d01a...`, equal to section 7, with commit `2149b64`. `reports/figures/` was restored from
  git.

No em dash or en dash appears in any file written in this step.

## 9. Next step

Step 0 is accepted (`2149b64`) and swept. Step 1, `lib.geopotential`, starts when the author says
go.
