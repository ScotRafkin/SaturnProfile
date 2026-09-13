# REPORT 02, Step 6. Anchoring the geoid on the equatorial radius

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.8, Step 6, with
the `SPEC_01_Lindal_Tool_Chain.md` v0.20 amendments to Steps 5 and 9 and
`SPEC_00_Architecture_and_Data_Files.md` v0.15 section 7.1. Status: **accepted
(REVIEW_02_step6, 13 September 2026). The Lindal default anchor rule is `equatorial_radius`**
(SPEC_02 v0.9 decision 9). SPEC_02 closes at this step's acceptance commit and sweep.

**Seven of seven acceptance checks pass**, with check 5 restated at SPEC_02 v0.9 as the `r³`
scaling of finding 1. As first run against v0.8 it read six of seven: the unscaled asymmetry
differs from Step 7 by -4.674 m.

Kinds T and D were verified and kept, with their hashes and commit `ece58d2` unchanged. The
manifest was rewritten under the new default. The product was rebuilt clean at the sweep after the
acceptance commit, and its hash is in section 7.

## 0. Before this step, per REVIEW_02_step5

| Commit | What |
|---|---|
| `2b2f347` | SPEC_02 Step 5 with the two figure changes: F4's fractional uncertainty on a top axis, F1's inset label. The inset was also raised clear of the panel's tick labels. SPEC_00 v0.15, SPEC_01 v0.20 and SPEC_02 v0.8 were committed with it. |
| `13b8bf4` | The sweep. `lindal_refractivity.nc` was rebuilt clean at `2b2f347`, SHA-256 `a91f75a8...`, recorded in REPORT_02_step5 section 7, superseding `327f8141...`. The attached figures were rendered from it. |

The status lines of SPEC_02 v0.8 and `STATE.md` were read before starting.

---

## 1. The comparison table

From `reports/step02_6/accept_step02_6.py` check 6. Each column is the whole anchor fixed point
under that rule, on the committed inputs. The uncertainty rows are the full Step 3 propagation,
run for the two rules the author is choosing between.

| quantity | `mean_polar_radius` | `north_pole` | `south_pole` | **`equatorial_radius`** | expected for `equatorial_radius` (v0.8) |
|---|---|---|---|---|---|
| `r0` (km) | 58,519.883 | 58,537.683 | 58,502.116 | **58,516.188** | 58,516.0 to 58,516.5, inside |
| `phi_c` (deg) | 30.804949 | 30.801965 | 30.807925 | **30.805568** | shift under 0.002 deg: +0.000619, inside |
| north polar radius (km) | 54,423.627 | 54,438.000 (anchor) | 54,409.277 | **54,420.643** | near 54,420.6, inside |
| south polar radius (km) | 54,452.373 | 54,466.768 | 54,438.000 (anchor) | **54,449.384** | near 54,449.3, inside |
| mean polar radius (km) | 54,438.000 (anchor) | 54,452.384 | 54,423.639 | **54,435.014** | 54,434.7 to 54,435.0, **14 m above**, from the `r³` scaling (finding 1) |
| polar asymmetry (km) | 28.7450 | 28.7676 | 28.7226 | **28.7404** | v0.9: 28.7403 (Step 7 scaled by `r³`) to 1 m, **0.056 m**, inside; v0.8's unscaled 28.745 was 4.7 m off (finding 1) |
| equatorial radius (km) | 60,370.999 | 60,390.265 | 60,351.770 | **60,367.000** (anchor) | |
| fixed point iterations | 5 | 5 | 5 | 5 | at most 8 |
| anchoring residual (m) | -1.7e-6 | 0 | -2.6e-5 | **-1.05e-4** | below 1e-3 |
| `dr0/dr_anchor`, total | 1.237328 | | | **0.923947** | 0.85 to 0.97, inside |
| anchor term on `r0` (km) | 12.373 | | | **3.696** | 3.4 to 3.9, inside |
| label term on `r0` (km) | 10.729 | | | **10.727** | about 10.73, inside |
| `radius_uncertainty_m` (km) | 16.377 | | | **11.346** | 11.2 to 11.4, inside |
| `phi_c` uncertainty (deg) | 0.109208 | | | **0.109191** | near 0.109, inside |

Lindal's Fig. 9 caption gives the mean polar radius as 54,438 +- 10 km. The equatorial anchor
predicts 54,435.01 km, 2.99 km below it and well inside his +-10 km. It predicts the south pole
28.74 km above the north; the caption says "of order 10 km".

## 2. The wind scaling diagnostic

From check 7. The whole reduction was rerun on in-memory copies of the inputs, with the kind W wind
scaled by 0.95 and by 1.05; all three components were scaled, so the poles stay at zero and the
sum identity holds. Report only: no file is written.

| rule | `u` x | `r0` (km) | shift | `phi_c` (deg) | shift | north (km) | south (km) | equatorial (km) |
|---|---|---|---|---|---|---|---|---|
| `equatorial_radius` | 0.95 | 58,518.789 | +2.601 | 30.805245 | -0.000324 | 54,426.198 | 54,453.492 | 60,367.000 |
| `equatorial_radius` | 1.05 | 58,513.582 | -2.606 | 30.805893 | +0.000325 | 54,415.082 | 54,445.270 | 60,367.000 |
| `mean_polar_radius` | 0.95 | 58,516.506 | -3.377 | 30.805627 | +0.000678 | 54,424.354 | 54,451.646 | 60,364.530 |
| `mean_polar_radius` | 1.05 | 58,523.265 | +3.381 | 30.804270 | -0.000679 | 54,422.900 | 54,453.100 | 60,377.481 |

**Half range of `r0` over `u` x 0.95 to 1.05: 2.60 km under `equatorial_radius`, 3.38 km under
`mean_polar_radius`.** The specification expected "a few kilometers" and "under one kilometer"
respectively. The first is met. The second is not, and the two move in **opposite senses**.
Finding 2.

---

## 3. What was built

1. **`lib.geoid.wind_geoid`, the `equatorial_radius` rule** (SPEC_01 v0.20 Step 5).
   - `r_anchor` is the radius at the equator, read at an exact node. The equator is now always
     a node, `0.0`, and the dense grid's floating point neighbor of zero (about 1e-14 rad) is
     dropped, so no rule reads a radius beside it.
   - `WindGeoidResult` gains `equator_radius_m` and `anchor_node_latitude_rad`: 0.0 for this
     rule, the pole for the pole rules, None for the two-node mean rule.
   - **The first secant correction now scales the north polar start instead of shifting it**
     (decision 2).
2. **`lib.geoid.reference_geoid` gains `anchor_latitude`**, default the pole. The no wind surface
   passes through the anchor radius where the rule anchors it. `refrac.anchor` and F3 pass the
   equator under `equatorial_radius`. Without this the no wind reference would have treated
   60,367 km as a polar radius (decision 3).
3. **`lib.control`** (SPEC_00 v0.15 section 7.1).
   - It accepts `equatorial_radius` and exports `ANCHOR_QUANTITY_FOR_RULE`, the one table of
     which kind D quantity each rule anchors on.
   - It refuses a rule with the other's quantity.
4. **`tools/lindal/build_inputs.py`** (SPEC_01 v0.20 Step 9).
   - It reads `geoid_anchor_rule` and `geoid_anchor_quantity` and refuses a pair that does not
     belong together.
   - It writes `[diagnostics]` from `diagnostics_figures`, `diagnostics_format` and
     `diagnostics_dpi`.
   - **It keeps kinds T and D when their content is unchanged** (decision 4). Its final line
     now says "ready" for each file, since a kept file was not written.
5. **`lindal_build.toml`**: `geoid_anchor_quantity = "radius_equatorial_m"`,
   `geoid_anchor_rule = "equatorial_radius"`, and the three diagnostics keys. The rewritten
   `lindal_reduction.toml` carries them.
6. **`refrac`**.
   - `FrozenAnchor` gains `equator_radius_m`.
   - The anchor rule spread is the full fixed point under **every other** manifest rule, each
     on its own quantity. Under `equatorial_radius` that is `mean_polar_radius`, `north_pole`
     and `south_pole`.
   - `reduction_record` adds `polar_radius_mean_m` and `equatorial_radius_marched_m`, with a
     `polar_radii_note` saying which of them are predictions under the rule in force.
   - `anchor_partials` needed no change. It differentiates with respect to whatever the manifest
     anchors on, and now takes `radius_equatorial_uncertainty_m`, 4 km `1sigma`.
7. **Earlier suites that the new default would have broken** (decision 5).
   - The SPEC_02 Step 2 and Step 3 suites pin `mean_polar_radius` in memory.
   - The Step 1 suite reads the new manifest hash from this report's section 7.
   - The Step 1 suite's anchor rule edit matches any rule.
   - The Step 5 suite adds `[diagnostics]` only to a manifest that lacks it.

## 4. Acceptance results

Run by `reports/step02_6/accept_step02_6.py`; full output in `reports/step02_6/output.txt`.

| Check | Measured |
|---|---|
| 1. `casspian-lindal-inputs` writes the new manifest, keeps T and D, refuses a mismatch | **Pass.** `[geoid]` reads `equatorial_radius` on `radius_equatorial_m`; `[diagnostics]` reads `{figures: true, format: png, dpi: 150}`. Manifest SHA-256 `97e542d6...`, was `d9080156...`. `lindal_thermo.nc` and `lindal_geodesy.nc` keep the hashes REPORT_01_step9 records and commit `ece58d2`. A control copy pairing `equatorial_radius` with `radius_polar_m` is refused with "do not belong together". |
| 2. `lib.control` accepts the pair and refuses the mismatches | **Pass**, four cases. `equatorial_radius` on `radius_equatorial_m` is accepted. Each rule on the other's quantity is refused. `mean_polar_radius` on `radius_polar_m` is accepted. |
| 3. The equator is a march node; residual below 1e-3 m | **Pass.** Anchor node latitude exactly 0.0. Marched radius there 60,366,999.999895 m against 60,367,000 m, residual -1.05e-4 m. The radius asked for at 0.0 equals the node value exactly. |
| 4. The product under `equatorial_radius` | **Pass.** `casspian-refrac` on the committed manifest converges in 5 iterations (30.852408, 30.8028834, 30.8057229, 30.8055596, 30.805569, 30.8055684). `anchor_isobar_radius_m` = 58,516,188.288884744 m equals an independent march at `phi_c` exactly. F1 to F4 and the combined PDF were rendered without being asked, from the manifest's `[diagnostics]`. The section 2.2 listing is complete. Product SHA-256 `53be01a7...` (working tree, `-dirty`). |
| 5. The asymmetry equals the Step 7 value scaled by the cube of the ratio of mean polar radii, to 1 m (restated at v0.9) | **Pass.** Expected 28,745.045 m x (54,435.014 / 54,438.000)^3 = **28,740.314 m**; measured **28,740.370 m**; difference **+0.056 m**. Against v0.8's unscaled wording the same measurement failed by -4.674 m (finding 1). |
| 6. The comparison table | **Pass**, reported in full in section 1. |
| 7. The wind scaling diagnostic | **Pass**, four runs, reported in full in section 2. |

## 5. Findings

1. **The polar asymmetry is not independent of the anchor; it scales with the size of the
   surface.** SPEC_02 v0.8 expects "shape unchanged by the anchor, the asymmetry stays 28.745 km
   to 1 m". Under `equatorial_radius` the whole surface sits lower. The mean polar radius is
   2.99 km, or 5.49e-5 of itself, below the `mean_polar_radius` build. The asymmetry is 4.67 m,
   or 1.63e-4 of itself, smaller: **2.96 times the relative change in radius.**

   That is what Eq. B3 requires. The slope `r0 G_phi / g` carries the centrifugal and wind
   accelerations, which grow as `r`, divided by gravity, which falls as `1/r^2`. A term built
   from them scales close to `r^3`. The pole rules show the same pattern: `north_pole` anchoring
   raises the surface and the asymmetry grows to 28.768 km, `south_pole` lowers it and it
   shrinks to 28.723 km.

   The same effect puts the predicted mean polar radius 14 m above the expected range, which
   assumed the asymmetry unchanged. Nothing here is a numerical error. The residual is 1e-4 m, and
   the scaling is measured across four independent fixed points. Proposed, and **adopted at
   SPEC_02 v0.9**: the acceptance says the asymmetry scales with the mean polar radius to the
   third power, to 1 m. That gives 28,745.045 x (54,435.014 / 54,438.000)^3 = 28,740.314 m
   against the measured 28,740.370 m, a difference of 0.056 m, and check 5 passes as restated.

2. **Polar anchoring is not better conditioned against the wind; it is worse.** The expectation
   was a few kilometers under `equatorial_radius`, because the march crosses the jet, and under
   one kilometer under `mean_polar_radius`. Measured: **2.60 km and 3.38 km, with opposite
   signs.**
   - **Equatorial anchor.** The equator is held and a stronger wind raises the bulge relative to
     the poles, so the mid latitudes come down: `r0` falls 2.6 km for +5 percent.
   - **Polar anchor.** The mean of the poles is held and the same stronger bulge pushes the
     equator up by 6.5 km and the mid latitudes up with it: `r0` rises 3.4 km.

   The profile latitude, 30.8 deg, sits nearer the equator in the dynamical height than the pole
   does. At the anchor the dynamical height is 65 km of the equator's 127 km, so a wind change
   moves it relative to the equator by less than relative to the poles. The spec's reasoning
   counted the wind the march integrates between the anchor and the profile, but under polar
   anchoring the march also crosses the jet, from the south pole to the north, to fix the mean.
   On this measure too the equatorial anchor is the better conditioned one.

3. **For the decision, the three numbers side by side.**

   | | `mean_polar_radius` | `equatorial_radius` |
   |---|---|---|
   | declared uncertainty propagated to `r0` | 16.38 km | **11.35 km** |
   | `r0` shift for a 5 percent wind error | 3.38 km | **2.60 km** |
   | spread of `r0` over the other three rules | -17.8 to +17.8 km | -13.9 to +21.5 km |

   The spread under `equatorial_radius` is relative to 58,516.188 km: `mean_polar_radius` +3.70
   km, `north_pole` +21.49 km, `south_pole` -14.07 km. Against Lindal's own numbers, the
   equatorial anchor predicts his mean polar radius to 3.0 km, inside his +-10 km. The polar
   anchor predicts his equatorial radius 4.0 km high, the Step 7 margin, inside his +-4 km only at
   its edge. On every measure in this table the equatorial anchor is the better conditioned. The
   caption's word "mean" supports reading his polar number as a mean; it does not make that
   number the better anchor. **The decision is the author's.**

## 6. Decisions

1. **The asymmetry check was left failing and reported, not loosened.** The failure is a
   property of the physics, and the report gave the scaling that would replace it (finding 1).
   The review adopted that scaling at SPEC_02 v0.9, and check 5 was restated and rerun.

2. **The first secant correction scales the start.** Shooting from the north polar start with an
   equatorial anchor, the old first step, `start - residual`, assumed the anchored radius moves
   one for one with the start, which holds at a pole and is about 10 percent wrong at the equator.
   With it, the secant stopped inside the 1 m `convergence_m` at a residual of -1.24e-3 m,
   above the 1e-3 m acceptance bound. The new first step, `start x r_anchor / outcome`, uses the
   near homogeneity of Eq. B3 in `r0` and gives -1.05e-4 m in the same five iterations.

   Under `mean_polar_radius` the frozen `r0` moves by **1.1e-6 m** (58,519,883.28568849 to
   ...28568964), because the secant takes a different path to the same tolerance. That is far
   below every recorded figure: Steps 2 to 5 quote `r0` to the millimetre or coarser. Their suites
   pass unchanged (section 7).

3. **The no wind reference is anchored where the rule anchors.** `reference_geoid` took its
   radius as polar. Under `equatorial_radius`, `refrac.anchor` would have passed 60,367 km as a
   polar radius, a no wind surface 5,900 km too large at the pole. The no wind reference is now
   anchored on the equator for this rule: `phi_c` 30.799334 deg, `r0` 58,567.166 km,
   `psi` 5.500666 deg. That is the no wind surface through the same equatorial radius, for the
   record. F3 draws the same surface.

4. **Kinds T and D are kept by content, not by hash.** The spec says they "are verified, not
   rewritten, unless their hashes differ". A rewrite would change the hash through `created_at`
   alone, so hashes cannot decide it. The tool builds each dataset in memory and compares it with
   the file on disk through `xarray.Dataset.identical`. It ignores the attributes that record the
   writing and not the content: the writer's globals, `history`, `control_file` and
   `input_hashes`, whose control file entry changes whenever a manifest choice does. Both files
   compare identical and are kept, so their SHA-256 and commit `ece58d2` stand. **Consequence:**
   their `control_file` attribute still names the `lindal_build.toml` hash they were written
   under, `f48930b1...`, not today's.

5. **Earlier suites follow the rule they accepted.** The SPEC_02 Step 2 and Step 3 suites test
   values accepted under `mean_polar_radius`, so they pin that rule and quantity on the manifest
   in memory with `dataclasses.replace`. The Step 4 and Step 5 suites run on the committed
   manifest. Step 4's expectations hold under either rule: the spread under the two pole rules
   is the same fixed point whichever rule is in force, and the product is compared with a fresh
   reduction under the same manifest. If the author reverts the default, no suite changes.

## 7. Hashes

The manifest and the two kept inputs, and the product. SPEC_02 Step 1's suite reads the manifest
row from here, superseding REPORT_01_step9 section 6. The product was rebuilt on a clean tree at
the sweep after the acceptance commit, by `casspian-refrac occul_data/lindal/lindal_reduction.toml`
with its figures, under `equatorial_radius`. This is the first hash of the product under the new
default; it supersedes `a91f75a8...` of REPORT_02_step5 section 7.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_reduction.toml` | manifest | `97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee` |
| `lindal_thermo.nc` | thermo | `ff9c12991399f94c08caaae80b1a6a91f45a45feaebc4d059ad563b45a720e85` |
| `lindal_geodesy.nc` | geodesy | `16e0c7d551c3026ce5623141b1a2f86187c013a2dbf929f3e330f79f410f787e` |
| `lindal_refractivity.nc` | refractivity | `888517ccdea0688ea9059f434f07b1384ee117eb57d3c58a5f45ea2094f82a12` |

## 8. Regression

The ten SPEC_01 suites and the SPEC_02 Step 1 to 5 suites were rerun after this step's changes (`reports/step02_6/regression.txt`). All pass: 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 6, 7, 9, 9, 7 for SPEC_02 Steps 1 to 5. The Step 2 and Step 3 suites ran with `mean_polar_radius` pinned in memory; the Step 1 suite matched the new manifest hash from section 7; the Step 4 and Step 5 suites ran on the committed manifest under `equatorial_radius`.

After the review only `reports/step02_6/accept_step02_6.py` changed, to restate check 5; no
library, tool or `refrac` code changed, so the regression above stands for the committed code. The
Step 6 acceptance was rerun in full: seven of seven.

The Step 5 suite's rerun during that regression had re-copied the seven figures attached to
REPORT_02_step5 from the product then on disk, built under `equatorial_radius` from the working
tree. The committed copies under `reports/figures/` were restored before the commit, so the Step 5
record still shows the clean `2b2f347` product it was accepted on.

No em dash or en dash appears in any file written in this step.

## 9. Next step

The author decided `equatorial_radius` (SPEC_02 v0.9 decision 9). The acceptance commit and the
sweep close SPEC_02. Nothing further starts until SPEC_03 is written.
