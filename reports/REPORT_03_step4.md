# REPORT 03, Step 4. `forward.production`, `casspian-forward`, the closure run, the product, F5 and F6

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 15 September 2026, refreshed 16 September 2026. Specification: `SPEC_03_Forward_Production.md`
v0.12, Step 4, with the v0.11 additions (the band clause as a count, the kind C `source_statement`
rule); refreshed on v0.13. Status: **accepted** (REVIEW_03_step4).

**Refreshed after the review (16 September 2026).** Findings 1 and 4 applied as ruled: the closure
statistics and check 5 on the v0.13 bins, and `date = "1981-08-26"` in the closure namelist. The
acceptance was rerun: **thirteen of thirteen pass.** Check 5, check 12 (`epoch` now 1981-08-26),
decision 5, decision 10 (new) and section 6 are refreshed; the text below the refresh notes is the
report as reviewed. Output of the rerun in `reports/step03_4/output.txt`; the reviewed run's output
kept as `output_v012.txt`.

As reviewed:

**Twelve of thirteen acceptance checks pass.** Check 5 fails, on the reading of "the top decade"
rather than on the physics (finding 1): read as the acceptance's partition, above 10 mbar, the largest
residual is 2.119e-3 against a bound of 1.5e-3; read as the specification's own expected value reads
it, above 1 mbar, it is 9.74e-4 and every bound passes. The check was left as it ran. As first run,
eleven of thirteen passed: check 12 compared a one-element array attribute as a list, and netCDF reads
it back as a scalar; the script now compares `atleast_1d`, and nothing in the product changed.

The closure reproduces the reviewing agent's expected values to the figures printed: the largest
residual above 1 mbar 9.74e-4 (9.8e-4), between 10 and 100 mbar 2.61e-3 (2.6e-3), below 100 mbar
1.71e-3 (1.7e-3), the bottom row −3.28e-3 (−3.3e-3), the mean below 10 mbar −4.90e-4 (−4.9e-4); and
both negative controls fail the mean bound as they must.

## 0. Before this step

| Commit | What |
|---|---|
| `1dfe483` | REVIEW_03_step3, SPEC_00 v0.19, SPEC_01 v0.28, SPEC_03 v0.12, `STATE.md`. |
| `6ae113e` | The Step 3 acceptance commit. |
| `5965d1e` | The Step 3 sweep record: every product and run input clean at `6ae113e`. |
| `145bd3c` | REVIEW_03_step3 addendum, the sweep verified; Step 3 closed. |

The status line of SPEC_03 v0.12 ("Step 4 proceeds after the sweep record when the author says go")
and `STATE.md` were read, and the author said go. The anchor and the run's four inputs are the clean
products of the Step 3 sweep, so no relaxation was needed; only the product this step writes carries
`-dirty`.

---

## 1. What was built

**`src/casspian/forward/production.py`**
- `produce(profile, inputs, gauge, p_b)`: pure orchestration on arrays and the loaded inputs, no file
  access. `u(phi_c)` from the wind's reference level by `refrac.anchor.wind_of_latitude`, the function
  the reduction uses (REVIEW_03_step1 decision 6); `|g_eff|`, `psi`, `Phi` by
  `lib.geopotential.geopotential_along_profile`; `R_bar` and `m_bar` from the run's kind C and its
  species group (`mean_properties`); `n = N / R_bar`, `rho` by `lib.hydrostatic.density`; the layer
  masses, `p` from the top and `T` by `lib.hydrostatic`. Returns a frozen `Production`.
- `profile_from_anchor(anchor)`: the `Profile` (radius, field-line height, refractivity, `phi_c`) a
  kind N anchor carries.
- `closure_statistics(...)`: the residual `p / p_tab − 1` at every level and its statistics by bin.
- `run(namelist_path)`: `read_run_namelist` and `load_run_inputs` (Step 3's checks, closure
  comparison included); `p_b` by `anchor_profile_top` (the top level, checked to be the least
  pressure); `produce` on the anchor's levels at the gauge level; kind `profile` written with the
  anchor under `anchors/lindal`, the four inputs under `inputs/<kind>`, `namelist` (`text`, `sha256`,
  `resolved`) and `production_record`; F5 and F6 rendered from the written file when the namelist
  asks. Returns a `RunResult`.
- `main`, the entry point `casspian-forward <namelist>`: takes the path, runs, prints what was
  written.

**`lib.schema`**: the kind C reader enforces the globals its `conditional_globals` declares for the
file's `composition_role` (`source_statement` for `reduction`), and refuses a `composition_role`
outside `reduction` and `forward` (SPEC_03 v0.11 Step 4, REPORT_03_step3 finding 3). No product
changes and nothing is rebuilt; every composition written carries `source_statement`.

**`pyproject.toml`**: the entry point `casspian-forward`; the package reinstalled.

## 2. Acceptance results

Run by `reports/step03_4/accept_step03_4.py`; full output in `reports/step03_4/output.txt`. The product
is `forward/lindal_closure/output/lindal_closure_profile.nc`, with F5, F6 and the combined PDF under
`output/figures/`.

| Check | Measured |
|---|---|
| 1. `casspian-forward` writes the product; it reads back as kind `profile` with the anchor and the four inputs embedded identical | **Pass.** Exit 0; F5, F6 and the combined PDF written by the driver. Reads back as a DataTree; `anchors/lindal`: 11 groups identical (`xarray.Dataset.identical`) to `lindal_refractivity.nc` on disk; `inputs/composition` 2 groups, `inputs/gravity`, `inputs/rotation`, `inputs/wind` one each, identical to the run's files. No `input_hashes` warning on read. |
| 2. `mean_refractivity_m3`, `mean_molar_mass_kg_mol`, `number_density_m3` equal the anchor's to 1e-14 | **Pass.** Largest relative differences 0.0, 0.0 and **2.2e-16**. |
| 3. `pressure_Pa[0]` equals `p_b` exactly and `p_b` is the anchor's top tabulated pressure | **Pass.** 19.952623149688797 Pa all three ways, and the global `boundary_pressure_Pa`. |
| 4. `Phi` is zero at the gauge level | **Pass.** Level 29 at 10,000.0 Pa, `Phi` 0.0; `Phi` from 2.852355e6 to −1.043743e6 m²/s². |
| 5. The residual bounds (refreshed, v0.13) | **Pass.** Above 2 mbar (the ten top levels): max **9.740e-4** (bound 1.5e-3). From 2 mbar down, bottom row excluded: max **2.607e-3** at 10 mbar (bound 3e-3). Bottom row (1298.48 mbar) **−3.281e-3** (bound 4e-3). Mean below 10 mbar **−4.895e-4** (bound 1.5e-3). Against the expected values: 2 to 10 mbar **2.119e-3** at 6.31 mbar (2.2e-3); 10 to 100 mbar 2.607e-3 (2.6e-3); below 100 mbar excluding the bottom row 1.714e-3 (1.7e-3). `production_record` carries the same numbers. As reviewed, on v0.12: **fail on the reading of "the top decade"** (finding 1). From 10 mbar down, bottom row excluded: max **2.607e-3** (bound 3e-3). "Top decade" read as above 10 mbar: max **2.119e-3** at 6.31 mbar (bound 1.5e-3), the failure; read as above 1 mbar: **9.74e-4**. Bottom row (1298.48 mbar) **−3.281e-3** (bound 4e-3). Mean below 10 mbar **−4.895e-4** (bound 1.5e-3), bottom row excluded −4.314e-4. `production_record` carries the same numbers. |
| 6. The temperature residual equals the pressure residual to 1e-12 | **Pass.** Largest difference **4.4e-16**. |
| 7. The negative controls fail the mean bound | **Pass.** One factor, `g_k (h_(k+1) − h_k)`: mean below 10 mbar **−5.072e-3** (about −5.1e-3), bottom row −7.842e-3 (about −7.8e-3). Fully radial, `g_k (r_(k+1) − r_k)`: mean **−9.632e-3** (about −9.6e-3), bottom row −1.240e-2 (about −1.2e-2). Both computed in the script only. |
| 8. F5 and F6 by the driver and by `casspian-plots` by hand, identical apart from the footer | **Pass.** Both 1650 x 750 px, identical above the 19 px time band. |
| 9. The residual against the envelope band, levels outside counted | **Pass** as a description. **21 of 66** levels lie outside the band (8.6e-4 to 1.57e-3): at 0.316 and 0.398 mbar; 5.0 to 15.9 mbar (five levels, up to +2.607e-3 at 10 mbar); 87.1 mbar (+2.03e-3); and twelve levels from 398 to 1298 mbar (−9.3e-4 to −3.28e-3, the bottom row). The full list with each residual and band value is in the output. |
| 10. The kind C `source_statement` rule | **Pass.** A reduction composition without `source_statement`: "composition_role is 'reduction', so the global attribute 'source_statement' is required: a reduction composition carries the source's own words (SPEC_00 section 6.2)". A forward composition without it: accepted. |
| 11. Beyond the specification: `produce` reproduces the product | **Pass.** `produce` on the loaded arrays gives every modeled variable bit for bit; `u(phi_c)` 2.16708095051181 m/s, as the record. |
| 12. Beyond the specification: season, companions, paths, record | **Pass** (after the script fix noted at the head). Refreshed: `epoch` **"1981-08-26"**, from the namelist's `date` (finding 4), and the resolved namelist carries the date (as reviewed, the SPEC_00 string "none: season declared as solar longitude only"); `solar_longitude_deg` 18.2; `anchor_solar_longitudes_deg` 18.2. Six companions, all NaN, `uncertainty_terms_unstated` "composition, anchor_radius, label_latitude". `input_hashes` relative (`../../../occul_data/lindal/lindal_refractivity.nc`, `../inputs/...`, `../lindal_closure.toml`). The closure comparisons "identical" for all four kinds; the projection notes copied from the anchor's record; the layer rules named. |
| 13. Beyond the specification: hash | **Pass**, provisional, section 7. |

F5 and F6 were rendered and looked at. F6 draws the residual against the band with the boundary and
the gauge marked, the gauge labeled 100 mbar from `gauge_isobar_Pa` (REVIEW_03_step3 ruling 6). In F5
the geopotential and height curves still lie on each other (finding 2). They carry the working-tree
footer; the figures attached for the author's acceptance by eye are rendered from the clean product at
the sweep (decision 8).

## 3. Findings

1. **"The top decade" has two readings, and only one of them agrees with the specification's own
   expected value.** Ruled: the specification's fault; the bins restated in pressure at SPEC_03 v0.13,
   check 5 passes on them (decision 10). The acceptance bounds the residual "within 3e-3 at every level from 10 mbar
   down, excluding the bottom row, within 1.5e-3 in the top decade". Read as a partition, the top
   decade is above 10 mbar, and there the largest residual is 2.119e-3 at 6.31 mbar, over the bound.
   The expected values say "within ±9.8e-4 in the top decade", which is the largest residual above
   **1 mbar**, measured 9.74e-4 (the rows Table I prints to two figures, the sense §0 gives the
   words). Under that reading every bound of the acceptance passes, but the levels from 1 to 10 mbar
   (largest residual 2.119e-3) fall under no bound. The check was run on the partition reading, as
   decided before the run (decision 5), and fails; it was not re-read after the result. Proposed: the
   top decade is above 1 mbar (bound 1.5e-3), and the 3e-3 bound runs from 1 mbar down, excluding the
   bottom row; measured, 9.74e-4 and 2.607e-3, both inside.
2. **F5: the geopotential and the height lie on each other.** On their two axes the curves coincide,
   since `Phi` is close to a constant `g` times `h`, so the height curve is not distinguishable. The
   gauge label now reads the declared 100 mbar. Deferred to this step by REVIEW_03_step3 ruling 6; left
   for the author's acceptance by eye rather than restyled unasked. Ruled: left for the author's
   figure pass. One option: draw the height as
   markers at the levels, or plot `Phi − g_anchor h` on the twin axis to show the departure.
3. **SPEC_00 §7.2 says `forward` writes "the fully resolved namelist (absolute paths, hashes,
   defaults filled)", which v0.18 §5 forbids.** The product records the resolved namelist in the
   `namelist` group's `resolved` attribute with every path relative to the product and hashed, and
   writes no separate file into `output/` (decision 6). Proposed: §7.2 says relative paths. Ruled:
   withdrawn at SPEC_00 v0.20; decision 6 is the rule.
4. **The product's `epoch` is the SPEC_00 string, while every input is dated 1981-08-26.** Deliverable
   3 dates the product by the namelist's `date`, or the string for a season declared without a date;
   the closure namelist declares no `date`. Proposed, for the author: `date = "1981-08-26"` in the
   closure namelist, so that the closure product carries the date of the observation it closes on, as
   SPEC_00 v0.19 does for every other file of the run. Adopted and applied (SPEC_03 v0.13); the
   product's `epoch` is 1981-08-26.
5. **21 of 66 levels lie outside the envelope band.** As v0.11 anticipated: the band is the leading-
   order budget (height rounding over the scale height), and the residual between 5 and 15 mbar and
   below 400 mbar exceeds it. Recorded, not a bound.

## 4. Decisions

1. **`produce` takes the loaded inputs as an object with `composition`, `gravity`, `rotation` and
   `wind`** (the `RunInputs` of Step 3 serves), and a `Profile` dataclass for the anchor's arrays, so
   that SPEC_04 passes a transferred profile without a file.
2. **`R_bar` and `m_bar` are formed in the order the species group lists the species**, with
   `lib.reduction.mean_over_species`, the arithmetic the reduction used; they equal the anchor's to
   round-off (check 2).
3. **`p_b` by `anchor_profile_top` is the anchor's level 0**, and the run refuses an anchor whose level 0
   is not the least pressure; `boundary_level_index` is 0.
4. **The uncertainty companions** carry `uncertainty_terms_unstated` from the anchor's own term lists
   (`refractivity_uncertainty` and `radius_uncertainty_m`), "composition, anchor_radius,
   label_latitude", with the reason in `uncertainty_note`.
5. **The closure statistics bins**, decided before the run: above 10 mbar `p < 1000 Pa`, above 1 mbar
   `p < 100 Pa`, 10 to 100 mbar `1000 <= p < 10000 Pa`, below 100 mbar `p >= 10000 Pa`; a level on an
   edge belongs to the deeper bin; the bottom row is the deepest level. Both readings of the top decade
   are recorded in `production_record`, and the acceptance was run on the partition reading (finding 1).
   Superseded by v0.13's bins (decision 10).
6. **The resolved namelist** is an attribute `resolved` of the `namelist` group: the run's name, mode,
   season and date, the anchor, inputs and namelist with relative paths and hashes, `p_b` rule and
   location, the gauge isobar and its level, and the diagnostics as filled (finding 3).
7. **`production_record`** carries: mode, run, anchor slug; gauge isobar and level; `p_b`, its rule,
   location and level; `u(phi_c)`; the geopotential and hydrostatic rules; the closure comparison rule
   and, per kind, its result and the attributes dropped per group; the anchor's projection rule,
   residual, drift and note (as `anchor_*`); the residuals of pressure and temperature at every level
   and their statistics; the CODATA release, version and commit.
8. **The figures attached to this report** for the author's acceptance by eye are rendered from the
   clean product at the sweep, as for Step 0; the working-tree figures carry `-dirty` in their footer.
9. **The product's global `mode`** is written beside the SPEC_00 globals, so a reader sees the run mode
   without opening the record.
10. **The v0.13 bins (refresh).** `closure_statistics` records the largest residual above 2 mbar,
    2 to 10 mbar, 10 to 100 mbar, below 100 mbar excluding the bottom row, and from 2 mbar down
    excluding the bottom row (the acceptance's bound), the bottom row's value and the mean below 10
    mbar, with the count of levels above 2 mbar; the two readings of the top decade are no longer
    recorded. **The 2 mbar edge is the grid level `10**2.3` Pa (199.526 Pa)**, not 200 Pa: that is the
    row Table I prints as 2.00 mbar, so it lies on the edge and belongs to the deeper bin, and the bin
    above 2 mbar holds the ten top levels the specification counts. With 200 Pa as the edge the row
    would fall above 2 mbar and the bin would hold eleven; its residual is −1.022e-4, so no bounded
    value changes under either reading. The 10 and 100 mbar edges are grid levels already.

## 5. Regression

`lib.schema` changed, so the full regression was run after the acceptance
(`reports/step03_4/run_regression.sh`, results in `reports/step03_4/regression.txt`) on the clean
products of the Step 3 sweep, no relaxation. **All pass:** 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01;
6, 7, 9, 9, 7, 7 for SPEC_02 Steps 1 to 6; 9, 8 and 16 for the SPEC_03 Step 1, 2 and 3 acceptances.

Steps 02_4 to 02_6 rewrite `lindal_refractivity.nc` and its figures and 02_5 recopies the committed
Step 5 report figures, so the swept product and figures were copied aside and restored: SHA-256
`d2886aed...` before and after, the value REPORT_03_step3 section 7 records; `reports/figures/`
restored from git. The kind C `source_statement` rule refused no product of the chain or of the run.

No em dash or en dash appears in any file written in this step.

## 6. Hashes

Provisional, from the working tree, before the acceptance commit; the clean product's hash is
recorded at the sweep (section 8).

| File | Kind | SHA-256 (provisional) |
|---|---|---|
| forward/lindal_closure/output/lindal_closure_profile.nc | profile | b3ea78219b772e859a3e08d5bfd26d233f1cc8c37f2970063ab1d0510cf1f0f6 (`f9935b7-dirty`, refreshed) |

As reviewed: `b162a836de80d2a4fccb86c3303ac1ef2cd1adc7afbd97dcb7ac17bf7157224e` (`145bd3c-dirty`).
The product changed with the namelist's `date` (its `epoch`, the resolved namelist, the namelist hash)
and the v0.13 statistics in `production_record`. The four inputs and the anchor carry the hashes the
reviewed run recorded, and every residual statistic reproduces the reviewed value to the figures printed.

The product is under `output/`, which `.gitignore` excludes; it is regenerable from the namelist.

## 7. Next step

The review. Nothing is committed until the review is in and the author says go. After it: the
acceptance commit, the closure rerun on the clean tree with F5 and F6 attached, and the product's hash
recorded.
