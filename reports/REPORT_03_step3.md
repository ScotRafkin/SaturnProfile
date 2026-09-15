# REPORT 03, Step 3. The season identifier, the run directory and its inputs, the run namelist, kind `profile`, and the reader warning

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 15 September 2026. Specification: `SPEC_03_Forward_Production.md` v0.10, Step 3 (deliverables
0 to 4), with `SPEC_00_Architecture_and_Data_Files.md` v0.18, `SPEC_01_Lindal_Tool_Chain.md` v0.26
and `SPEC_02_Refrac_and_Diagnostics.md` v0.11. Status: **accepted (REVIEW_03_step3, 15 September
2026)**, with the author's decision on finding 1 applied before the acceptance commit (SPEC_00
v0.19, SPEC_01 v0.28, SPEC_03 v0.12).

**Refreshed after the review.** Finding 1 was ruled: every file of the chain and of the closure run
carries `epoch = "1981-08-26"`, the occultation date, with the source's own dating in `epoch_note`;
the G and R build-file sections carry `epoch`, and the reader parses it. The rule was applied, the
chain rebuilt once more by `rebuild.sh`, and the acceptance rerun: sixteen of sixteen. Checks 1, 2,
15 and 16 and section 7 below are from that rerun; the F5 gauge marker is now labeled with the
declared `gauge_isobar_Pa` (ruling 6). The first run's results are superseded.

**Sixteen of sixteen acceptance checks pass**: the specification's checks (1 to 4, 6 to 9, 11 to
14) and four beyond it (5, 10, 15, 16). The whole reduction chain was rebuilt in the working tree,
gravity and rotation included; **no value in any product changed** (check 15), and kind N and its
`reduction_record` are unchanged bit for bit apart from the season attributes it gains.

The acceptance ran on a candidate kind N built in memory and on run inputs that carry `-dirty`, by
the SPEC_03 v0.6 section 0 procedure for a step that changes an input file. The registered products
are rebuilt clean at the sweep.

## 0. Before this step

| Commit | What |
|---|---|
| `3ba6d9c` | The Step 2 record: no product changed, `STATE.md` accepted. |
| `fbb0861` | REVIEW_03_step1 with order of work item 3 corrected. |
| `feef1ea` | The author's documents: SPEC_00 v0.18, SPEC_01 v0.26, SPEC_02 v0.11, SPEC_03 v0.10, `STATE.md`, `docs/CASSPIAN_Seasonal_Design_Note.md` v0.3. |

The status line of SPEC_03 v0.10 ("Step 3 proceeds") and `STATE.md` were read, and Step 3 in full
with deliverable 0 and the two ruling blocks of section 8 dated 15 September, before building.

---

## 1. What was built

**`lib`**
- `lib.io`: `recorded_path(path, product)`, the one helper every writer records a path through,
  relative to the product's directory with `/`, refusing a path on another drive;
  `input_hash_entry` and `input_hashes` now take the product as a required argument;
  `hash_entries`; `warn_changed_inputs`, called by `read` for every kind that carries
  `input_hashes`, resolving each entry against the file's own directory and silent when it does
  not resolve; kind N keeps its refusal on a missing recorded hash; kind `profile` must carry an
  `anchors/<slug>` group, checked on write and on read.
- `lib.schema`: `epoch` required on every kind; `_check_season`: exactly one of
  `solar_longitude_deg` (in [0, 360), with `solar_longitude_source`) and
  `season_absent_meaning = "uniform"`; kind `profile` registered (deliverable 3's variables,
  scalars, globals `boundary_pressure_Pa`, `anchor_solar_longitudes_deg`, `input_hashes`, and the
  groups), with `_check_profile`: the tabulated pair together, the geopotential coordinate strictly
  monotonic with `positive = "up"`, every `modeled` variable with its companion present and all NaN,
  a season (not `uniform`), role `forward`.
- `lib.control`: `build_role`; `read_run_namelist` (`RunNamelist`, `RunAnchor`); `load_run_inputs`
  (`RunInputs`); `check_closure_inputs` (`ClosureComparison`) with `CLOSURE_DROPPED_ATTRIBUTES` as
  ruled in section 8 ruling 2.

**Writers (every path relative to the product, SPEC_00 v0.18)**
- `casspian-lindal-raw`: `epoch` (the observation date), `solar_longitude_deg`,
  `subsolar_latitude_deg`, `solar_longitude_source` from `[source]`, refused when absent.
- `casspian-gravity-file`, `casspian-rotation-file`: required `role`; `season_absent_meaning =
  "uniform"`; the prose `epoch` kept (finding 1).
- `casspian-wind-from-curve`: required `role`; `epoch`, `epoch_note` and the season from the data
  properties file's `[epoch]` table, refused when absent; `latitude_conversion_inputs` relative.
- `casspian-composition-lindal`: required `role`, written as `role` and `composition_role`;
  `epoch` the observation date; `season_absent_meaning = "uniform"`; `raw_bundle`,
  `master_table`, `master_table_hash` relative.
- `casspian-lindal-inputs`: required `role` in `[stage_two]`; kinds T and D carry the season from
  the raw bundle; `raw_bundle`, `control_file`, `input_hashes` relative to each product.
- `casspian-refrac`: kind N carries `epoch`, `solar_longitude_deg`, `solar_longitude_source` from
  kind T; `input_hashes` relative to the product.

**Transcriptions and control files**
- `occul_data/lindal/raw/lindal_scalars.toml` `[source]`: `solar_longitude_deg = 18.2`,
  `subsolar_latitude_deg = 8.06`, `solar_longitude_source` in SPEC_01 v0.25's words.
- `data_static/winds/smith1982_fig4.toml` `[epoch]`: `value = "1981-08-25"`, `epoch_note`,
  `solar_longitude_deg = 18.2`, `solar_longitude_source` (finding 4).
- `occul_data/lindal/lindal_build.toml`: `role = "reduction"` in every section.

**The run directory (deliverable 1)**
- `forward/lindal_closure/lindal_closure.toml`, the namelist as deliverable 2 prints it.
- `forward/lindal_closure/lindal_closure_build.toml`: the four sections of `lindal_build.toml` with
  `role = "forward"`, the run prefix, outputs under `inputs/`, the wind pointed at the run's own
  gravity and rotation files and at the raw bundle (section 8 ruling 4).
- `src/casspian/tools/run/run_inputs.py`, entry point `casspian-run-inputs`.
- `forward/README.md` rewritten; `.gitignore` gains `forward/*/inputs/*.nc`; `pyproject.toml`
  gains the entry point (package reinstalled).

**Plots**
- `tools/plots/figures_profile.py`: F5 (geopotential against the hydrostatic pressure, height on a
  twin axis, boundary and gauge levels marked) and F6 (`pressure_Pa / pressure_tabulated_Pa - 1`,
  the envelope band, boundary and gauge marked; rendered only with the tabulated pair).
- `figures_product.py`: F5 and F6 removed; kind N renders F1 to F4 and reports nothing skipped.
- `render.py`: kind `profile` dispatched, with a combined PDF.

**Suites changed** (each other check unchanged)
- `step1/accept_step1.py`, `step1/verify_review_changes.py`: the synthetic files carry `epoch` and
  `season_absent_meaning`, now required.
- `step3/accept_step3.py`: the generated iess2019 control file declares `role`.
- `step02_5/accept_step02_5.py`: check 1 now requires nothing skipped; check 7 renders kind N and
  requires F1 to F4 with nothing skipped, the synthetic kind N with forward fields removed (the
  "F5 and F6 skipped" check retired, as the specification says).

## 2. Acceptance results

Run by `reports/step03_3/accept_step03_3.py` after `rebuild.sh` (the chain, the candidate kind N by
`build_candidate.py`, the run's inputs by `casspian-run-inputs`); full output in
`reports/step03_3/output.txt`, the rebuild log in `rebuild.txt`.

| Check | Measured |
|---|---|
| 1. Every product of the rebuilt chain and every run input validates with the season attributes, every one dated 1981-08-26 (v0.12) | **Pass**, twelve files, every `epoch` 1981-08-26. Raw bundle, T, D, candidate N: `Ls` 18.2. W (reduction and run): `Ls` 18.2, `epoch_note` "The curve combines Voyager 1 imaging (November 1980, Ls about 8 deg) and Voyager 2 imaging (late August 1981); it is assigned the occultation date". G and R (reduction and run): `uniform`, `epoch_note` "Pioneer 11 encounter, 1979" and "Voyager 1 and 2, 1980 to 1981". C (reduction and run): `uniform`. |
| 2. Neither, both, and a non-ISO epoch are refused | **Pass.** "carries neither solar_longitude_deg nor season_absent_meaning; SPEC_00 v0.17 section 5 requires exactly one, so that a missing season is never silent"; "carries both solar_longitude_deg and season_absent_meaning; SPEC_00 v0.17 section 5 allows exactly one"; (v0.11) "epoch is 'Voyager 1 and 2, 1980 to 1981'; SPEC_00 v0.19 section 5 requires an ISO 8601 date. A source's own dating goes in epoch_note." |
| 3. Kind N's season equals kind T's | **Pass.** `epoch`, `solar_longitude_deg` and `solar_longitude_source` identical. |
| 4. `sin δ_s = sin ε sin Ls` to 0.05° | **Pass.** `asin(sin 26.73° sin 18.2°)` = **8.0759°** against the transcribed 8.06°, difference **+0.016°**. |
| 5. Beyond the specification: `role` from the build file | **Pass.** Every reduction product `reduction`, every run input `forward`, both compositions' `composition_role` matching; a `[gravity]` section without `role`: "[gravity] is missing required key(s) ['role']". |
| 6. Only relative paths, scanned in every group | **Pass.** No drive letter and no leading `/` in any of 5, 4, 4, 2, 2, 9, 6 recorded paths of the reduction files, 34 in the candidate N (with its embedded groups), 6, 2, 2, 9 in the run inputs. Kind T records `raw/lindal_raw.nc` and `lindal_build.toml`; the run wind records `lindal_closure_gravity.nc` and `lindal_closure_rotation.nc`. |
| 7. The `input_hashes` warning | **Pass.** The edited kind T copy with `raw/lindal_raw.nc` beside it warns once ("input_hashes entry raw/lindal_raw.nc does not match the file now at that path beside it"); the same copy alone reads silently; the unedited kind T in place reads silently. |
| 8. The four run inputs | **Pass.** Each validates, `profile_or_run = "lindal_closure"`, role `forward`, prefix `lindal_closure_`; hashes differ from the reduction's (composition 2ae90640... against aa75adcd..., gravity 3492ef30... against f2b2ec29..., rotation fd0da383... against cbc8d1cd..., wind 9b567856... against dfa33001...). |
| 9. Content-identical to the embedded copies | **Pass**, all four. Dropped from the root of composition: `casspian_git_commit`, `composition_role`, `created_at`, `created_by`, `history`, `input_hashes`, `profile_or_run`, `raw_bundle`, `role`, `title`; from its `species` group: `master_table`, `master_table_hash`. Gravity and rotation: `casspian_git_commit`, `created_at`, `created_by`, `history`, `input_hashes`, `profile_or_run`, `role`, `title`. Wind: those, and `decomposition_geometry` and `latitude_conversion_inputs`. |
| 10. Beyond the specification: the driver refuses a build file that is not a run's | **Pass.** A `[stage_two]` section; a section with `role = "reduction"`; an output outside `inputs/`. |
| 11. The namelist loads and resolves | **Pass.** The committed namelist parses (mode closure, `Ls` 18.2, `p_b_rule` anchor_profile_top, gauge 10,000 Pa, product under `output/`). A case directory with the anchor at the candidate loads: gauge level **29**, all four closure comparisons identical. |
| 12. Thirteen refusal cases | **Pass**, each message below. |
| 13. Synthetic kind `profile` | **Pass.** Written and read back as a DataTree with 20 groups; refused without `temperature_tabulated_K`, with the coordinate not monotonic, and without `pressure_uncertainty_Pa`. |
| 14. `casspian-plots` | **Pass.** On the synthetic profile F5 and F6 and the combined PDF, nothing skipped; the envelope's pressure term zero at 65 of 66 levels and **3.85e-6** at the bottom row, the height term **8.6e-4 to 1.57e-3** (50 m over a scale height of 31.9 to 57.9 km). On the candidate kind N, F1 to F4 and the PDF, nothing skipped (CLI and `render`). |
| 15. Beyond the specification: the rebuild changed no value | **Pass.** Kind N: 24 variables and 61 `reduction_record` attributes unchanged bit for bit (commit and input hashes not compared). Every reduction file: no value differs; attributes added only (raw: `epoch`, `solar_longitude_deg`, `solar_longitude_source`, `subsolar_latitude_deg`; T, D: `epoch` and the season pair; G, R: `epoch_note` and `season_absent_meaning`, `epoch` changed from prose to 1981-08-26; W: `epoch_note` and the season pair, `epoch` changed; C: `epoch`, `season_absent_meaning`). Compared against the Step 0 products copied before the first rebuild, which the rerun kept. |
| 16. Beyond the specification: hashes | **Pass**, provisional, section 7. |

**The thirteen refusals (check 12), quoted:**
1. `[inputs]` missing: "required section [inputs] is missing."
2. `[target]`: "[target] is not accepted in closure mode: the target of a closure is the anchor's own latitude phi_c (SPEC_03 Step 3 deliverable 2)."
3. `[grid]`: "[grid] is not implemented in this specification (SPEC_03 Step 3 implements the closure subset of SPEC_00 section 7.2; SPEC_04 adds it)."
4. An unknown key: "[run] carries the unknown key 'flavor'. SPEC_00 section 7 makes an unknown key an error."
5. A second `[[anchors]]`: "closure mode takes exactly one [[anchors]] entry; this namelist has 2."
6. A slug that does not match: "[[anchors]] slug = 'cassini' but lindal_refractivity.nc carries profile_or_run = 'lindal'."
7. `p_b_Pa` with `p_b_rule`: "[hydrostatic_boundary] declares both p_b_rule and p_b_Pa; SPEC_00 section 7.2 v0.16 requires exactly one."
8. Gauge 9,000 Pa: "gauge_isobar_Pa = 9000.0 is not a tabulated level of the anchor; the nearest is level 28 at 8709.635899560806 Pa."
9. A `-dirty` anchor (the relaxation turned off for the case): "lindal_refractivity.nc (anchor) carries casspian_git_commit = 'feef1ea...-dirty'. SPEC_00 section 8: a file offered as an input to forward may not carry -dirty."
10. An input without the run prefix: "[inputs] wind = 'inputs/lindal_wind.nc' does not carry the run prefix 'lindal_closure_'."
11. A geodesy file: "[inputs] names a geodesy file. forward refuses one."
12. The composition with `x_He` +0.01 and `x_H2` −0.01: "the run's composition input is not content-identical to the anchor's embedded copy, so this is not a closure. First difference: composition /: variable x_H2 differs in its values. All differences: ... x_H2 ...; ... x_He ..." Only the two mole fraction edits were needed: the kind C reader carries no closure-declaration check (that check is `refrac.reduce.composition_closure`), so no declared share was edited (finding 3).
13. `solar_longitude_deg = 18.3`: "[run] solar_longitude_deg = 18.3 but the anchor lindal_refractivity.nc is at 18.2; not a closure."

F5 and F6 of the synthetic profile were rendered and looked at (`reports/step03_3/synthetic/figures/`). They carry a synthetic residual and a `-dirty` footer and are not attached; the attached figures are Step 4's, from a real production.

## 3. Findings

REVIEW_03_step3 ruled on each: finding 1 by the author's decision (one date for the chain, applied
and rerun, see the head of this report); finding 2 as intended; finding 3 to be enforced at Step 4;
finding 4 accepted; finding 5 restated in SPEC_03 v0.11; finding 6 deferred to Step 4's figures, the
gauge label changed now.

1. **SPEC_00 v0.17 asks `epoch` to be an ISO date; SPEC_01 v0.25 keeps the prose epochs of G and R.**
   The harmonic set's `epoch` is "Pioneer 11 encounter, 1979" and the rotation system's "Voyager 1
   and 2, 1980 to 1981"; SPEC_01 v0.25 Step 3 says their epoch is "as already recorded". The reader
   therefore requires `epoch` and does not parse it. Proposed: ISO dates transcribed for both in
   `data_static` with the present text as `epoch_note`, and the reader then parses `epoch` as a date
   or the single SPEC_00 string for a season declared without a date. The dates are the author's to
   transcribe; none was typed here.
2. **The registered kind N is refused by the new reader until the sweep rebuilds it.**
   `occul_data/lindal/lindal_refractivity.nc` is still the Step 0 product, which has no `epoch` and
   no season. The committed namelist points at it, so `load_run_inputs` on the committed namelist
   succeeds only after the sweep; the acceptance loaded a case directory pointed at the candidate
   (check 11). The same holds for any old product offered to the new reader, which is the intent of
   the rule.
3. **The kind C `source_statement` rule is declared and never enforced.** `lib.schema` carries
   `conditional_globals={"reduction": ("source_statement",)}` on kind C, but no code reads
   `conditional_globals` for kind C. The composition tool writes `source_statement` in both roles,
   so nothing is missing today. Not changed in this step; recorded for a later clean-up.
4. **The wind's season source is written by the coding agent.** SPEC_01 v0.25 names the values
   (1981-08-25, `Ls` 18.2) and the `epoch_note`, but not the wording of `solar_longitude_source` for
   the wind. It reads: computed by the reviewing agent from the same elements and pole, "the method
   that gives 18.2 deg for the Voyager 2 ingress of 1981-08-26; not stated by the source; to be
   replaced by the ephemeris value when the season tool exists". `Ls` changes by about 0.03° a day,
   so the one day between the two dates does not show at the precision given.
5. **F6's envelope needed a column-wise precision, found while looking at the first render.** Read
   value by value, a height printed as 90.0 km has no decimals once it is a float, and the first run
   drew a height term of up to 1.54e-2 (half a unit of 500 m). The precision is now the finest any
   value of the column shows, 0.01 mbar and 0.1 km for Table I; the envelope reads as in check 14.
6. **F5 presentation.** The geopotential and the height curves nearly coincide on their two axes, so
   the height curve is hard to see, and the gauge label reads the hydrostatic pressure at the gauge
   level (99.9668 mbar on the synthetic file) rather than the tabulated 100 mbar. For the author's
   appearance pass.

## 4. Decisions

1. **The acceptance ran by the SPEC_03 v0.6 section 0 procedure.** The candidate kind N was built
   in memory (`build_candidate.py`) under `relaxed.py`; the relaxation is named in the output and
   turned off for the `-dirty` anchor case.
2. **`recorded_path` refuses a path it cannot express relatively** (another drive) rather than
   recording it absolute.
3. **`epoch` is required and not parsed** (finding 1). **`solar_longitude_deg` must lie in
   [0, 360)**; kind `profile` must carry a season, not `uniform`.
4. **The raw bundle's `role` stays `reduction` in code.** `casspian-lindal-raw` has no build-file
   section to carry a `role` key; the bundle is the transcription of a source.
5. **Kind `profile` rules beyond the table**: every `modeled` variable's companion must be all NaN,
   not only present; role `forward`; an `anchors/<slug>` group checked on write and read. The
   latitude marker is `latitude_planetocentric_absent_meaning = "point"`, the SPEC_00 section 5
   form, where the specification wrote `latitude_absent_meaning`.
6. **The namelist**: the file must be named `<run>.toml` for its `[run] name`; `date`, when given,
   must be an ISO date; the output directory must be a directory of the run directory and the
   product carry the run prefix; `p_b_Pa` alone is refused in closure mode, where the rule is the
   only form accepted; `transfer` is refused as not implemented.
7. **The loader** refuses in this order: anchor missing, anchor `-dirty`, slug; each input's
   existence, `-dirty`, run prefix, role, composition role; wind against rotation; components and
   poles; coverage in pressure and latitude; season; composition levels; gauge level (exact, nearest
   named); the closure comparison. The anchor's season is compared with `float` equality.
8. **The closure comparison** drops the listed attributes and any attribute holding `sha256:` from
   variables as well as groups, and its refusal names every difference, not only the first (in case
   12 the first alphabetically is `x_H2`, and `x_He` is named in the list).
9. **The driver** checks that the file is `<run>_build.toml`, carries exactly the four sections,
   and that every section is `role = "forward"`, carries the run prefix, and writes under `inputs/`.
10. **The F6 envelope** is formed from the anchor's embedded kind T: the pressure term is zero on
    levels that lie on the grid parsed from `pressure_grid_rule`, and half a unit of the printed
    pressure column elsewhere; the height term uses the scale height `R T / (M g)` with
    `g = |dPhi/dh|` from the profile's own geopotential and height.
11. **The synthetic profile** uses the Step 1 geopotential on the candidate and the tabulated
    pressure with a synthetic residual of `5e-4 sin(k)`, not the Step 4 production.

## 5. Suites not run before the commit

Steps 02_1 (it compares input and manifest hashes with the report tables, which this report fills
only at the sweep), 02_4, 02_5 and 02_6 (they run `casspian-refrac` on the on-disk inputs, which
refuses them until they are clean, and 02_4 to 02_6 rewrite the registered product), and the SPEC_03
Step 1 and Step 2 acceptances (they read the registered kind N, which the new reader refuses until
the sweep rebuilds it, finding 2). They run after the acceptance commit and the sweep, and their
results will be added to section 8.

## 6. Manuscript and specification notes

The seasonal design note's decisions are carried by the files, not the manuscript, at this step.
For the review: finding 1 (ISO epochs for G and R), finding 3 (kind C `source_statement`).

## 7. Hashes

Provisional, from the working tree of the rerun (`1dfe483...-dirty`), for the record only; not in
the row format the Step 02_1 suite reads. The registered hashes are recorded here after the sweep.
The first run's provisional hashes are superseded.

| File | Kind | SHA-256 (provisional) |
|---|---|---|
| occul_data/lindal/raw/lindal_raw.nc | raw | f74df3a51fbb3b30468c0d96ccef02818f283c7b8f2fff923f3bbed929367414 |
| occul_data/lindal/lindal_thermo.nc | thermo | 9bc052d0e53327d1db0a5fcbb1d1c4b78e862cb06a364ef0e1c2014840c8df5f |
| occul_data/lindal/lindal_geodesy.nc | geodesy | 6103b4c9746fcf27fed3f79b1846339e16784d55e427b7f01456b6829cd3595d |
| occul_data/lindal/lindal_gravity.nc | gravity | 4948bf1b8408406268bbcae38dd0afabb10d8099ce6a6a5fb8d84d17d601a717 |
| occul_data/lindal/lindal_rotation.nc | rotation | eba8612d3aeb494400f10878006dbe3c6a638aad43c05d631ae491084bf4a64d |
| occul_data/lindal/lindal_wind.nc | wind | 7d445c9812e25e984c7af05789aba5ab06fc4084c43157671ce30cb60244d7c2 |
| occul_data/lindal/lindal_composition.nc | composition | 45e4aad47714304a48a59f3042f9fa3d0c67e61c2bc03553daa514a79cde6298 |
| occul_data/lindal/lindal_reduction.toml | manifest | 97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee (unchanged) |
| candidate lindal_refractivity.nc | refractivity | c17a67a1c0f861226a70134ab0e85b9d10eca54aadd1dfec12e5edb3e4a01d5a |
| forward/lindal_closure/inputs/lindal_closure_composition.nc | composition | d76db689b350b647b6a68f715a614c655fc771a63f36671a37175b876f99e2b7 |
| forward/lindal_closure/inputs/lindal_closure_gravity.nc | gravity | 12f0fe891b419c073bbae6c579d855f3f3bd61d8ecda6630265c341e1d90f8c7 |
| forward/lindal_closure/inputs/lindal_closure_rotation.nc | rotation | 55c66e3111e30698cf0192713fa4d4be5a27493c9ada9a77aafc701632d8b601 |
| forward/lindal_closure/inputs/lindal_closure_wind.nc | wind | 11e6e8614548ac2b182488c8425a059b07e2a8b697b869e7a5021034a2910534 |

## 8. Regression

Run by `reports/step03_3/run_regression.sh` on the rebuilt working tree products; results in
`reports/step03_3/regression.txt`, each suite's output beside it. **All pass:** 6, 14, 8, 5, 7, 7,
6, 8, 9, 6 for SPEC_01, and 7 and 9 for SPEC_02 Steps 2 and 3, which load the reduction's `-dirty`
inputs and ran through `relaxed.py`. Steps 02_1, 02_4, 02_5 and 02_6 and the SPEC_03 Step 1 and
Step 2 acceptances run after the sweep (section 5). Rerun after the review's finding 1 was applied
(the reader now parses `epoch`; the gravity and rotation tools read it from the build file; the
Step 1 verification fixtures and the Step 3 suite's control file carry ISO epochs): the same counts,
all pass.

No suite wrote into `occul_data/` or `forward/`.

No em dash or en dash appears in any file written in this step.

## 9. Next step

The review. Nothing is committed and Step 4 does not start until the review is in and the author
says go. After it: the acceptance commit, the clean rebuild of the chain, the product and the run's
inputs, the sweep with hashes recorded in section 7, then the suites of section 5.
