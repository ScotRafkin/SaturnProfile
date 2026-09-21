# REPORT 04, Step 0. Kind W in three parts, kind C as one structure, the transfer run, the namelist, the loader, the anchor object and the propagation hook

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 21 September 2026. Specification: `SPEC_04_Transfer.md` v0.9, Step 0. Status: **reported**.

**Refreshed under REVIEW_04_step0 and SPEC_04 v0.9 (21 September 2026).** The review accepted
this step with changes, and the changes are the author's decisions N and O and the rulings of
SPEC_04 §11. All of them are applied here and the whole step was rebuilt and rerun under them:
kind C is now one structure for every use, a field on `(level, latitude)`; the loader's
coverage, season, point-composition and retrieval refusals are gone, the season recorded and
warned instead; the mesh's extent is no longer declared in the namelist; the `step03_3`
baseline is refreshed. Sections 1 to 4 and 6 are the refreshed text; section 5 is the rerun
regression. The findings of the first filing are kept in section 4 with the ruling on each,
since the report is the record of how the step was worked.

**Fifteen of fifteen acceptance checks pass**, including check 14, which reruns the SPEC_03
Step 3 acceptance and now reports 16 of 16. The one failure of the first filing was that suite's
own check 15 against a baseline written before this step was in prospect; ruling 3 refreshed the
baseline, and the check now means "no value changed since SPEC_04 Step 0" and passes on the
numbers already measured.

**No value changed anywhere in the cascade**, measured twice over: check 6 compares thirteen
products against the pre-step copies group by group and finds only the intended changes of
shape, and the closure product's 22 variables and kind N's 23 are bit-identical to the pre-step
build. The closure product's produced pressure at the gauge level is 9998.465450581718 Pa, the
accepted value.

## 0. Before this step

| Commit | What |
|---|---|
| `545b385` | The SPEC_03 Step 4 sweep record, the last step record. |
| `5f23a8b` | The REVIEW_03_step4 addendum; SPEC_03 closed. |
| `6bc1f9f` | The author's documents: SPEC_04 v0.7 accepted, the 17 September handoff, `STATE.md`. |

The status line of SPEC_04 and `STATE.md` were read before the first filing and again before
this refresh. `reports/REPORT_04_preexecution.md` was filed first, as the status line requires.

## 1. What was built

**Deliverable 1, kind W in three parts.** `lib.schema`: the kind W spec carries
`u_reference_ms(latitude_planetocentric)` in place of `u_cylindrical_ms`, and
`decomposition_geometry` leaves `globals_required`; `check_wind_components` checks
`u_total = u_reference + u_shear` with the reference level wind broadcast along the pressure
axis, and refuses a file whose `u_total_ms` or `u_shear_ms` is not on
`(latitude_planetocentric, pressure)`. `tools/wind/build_wind.py` writes the reference level
column as `u_reference_ms` and `u_total - u_reference` as `u_shear_ms`, and writes neither
`decomposition` nor `decomposition_geometry`.

**Deliverable 2, kind C is one structure (decision O).** `lib.schema` makes
`latitude_planetocentric` a required dimension of kind C and `latitude_planetocentric_deg` a
required variable; the one-latitude form and its `latitude_planetocentric_absent_meaning`
marker are retired from the kind. `tools/composition/build_composition.py` takes
`latitude_grid_deg` as a required key and always writes every level variable on
`(level, latitude_planetocentric)` with `pressure_Pa` one dimensional, recording what it did in
`latitude_grid_rule`. A new `lib/composition.py` holds the one reading rule: `column_at(field,
latitude_deg)` returns the column by the interpolant of decision L, linear between the file's
nodes, refusing a latitude outside the grid, with integer and flag variables taken from the
nearer node rather than averaged. The readers use it: `refrac.reduce.reduce_profile` at the
anchor's own latitude, `forward.production.mean_properties` at the profile's, and the
composition panels of `tools/plots/figures_inputs.py` at the latitude the caller names, which
the panel title then carries. The reduction's and the closure's control files declare the grid
`[-90, 90, 1.0]`, so the change rides the same cascade as deliverable 1.

**The cascade.** `reports/step04_0/rebuild_chain.py` rebuilt, in order, the raw bundle, the six
reduction inputs, kinds T and D with the manifest, the kind N product, the closure run's four
inputs, the closure product, and the transfer run's four inputs. Hashes in section 6.

**Deliverable 3, the transfer run.** `forward/lindal_transfer/` with `lindal_transfer.toml` as
the specification prints it at v0.9 (the two spacings, no declared extent),
`lindal_transfer_build.toml` carrying the four sections under the run prefix with the
composition on `[-90, 90, 1.0]`, `inputs/` written by `casspian-run-inputs`, and `output/`.

**Deliverable 4, the namelist in transfer mode.** `lib.control`: `transfer` joins the
implemented modes; the vocabulary is chosen by the declared mode, so closure mode refuses
`[grid]`, `[numerics]`, `[estimation]`, `datum_isobar_Pa` and `[target]` exactly as before.
`[grid]` now carries `geopotential_spacing_m2s2` and `latitude_spacing_deg` only:
`geopotential_range_m2s2` is gone with no replacement, and with it the range check and the
`_LIST` type it needed, since Step 2 builds the extent from the anchors and grows it when a
traced curve needs more. `_check_numerics` checks the six `[numerics]` tables, refusing an
unknown table, a missing one, a scheme not implemented, and `[numerics.geopotential]` or
`[numerics.hydrostatic]` as rules SPEC_03 closed and not namelist keys. `weight` is the anchor's
role and only 1 or 0 is accepted.

**Deliverable 5, `load_run_inputs` in transfer mode.** `_load_transfer_inputs` loads every anchor
and the four inputs with the SPEC_03 checks and no closure comparison. Under decision N it makes
no coverage check (that is the interpolants' business, and `lib.composition.column_at` already
refuses a latitude outside its grid), no point-composition check (the form is retired) and no
retrieval-instance check (the schema refuses that file first). The seasons of W, C and every
anchor are recorded beside the run's in `RunInputs.seasons` and a difference is warned, never
refused.

**Deliverable 6, the anchor as it arrives.** `lib.control.LoadedAnchor`, one per `[[anchors]]`
entry, carrying on the anchor's own levels its `ln N`, `radius_m`, `label_pressure_Pa`,
`sigma_ln_N_measurement` (`scale * sigma_N / N`, an uncertainty of `ln N`, as ruling 1 fixed),
`sigma_ln_N_season` with `season_term`, its latitude, its gauge level index, its weight and
scale, its own season and the run's, and the hook's record.

**Deliverable 7, the propagation hook.** `forward/propagate.py`: `propagate(anchor, season)`
returns the anchor with `N`, its levels and its season term untouched, and a record carrying the
anchor's season, the run's, the seasonal distance the shorter way round the seasonal circle, and
`propagation = "none: propagator not implemented"`.

**One refactor, to make the relaxation possible.** `refrac` and `forward` each refused a `-dirty`
input at their own site with their own message. This step rebuilds input files on a working tree,
so the two were merged into `_refuse_dirty_commit(path, attrs, consumer, what=None)`, which
reproduces both message forms byte for byte (check 0 quotes both), and `_refuse_dirty` delegates
to it. The acceptance relaxes that one function and nothing else.

## 2. Decisions

Decisions 1, 3, 4 and 6 to 13 of the first filing were accepted as reported; 2 and 5 were
superseded. The list below is the refreshed one.

1. **`sigma_ln_N_measurement` is the uncertainty of `ln N`**, `measurement_uncertainty_scale *
   refractivity_uncertainty / refractivity`. Finding 1, ruled correct; deliverable 6 now says so.
2. **Superseded by ruling 2.** The `[grid]` range check is removed, not relocated:
   `geopotential_range_m2s2` leaves the namelist and `lib.control` entirely.
3. **`model_error_correlation_length_deg` other than 0.0 is refused**, naming Eq. A35 and the
   combination specification. Kept under decision N: a recorded-and-ignored value would be the
   silent kind.
4. **`measurement_uncertainty_scale` must be positive.** It multiplies an uncertainty.
5. **Superseded by decision O.** The point-composition rule and its planetographic pairing go
   with the form they policed.
6. **`latitude_grid_deg` requires the step to divide the interval exactly**, so both ends are
   nodes and the grid a reader sees is the grid the control file declares.
7. **The gridded composition records `latitude_grid_rule`**, prose the model does not read
   saying that the source's one column was written at every node.
8. **`decomposition_geometry` stays in `CLOSURE_DROPPED_ATTRIBUTES`.** The attribute is retired,
   so the entry is inert; removing it would change SPEC_03 closed code to no effect.
9. **The swept copies were rebuilt in a git worktree at `5f23a8b`.** The pre-step products were
   overwritten by the cascade before they were copied aside. Recorded so the sweep's own copies
   replace them.
10. **`reports/step03_3/candidate/` was rebuilt** from the new-form inputs with its own
    `build_candidate.py`, twice: once for kind W and again for kind C. No acceptance script was
    edited to make it pass.
11. **`reports/step02_1/accept_step02_1.py` reads one more report**, `REPORT_04_step0.md`, so
    that section 6 below is read. Section 7 of the handoff allows this.
12. **`reports/step7/accept_step7.py` checks 2 and 6 were moved to the three-part form.**
    Finding 5, accepted by the review.
13. **The regression is structured on `step03_3`'s, not `step03_4`'s.** Section 5.
14. **`lib/composition.py` is a new module.** Decision O says every reader takes kind C as a
    field, so the one reading rule lives in one place rather than in each reader. The
    specification names `lib.windfield`, `lib.mesh` and `lib.kernel` for later steps and does
    not name this one; it is the smallest home for a rule three readers share.
15. **Flags are not interpolated in latitude.** `column_at` interpolates floating point
    variables and takes an integer or flag variable, such as `nh3_provenance`, from the nearer
    node: the average of two provenance codes means nothing. For the uniform field in hand every
    node carries the same code, so this decides nothing today and is recorded for the files that
    will vary.
16. **A composition panel draws one column and names it.** `panel_mole_fractions` and
    `panel_mean_properties` take a latitude; the product figures pass the anchor's, and a kind C
    file offered on its own has no latitude in its context, so the panel draws the node nearest
    the equator and the panel title says which latitude is drawn.
17. **`reports/step03_3/before/` was refreshed** at this step's products, as ruling 3 directs.
18. **`reports/step8/accept_step8.py` and `reports/step02_3/accept_step02_3.py` were moved
    to the field form of kind C.** Finding 6, the same class of change as decision 12.

## 3. Acceptance results

Script `reports/step04_0/accept_step04_0.py`, output `reports/step04_0/output.txt`. The products
were rebuilt on this working tree and carry `-dirty`, so the one refusal point is relaxed for the
script's duration; check 0 shows the refusal firing first and quotes it.

| Check | Result |
|---|---|
| 0. The `-dirty` refusal fires, then is relaxed, named | **Pass.** Refused: `lindal_refractivity.nc (anchor lindal) carries casspian_git_commit = '6bc1f9f...-dirty'.` |
| 1. Kind W in three parts, the retired items gone, all three files | **Pass.** `u_reference_ms(latitude_planetocentric)` shape (361,); no `u_cylindrical_ms`, no `decomposition`, no `decomposition_geometry`. |
| 2. `u_total` unchanged to the bit against the pre-step copies | **Pass.** `u_total_ms`, `u_total_uncertainty_ms`, `value_provenance` and both coordinates identical in the reduction and closure files; `u_reference_ms` equals the pre-step `u_cylindrical_ms` at column 50, the 1.0e5 Pa reference level; the transfer wind equals the closure wind bit for bit. |
| 3. The sum identity, and the shear at the reference level | **Pass.** Worst `\|u_total - (u_reference + u_shear)\|` 0.000e+00 m/s in all three; max `\|u_shear\|` at the reference level 0.000e+00 m/s. |
| 4. A copy with the identity broken is refused | **Pass.** `u_total does not equal u_reference + u_shear; worst departure 1.000e-03 m/s against a field scale of 4.905e+02 m/s`. |
| 5. Kind C is one structure carrying the source's column at every node | **Pass.** The transfer composition on 181 nodes, -90 to 90, one unique step of 1.0; `x_H2`, `x_He`, `x_NH3`, `x_H2_uncertainty` and `nh3_provenance` equal across columns and equal to the pre-step point file bit for bit; `pressure_Pa` still one dimensional and unchanged; no absence marker. The reduction's and the closure's compositions are fields on `{'level': 66, 'latitude_planetocentric': 181}`, and the column each reader takes, at `phi_c` = 30.805568 deg, equals the pre-step point file bit for bit in every species variable. A latitude outside the grid is refused by the interpolant. |
| 6. The cascade changed no value, group by group | **Pass.** Thirteen products compared. The only differences anywhere: `u_reference_ms` added and `u_cylindrical_ms` removed with `decomposition` and `decomposition_geometry` gone, in the three wind files and in the wind group embedded in kind N and in the closure product; `latitude_planetocentric_deg` added and `latitude_planetocentric_absent_meaning` gone, in the three compositions and in the composition groups embedded in kind N and in the closure product. Every variable present in both carries the same values, a variable that gained the latitude axis checked column by column. |
| 7. `casspian-run-inputs` wrote the four transfer inputs | **Pass.** Each validates under its kind with `role = 'forward'` and `profile_or_run = 'lindal_transfer'`. |
| 8. The transfer namelist loads and resolves | **Pass.** mode `transfer`, target 10.0, gauge 1.0e4 Pa, datum 1.0e5 Pa, grid spacing 5.0e3 and 0.05 deg with no declared extent, six numerics tables, estimation as declared. |
| 9. The anchor as it arrives, with its two uncertainty columns | **Pass.** `lindal` at 30.805568 deg, 66 levels, gauge level 29 at 10000 Pa, weight 1.0; `sigma_ln_N_measurement` 2.331845e-02 at every level, equal to `scale * sigma_N / N`; `sigma_ln_N_season` all zero, `season_term = 'absent'`. |
| 10. The hook, the identity, and its record | **Pass.** `ln N` unchanged bit for bit, season term still zero, `propagation = 'none: propagator not implemented'`; the same anchor moved to season 200.0 gives `seasonal_distance_deg` -178.2 and `season_matches_run` False, with `N` untouched. |
| 11. Four refusal cases, each message quoted (decision N: these and no others) | **Pass.** A missing `[target]`; an unknown `[numerics]` table; a scheme not implemented; a `weight` of 0.5. |
| 12. A season that differs is recorded and warned, never refused | **Pass.** A wind and a composition edited to `solar_longitude_deg = 200.0` each load, each raise the warning, and each are recorded as 200.0 beside the run's 18.2. The unedited run records `{'run': 18.2, 'wind': 18.2, 'composition': 'uniform', 'anchor:lindal': 18.2}`. |
| 13. M > 1, the validation anchor, and the edited composition | **Pass.** Two `[[anchors]]` entries pointing at one file load as `[('lindal', 1.0, True, 1.0), ('lindal', 0.0, False, 2.0)]`; a scale of 2.0 doubles `sigma_ln_N_measurement`; a composition with the node at 45 deg edited loads, refused nowhere. |
| 14. The closure namelist and the SPEC_03 Step 3 acceptance unchanged | **Pass.** The closure namelist loads with all four closure comparisons identical and every transfer-only field empty; `reports/step03_3/accept_step03_3.py` reports **16 of 16** (exit 0), its check 12 being the thirteen closure refusals. |

## 4. Findings, and the ruling on each

Findings 1 to 5 are the first filing's, each with the ruling that settled it. Finding 6
is new, found by the rerun regression.

**1. Deliverable 6 defined `sigma_ln_N_measurement` as a quantity that is not an uncertainty of
`ln N`, and the error was invisible at M = 1.** Measured on the rebuilt kind N:
`refractivity_uncertainty` runs from 1.1750e-9 at the top level to 7.2541e-6 at the bottom, a
factor of 6.2e3 ordered by level, while `sigma_N / N` is 2.331845e-2 at every one of the 66
levels. Taken literally the weights of Eq. A30 would fall by 3.8e7 from bottom to top, and at
M = 1 the estimate is the identity so nothing would show until a second anchor arrived.
**Ruled (§11.1):** correct, the specification's parenthetical was dimensionally wrong;
deliverable 6 now reads `sigma_N / N` times the scale. Implemented as ruled.

**2. The acceptance's `[grid]` refusal could not be made at Step 0**, since an anchor's `Phi_k`
is the field line integral Step 1 forms. **Ruled (§11.2):** the check is removed rather than
moved. The extent is not something a user declares and the code polices; Step 2 builds it from
the anchors' levels and extends it when a traced curve needs more. `geopotential_range_m2s2`,
its checks and the `_LIST` type are gone from the namelist and `lib.control`. Decision 2 of the
first filing is superseded.

**3. The SPEC_03 Step 3 acceptance failed its own check 15, on the change this step makes.**
It compared against `reports/step03_3/before/`, copies taken before SPEC_03 Step 3, and reported
for `lindal_wind.nc` alone `values differing ['/:u_cylindrical_ms']` with `decomposition` and
`decomposition_geometry` gone. **Ruled (§11.3):** refresh the baseline at the Step 0 acceptance
commit, so the check means "no value changed since SPEC_04 Step 0". Done (decision 17); the suite
is now 16 of 16 and check 14 passes on the numbers already measured. Nothing was loosened: check
6 of this step is the same comparison against the pre-step copies and passes.

**4. The retrieval-instance refusal was unreachable through `lib.io.read`**, which refuses a
missing `radius_m` as a schema fault first. **Ruled (§11.4):** removed from the loader as dead
code by the schema's prior claim; the retrieval leg adds its own check when the instance exists.
Done.

**5. The SPEC_01 Step 7 acceptance read two items the Appendix retires**, raising
`KeyError: 'u_cylindrical_ms'`. The two checks were moved to the three-part form: check 2 now
sets the whole polar column of `u_total_ms` and the one polar value of `u_reference_ms`, since
the reference level wind is one value per latitude and moving a single column would break the
identity elsewhere; check 6 forms `u_total - (u_reference + u_shear)`. **Ruled (§11.5):**
accepted as the Appendix's amendment applied to a closed step's acceptance record. `step7` is
8 of 8 with its sum identity departure and maximum shear both zero, as before.

**6. Two more accepted acceptances read kind C in the retired one-latitude form.** The regression
found them, which is what it is for. `reports/step8/accept_step8.py` raised
`IndexError: index 9955 is out of bounds for axis 0 with size 66` at its ammonia clamp check,
indexing a (66, 181) field with indices meant for 66 levels;
`reports/step02_3/accept_step02_3.py` raised `ValueError: The truth value of an array with more
than one element is ambiguous` comparing `x_NH3[k2] == 0.0`. Both were moved to the field form
the way finding 5 moved `step7`: each takes the column by `lib.composition.column_at`, `step8`
at the file's own first node and `step02_3` at the anchor's latitude, which is the column the
reduction reads. Each first asserts that no variable varies with latitude, so that taking one
column of the field hides nothing and the suite would fail loudly rather than quietly if a
later composition did vary. Measured after the change: `step8` is 9 of 9 and `step02_3` is 9 of
9, both their reference counts. This is the Appendix's "SPEC_01 Step 8 amended the same way"
applied to two accepted steps' acceptance records, the same class of change the review accepted
as finding 5, and it is recorded here for the same reason: it edits closed steps' scripts.

Nothing else arose from the refresh. The work of decisions N and O went in without a surprise:
the one thing worth recording is that the column the reduction now reads is bit-identical to the
pre-step point file, which is what makes "no value changes" true rather than approximately true.

## 5. Regression

`reports/step04_0/run_regression.sh`, output `reports/step04_0/regression.txt`. The template is
`reports/step03_3/run_regression.sh`, not `step03_4`'s, because this step rebuilds input files
and `step03_3` is the last step that did. The suites fall in three groups, as section 6 item 4
of the handoff of 17 September and SPEC_03 section 0 direct: those that read the files run
directly; `step02_2` and `step02_3`, which run `refrac` on them, run through
`reports/step03_3/relaxed.py`, which relaxes the `-dirty` refusal in memory and marks each
relaxed commit visibly; and those that run `refrac` or `forward` on the registered products, or
read the registered kind N, run after the sweep.

| Suite | Result | Reference |
|---|---|---|
| `step1/accept_step1` | 6 of 6 | 6 |
| `step1/verify_review_changes` | 14 of 14 | 14 |
| `step2/accept_step2` | 8 of 8 | 8 |
| `step3/accept_step3` | 5 of 5 | 5 |
| `step4/accept_step4` | 7 of 7 | 7 |
| `step5/accept_step5` | 7 of 7 | 7 |
| `step6/accept_step6` | 6 of 6 | 6 |
| `step7/accept_step7` | 8 of 8, after the change of finding 5 | 8 |
| `step8/accept_step8` | 9 of 9, after the change of finding 6 | 9 |
| `step9/accept_step9` | 6 of 6 | 6 |
| `step03_3/accept_step03_3` | 16 of 16, after the baseline refresh of ruling 3 | not in the reference list |
| `step04_0/accept_step04_0` | 15 of 15 | this step |
| `step02_2/accept_step02_2` (relaxed) | 7 of 7 | 7 |
| `step02_3/accept_step02_3` (relaxed) | 9 of 9, after the change of finding 6 | 9 |

Every SPEC_01 suite returns its reference count, 6, 14, 8, 5, 7, 7, 6, 8, 9, 6, and the two
relaxed SPEC_02 suites return 7 and 9, theirs.

Deferred to the sweep, each needing clean on-disk products: `step02_1` (whose recorded hashes
the sweep replaces), `step02_4`, `step02_5`, `step02_6`, `step03_1`, `step03_2`, `step03_4`.
They are rerun on the swept products and recorded in the sweep record, as REPORT_03_step3 did.

Restored after the run, measured: the registered kind N SHA-256 is
`c3662500e5db36f0bc887821447911aefbbb9988bb11de09b54202690ede7813` before and after, equal; the
closure run's inputs and product restored equal; the transfer run's inputs restored equal.
Porcelain afterwards holds only this step's own source and documents and the author's modified
design note.

## 6. Hashes, provisional until the sweep

Every product below was built on this working tree and carries `casspian_git_commit`
`6bc1f9ffdd96799672b7c8e9922058f02456449b-dirty`. They are rebuilt on the clean tree at the sweep
and this table is replaced by the sweep record then.

| `file` | kind | `64 hex digits` |
|---|---|---|
| `occul_data/lindal/raw/lindal_raw.nc` | raw | `b4a093c8aa4154698afa2d79f5bb60fbbf9124a27d54d590c82d1554286f49ce` |
| `occul_data/lindal/lindal_thermo.nc` | thermo | `85e9c3ca0228e815694928eb2dda8033fd83f59d91587e5479a5050bda50151e` |
| `occul_data/lindal/lindal_geodesy.nc` | geodesy | `3ebd61ce4a1148469f1e378fc02a2d912839dc43cd7d3a3c2cd63eed83db0bc1` |
| `occul_data/lindal/lindal_gravity.nc` | gravity | `8934c1425f8844b86e197df5a7274725e886675ae5f2d4d3f66dce68561ba03b` |
| `occul_data/lindal/lindal_rotation.nc` | rotation | `0f7a643d6ae989fb7bf91eddc07e5ddf807b60eb7ef8e2a52183f26c4a6ed0a7` |
| `occul_data/lindal/lindal_wind.nc` | wind | `dd8d2b15c6c5a0576050ec69f922e5ef97950d19b8d422f23e63019898bf83ed` |
| `occul_data/lindal/lindal_composition.nc` | composition | `35b30ff89faa33f606f87a84df7e4bd8cad9e0400282aa8ee21792307d3d5e9d` |
| `occul_data/lindal/lindal_reduction.toml` | manifest | `97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee` |
| `occul_data/lindal/lindal_refractivity.nc` | refractivity | `c3662500e5db36f0bc887821447911aefbbb9988bb11de09b54202690ede7813` |
| `forward/lindal_closure/inputs/lindal_closure_gravity.nc` | gravity | `52f55b782724d0036e006b6acdb447d2316b2445fe14e1d90d9ada4ea64c9678` |
| `forward/lindal_closure/inputs/lindal_closure_rotation.nc` | rotation | `8a876cf4144115d94e208696baddaf6196aece78b0a7cd3e88c70b4c3d35b89b` |
| `forward/lindal_closure/inputs/lindal_closure_wind.nc` | wind | `6fa3d9ce908a43f5fc12374aa10e91bf128763e4094c3f2e87c4722e00fe40ec` |
| `forward/lindal_closure/inputs/lindal_closure_composition.nc` | composition | `48b207375bebbe7431d572c4102c345d910731178499df58077b4e7b2f66bda9` |
| `forward/lindal_closure/output/lindal_closure_profile.nc` | profile | `7e71ace145c4c78a0ddaae90fa632f81e92d9fbb3103ebdc364adc91051f5c03` |
| `forward/lindal_transfer/inputs/lindal_transfer_gravity.nc` | gravity | `f2f4a641564ad57464f568a646e1af144962b39e601f39e65570163a61ddd3a5` |
| `forward/lindal_transfer/inputs/lindal_transfer_rotation.nc` | rotation | `9cf9f57f47c7358ddae176719a5f506969a008e9167e6a680daa62339dc49cef` |
| `forward/lindal_transfer/inputs/lindal_transfer_wind.nc` | wind | `d8f46ab1939b09141c2760c4a39f7d8a57c58bee0afff9f61e8e02c198c2f0f3` |
| `forward/lindal_transfer/inputs/lindal_transfer_composition.nc` | composition | `7e9157f9e8de5714890f395294c4cdd7122c3b9d36a56ab66c3e7666bd45085a` |

REPORT_03_step4's product hash `e6693173...` is superseded by the closure product row above, and
that row is superseded in turn at the sweep.

## 7. Next step

The review's order of work: commit the author's documents, then the acceptance commit for Step 0,
then the sweep with the deferred suites and the swept hashes, then `STATE.md` to accepted. Step 1
follows: `lib.geoid.through_anchor`, `lib.windfield`, and the anchors' geopotential under the
run's wind with `produce` taking the wind along the column.
