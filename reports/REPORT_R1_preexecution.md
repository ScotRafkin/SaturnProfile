# REPORT R1, pre-execution: the classification under the author's rule

SPEC_R1 v0.3, Step 1, before any code changes. 22 October 2026. For the author and the reviewing
agent. `main` is at `8adc7b5`, pushed; the only changes in the tree are the specification files.

The rule (SPEC_R1 §0): a refusal, test check or required key stays only if without it the code
would compute a wrong number without saying so, or fail with no clear message. Every refusal in
`lib.control`, `lib.schema`, `lib.io`, the tools' control readers and the suites was read, and every
attribute that a writer writes was traced to its readers. The classification below is what Step 1
would implement. Two questions in §4 need a ruling first, because the answer changes how much of
every product is rewritten.

## 1. Removed

### 1.1 Named by SPEC_R1 §1 item 1

| Item | Where | Reason |
|---|---|---|
| The closure content comparison: `check_closure_inputs`, `CLOSURE_DROPPED_ATTRIBUTES`, `ClosureComparison`, `_strip`, `_dropped_names`, `_differences`, `_same_value`, `RunInputs.closure`; the `closure_comparison_*` and `closure_dropped_*` keys of `production_record` | `lib/control.py`, `forward/production.py` | Named. No number depends on it. |
| The closure season equality | `control.py` line 1374 | Named. |
| `profile_or_run` against the `[[anchors]]` slug, closure and transfer | `control.py` lines 1309, 1509 | Named. |
| `profile_or_run` against the run name, `role = "forward"`, `composition_role = "forward"` on run inputs, closure and transfer | `control.py` lines 1327 to 1341, 1528 to 1542 | Named. |
| The run prefix on input files and on the product | `control.py` `run_input`, line 1047 | Named. |
| The `p_b_location` value refusal | `control.py` line 957 | Named. One value, read by nothing but the record. |
| `step9` check 6, the directory listing | `tests/step9/accept_step9.py` | Named. |
| The `-dirty` refusal, `_refuse_dirty_commit` and `_refuse_dirty` | `control.py` lines 1253 to 1274 | Named. `casspian_git_commit` is still written, `-dirty` and all, by `lib.io.git_commit`. |

### 1.2 By the rule

| Item | Where | Reason |
|---|---|---|
| The manifest's slug prefix on inputs and product | `control.py` `inside()` | A name. |
| The namelist file named for its run (`path.stem != name`) | `control.py` line 917 | A name. |
| `casspian-run-inputs`: `role = "forward"`, `prefix` equal to the run name, `output` carrying the run prefix | `tools/run/run_inputs.py` lines 77 to 91 | Labels and names. |
| `build_role` value check | `control.py` line 123 | A label that no code reads after SPEC_R1 (§4, question 2). |
| The prefix check on manifest file names | `tools/lindal/build_inputs.py` line 380 | A name. |
| "Mode not implemented" | `control.py` line 843, `production.py` line 386 | Dead code: every mode is implemented. |
| Bare names `latitude`, `r`, `z`, `p` refused | `schema.check_names` | Nothing reads a bare name, so nothing reads one wrongly. |
| `role` value, kind N `role = "reduction"`, profile `role = "forward"` | `schema._check_globals` | Labels. |
| `thermo_instance` value, the 13 source profile globals required, and forbidden on other instances; a stated `latitude_*_deg` global required | `schema._check_globals` | Descriptive. The one that is read, `latitude_planetographic_deg`, fails with a clear `KeyError` naming it in `refrac`. **Not placed instead** if the author wants the clear message kept (§3 item 7). |
| `pressure_grid_rule` required beside `pressure_printed_Pa` | `schema._check_globals` | The figure reads it with a default and falls back to no grid term. |
| `composition_role` value and `source_statement` for role reduction | `schema._check_globals` | Labels and prose. |
| `vertical_coordinate` value, and the named variable present | `schema._check_globals`, `_check_variables` | Read by nothing. |
| `epoch` an ISO date | `schema._check_season` | Read by nothing. `epoch` on kind T is copied into kind N and read by nothing there. |
| `solar_longitude_source` required | `schema._check_season` | Prose. |
| `units`, `long_name`, `provenance` required on every variable; `provenance` in the vocabulary | `schema._check_variable_attributes` | Read by nothing but the schema. `provenance` is copied into kind N and read by nothing there. |
| `uncertainty_kind` required and in the vocabulary on every companion | `schema._check_variable_attributes` | `refrac` refuses a finite companion without a convertible kind itself (`reduce.py` line 132, `reduction.py` line 63); those stay. Elsewhere nothing reads it. |
| `positive` and `direction` required on vertical coordinates; profile `positive = "up"` | `schema._check_coordinates`, `_check_profile` | Read by nothing. Monotonicity stays (§2). |
| Profile companions present and NaN | `schema._check_profile` | A policy, read by nothing. |
| `<dimension>_absent_meaning` markers and their values, and the stated latitude for `point` | `schema._check_dimensions` | Read by nothing. Required dimensions stay (§2). |
| Required globals read by nothing: `title`, `source`, `epoch`, `role` on every kind; `profile_or_run` on C, D, G, R, T, W, raw; D `citation`, `fit_inputs`, `fit_latitude_convention`; G `GM_scope`; R `citation`; W `observation_level_Pa`, `source_latitude_convention`; profile `boundary_pressure_Pa`, `anchor_solar_longitudes_deg`, `input_hashes` | `schema` `AUTHOR_REQUIRED_GLOBALS`, `globals_required` | Read by nothing. The requirement goes; whether the writers stop writing them is §4, question 1. |
| Writer globals required on read (`Conventions`, `created_by`, `created_at`, `casspian_git_commit`, `history`) | `io.read` line 375, `schema` `writer_filled` | Read by nothing that fails without them. `casspian_kind` and `casspian_schema_version` stay. |
| Flag variables int8 with `flag_values` | `io._check_types` | Read by nothing. `flag_meanings` stays where the figures read it (§2). |
| `observation_level_justification` still written | `tools/lindal/lindal_wind.py` line 317 | SPEC_00 v0.24 dropped it; this writer was missed. |

Suite checks that only exercised a removed item are retired; multi-case checks lose the cases.
The ones found: `step1` checks 5 and some cases of `verify_review_changes`; `step02_1` 2 and 5 (in
part); `step02_4` 6; `step03_3` 2, 5, 10, 12 (in part); `step03_4` 10, 12 (in part); `step04_0` 0,
11 (in part), 14 (in part); `step05_3` 3; `step9` 6. Each is read in full at implementation and the
counts go into `tests/run_regression.sh` with old and new values. The relaxation helpers
(`tests/step03_3/relaxed.py`, `step03_0/relaxed.py`, and the in-script replacements in `step04_0`,
`step04_4`) go with the refusal.

## 2. Kept

| Item | What it prevents |
|---|---|
| Unknown keys and sections refused, in every control reader and `load_section` | A typo silently dropping a value, the run proceeding on a default. |
| Missing required keys, wrong types, non-positive values where a positive is used | An obscure failure deep in the arithmetic, or a wrong number. |
| Manifest anchor rule and quantity belonging together (`control.py` line 373, `build_inputs.py` line 131) | The geoid anchored on the wrong radius: a wrong number. |
| The reduction's cross-checks: rotation system and rate of W against R; reference level a node; anchor isobar a geodesy surface and a thermo level; geoid anchor surface a geodesy surface; anchor quantity present; composition levels equal to thermo levels | Each a wrong number or an interpolation the reduction does not do. |
| `_admit_wind`: `u_total_ms` present, both poles nodes with zero wind | An obscure failure, and a singular meridional gravity at the pole. |
| Closure: exactly one anchor; `[target]` refused; `p_b_rule` required; composition levels equal to the anchor's; gauge isobar a tabulated level; wind coverage spans the anchor | A value silently ignored (wrong number), or an extrapolation. |
| Transfer and namelist: weight 0 or 1; `measurement_uncertainty_scale` positive; exactly one of `p_b_rule` and `p_b_Pa`; known `p_b_rule`; schemes and gauge rule implemented; correlation length 0; target latitude a latitude; solar longitude in [0, 360); date ISO | A choice silently replaced by another, or a wrong number. |
| `harmonic_convention` a known code | Coefficients read under the wrong convention: a wrong number. |
| Required variables and dimensions of every kind; uncertainty companions present where `refrac` reads them | An obscure `KeyError` in the arithmetic. |
| Coordinates strictly monotonic; profile geopotential monotonic | Wrong interpolation. |
| Mole fractions sum to one | A wrong refractivity. |
| float64 for every float variable | Silent loss of precision. |
| `casspian_kind` present and as asked; schema version not newer | A file read as the wrong kind. |
| Required groups of N, C, profile; profile carries `anchors/<slug>` and one record group | An obscure `KeyError`. |
| `profile_or_run` on kinds N and profile | `tools/plots/render.py` line 77 reads it to name figures and fails with a bare `KeyError`. |
| `vertical_structure`, `flag_meanings`, `system_name`, `rotation_system_name`, `rotation_rate_rad_s`, `coverage_*`, `closure_rule`, `closure_species` | Read by the code or the figures. |
| `recorded_path` across drives; attribute types netCDF can store | An obscure failure in `os.path.relpath` or the netCDF library. |
| The numerical refusals of `lib` and the tools (extrapolation, crossing isobars, non-convergence, shapes, PDF tracing) | Named numerical failures and extrapolation, as the 10-03 handoff §6 keeps them. |

## 3. Not placed, kept until ruled

1. **Path containment:** manifest inputs inside the manifest's directory; run inputs in `inputs/`;
   output directory inside the run directory. No wrong number, but they keep a run from writing
   over registered data (SPEC_00 §2.2, "nothing under `forward/` may write here").
2. **The season present on every file** (exactly one of `solar_longitude_deg` and
   `season_absent_meaning`) and its range. The transfer loader records an absent season as
   `unstated`, so no number is wrong; but SPEC_00 v0.17 made it required so that a missing season
   is never silent.
3. **`package_version`** refuses when the installed and literal versions disagree. Only the
   recorded version would be wrong.
4. **Kind N's hash cross-check** (`io._check_refractivity_hashes`): detects an edited
   `input_hashes`, not edited values.
5. **`casspian-new-run` over an existing directory:** `mkdir` would fail anyway, with a less clear
   message.
6. **`reject_physical_values`:** a number in a section whose keys are paths and names. Unknown keys
   are refused anyway; this catches a number given for a known key.
7. **The source profile's stated latitude global** (§1.2): if the clear message is wanted, it moves
   here.

## 4. For the author, before implementation

1. **What "writers stop writing keys and attributes that no code reads" covers.** Read literally it
   removes almost every record the products carry: `production_record` and `transfer_record` in
   full (the pressure identity and the outer loop passes, which the reports quote), most of
   `reduction_record` (the fixed point iterates, the anchor rule spread, the propagation terms),
   every citation and `source`, `title`, `epoch`, the uncertainty notes, and the embedded
   namelist. None of it is read by code; all of it is read by people and by the suites.
   Proposed: writers stop writing only what existed to satisfy a removed check: `role`,
   `composition_role`, `profile_or_run` on the input kinds, `thermo_instance` and the source profile
   set it governed, `p_b_location`, and the closure comparison keys. Records, citations and
   provenance stay.
2. **The control keys that go with them.** `role` in every build-file section and `p_b_location`
   in every namelist would become keys no code reads, which the unknown-key rule then refuses.
   Proposed: remove them from the readers' tables and from every committed build file and namelist
   in the same step (`occul_data/lindal/lindal_build.toml`, `forward/*/`), and from
   `casspian-new-run`'s templates if any. The uncommitted `forward/lindal_iris_*` runs are rebuilt by
   their suites from build files that would be edited too.

The full regression, in chunks, follows implementation as SPEC_R1's acceptance asks.
