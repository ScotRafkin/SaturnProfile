# REPORT R1, Step 1. Before acceptance

Coding agent, 22 October 2026. Specification: `docs/specs/SPEC_R1_Regression_Fixes.md` v0.4, with
§2's rulings on `reports/REPORT_R1_preexecution.md`. Working tree on `main` at `8adc7b5`,
uncommitted.

**In short.**
- Every refusal, required key and suite check was judged by the author's rule. The removed, kept
  and not placed lists are §1 to §3, one line each.
- The closure content comparison, the closure season equality, the name and label refusals, the
  `-dirty` refusal and `step9` check 6 are gone, as named. About twenty-five more requirements went
  by the rule, most of them in `lib.schema`.
- Writers stop writing `role`, `composition_role`, `profile_or_run` on the input kinds,
  `thermo_instance` and the source profile set, `p_b_location`, the closure comparison keys and
  `observation_level_justification` (§2 ruling 1). `role` and `p_b_location` left every reader
  and every committed control file (ruling 2).
- One change beyond the lists was needed and is reported in §4: `casspian-lindal-inputs` compares
  only the globals it still writes when it decides to keep kinds T and D.
- Thirteen suites changed count; seven were added to the driver. The full regression is in §6.
- No computed value changed (§6).

## 1. Removed

### 1.1 Named by SPEC_R1 §1 item 1

| Item | Where | Reason |
|---|---|---|
| The closure content comparison: `check_closure_inputs`, `CLOSURE_DROPPED_ATTRIBUTES`, `ClosureComparison`, `_nodes`, `_strip`, `_dropped_names`, `_differences`, `_same_value`, `RunInputs.closure`; the `closure_comparison_rule`, `closure_comparison_<kind>` and `closure_dropped_<kind>` keys of `production_record` | `lib/control.py`, `forward/production.py` | Named. Decides only whether a closure is refused. |
| The closure season equality | `lib/control.py` `load_run_inputs` | Named. |
| `profile_or_run` against the `[[anchors]]` slug, closure and transfer | `lib/control.py` | Named. |
| `profile_or_run` against the run name, `role = "forward"`, `composition_role = "forward"` on run inputs, closure and transfer | `lib/control.py` | Named. |
| The run prefix on input files and on the product; run inputs required under `inputs/` | `lib/control.py` `read_run_namelist` | Named; §2 ruling 3.1. |
| The `p_b_location` value refusal, and the key | `lib/control.py`, `forward/production.py` | Named; §2 ruling 2. |
| `step9` check 6, the directory listing | `tests/step9/accept_step9.py` | Named. |
| The `-dirty` refusal, `_refuse_dirty_commit`, `_refuse_dirty` | `lib/control.py` | Named. `casspian_git_commit` is still written with its `-dirty` mark by `lib.io.git_commit`. |

### 1.2 By the rule

| Item | Where | Reason |
|---|---|---|
| The manifest's inputs inside its own directory, and the slug prefix on inputs and product | `lib/control.py` `read_reduction_manifest` | A location and a name (§2 ruling 3.1: data may come from anywhere). |
| The namelist file named for its run | `lib/control.py` | A name. |
| "Mode not implemented" | `lib/control.py` | Dead: every mode is implemented. `forward/production.py`'s guard against a transfer namelist stays, its message corrected to name the transfer driver. |
| `build_role`, `BUILD_ROLES` and the `role` key of every build section | `lib/control.py`, the five build tools | A label no code read. |
| `build_epoch`'s ISO check | `lib/control.py`, `tools/gravity/` | `epoch` is read by nothing; it is recorded as given. |
| `reject_physical_values` | `lib/control.py`, `tools/gravity/` | §2 ruling 3.6. |
| `casspian-run-inputs`: the file named `<run>_build.toml`, `role = "forward"`, `prefix` equal to the run, outputs under `inputs/` with the prefix | `tools/run/run_inputs.py` | Names and labels. The output must still lie inside the run directory (§2). |
| `casspian-new-run`'s own existing-directory refusal | `tools/run/new_run.py` | §2 ruling 3.5. `target.mkdir()` has no `exist_ok`, so it raises `FileExistsError` on an existing run; a comment says so. |
| The prefix on the files a manifest names | `tools/lindal/build_inputs.py` | A name. |
| Bare names `latitude`, `r`, `z`, `p` | `lib/schema.check_names` | Nothing reads a bare name. |
| `role` values; kind N role reduction; profile role forward | `lib/schema.py` | Labels. |
| `thermo_instance` and its values; the source profile globals required, and forbidden on other instances; a stated `latitude_*_deg` global | `lib/schema.py` | Descriptive. `refrac` reads `latitude_planetographic_deg` and its `KeyError` names it (§2 ruling 3.7). |
| `pressure_grid_rule` beside `pressure_printed_Pa` | `lib/schema.py` | The figure reads it with a default. |
| `composition_role` and `source_statement` | `lib/schema.py` | Labels and prose. |
| `vertical_coordinate` value and the named variable present | `lib/schema.py` | Read by nothing. |
| The season identifier: `epoch` an ISO date, exactly one of `solar_longitude_deg` and `season_absent_meaning`, the range on files, `solar_longitude_source` | `lib/schema._check_season` | §2 ruling 3.2. The namelist's season and its range stay. |
| `units`, `long_name`, `provenance` on every variable; the provenance vocabulary; `uncertainty_kind` and its vocabulary on every companion | `lib/schema._check_variable_attributes` | Read by nothing but the schema. `refrac` keeps its own refusal of a finite companion with no convertible kind (`reduce.py`, `reduction.py`). |
| `positive` and `direction` on vertical coordinates; profile `positive = "up"` | `lib/schema.py` | Read by nothing. Monotonicity stays. |
| Kind profile's companions present and NaN | `lib/schema._check_profile`, `_PROFILE` | A policy, read by nothing. The writers still write them. |
| `<dimension>_absent_meaning` and its values, and the stated latitude for `point` | `lib/schema._check_dimensions` | Read by nothing. |
| Required globals read by nothing: `title`, `source`, `epoch`, `role` on every kind, `profile_or_run` on C, D, G, R, T, W and raw; D `citation`, `fit_inputs`, `fit_latitude_convention`; G `GM_scope`, `epoch`; R `epoch`, `citation`; W `epoch`, `observation_level_Pa`, `source_latitude_convention`; profile `boundary_pressure_Pa`, `anchor_solar_longitudes_deg`, `input_hashes`; T `vertical_coordinate`; C `vertical_coordinate` | `lib/schema.py` | The requirement goes. The writers keep writing the descriptive ones (§2 ruling 1). |
| The vocabularies only those checks used: `PROVENANCE_VALUES`, `UNCERTAINTY_KINDS`, `THERMO_INSTANCES`, `SOURCE_PROFILE_GLOBALS`, `ABSENT_MEANINGS`, `ROLES`, `VERTICAL_COORDINATES`, `POSITIVE_VALUES`, `DIRECTION_VALUES`, `VERTICAL_COORDINATE_HINTS`, `FORBIDDEN_BARE_NAMES`, `AUTHOR_REQUIRED_GLOBALS`, `SEASON_ABSENT_MEANINGS`, `PROFILE_MODELED_VARIABLES`; `KindSpec.conditional_globals`; `validate`'s `writer_filled` | `lib/schema.py` | Nothing else used them. |
| Writer globals required on read (`Conventions`, `created_by`, `created_at`, `casspian_git_commit`, `history`) | `lib/io.read` | Read by nothing that fails without them. |
| Flag variables int8 with `flag_values` and `flag_meanings` | `lib/io._check_types` | Read by nothing that fails without them; the figures read `flag_meanings` where present. |
| The group naming rule applied to groups on write | `lib/io.write` | Went with `check_names`. |
| `package_version`'s agreement with `__version__` | `lib/io.py` | §2 ruling 3.3. The installed version is recorded. |
| Kind N's hash cross-check, `_check_refractivity_hashes` | `lib/io.py` | §2 ruling 3.4. |

### 1.3 Writers (§2 ruling 1)

| Writer | Stops writing |
|---|---|
| `build_raw.py` | `profile_or_run`, `role` |
| `build_inputs.py` (T, D) | `profile_or_run`, `role`; on T `thermo_instance` and the source profile set: `spacecraft`, `event`, `observation_date`, `frequency_bands`, `longitude_deg`, `longitude_system`, `height_datum`, `source_top_boundary`, `source_gravity_citation`, `source_rotation_system`, `source_wind_citation`. `latitude_definition` (copied by `refrac`) and `raw_bundle` (provenance, on every kind) stay. The raw bundle keeps the source's statements verbatim. |
| `build_gravity.py`, `build_rotation.py`, `build_wind.py` | `profile_or_run`, `role` |
| `build_composition.py` | `profile_or_run`, `role`, `composition_role`; the default title loses the role |
| `lindal_wind.py` | `observation_level_justification` (missed at SPEC_13 Step 3) |
| `refrac/product.py` (N), `forward/production.py` and `forward/transfer.py` (profile) | `role`. `profile_or_run` stays on N and profile: the figures read it. The transfer record's per-anchor `role` (construction or validation) is a result and stays. |
| `forward/production.py` | `p_b_location` in `production_record` and the resolved namelist; the three closure comparison keys |

### 1.4 Control files (§2 ruling 2)

`role` left every section of `occul_data/lindal/lindal_build.toml`, the 17 committed run build
files under `forward/` and `tests/step13_1/sl_wind_build.toml`. `p_b_location` left the 17
committed namelists. One header comment of `forward/lindal_closure/lindal_closure_build.toml` that
named `role` and "forward checks that they do" was reworded. The suites' inline namelists
(`step04_4`, `step04_5`) lost the key too.

## 2. Kept

| Item | What it prevents |
|---|---|
| Unknown keys and sections refused, in `load_section`, the manifest and namelist parsers, the shear cases and the curve formats | A typo silently dropping a value. |
| Missing required keys, wrong types, non-positive values where a positive is used | An obscure failure in the arithmetic, or a wrong number. |
| Manifest anchor rule and quantity belonging together (`control.py`, `build_inputs.py`) | The geoid anchored on the wrong radius. |
| The reduction's cross-checks: rotation system and rate of W against R; reference level a node; anchor isobar a geodesy surface and a thermo level; geoid anchor surface; anchor quantity present; composition levels equal to the thermo levels | Each a wrong number, or an interpolation the reduction does not do. |
| `_admit_wind`: `u_total_ms` present, both poles nodes with zero wind | An obscure failure, and a singular meridional gravity. |
| Closure: one anchor; `[target]` and the transfer sections refused; `p_b_rule` required; composition levels the anchor's; gauge a tabulated level; wind coverage spans the anchor | A value silently ignored, or an extrapolation. |
| Namelist: weight 0 or 1; scale positive; one of `p_b_rule` and `p_b_Pa`; known rule; schemes, gauge rule, correlation length 0; target a latitude; season in [0, 360); date ISO | A choice silently replaced by another, or a wrong number. |
| A run writes only inside its own directory: the namelist's output directory and product, `casspian-run-inputs`'s outputs; `casspian-new-run` never writes over a run | §2 ruling 3.1 and 3.5: a run cannot overwrite registered data or a kept hypothesis. |
| `harmonic_convention` a known code | Coefficients read under the wrong convention. |
| Required variables, dimensions and groups of every kind; uncertainty companions where `refrac`, `forward` and the figures read them (T, D, G, W, N) | An obscure `KeyError`. |
| Required globals some code reads: G `harmonic_convention`; R `system_name`; W `rotation_system_name`, `rotation_rate_rad_s`, `vertical_structure`, `coverage_pressure_Pa`, `coverage_latitude_planetocentric_deg`; N and profile `profile_or_run` | An obscure `KeyError` in the loaders or the figures. |
| Coordinates strictly monotonic; profile geopotential monotonic; the tabulated pair both or neither | Wrong interpolation; an obscure failure in the figures. |
| Mole fractions sum to one | A wrong refractivity. |
| float64 for every float variable; attribute types netCDF can store | Silent loss of precision; an obscure netCDF failure. |
| `casspian_kind` present and as asked; schema version not newer | A file read as the wrong kind. |
| Profile carries `anchors/<slug>` and one record group | An obscure `KeyError`. |
| `recorded_path` across drives | An obscure failure in `os.path.relpath`. |
| The numerical refusals of `lib` and the tools: extrapolation, crossing isobars, non-convergence, shapes, the PDF tracing | Named numerical failures and extrapolation. |

## 3. Not placed

None. §2 of SPEC_R1 v0.4 ruled the seven of the pre-execution report.

## 4. Beyond the lists: kinds T and D kept by `casspian-lindal-inputs`

`build_inputs._write_unless_unchanged` keeps the kind T and D files on disk when they hold exactly
what the run would write, so that a manifest choice does not rewrite them (SPEC_01 v0.20 Step 9).
It compared every global. Once the writer stops writing `role`, `thermo_instance` and the source
profile set, the registered T and D, written before, no longer compare equal, and
`step02_6` (which runs the tool on the registered directory) rewrote them in place. Measured on the
first run of this step: `step02_6` check 1 failed on "keeps kinds T and D", and
`occul_data/lindal/lindal_thermo.nc` and `lindal_geodesy.nc` were left modified. They were
restored from git.

The comparison now drops from both sides every global the writer no longer writes (`_content`
takes the set of globals the new dataset carries). A file written before SPEC_R1 is kept, as
"old files keep theirs and still read" asks; a global the writer still writes must still be equal.
No other in-place writer was found among the suites: `refrac`'s kind N is restored by the driver,
and every other rebuild runs in a copy.

## 5. Suites

| Suite | Old | New | What changed |
|---|---|---|---|
| `step1/accept_step1` | 6 | 5 | check 5, a bare `latitude`, retired |
| `step1/verify_review_changes` | 14 | 4 | 2a to 2c, 3, 4a to 4e and 6 retired (vocabularies, `thermo_instance`, `positive`) |
| `step3/accept_step3` | 5 | 5 | its generated build file lost `role` |
| `step9/accept_step9` | 6 | 5 | check 6 retired; check 5 no longer asks for the prefix or the same directory |
| `step02_1/accept_step02_1` | 6 | 5 | check 2 (`-dirty`) retired; check 5 lost the outside-path case |
| `step02_4/accept_step02_4` | 9 | 8 | check 3 (the hash cross-check) retired |
| `step03_3/accept_step03_3` | 16 | 13 | checks 2 (season, epoch), 5 (role) and 9 (content comparison) retired; check 8 lost the label conditions; check 10 lost the role case and its outside case is now outside the run directory; check 12 lost five cases (slug, `-dirty`, prefix, the composition edit only the comparison caught, the season) and is 8 cases; check 13 lost the companion case; `relaxed.py` deleted |
| `step03_4/accept_step03_4` | 13 | 12 | check 10 (`source_statement`) retired; check 12 no longer asks for the comparison keys |
| `step04_0/accept_step04_0` | 14 | 13 | check 0 (`-dirty`) retired; check 6 keeps its own copy of the provenance list it sets aside; check 14 no longer asks for the comparison |
| `step05_1/accept_step05_1` | 11 | 11 | check 1 compares values only |
| `step05_3/accept_step05_3` | 7 | 6 | check 3 (the content comparison) retired |
| `step13_1/accept_step13_1` | 3 | 3 | `profile_or_run` and `role` added to the attributes a wind may have dropped |
| `step04_4`, `step04_5`, `step05_2`, `step06_1`, `step11_2`, `step12_1`, `step13_2` and the helper scripts of `step04_0`, `step04_4`, `step05_4`, `step06_1` | | | the in-process `-dirty` relaxation removed, and its prose |

`tests/run_regression.sh` lists the seven newer suites after `step08_1` with 3, 2, 5, 2, 4, 3 and
2, and the amended counts above.

Two expectations moved during the regression, with the counts unchanged:
- **`step03_3` check 7**, beyond the specification, asked that the registered kind T read silently
  in place. It now warns once, for its `input_hashes` entry `lindal_build.toml`, which this step
  edited (§1.4). That is the warning working. The check now asks for silence apart from that entry,
  and prints both.
- **`step05_3` check 6**: the second `casspian-new-run` is refused by `mkdir`'s `FileExistsError`,
  not by a `ControlFileError` (§2 ruling 3.5). The suite catches that exception.

## 6. The full regression

`tests/run_regression.sh`, 38 suites, run in six chunks so as to stay inside the background-job
limit (`reports/regression/regression_R1_chunk{A,B,C1,C2,D,E}.txt`), each chunk restoring the
registered kind N and the closure and transfer files and reporting them equal. Per-suite logs:
`reports/regression/regression_accept_*.txt`.

- **First pass: 35 of 38 at their counts.** `step03_3` 12 of 13 (check 7), `step04_0` 12 of 13
  (check 14 reruns `step03_3`), `step05_3` stopped at check 6 on the `FileExistsError`. Each is
  explained in §5.
- **After the two expectation changes**, the three rerun through the driver
  (`reports/regression/regression_R1_rerun.txt`): `step03_3` 13 of 13, `step04_0` 13 of 13,
  `step05_3` 6 of 6. **38 of 38 at their reference counts.**
- **`step05_3` checks 4 to 7 ran** for the first time since SPEC_11 Step 2: the closure product
  equals the registered one in all nine computed variables; the transfer product's N, p, T,
  altitude, r0 and pressure identity equal the registered ones; the new run and F9 as before.
- **No value changed.** Every check that compares rebuilt values with registered ones passed:
  `step02_6` (T and D kept), `step04_0` check 6 (every product, group by group), `step04_5` check
  11, `step05_2`, `step05_3` checks 4 to 6, `step06_1`'s rebuilt closure. The registered kind N
  hash is `b79c60cb...` before and after every chunk.
- The seven newer suites, run through the driver for the first time: `step09_1` 3, `step10_1` 2,
  `step11_1` 5, `step11_2` 2, `step12_1` 4, `step13_1` 3, `step13_2` 2.
- Suite time: 9775 s (2 h 43 min) in the first pass, plus 529 s for the three reruns.

## 7. Not done here

- The untracked reference documents in `docs/` and the specification files go into the
  acceptance commit (§2 ruling 4). Until then every product built carries `-dirty`, which no
  longer blocks a run.
- `tests/step03_0/` still carries its own `relaxed.py`. The suite is not in the driver (it checks
  a one-time rebuild against `before/` products today's reader no longer reads), and its relaxation
  patches `_load_into_memory`, which still exists, so it is left as it is.
- SPEC_00 v0.25's cut of §7.2 and §8 to the kept list is the reviewing agent's.
