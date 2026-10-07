# REPORT 08, pre-execution. The coding agent's reading of SPEC_08 v0.2

Coding agent, 7 October 2026. Specification: `docs/specs/SPEC_08_Kind_W_Parts.md` v0.2, §4.
Working tree on `main` at `104f031` (the report-only commit after `84c870a`; nothing under `src` or
`tests` differs). Measurements by `tests/step08_1/preexecution_measure.py`, output
`reports/step08_1/preexecution.txt`. It simulates deliverable 1 in memory (the four parts made
optional in the kind W table, the two wind checks removed from `validate`) and leaves `lib.control`
and the shear tool as they are on disk, so it shows what the loaders and the tool do today with the
inputs SPEC_08 lets through. It writes only under `reports/step08_1/preexecution/`.

**In short.** §1's list of model loaders is complete and its description of the code is right in
all but two details (item 1). No accepted suite beyond those §2 names asserts the sum identity
or the polar rule through `lib.io.read` or `validate`. Three more suites carry a comment that
becomes inaccurate, and their conditions do not change (item 2). Nothing in `lib.io` relies on the
two checks (item 3). No registered file or fixture changes, and no fixture becomes refused (item 4).
Two deliverables cannot be written as stated, both in the shear tool (R1, R2). One check is lost
with `check_wind_components` and should be kept (R3). The rest are construction choices for checks
3, 4 and 6 and docstrings that become false (R4 to R9).

## 1. Every reader of the four variables in `src`

| Reader | Reads | Path |
|---|---|---|
| `lib.windfield.WindField.__init__` | `u_total_ms`; `reference_level_pressure_Pa`, eagerly | transfer (`forward.transfer`, three sites); shear tool; F9 |
| `refrac.anchor.wind_of_latitude` | `u_total_ms` at the column where `pressure_Pa == reference_level_pressure_Pa` (index `[0]`, so the level must be a node) | reduction (`geoid_setup`); closure (`forward.production.produce`) |
| `lib.control.load_reduction_inputs` | `reference_level_pressure_Pa`, refused when not a node | reduction |
| `lib.control.load_run_inputs`, `_load_transfer_inputs` | all three parts through `check_wind_components`; `u_total_ms` through `check_wind_poles` | closure, transfer |
| `lib.schema.check_wind_components`, `check_wind_poles` | all four | every `read` and `write` of kind W |
| `tools.wind.shear` | `u_total_ms` and `reference_level_pressure_Pa` through `WindField`; `u_reference_ms` and `u_shear_ms` in `assemble` | the `[shear]` build step |
| `tools.plots.figures_profile.figure_9` (F9) | `u_reference_ms`, `u_total_ms`, `reference_level_pressure_Pa`, and `WindField` | transfer products |
| `tools.plots.figures_inputs.wind_reference_column`, `wind_profile_at` | `u_total_ms`, `u_total_uncertainty_ms`, `value_provenance` at the reference level | kind N figures F1 and F3 (`figures_product`); `casspian-render` of a kind W file (`_figure_wind`) |
| `tools.wind.build_wind` | writes all three parts and the level | `casspian-wind-from-curve` |

§1's three loaders are the only places the model takes a wind; `forward.estimate` and
`forward.propagate` read none of the four. `u_shear_ms` has no reader outside the schema and the
shear tool. Two details of §1 differ from the code:

- **`WindField` reads the reference level too,** in its constructor, not only `u_total_ms`. That is
  why a file without the level cannot run the transfer, and it is the obstacle in R1.
- **The transfer production does not record the reference level.** `forward.transfer` reads it only
  through `WindField.__init__`, which stores it for `reference_wind()`. No code in `src` calls
  `reference_wind()`; `step04_1` check 7 does, in its detail line only. The level reaches the product only inside the
  verbatim `inputs/wind` copy. Decision 1 stands, because the transfer cannot build its field
  without the level. Its last sentence would read more accurately as "a change to
  `lib.windfield`", since making the read lazy is one line there and nothing in `forward` changes.

## 2. Suites that assert the sum identity or the polar rule through `read` or `validate`

Every `tests/` script was searched for the four names, `check_wind`, "sum identity", "pole", and
every `CasspianSchemaError` or `cio.read(..., "wind")` refusal. Beyond §2's list:

- **None asserts either rule.** The other refusal checks that touch kind W are listed here.
  - `step1` check 3: a raw file relabeled `wind`. It is refused by the missing global
    `rotation_system_name` today and under the simulated schema; the message is identical in both
    runs.
  - `step1/verify_review_changes` 2b and 2c: attribute checks on the companion, which do not change.
- **No diagnose script asserts either rule.** `step04_4/diagnose_anchor_placement`, the two
  `step04_5` diagnose scripts, `step04_5/fields.py`, `step05_4/run_experiments` and
  `step06_1/cycled_runs` write sheared winds with all three parts consistent and only rely on
  `write` admitting them.
- **Comments that become inaccurate, conditions unchanged** (named so that none is edited silently,
  see R9):
  - `step04_3`, the module docstring (line 14) and `synthetic_field`'s docstring (line 230): "a field
    that is nonzero at the poles is not a kind W file". After SPEC_08 it can be a kind W file. It
    is not one the model runs.
  - `step04_1` check 13, detail text: the synthetic field is "never read through the schema". That
    stays true.
  - `step02_6` check 7, comment (line 270): "scaling all three keeps the sum identity". That stays
    true, though nothing requires it any more.
- **`step04_0` checks 1 and 3** assert three parts and an exact sum on the three registered winds.
  They compute this from the file, and it stays true of every file the present tools write. They
  are of the same kind as the §2 wording items, and no wording there mentions the reader.

## 3. `lib.io` and the two checks

Nothing in `lib.io` relies on them.

- **`write`** calls `validate` on the root dataset only. Groups get `check_names` and the attribute
  type check, never a kind's `validate`.
- **`read`** validates the root under the asked kind. A kind `profile` or `refractivity` file's
  embedded `inputs/wind` group is never validated as kind W, on write or on read. The closure
  loader's comparison with the anchor's copy (`check_closure_inputs`) compares content and needs
  neither rule.
- **The only dependence is wording:** `read`'s docstring lists "wind components that do not sum to
  their total" among its refusals (R8).

## 4. Whether a registered file or fixture could change

**None can.** No arithmetic changes. `build_wind` and `casspian-wind-from-curve` are unchanged. The
shear tool's output is unchanged for an input carrying all four variables, and every registered
input does: measured on all five registered winds, each 361 by 61, reference level 1e5 Pa a node,
`u_total` exactly 0 at both poles, sum departure exactly 0. No product embeds a wind that is
re-derived. No suite writes its `output.txt` into git.

**One fixture changes behavior, not content.** `tests/step04_0/fixtures/swept/` holds two winds in
the pre-SPEC_04 form (with `u_cylindrical_ms`). `read` refuses them today only because
`u_reference_ms` is required (measured). After SPEC_08 they read; since they carry `u_total_ms`,
the level and zero poles, the loaders would admit them. No suite reads them through `lib.io`:
`step04_0` opens them with `open_any`, which applies no schema check. This agrees with §1, which
checks no extra variable. The `step03_3` "before" wind already reads today and is unchanged.

## 5. For the ruling

**R1. The shear tool cannot carry an input without the reference level (deliverable 5 against
deliverable 3).** `u_at_pressure` builds a `WindField`, whose constructor reads
`reference_level_pressure_Pa`. An input without it fails with a `KeyError` before any case runs
(measured). "Carried when present" therefore needs one of two things:

- (a) The tool also requires `reference_level_pressure_Pa` and refuses an input without it by
  name, as it does `u_total_ms`. Its output could not run the model anyway (§1). **Recommended.**
- (b) `WindField` reads the level lazily. That is a `lib` change, outside deliverable 3.

**R2. The shear part when the input's parts are incomplete.** Today `assemble` writes into the
input's own `u_shear_ms` (measured `KeyError` when it is absent). Two cases need a rule:

- **`u_reference_ms` without `u_shear_ms`.** Deliverable 5 writes a shear, but the variable must be
  created with `units`, `long_name` and `provenance`, which the input does not give. Two options:
  - (a) write the shear only when the input carries both parts (replace, never add);
    **recommended**;
  - (b) create the variable with stated attributes (`provenance = "derived"`).
- **`u_shear_ms` without `u_reference_ms`.** The deep copy would carry the input's shear, which is
  false of the new total. I read "the output has no shear part" as dropping it, and ask that the
  ruling confirm this.

**R3. A dimension check is lost.** `check_wind_components` holds the only check, for any kind, that
`u_total_ms` and `u_shear_ms` lie on `(latitude_planetocentric, pressure)`. `VarSpec.dims` is
declared in the table but no code in `lib.schema` enforces it. After deletion, a transposed
`u_total_ms` on a square grid would be read wrong by `wind_of_latitude` (`[:, column]`) with no
refusal. `WindField` checks the shape only. Under decision N this is a schema fault. Options:

- (a) keep the dimension half of `check_wind_components` in `validate`, applied to whichever of
  `u_total_ms`, `u_shear_ms` (latitude and pressure) and `u_reference_ms` (latitude) is present;
  **recommended**;
- (b) check `u_total_ms` only, in the loaders;
- (c) drop it.

**R4. Where the loaders refuse (a statement of how I would write it, for confirmation).**
`check_wind_poles` returns silently when `u_total_ms` is absent. In closure mode, the comparison
with the anchor's copy fires before anything that names the variable. The measurements on the
simulated schema show the gap:

| Copy | Reduction loader | Closure loader | Transfer loader |
|---|---|---|---|
| without `u_total_ms` | admitted | refused, not a closure | admitted |
| without the level | bare `KeyError` | refused, not a closure | admitted |
| north pole at 1e-6 m/s | admitted | refused by the sum identity, which fires before the polar rule | refused by the sum identity, which fires before the polar rule |

I would place the presence check, then the polar rule, directly after the rotation check in all
three loaders, so that check 3's message names the variable and check 4's comes from the polar
rule.

**R5. How checks 3 and 4 hand a copy to a loader.** The loaders accept inputs only from the
manifest's own directory with the slug prefix, or from the run's `inputs/` with the run prefix. A
copy under `reports/step08_1/` cannot be named in a manifest or namelist without copying the whole
input set and the namelist. Two parts to the proposal:

- **How the loader finds the copy.** Read the registered manifest and namelists with the real
  readers, then point the frozen object's `inputs["wind"]` at the copy with `dataclasses.replace`,
  as `step02_6` does with its inputs. The measurement script does this. The alternative is to copy
  each directory under `reports/step08_1/`.
- **How the copy is written.** With `xarray.to_netcdf`, attributes kept, not `cio.write`.
  `cio.write` would stamp `-dirty` from the working tree, and the dirty refusal would fire first.
  Each loader takes its own registered wind, since the prefix rule ties the file to the run.
  Measured cost: the three loaders on the registered inputs take 3.0, 7.4 and 7.9 s, so checks 3
  and 4 take about 2 minutes.

**R6. How check 6 gets a run without `u_reference_ms`.** No present tool builds one: `build_wind`
writes all three parts and the shear tool carries what its source has. Building one means a
hand-made source and a transfer run. Proposal: draw F9 from the registered transfer product held
in memory with `u_reference_ms` dropped from its `inputs/wind` group. F9 recomputes nothing from
the run (its docstring), so the figure is the one such a run would give. The figure's data
dictionary omits `u_reference_ms` when the file lacks it. `step05_3` check 7 reads it for the
transfer run, where it is present, and is unchanged.

**R7. The restated `step7` check 2.** The copy perturbed `u_total_ms` and `u_reference_ms` together
to get past the sum identity. Restated, I would perturb `u_total_ms` only, the variable the rule
reads, and hand the copy to `load_reduction_inputs` as in R5.

**R8. Docstrings in `src` that become false.** Deliverable 3 forbids any change to `refrac` and to
the rest of `lib`. I ask whether docstrings count; the recommendation is to restate them (text
only, no code):

- `lib.io.read` (the sum identity among its refusals);
- `lib.control.load_run_inputs` (line 1264, "the wind's parts do not sum");
- `refrac.anchor.wind_of_latitude` (the poles "checked by `read`"; after this step, checked by the
  loader);
- `lib.windfield.WindField.reference_wind` ("a pressure node, which kind W requires"; the
  reduction's loader requires it, not kind W);
- the `tools.wind.shear` module docstring and `assemble` (the sum identity, and "a non-finite value
  reaches the poles, where the schema refuses the file"; after this step the tool writes such a
  file and the run's loader refuses it, since `max|NaN| != 0`). No present case produces a
  non-finite value.

**R9. The suite comments of item 2.** The proposal is to restate the two `step04_3` docstrings and
leave the `step04_1` and `step02_6` texts, which remain true. Each would be named in the report.

**Noted, no change proposed.**

- **`casspian-render` on a data-only kind W file.** `_figure_wind` reads `u_total_ms` and the level,
  so a shear-only or reference-only file fails with a `KeyError`. SPEC_09 will want to view its
  digitized file. I would leave this to SPEC_09.
- **No companion for the reference or shear part.** The table gives none. Measured: a
  reference-only file carrying `u_reference_uncertainty_ms` with `uncertainty_kind` writes and
  reads, since extra variables are admitted. This is for SPEC_11 to name.
- **A one-level pressure coordinate.** `_check_coordinates` skips any coordinate with fewer than 2
  values, so a one-level file is not asked for `positive` or `direction`. Check 2's total-only
  file passes in part because of this. The behavior is existing.
- **The reference level as a pressure node.** The forward loaders do not require it, and only the
  closure production needs it. That run's wind is content-identical to the anchor's copy, which
  passed the reduction's loader, and the transfer interpolates.
- **Check 2 measured on the simulated schema.** Its three files (shear only on three levels, total
  only on one, reference only), on latitudes inside plus or minus 85 deg, write and read back
  array-equal.
- **The regression of record.** It took 6923 s (about 115 minutes) at `7652852`, not the 100
  minutes §2 states.
