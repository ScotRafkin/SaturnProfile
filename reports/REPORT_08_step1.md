# REPORT 08, Step 1. Kind W holds what the data hold

Coding agent, 7 October 2026. Specification: `docs/specs/SPEC_08_Kind_W_Parts.md` v0.3, §2 as §6
rules on it. Working tree on `main` at `104f031`, uncommitted. Acceptance:
`tests/step08_1/accept_step08_1.py`, **6 of 6 checks pass**; output in
`reports/step08_1/output.txt`, figures in `reports/step08_1/figures/`. Regression: §4.

**In short.** Kind W now requires the reference wind and its level, and nothing else of the
wind. The sum identity, its dimension check and the polar rule are gone from `validate`. The three
loaders refuse a wind without `u_total_ms` by name, then apply the polar rule. The shear tool
refuses a source without the total, and writes a shear only where its source carries one. F9 and
`casspian-render` draw the parts a file carries. No arithmetic changed and no registered file moved.
`step04_0` drops from 15 checks to 14, `step7` check 2 is restated, and seven suites change wording
only.

## 1. The changes

| File | Change |
|---|---|
| `lib/schema.py` | Kind W table: `u_total_ms` and `u_shear_ms` optional; `u_reference_ms` and `reference_level_pressure_Pa` required, as before (ruling on what kind W is). `check_wind_components` deleted (R3 dropped with it). `validate` no longer calls `check_wind_poles`, which stays for the loaders; its docstring says where it is applied. The kind's `notes` restated. |
| `lib/control.py` | `_admit_wind(wind, where)`: refuses a wind without `u_total_ms`, naming it, then applies `check_wind_poles` as a `ControlFileError`. It is called directly after the rotation check in `load_reduction_inputs` (new), `load_run_inputs` and `_load_transfer_inputs` (replacing their `check_wind_components` and `check_wind_poles` block). Docstrings of both public loaders restated. |
| `lib/io.py`, `lib/windfield.py`, `refrac/anchor.py` | Docstrings only (R8): `read` no longer lists the sum identity; `reference_wind` says the reduction's loader, not kind W, requires the level as a node; `wind_of_latitude` says the poles are checked by the loaders, not by `read`. |
| `tools/wind/shear.py` | `construct` refuses a source without `u_total_ms`, by name, for every case. `assemble` replaces `u_shear_ms` only where the source carries it, and never creates it (ruling 2). Module and `assemble` docstrings restated. |
| `tools/plots/figures_profile.py` | F9 draws the source line only when `inputs/wind` carries `u_reference_ms`; its data dictionary then omits the key (R6). |
| `tools/plots/figures_inputs.py` | `_figure_wind` (`casspian-render` of kind W, ruling 10). Left panel: the total at the reference level, with its band, where the file carries the total on that level, as before; otherwise the file's `u_reference_ms`, with no band, and provenance shading only where the level is a pressure node. Right panel: `u(p)` of the total where present, `u_shear(p)` where only the shear is, a note and no curves where neither. A one-level file is drawn with markers. `wind_profile_at` takes the variable name, defaulting to the total. |

Nothing else in `lib`, `refrac` or `forward` changed (deliverable 3). `casspian-wind-from-curve`
is unchanged (deliverable 6). A registered wind takes exactly the paths it took before, because it
carries all four variables, its total on the reference level, and zero poles.

## 2. Acceptance, `tests/step08_1/accept_step08_1.py`

1. **Nothing moves.**
   - `git status` shows the five registered winds, the two run products and the kind N file
     unmodified.
   - Each registered wind reads as kind W.
   - Each of the three loaders admits the registered inputs.
2. **A file holds what the data hold** (ruling 11). Three files are cut from the reduction's wind
   on its 339 latitudes inside plus or minus 85 deg (`-84.5` to `84.5`). Each is written through
   `lib.io.write`, read back, and every variable is array-equal:
   - the reference wind on one level;
   - the reference wind with a shear on three levels (`u_reference x 0.1 ln(p / p_ref)`, nonzero
     so equality means something), no total;
   - the reference wind with the total on one level.

   The shear and the level without the reference wind are refused on read: "required variable
   'u_reference_ms' is missing (kind wind)".
3. **The model takes only what it can run on** (ruling 4). All 9 cases are refused, each message
   naming its variable:
   - each loader, given a copy without `u_total_ms` and its companion, refuses it with "the wind
     carries no u_total_ms";
   - a copy without `reference_level_pressure_Pa` or without `u_reference_ms` is refused on read by
     the schema, inside the loader.
4. **The polar rule at the model.** For each loader, a copy of its registered wind with the north
   pole at 1e-6 m/s reads through `lib.io`. The loader refuses it with "u_total_ms at +90 degrees
   is 1.000e-06 m/s, not exactly zero". In closure and transfer mode it is now the polar rule that
   fires, not the sum identity (pre-execution R4).
5. **The shear tool.** The source is the transfer run's wind with `u_shear_ms` removed. Under the
   case of SPEC_05 Step 1 check 5, the output's total departs from `u_s F(p)` by 0.000e+00 relative
   (bound 1e-12), and the output has no `u_shear_ms`. A source without `u_total_ms` is refused by
   name, and nothing is written.
6. **Figures** (rulings 6, 10, 12). For the author to view, all under `reports/step08_1/figures/`:
   - `F9_transfer.png`: both lines;
   - `F9_without_reference.png`: the total line only, drawn from the registered transfer product in
     memory with `u_reference_ms` dropped from `inputs/wind`;
   - `casspian-render` of check 2's three files, in `reference_only/` (right panel blank, with its
     note), `reference_and_shear/` (`u_shear(p)`) and `reference_and_total/` (`u(p)`, one level, as
     markers).

**One reading of the specification, stated.** §2 check 5 asks for "an input without
`u_reference_ms`". Since ruling 1, such an input is not a kind W file. By ruling 2 the case that
remains is an input without `u_shear_ms`, and check 5 tests that.

## 3. Accepted suites

| Suite | Change | Old | New |
|---|---|---|---|
| `step04_0` | check 4 retired (decision 3); a one-line marker keeps the numbering | 15 checks | 14 |
| `step7` | check 2 restated (R7): `u_total_ms` alone set to 1e-6 m/s at the north pole; the copy reads, and `load_reduction_inputs`, pointed at it through the registered manifest, refuses it. Check 6: description and one detail line no longer claim the polar check on read | 8 | 8 |
| `step05_1` | checks 3 and 8: "sum identity included" and "sum identity on read" removed from the description and the detail | 11 | 11 |
| `step04_2` | check 8: detail says the sum identity is computed here from the file, and that `write` validates the schema | 13 | 13 |
| `step04_4`, `step04_5` (and `step04_5/fields.py`, its docstring) | the report line no longer says the reader checks the sum identity; it gives the value computed from the file read back | 11, 11 | 11, 11 |
| `step04_3` | the two docstrings of R9: "not a kind W file" becomes "not a wind the model runs on" | 9 | 9 |
| `step04_1`, `step02_6` | left, as R9 rules: their texts stay true | 13, 7 | 13, 7 |
| `step08_1` | added to `tests/run_regression.sh` at 6 | | 6 |

All conditions are unchanged except `step04_0` check 4 (retired) and `step7` check 2 (restated).

## 4. Regression

The full set, as §2 requires for a change to `lib`, by `tests/run_regression.sh` on the registered
products as committed, with this step's working tree: **31 of 31 suites at their reference
counts**. That is the 30 of the record, with `step04_0` at its new 14, and `step08_1` at 6. The run
took 7046 s in all; the longest suites were `step04_3` (1476 s), `step05_4` (1040 s), `step04_4`
(962 s) and `step04_5` (947 s). The registered kind N file and the closure and transfer inputs and
products were restored equal; the kind N sha256 is `b79c60cb...` before and after. No suite left a
path dirty.
- The run began with 21 paths of the step's own.
- From `step6` on it shows 22: the extra path is this report, written while the run was in
  progress.

The output is in `reports/regression/regression.txt`.

The author asked during the run whether the full set was needed. On that question, this step
reached nearly the whole set anyway. 26 of the 31 suites read a kind W file or go through a loader.
The ones that exercise what changed beyond the read check are:
- `step04_0` and `step7`, whose refusal tests were retired or restated;
- `step04_4` and `step04_5`, which write and run their own sheared winds;
- `step05_1` to `step05_4`, whose winds go through the shear tool, with F9 in `step05_3`;
- `step08_1`.

No number could move, since neither the data nor the arithmetic changed. The rerun shows that no
suite's inputs are now refused, and that no refusal test depended on the checks removed from
`validate`.

## 5. Notes

- **`tests/step08_1/preexecution_measure.py`** is the reading's record at `104f031`. It simulates
  the v0.2 schema, so it no longer describes the tree, and it is not in the regression, as
  `tests/step07_1/preexecution_measure.py` was not.
- **The line endings.** `.gitattributes` sets `eol=lf`; every edited file was checked LF in the
  working tree (`git ls-files --eol`).
- **`wind_profile_at` interpolates in latitude with `np.interp`,** which holds the end value
  outside the grid. For a data file narrower than plus or minus 60 deg the right panel would show
  held values at the outer latitudes. This is the existing behavior, met in no file here (check 2's
  files reach 84.5 deg), and I have not changed it.
