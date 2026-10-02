# REPORT 05, Step 3. The shear step in every forward run, and the new-run tool

Coding agent, 2 October 2026. Specification: `docs/specs/SPEC_05_Shear_Experiments.md` v0.6, §3.
Working tree on `main` at `32278b0`; nothing of this step is committed.

**Result.** `tests/step05_3/accept_step05_3.py` passes 7 of 7. The migration moves no computed
value: the closure product equals the registered one under `step04_5` check 11's comparison, and
the transfer product, and a run made by `casspian-new-run`, equal the registered transfer product
in delivered N, p, T, altitude, `r0` and pressure identity. The full regression, all 29 suites, passes at every reference count. The registered run files are rebuilt by the sweep after the acceptance commit
(section 4).

## 1. What changed

| File | Change |
|---|---|
| `src/casspian/tools/run/run_inputs.py` | Five sections in order, gravity, rotation, wind, shear, composition; `[shear]` required; its `role` and `prefix` checks skipped (it carries neither), its output checked like every section's. |
| `forward/lindal_closure/lindal_closure_build.toml`, `forward/lindal_transfer/lindal_transfer_build.toml` | `[wind] output` is `inputs/<run>_wind_source.nc`; `[shear]` added, `source` that file, `output` `inputs/<run>_wind.nc`, `case = "identity"`. No namelist changes. |
| `src/casspian/tools/run/new_run.py`, `pyproject.toml` | `casspian-new-run <name> --from <run directory>`: the two control files copied beside the source with every whole-token occurrence of the source run's name replaced (not preceded by a letter, digit or underscore, not followed by a letter or digit, so `<run>_` prefixes follow), nothing else changed; refuses an existing directory and a source without its two files. |
| `src/casspian/tools/plots/figures_profile.py` | F9's left panel: the source wind labelled `source wind, assigned to <p> mbar`, and `u_total at <p> mbar` drawn beside it, read at the file's `reference_level_pressure_Pa` by `WindField`, the model's rule. The right panel draws a constant field as one colour with its value (section 3, at the author's request). |
| `docs/RUNBOOK.md` | v0.4: the five sections, `casspian-new-run`, and `pip install -e . --no-deps` after a pull that adds a console entry. |
| `tests/step04_0/sweep.py` | `--runs`: rebuild only the two runs' inputs and products (section 4). The kind `_wind_source.nc` is known to its hash table. |
| `tests/step05_3/accept_step05_3.py`, `tests/run_regression.sh` | The acceptance, and its row `step05_3/accept_step05_3 7`. |

## 2. The acceptance

Everything runs in a copy, `reports/step05_3/tree/`, in the repository's layout: `data_static/`,
the registered Lindal reduction as it is (this step does not reach it), and the two runs' migrated
control files. The runs' inputs and products are rebuilt there from the working tree with the
`-dirty` refusal relaxed for the script only, and compared by value with the registered files.

| Check | Measured | Result |
|---|---|---|
| 1. no `[shear]`: refused naming the five sections | `... carries exactly the sections ['gravity', 'rotation', 'wind', 'shear', 'composition']; unknown none, missing ['shear']` | pass |
| 2. both build files write both wind files; `_wind.nc` array-equal to `_wind_source.nc` | each run writes five files; 12 of 12 variables equal in both runs | pass |
| 3. the closure content comparison against the anchor's embedded copies | composition, gravity, rotation and wind identical; the identity-sheared wind's new `history` and `input_hashes` are among the dropped attributes | pass |
| 4. closure product against the registered one (`step04_5` check 11's nine variables) | 9 of 9 array-equal | pass |
| 5. transfer product: N, p, T, altitude, `r0`, pressure identity | 6 of 6 array-equal | pass |
| 6. `casspian-new-run` | 18 lines of the build file and 7 of the namelist differ, every one the name (renaming back gives the source exactly); a second call refused (`... exists; casspian-new-run makes a new run and never writes over one`); the new run's product 6 of 6 array-equal to the registered transfer product | pass |
| 7. F9 carries both lines | migrated transfer: both labelled, `u_total` at 1000 mbar array-equal to the source wind (coincident); `uniform, c = 0`: both labelled, `u_total` at 1000 mbar exactly zero under a source wind up to 490.5 m/s | pass |

Renders for viewing, written by check 7: `reports/step05_3/F9_migrated_transfer.png` and
`reports/step05_3/F9_uniform.png`.

## 3. F9's constant field (the author's request)

The first render of the `c = 0` run showed the right panel's colour bar running from -1.5e-14 to
1.5e-14: with nothing to span, matplotlib's automatic contour levels were set at round-off, a scale
that reads as structure. A field whose largest and smallest values are equal is now drawn as one
colour, the colour bar carries its one value, and the panel says `u_total = <value> m/s
everywhere`. A field with any range is drawn as before. The author confirmed the case: `uniform`,
`c = 0`, `u_total` exactly zero at every node.

## 4. The registered files

The migration changes what the two runs' build files write, so the registered run files change:
`<run>_wind_source.nc` is new in each run's `inputs/` (committed, under the existing `.gitignore`
exceptions), `<run>_wind.nc` is now the shear step's output (an identity copy, with its own
`input_hashes` and writer's stamps), and both products embed that file. The computed columns do not
move (checks 4 and 5). They are rebuilt on the clean tree after the acceptance commit by
`python tests/step04_0/sweep.py --runs`, which rebuilds only the two runs' inputs and products: the
reduction chain is not reached by this step, and rebuilding it would change every chain file's
bytes and so the hashes `step02_1` and `step02_6` read, for no change of content. The rebuilt
files are then compared with the committed ones by value and committed, and the comparison is
added to this report.

## 5. The regression

`tools/run`, the build files and F9 reach every suite that builds forward inputs, reads a forward
product or draws a figure; with the regression now about 100 minutes, the full set was run.

| Suite | Checks | Time (s) | | Suite | Checks | Time (s) |
|---|---|---|---|---|---|---|
| `step1/accept_step1` | 6 of 6 | 3 | | `step02_6/accept_step02_6` | 7 of 7 | 297 |
| `step1/verify_review_changes` | 14 of 14 | 3 | | `step03_1/accept_step03_1` | 9 of 9 | 9 |
| `step2/accept_step2` | 8 of 8 | 2 | | `step03_2/accept_step03_2` | 8 of 8 | 15 |
| `step3/accept_step3` | 5 of 5 | 3 | | `step03_3/accept_step03_3` | 16 of 16 | 157 |
| `step4/accept_step4` | 7 of 7 | 0 | | `step03_4/accept_step03_4` | 13 of 13 | 70 |
| `step5/accept_step5` | 7 of 7 | 43 | | `step04_0/accept_step04_0` | 15 of 15 | 288 |
| `step6/accept_step6` | 6 of 6 | 29 | | `step04_1/accept_step04_1` | 13 of 13 | 27 |
| `step7/accept_step7` | 8 of 8 | 89 | | `step04_2/accept_step04_2` | 13 of 13 | 533 |
| `step8/accept_step8` | 9 of 9 | 3 | | `step04_3/accept_step04_3` | 9 of 9 | 1488 |
| `step9/accept_step9` | 6 of 6 | 5 | | `step04_4/accept_step04_4` | 11 of 11 | 867 |
| `step02_1/accept_step02_1` | 6 of 6 | 11 | | `step04_5/accept_step04_5` | 11 of 11 | 895 |
| `step02_2/accept_step02_2` | 7 of 7 | 43 | | `step05_1/accept_step05_1` | 11 of 11 | 21 |
| `step02_3/accept_step02_3` | 9 of 9 | 84 | | `step05_2/accept_step05_2` | 3 of 3 | 194 |
| `step02_4/accept_step02_4` | 9 of 9 | 164 | | `step05_3/accept_step05_3` | 7 of 7 | 185 |
| `step02_5/accept_step02_5` | 7 of 7 | 180 | | | | |

29 of 29 at their reference counts, 95 minutes of suite time; the registered kind N
(`b79c60cb...`), the closure inputs and product, and the transfer inputs and product restored
exactly, and afterwards no registered file differs from the commit.

**A run that was discarded.** A first attempt at this regression was started in a way that kept it
running out of sight while a second was started; the two drivers set aside and restored the same
registered files, several suites crashed and one file was locked. Both were stopped, every
registered file was verified identical to the commit (`git diff -- '*.nc'` empty, no
`skip-worktree` mark), and the regression above was run once, alone.

## 6. Accepted suites whose expectations change

None. `step03_3` (the fixture comparison of its check 9 and the run input list of check 8) and
`step04_0` check 6 (the closure wind against the swept fixture under its list of allowed
differences) pass unchanged on the migrated build files, as do `step04_4` and `step04_5`, whose
namelists name `<run>_wind.nc`, still the file the namelist reads.

## 7. Findings

1. **The runbook has statements that are now out of date**, outside this step's deliverable and
   left for the author: section 1 says the repository holds no netCDF (the registered products are
   committed since `4d32efe`); section 4 quotes the closure product's SHA-256 at `e71d59e`; section
   6 gives about seven hours for the regression (now about 100 minutes) and `26 of 26` (now 29).
   Corrected in runbook v0.4 under REVIEW_05_step3: the committed products stated, a rebuild compared
   by value rather than by hash, about 100 minutes and 29 of 29.
2. **`casspian-new-run`'s token rule** replaces the source name when it is followed by `_`, so that
   `<run>_` prefixes follow. A source name that is a prefix of another name in the same files (a
   run `a` in files that also mention `a_b` as a different run) would be renamed there too; no
   current run's files contain one.
3. **F9 of an identity run keeps the source's `vertical_structure`** (`altitude independent`) in its
   right panel title, as §1.7 says identity leaves every attribute as it is.

## 8. What is committed on acceptance

The source, control, test and runbook files of section 1 and this report; then, after the sweep,
the rebuilt registered run files.
