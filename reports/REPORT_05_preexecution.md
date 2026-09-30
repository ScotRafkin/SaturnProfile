# REPORT_05_preexecution: SPEC_05 v0.3 read against the repository

Coding agent, 30 September 2026. SPEC_05 v0.3 read against `main` at `4d32efe`, as the handoff of
29 September section 3 asks. Nothing was changed. The author's three rulings given in chat on 30
September (the registered transfer product, the Step 2 recommit of both registered products, the
Step 1 line that the input is not written over) are taken as given and not repeated here, except
where a finding below touches them.

Each finding says what the specification states, what is on disk, and what I propose. Findings
1 to 4 are places where the specification's description of the code does not match the code.
Findings 5 to 8 are acceptances I cannot write as stated. Findings 9 to 15 are choices the
specification leaves open that I think the author should make.

## 1. Where the specification does not match the code

**1. The closure bit-identity check is in `step04_5`, not `step03_4`.** §3 deliverable 3 and Step
2 check 4 cite "the bit-identity check `step03_4` already makes". `step03_4` check 1 compares the
closure product's embedded inputs with the files on disk, group by group; it does not compare the
product with the registered one. The check that does is `step04_5` check 11
(`tests/step04_5/accept_step04_5.py:953`): it reruns the closure namelist in a copy and compares
nine computed variables with the registered product by `array_equal`, never the embedded copies or
the attributes. Proposed: cite `step04_5` check 11. That check is also the right one for ruling 2,
since it compares exactly the computed columns.

**2. F9 is drawn in transfer mode only.** `tools/plots/figures_profile.py:484` draws F9 only for a
transfer product. Step 2 check 7 uses two transfer runs, so the check can be written; the
specification should say that the closure run has no F9, so nobody expects one.

**3. A namelist has no `title`.** §3 deliverable 4 has `casspian-new-run` change "every prefix,
output name, `[run] name` and title". The build file sections carry `title`; the namelist carries
`[run] description`, which does not name the run in `lindal_transfer.toml`, and `[output] product`,
which does. Both control files also name the run in their header comments. See finding 12 for the
rule I propose.

**4. `step03_4` is not the only accepted suite that names the runs' wind files.** §3 names
`step03_3` as one a migration may reach. Four suites name `inputs/lindal_closure_wind.nc` or
`inputs/lindal_transfer_wind.nc`: `step03_3` (its run input list and its namelist cases), `step04_0`
(check 0's `-dirty` copy and check 6, which compares the closure wind with the swept fixture under
a list of allowed differences), `step04_4` and `step04_5` (namelists they write). The names stay
valid after the migration, since `[shear]` writes `<run>_wind.nc`. What can move is `step04_0`
check 6: the identity-sheared wind carries a new `history` line, new `input_hashes` and new writer
globals. Whether its allowed-differences list already admits them is a measurement for Step 2, and
the report will name the suite either way.

## 2. Acceptances I cannot write as stated

**5. Step 3 at two geopotential spacings.** Every run from 2 to 9 is to be made at 5e4 and 2.5e4.
The spacing is `[grid] geopotential_spacing_m2s2` in the run's committed namelist, and
`read_run_namelist` requires the namelist to be named for its run and its output directory to be
inside the run directory. So a second spacing cannot be run in place without editing a committed
file. Proposed: the run of record is the run's own directory at 5e4, built and run as the runbook
says; `run_experiments.py` makes the 2.5e4 run in a copy under `reports/step05_3/spacing_2p5e4/`,
with the anchor path rewritten one level deeper, which is how `step04_5` check 11 reruns the
closure namelist. Nothing committed is edited.

**6. Run 7f is not "only its `[shear]` section" different.** §4's opening sentence says the runs
differ only in `[shear]` and, for 3b, the target; run 7f differs from run 7 in `[wind]
pressure_grid_Pa` (the table and the paragraph after it say so). Proposed: add 7f to the
exception in the opening sentence.

**7. Step 2 check 5 and Step 3 check 1 compare with a product that is not registered.** Covered by
ruling 1. I note only the order it implies: the transfer product has to be rebuilt on a clean tree
and committed before Step 1, which needs `git status --porcelain` empty at the time, including
untracked files such as this report.

**8. Step 3 check 3 needs the name of `r0` in the product.** "Run 3a's `r0(10 N)`" is read from the
transfer product. I will find the variable before Step 3 and cite it; I record it here only so
that the check's wording can name it at v0.4.

## 3. Choices I think the author should make

**9. `u_s` from the model's own interpolant, or from the tool's own part.** §1.3 asks for `u_s`
linear in ln p "(the model's own vertical rule)", and Step 1 deliverable 1 lists the ln p
interpolation as a part tested on its own. `lib.windfield.WindField` already is that rule, and
refuses outside the grid. Calling it guarantees that the tool's `u_s` is the value the model would
read at `p_s`; a separate part could differ from it by rounding. Proposed: the part is a thin call
to `WindField` at the latitude nodes, tested on its own against the closed form. Step 1 check 2
needs `u_total` array-equal at a node, so the check will show whether the interpolant returns a
node value exactly.

**10. Values the specification does not list as malformed.** §1.4 lists the refusals. Not listed:
`p_s` zero or negative (it cannot be inside the grid, so the grid refusal catches it), `p_stop`
negative, `stop_fraction` or `scale` not finite, and `shape` not one of the two names. Decision N
refuses a malformed value, so I propose to refuse each of these with a message naming the key, and
to list them in Step 1 check 10. Please confirm, since check 10 says "and only those".

**11. The role of a `[shear]` section.** §1.2 says a reduction's wind is made without a shear case.
The tool reads `role` with `build_role`, which accepts every build role. Proposed: accept only
`forward`, since the specification gives no use for any other, and refuse the rest as a
malformed value. The alternative is to accept any role and leave the rule to the driver, which
already requires `forward`.

**12. What `casspian-new-run` renames.** Proposed rule: every occurrence of the source run's name
as a whole token in both files, comments included, is replaced by the new name; nothing else is
touched. That covers the prefixes, output names, `[inputs]` paths, `[output] product`, `[run]
name`, titles and comments, and makes check 6 a plain statement: the files differ from the source
only in those substitutions. It leaves `[run] description` as it is when it does not name the run,
so an experiment's description is the person's to edit.

**13. The output's `title` when the control file gives none.** §1.7 says every other global is
carried unchanged; `title` is an optional key. Proposed: without a `title` key, the input's title is
carried, with `, shear case <name>` appended, so the file names its case without being opened.

**14. The uncertainty's `long_name` without `uncertainty_ms`.** §1.7 says the field passes through
unchanged, its `long_name` stating that it is the source's uncertainty. That changes one
attribute of a field otherwise untouched. Proposed: yes, change it, for every case except
`identity`; check 9a compares values only, so it is unaffected. Please confirm.

**15. F9's second line.** Step 2 deliverable 5 draws `u_total at <p> mbar` beside `source wind,
assigned to <p> mbar`. `p_s` is not written to the file (§1.2), so the one pressure F9 can know is
the source's `reference_level_pressure_Pa`. Proposed: both lines at that pressure, `u_total`
interpolated there by the model's rule. Please confirm that this, and not `p_s`, is the pressure
meant.

## 4. Two things for the runbook

- Step 1 and Step 2 each add a console entry to `pyproject.toml`. A clone that pulls them has to
  run `pip install -e . --no-deps` again before `casspian-wind-shear` or `casspian-new-run` exists.
  The runbook v0.4 of Step 2 deliverable 6 is the place to say so.
- The acceptance input of Step 1, `forward/lindal_transfer/inputs/lindal_transfer_wind.nc`, is a
  committed file since `4d32efe`. Ruling 3's line in Step 1 covers it.
