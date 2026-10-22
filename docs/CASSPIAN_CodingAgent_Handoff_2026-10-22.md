# CASSPIAN coding agent handoff, 22 October 2026

For the incoming coding agent. Author of record: S. Rafkin (SwRI), who approves everything. This
follows `docs/CASSPIAN_CodingAgent_Handoff_2026-10-03.md`; its sections 4 and 5 (how the work is
judged, and practice) still hold, and are restated in short in section 6 here.

## 1. Before you start

- **Wait for the author.** Do not start work without a specification whose status line says it is
  accepted. A "ready" in chat does not accept a draft; the specification's status line and
  `docs/specs/STATE.md` gate the work.
- **Start from a clean `main`** at `7f03d9c` or later, pushed to `origin`.
- **The next work** is a small fix step for the three open items of section 4, with no new
  physics, ahead of SPEC_14. SPEC_14 then starts from the run `lindal_iris_v_sl_ii_110`
  (REVIEW_13_step3).

## 2. Where things stand

- **Specifications.** SPEC_00 to SPEC_13 are closed. SPEC_00 is at v0.24.
  - **SPEC_08:** kind W holds the three parts of a wind (`u_total_ms`, `u_reference_ms`,
    `u_shear_ms`) with the reference level required. The run's loader admits a wind by the presence
    of `u_total_ms` and the polar rule.
  - **SPEC_09:** the Voyager IRIS temperatures at 110, 290 and 730 mbar, digitized from Conrath
    and Pirraglia (1983) Fig. 1 (`occul_data/lindal/iris_temperatures.csv`).
  - **SPEC_10:** `iris_temperatures.fit`, a local linear fit with a Gaussian window that widens to
    3 effective points where the data thin.
  - **SPEC_11:** the Lindal-only wind (`tools/lindal/lindal_wind.py`): the cloud wind with the
    shear the IRIS gradients imply by the model's own balance. It is built by the case
    `lindal_iris` of `casspian-wind-shear`, and run as `lindal_iris_ii_110`.
  - **SPEC_12:** the Voyager 2 egress profile, Lindal et al. (1985) Fig. 4, digitized and converted
    to pressure (`tools/lindal/egress_profile.py`, `occul_data/lindal/voyager2_egress.csv`).
    `lindal_iris_ii_110` is targeted at the egress label, 31.2 S planetographic
    (-26.06095353762054 deg planetocentric).
  - **SPEC_13:**
    - the case `lindal_iris` takes `gradient` (`"slope"` or `"value"`) and `window`
      (`"adaptive"` or `"smooth"`);
    - the wind tool's curve path takes a table (`curve_format = "table"`,
      `gap_rule = "pchip_bridge"`, `uncertainty_source = "table"`, `ring_gap_deg`) for the
      Sanchez-Lavega et al. (2000) cloud wind;
    - the data property files carry only the values code reads, and kind W no longer requires
      `method`, `observation_level_justification` or `solar_longitude_source` (SPEC_00 v0.24);
    - two runs, `lindal_iris_v_ii_110` and `lindal_iris_v_sl_ii_110`.
- **The runs of the Lindal line.** All three target 31.2 S. Their inputs and products are not
  committed; each acceptance rebuilds them.

  | Run | Gradient | Window | Cloud wind | Minus the egress at 110, 290, 730 mbar |
  |---|---|---|---|---|
  | `lindal_iris_ii_110` | slope | adaptive | Ingersoll and Pollard | +3.91, +2.37, -0.95 K |
  | `lindal_iris_v_ii_110` | value | smooth | Ingersoll and Pollard | +4.95, +2.25, +0.26 K |
  | `lindal_iris_v_sl_ii_110` | value | smooth | Sanchez-Lavega | +5.00, +2.32, +0.40 K |

  The egress's sigma_T there is 2.27, 2.29 and 2.73 K. At 110 mbar the radio occultations and the
  IRIS retrievals disagree on the change from 36.3 N to 31.2 S: the radio shows +0.03 K, and the
  IRIS fit +4.8 to +5.2 K, which the runs carry (REPORT_12_step1, REPORT_13_step3).
- **Registered products** are unchanged since `7652852` and are committed. Their embedded copies
  of the inputs predate SPEC_00 v0.24 (section 4, item 1).
- **Unchanged from the 10-03 handoff:** decision L (the linear wind interpolant), the branch
  `spec07-l2-parked`, and the tabled redesign of the vertical.

## 3. The regression at `cf33cd6`

Run on 8 October 2026 at `7f03d9c`, which is `cf33cd6` plus the STATE.md hash. The tree was clean:
the five untracked reference documents in `docs/` were listed for the run in `.git/info/exclude`
and the listing was removed afterward. An untracked file makes `lib.io` stamp every product
`-dirty`, and the readers then refuse it.

- **`tests/run_regression.sh`, 31 suites: 29 at their reference counts**, 7621 s of suite time.
  The driver was run in four chunks (`reports/regression/regression_chunk{A,B,C1,C2}.txt`) to stay
  inside the background-job limit; each restored the registered kind N and the closure and
  transfer files, and each reports them equal.
  - `step9`: 5 of 6 (section 4, item 2).
  - `step05_3`: stopped at check 3, after 2 of 3 passed (section 4, item 1).
- **The seven suites the driver does not list, run by hand in order: all at their counts**, 2956
  s (`reports/regression/extra/summary.txt`).

  | Suite | Checks | Time |
  |---|---|---|
  | `step09_1` | 3 of 3 | 6 s |
  | `step10_1` | 2 of 2 | 3 s |
  | `step11_1` | 5 of 5 | 1607 s |
  | `step11_2` | 2 of 2 | 200 s |
  | `step12_1` | 4 of 4 | 201 s |
  | `step13_1` | 3 of 3 | 441 s |
  | `step13_2` | 2 of 2 | 498 s |

  `step11_1`, edited at SPEC_11 Step 2 and not run since, passes 5 of 5. The reruns reproduce
  their accepted numbers: for example, the round trips at 110 mbar are +3.67, +4.71 and +4.76 K,
  and the comparison with the egress is as in the table above.
- **No value changed.** Every check that compares rebuilt values passed. `step03_3` check 15
  reports that the rebuild changed no value in any product and that kind N is unchanged bit for
  bit. The only difference found is in attributes (item 1).
- **Per-suite logs** are in `reports/regression/regression_accept_*.txt` and
  `reports/regression/extra/`.

## 4. Open items, for the fix step

### Item 1. The closure run is refused (`step05_3`)

- **What happens.** In closure mode `lib.control.load_run_inputs` calls `check_closure_inputs`. It
  compares each run input with the anchor's embedded copy, group by group, with
  `xarray.Dataset.identical` after `_strip` removes `CLOSURE_DROPPED_ATTRIBUTES` (control.py
  line 660) and any attribute holding a `sha256:` hash, on the group and its variables.
- **Why it fails now.** The registered `occul_data/lindal/lindal_refractivity.nc` embeds a wind
  built before SPEC_00 v0.24. A wind rebuilt today lacks five attributes: four on the file
  (`epoch_note`, `observation_level_justification`, `solar_longitude_source`,
  `source_latitude_convention_source`) and `value_source` on the variable
  `reference_level_pressure_Pa`. Every value is equal. The refusal reads: "the run's wind input
  is not content-identical to the anchor's embedded copy, so this is not a closure". It raises in
  `step05_3` at line 158, `forward_production.run(CLOSURE / "lindal_closure.toml")`.
- **What it reaches.** Any closure run built from freshly rebuilt inputs: `casspian-run-inputs`
  followed by `casspian-forward` in closure mode. Transfer runs make no closure comparison and are
  unaffected.
  - `step05_3` checks 4 to 7 did not run. They compare a rebuilt closure and transfer product with
    the registered ones, and run `casspian-new-run` and F9's figure. Run them after the fix.
  - `step04_5` passes because it runs the closure from copies of the committed inputs, which
    still carry the attributes.
  - The other committed winds also carry them: `occul_data/lindal/lindal_wind.nc` (the
    reduction's input) and `forward/lindal_closure/inputs/` and `forward/lindal_transfer/inputs/`.
- **The reviewing agent's recommendation.** Compare values only and ignore descriptive
  attributes, as v0.24 intends; do not rebuild the registered product. The narrowest form is to
  add the five names to `CLOSURE_DROPPED_ATTRIBUTES`. `value_source` is a variable attribute and
  `_strip` already drops listed names from variables. The specification should say which form is
  wanted:
  - the five names;
  - every descriptive attribute;
  - values only.

  The first touches nothing else; the last changes what closure means for every kind.

### Item 2. The directory listing (`step9` check 6)

`occul_data/lindal/` holds `iris_temperatures.csv` (SPEC_09) and `voyager2_egress.csv` (SPEC_12),
which SPEC_00 §2.2's listing does not name. The check reports them as unexpected. Nothing else
differs, and the check has failed since SPEC_09, not since SPEC_13. Adding the two files to SPEC_00
§2.2 settles it. The suite reads its expected list from the specification's listing, so check
whether the suite or the specification is the place to change.

### Item 3. The driver does not list the newer suites

`tests/run_regression.sh` ends at `step08_1`. Add the seven suites of section 3 with their
reference counts. The counts are 3, 2, 5, 2, 4, 3 and 2.

Notes for the edit:
- **Order.** `step11_2`, `step12_1` and `step13_2` rebuild and rerun the three Lindal-line runs,
  and `step13_1` check 0 compares against `lindal_iris_ii_110`'s inputs as found. Keep them in
  numeric order.
- **What the driver restores.** It does not save or restore `forward/lindal_iris_*`, which are
  uncommitted and rebuilt by those suites.
- **Run time.** The newer suites add about 50 minutes, bringing the full set to about 2 h 56 min of
  suite time. That is beyond a single background job's limit in this environment: run it in
  chunks, as here, or detached.
- **Correction to the reviewing agent's note.** The newer suites were run, by hand, and all
  passed. What the item needs is the driver's list, not a first run.

## 5. Things the next agent should know

- **The clean tree.** The five untracked documents in `docs/` (the 10-03 handoff, the
  Sanchez-Lavega HTML and PDFs, and the Conrath et al. IRIS PDF) are not to be committed unless
  the author asks. While they are untracked, every product built is `-dirty`. List them in
  `.git/info/exclude` for a regression, or have the author commit or move them.
- **Files locked by viewers.** On this machine an image open in the Photos app, or a PDF open in
  a viewer, can make a write or a move fail with "Invalid argument" or "Device or resource busy".
  Ask the author to close the viewer; do not work around it.
- **Shell heredocs** mangle backslash escapes in inline Python. Write edit scripts to files in the
  scratchpad.
- **Latitudes.** In figures and outputs, give every latitude planetographic, with planetocentric
  in parentheses (REVIEW_12_step1). The wind file's own planetographic conversion (no-wind geoid)
  puts the target at 31.11 S, and the anchor's fixed point at 31.2 S. The round trips use the
  former, as SPEC_11 Step 2 did.
- **The round-trip diagnosis** (`tests/step13_2/accept_step13_2.py`) splits each difference into
  (A) to (D) as REPORT_11_step2 §3 did. It integrates to the exact end latitudes. The SPEC_11
  version used grid nodes only; SPEC_11 v0.5 and REVIEW_11_step2 record the correction.
- **The smooth window** is a running maximum of the adaptive width over ±4 deg, averaged with a
  Hann kernel of the same half-width. It is never narrower than the adaptive width by
  construction. It removes the derivative's spikes, but it also changes the fitted temperatures
  in the IRIS gap by up to 1.3 K (REPORT_13_step3).

## 6. How the work is judged, and practice (from the 10-03 handoff, in short)

- **General code.** The code must work on any profile, real or invented. An unexplained behavior
  is explained with a measurement before acceptance.
- **Bounds.** Tests are general, with physical bounds and none set for one case. Report the
  pressure identity on every run.
- **A missed bound.** Stop and file an interim report with every other check measured. Do not tune
  numerics, spacings or tolerances to pass.
- **Fixes.** One fix per specification, each revertible by git, with no switch between old and new
  paths. Replace, never add. Refuse only typos, extrapolation and named numerical failures.
- **Regression.** A full set only when a computed value can change or at a named milestone;
  otherwise the step's own suites. Name what was rerun and why.
- **Accepted suites.** A suite whose expectation moves is updated in the same step, with its old
  and new value in the report.
- **Dirty products.** Products from a dirty tree end `-dirty` and are refused. Acceptance scripts
  relax the refusal for their own process only.
- **Files and scripts.** Close every netCDF read before opening another. Do not edit a script
  while a shell is running it.
- **Commits.** Commit only after a step is reviewed and accepted, then push `main` without asking.
  The acceptance commit holds the code, suites, reports, reviews and specifications the review
  lists. A second commit adds that hash to the step rows of `STATE.md`. The step rows are the
  coding agent's; the header lines are the reviewing agent's.
- **Style:**
  - no em or en dashes in prose;
  - American spelling;
  - every measured number with its bound;
  - no attribution of any kind to an AI in commits or code;
  - the work is SwRI funded and proprietary, and nothing leaves the repository.
