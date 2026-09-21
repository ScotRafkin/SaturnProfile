# CASSPIAN coding agent handoff, 17 September 2026

For a coding agent starting fresh in the repository. This document supersedes
`docs/CASSPIAN_CodingAgent_Handoff_2026-09-14.md`, whose state section predates SPEC_03; what
still holds from it is carried here, so nothing else has to be read first. Where this document
and the repository disagree, the repository is right, and the disagreement is worth a sentence
in your first report.

---

## 1. What this is, and who does what

CASSPIAN is a one-dimensional Saturn atmospheric structure model for probe mission design: a
reference profile of temperature, pressure and density against altitude at a chosen latitude
and season, with quantified uncertainty. It is built from radio occultation profiles reduced
to refractivity (the transportable quantity), a gravity field, a rotation rate, a wind
hypothesis and a composition, and it moves refractivity between latitudes along isobars under
gradient wind balance. The physical specification is the author's manuscript,
`CASSPIAN_AtmosphericModel_Draft9_2.docx` (the draft of record for the labels), whose
equations the specifications cite by label (A1 to A40, B1 to B8). The manuscript and its notes
are the author's and stay outside the repository while the draft is in flux; the author supplies
the file. SPEC_04 cites A13
to A34, A39 and B3 to B7, and does not restate them.

Three parties. **The author**, S. Rafkin, makes every decision of substance and accepts every
step. **The reviewing agent** writes the specifications (`docs/specs/SPEC_*.md`), reviews every
report against its specification (`reports/REVIEW_0N_stepM.md`), rules on questions, and
measures expected values with code independent of the repository. **The coding agent**, you,
codes one step at a time against an accepted specification, runs the acceptance, writes the
report, and commits only after the review and the author's go. You never decide physics: a
choice the specification did not make is a decision stated in the report or a question to the
author, never a default in the code.

The work is Southwest Research Institute funded and proprietary. Implementation details do not
go to external collaborators.

---

## 2. Exact state at handoff

- **Repository:** `C:\Users\srafkin\SaturnProfile`, branch `main`, remote
  `https://github.com/ScotRafkin/SaturnProfile.git`.
- **Last accepted step record:** `545b385`, the SPEC_03 Step 4 record. `HEAD` is `5f23a8b`:
  the three commits after the record (`2cd41bd` the SPEC_03 closure, `e2e7eb4` the STATE
  update, `5f23a8b` the REVIEW_03_step4 addendum) are the author's documents, none a step.
  Confirm with `git log` and `git status` before anything else; the author's uncommitted
  documents (`docs/specs/SPEC_04_Transfer.md`, this file, and possibly a newer
  `docs/CASSPIAN_Seasonal_Design_Note.md`) may be untracked or modified. They are committed in
  their own commit, "the author's documents", when the author says so, never mixed into a step.
- **Specifications:** SPEC_00 v0.20 (architecture and data files, accepted); SPEC_01 v0.28
  (the Lindal tool chain, closed); SPEC_02 v0.11 (the reduction and diagnostics, closed);
  SPEC_03 v0.14 (the forward production and its hydrostatic closure, closed at `545b385`);
  **SPEC_04 v0.7 (the transfer), accepted by the author at v0.6 on 21 September 2026: Step 0
  proceeds once `reports/REPORT_04_preexecution.md` is filed.**
- **`docs/specs/STATE.md`:** every step of SPEC_01, SPEC_02 and SPEC_03 is `accepted`; the
  SPEC_04 section lists Steps 0 to 5 as `not started`.
- **Products on disk, all clean at `545b385` or earlier:** the reduction chain under
  `occul_data/lindal/` (raw bundle, kinds T, D, G, R, W, C and the kind N product
  `lindal_refractivity.nc`, SHA-256 `d2886aed...`, every file dated `epoch = "1981-08-26"`),
  and the closure run under `forward/lindal_closure/` (namelist, build control file, `inputs/`
  by `casspian-run-inputs`, `output/lindal_closure_profile.nc`, SHA-256 `e6693173...`, commit
  `1504056`). Hashes are recorded in REPORT_03_step3 section 7 and REPORT_03_step4 section 8.
  Git holds no netCDF product (see `.gitignore`); the reports hold the hashes.
- **Nothing is in progress.** The first act is `reports/REPORT_04_preexecution.md` (section
  4); Step 0 follows. A chat message saying "ready" never accepts a draft; the status line and
  `STATE.md` do.

---

## 3. Reading order

1. `docs/specs/SPEC_00_Architecture_and_Data_Files.md`: the architecture, the file kinds, the
   run namelist (§7.2), the sweep (§8).
2. `docs/specs/SPEC_03_Forward_Production.md` §0 (the step protocol, adopted by SPEC_04), §7
   and §8 (decisions and the rulings on the last pre-execution review, which show what a
   useful finding looks like), and Step 4 (the production interface SPEC_04 calls).
3. `docs/specs/STATE.md`.
4. `reports/REPORT_03_step4.md` and `reports/REVIEW_03_step4.md`: the last report and review,
   the model for yours.
5. `docs/specs/SPEC_04_Transfer.md` v0.6, all of it, with the manuscript's Appendix A and B
   beside it.
6. `docs/CASSPIAN_Seasonal_Design_Note.md` v0.5: the philosophy behind SPEC_04's decisions on
   the wind, the anchors and the seasons. Not a specification; nothing in it is built until a
   specification says so.
7. SPEC_01 and SPEC_02 as needed for the tools and the reduction you will touch at Step 0.

---

## 4. The task now: the record of the pre-execution review, then Step 0

The pre-execution review of SPEC_04 v0.5 was made on 17 September and ruled on in SPEC_04 §10
(v0.6); the author accepted the specification at v0.6 on 21 September. The first act is its
record: `reports/REPORT_04_preexecution.md`, the findings as sent and, against each, the ruling
of §10 marked accepted or contested with the reason. Then Step 0 by the procedure of section 6.
The standard that review set (a self-contradiction found by measurement, a wrong count, a
refusal case that would fail for the wrong reason, each with the rule proposed) is the standard
for every report.

Things to look at hard: the three-part kind W and what `casspian-wind-from-curve` has to
change; the loader's anchor object and the propagation hook (Step 0 deliverables 6 and 7);
the outer loop of the isobar map; the M = 2 identity test and what the synthetic anchor needs
to pass `lib.io.read` as kind N; the estimate's union levels and interpolation; the
`transfer_record` and `estimate` groups against `lib.schema`; the expected values, which are
instances for the Lindal anchor and are to be reproduced, not tuned to.

---

## 5. Environment

- **Windows 11.** Python 3.14 (`pythoncore-3.14-64`), numpy, scipy, xarray, netCDF4,
  matplotlib 3.10.
- **Install:** the package is installed editable. After adding an entry point to
  `pyproject.toml`, reinstall so the console script exists:

  ```
  pip install -q --no-deps --no-build-isolation -e .
  ```

- **Entry points:** `casspian-lindal-raw`, `casspian-gravity-file`, `casspian-rotation-file`,
  `casspian-wind-from-curve`, `casspian-composition-lindal`, `casspian-lindal-inputs` (the
  SPEC_01 tools); `casspian-refrac <manifest>` (the reduction, rendering figures when the
  manifest asks); `casspian-plots <file.nc> [--out DIR] [--format png|pdf] [--dpi N]` (the
  diagnostics); `casspian-run-inputs <build control file>` (a forward run's inputs);
  `casspian-forward <namelist>` (the forward run; `mode = "closure"` today, `"transfer"` after
  SPEC_04).
- **Shell:** Git Bash for POSIX commands; PowerShell also available. Run Python from the
  repository root; the acceptance scripts use relative paths.
- **Unicode:** set `PYTHONIOENCODING=utf-8` whenever a script prints a DataTree or any
  non-ASCII text; the console codec (cp1252) otherwise raises `UnicodeEncodeError`.

---

## 6. The procedure for every step

1. **Gate.** Read the status line at the top of the specification and `STATE.md`. A draft is
   not started.
2. **Build only the step.** Acceptance script under `reports/step04_<N>/`, for example
   `reports/step04_0/accept_step04_0.py` (ignored by git). Do not implement ahead.
3. **Every check prints its measured value.** A check the specification did not ask for is
   labeled "beyond the specification". Output to `output.txt` beside the script. Run the checks
   as decided before the run; when one fails, report the failure and the physics, and never
   loosen a check to pass it or re-read the words afterward to make it pass (REPORT_03_step4
   finding 1 is the model: the failure was the specification's wording, and it was reported as
   a failure and ruled on).
4. **A step that changes an input file** (SPEC_03 §0, author ruling of 14 September):
   `refrac` and `forward` refuse `-dirty` inputs, so such a step is accepted on a candidate
   product built in memory by the acceptance script, with the refusal relaxed inside that
   script only and the relaxation named in the output; the registered product is rebuilt on
   the clean tree at the sweep, and the suites that run the command-line tools on the on-disk
   products run after the sweep and are recorded in the report then.
5. **Write `reports/REPORT_04_step<N>.md`:** what was built; decisions (every choice the
   specification did not make); the acceptance table; findings (wherever a specification number
   or expectation is wrong, with the measured value and the reason); regression; hashes; next
   step. State failures plainly.
6. **Set the step to `reported` in `STATE.md`. Do not commit. Stop.** Tell the author the
   report is ready.
7. **After the review** (`reports/REVIEW_04_step<N>.md`), follow its order of work exactly:
   make the changes; rerun; commit the author's documents in their own commit; the acceptance
   commit for the step; push; the sweep (section 8) over `occul_data/` and `forward/`; set
   `STATE.md` to accepted with the commit; the record commit when the review asks for one.
8. **Never overwrite or delete an author file.** Look first.

**Standing rules for everything written:**

- no em dashes or en dashes anywhere, in code, comments, documents, figures or commit
  messages (a minus sign in a number is fine); check every file you write;
- American spellings;
- **no AI attribution in any commit message, document or code.** The author's rule overrides
  any tool reminder to add a `Co-Authored-By` or similar line; no commit in this repository
  carries one;
- no silent choices;
- do not explain elementary physics in reports; state what was measured;
- do not refer to text by paragraph or line number; cite sections, deliverables, equation
  labels and check numbers;
- every angle inside `lib` is in radians; degrees are converted at the boundaries;
- a reasonably documented, serviceable team tool is the goal. Do not get carried away with
  abstraction, configuration or defensive code the specification did not ask for.

---

## 7. The regression

Run every earlier suite after any change to `lib`, `refrac`, `forward` or `tools`, and record
the result lines in the step's `regression.txt` and the report.
`reports/step03_4/run_regression.sh` is the current template: it copies the swept product, its
figures and the closure run aside, runs every SPEC_01 and SPEC_02 suite and the SPEC_03 Step 1
to 3 acceptances, restores everything byte for byte, restores `reports/figures/` from git, and
prints the porcelain state. Add `step03_4/accept_step03_4` to its list, and each SPEC_04 suite
as it is accepted.

**Last full result (at `1504056`, REPORT_03_step4 section 5):** 6, 14, 8, 5, 7, 7, 6, 8, 9, 6
for SPEC_01; 6, 7, 9, 9, 7, 7 for SPEC_02; 9, 8, 16 for SPEC_03 Steps 1 to 3; all pass, the
product restored with its hash unchanged.

**Timing.** Most suites take seconds. `step7` takes 3 to 4 minutes (wind geoid convergence),
`step02_3` about 3.5, `step02_4` about 4, `step02_5` about 3, `step02_6` about 12 (eight extra
fixed points). The loop takes about 20 minutes; run it in the background and wait for it.
Never run two suites that write the product at once: `step02_4`, `step02_5` and `step02_6` all
rewrite `lindal_refractivity.nc` and `figures/`.

**Dependencies to respect:** `step02_1` checks input and manifest hashes against report tables
in the row format below, later rows overriding earlier; a new hash must be recorded in that
format or the script's report list extended. `step02_2` and `step02_3` pin the
`mean_polar_radius` anchor rule in memory. `step9` checks the directory listing.

```
| `file` | kind | `64 hex digits` |
```

**Side effects.** `step02_4` to `step02_6` rewrite the registered product from the working
tree, which then carries `<HEAD>-dirty`; `step02_5` overwrites the committed figures attached to
REPORT_02_step5. The template script restores both; if you run a suite by hand, restore by hand
(`git checkout -- reports/figures/`) before any commit.

---

## 8. The sweep, exactly

`lib.io.git_commit` marks a file `-dirty` when `git status --porcelain` shows anything,
untracked files included. A clean build needs a completely clean tree.

1. Commit and push everything for the step.
2. Make the tree clean: move or delete untracked scratch; `git stash push -u -- <paths>` for
   uncommitted later work, popped after the build. Confirm the porcelain output is empty.
3. Rebuild each product that carries `-dirty` with its own tool: the chain with the SPEC_01
   tools and `casspian-refrac occul_data/lindal/lindal_reduction.toml`; a run's inputs with
   `casspian-run-inputs` and its product with `casspian-forward`.
4. Read back `casspian_git_commit` from every rebuilt file and refuse if it ends in `-dirty`.
5. Compute each SHA-256 with `casspian.lib.io.sha256` and record it in the step report's hash
   table in the row format of section 7.
6. Commit that record ("Record the Step N sweep: ...") and push. Figures a report cites are
   rendered from the clean product, since the footer shows the commit.

---

## 9. Pitfalls that have cost time here

1. **Bash heredocs embedding Python with backslashes or triple quotes break**, and a broken
   chain once reported nothing while nothing ran. Write patch scripts to a file with the Write
   tool and run them, or use the Edit tool on files already read; build awkward characters
   with `chr()`; assert each replacement matched exactly once; `python -m py_compile` after
   every patch.
2. **The Edit tool refuses a file not read with the Read tool in the session.** `cat` does not
   count.
3. **`set -e` is not reliable in a long chained command.** Use explicit `|| exit 1` after every
   step that matters and read results back.
4. **Anything claimed in a report must be measured.** Two unmeasured numbers have reached
   reports before. Check every derived figure against the output before writing it.
5. **`STATE.md` rows:** several specifications have a Step 0 and a Step 4; match on the
   deliverable text, never on the step number.
6. **Figures need looking at.** Acceptance checks do not see overlapping labels or a curve
   lying on zero. Render, view the PNG, fix, then run the final acceptance. Attach figures
   from a clean product.
7. **xarray DataTree:** set a root dataset with `tree.dataset = ...`; `tree["/"] = ...` raises;
   `node.to_dataset(inherit=False)` for a group's own variables.
8. **The kind N reader validates hashes.** `lib.io.read(path, "refractivity")` refuses a file
   whose `input_hashes` no longer lists a hash recorded in `reduction_record`, and every
   derived kind warns when an input on disk differs. Copies made for refusal tests, and the
   synthetic anchor SPEC_04 Step 4 writes, must be built with that in mind.
9. **Anchor rules and quantities are paired** in `casspian.lib.control.ANCHOR_QUANTITY_FOR_RULE`;
   `lib.control` and `casspian-lindal-inputs` both refuse a mismatch.
10. **Long jobs** (the anchored fixed point takes 30 to 65 s; the wind geoid march and the
    outer loop of SPEC_04 will be comparable) belong in the background, and the report waits
    for the result.
11. **When the specification describes one thing two ways, raise it; do not choose.** The "2
    mbar edge" of REPORT_03_step4 came from a bin described as "above 2 mbar" in one place and
    "the ten top levels" in another; the right response was the finding, which it was.
12. **Read back what you write.** Files written to disk have silently failed to land more than
    once in this project. After writing a file that matters, read it back and compare.

---

## 10. Code map at `545b385`

| Module | What it holds |
|---|---|
| `lib/constants.py` | CODATA 2018 |
| `lib/schema.py` | kind registry (raw, thermo, geodesy, gravity, rotation, wind, composition, refractivity, profile), `validate`, uncertainty companions, writer globals, the season identifier rules |
| `lib/io.py` | `write`, `read` (a DataTree for the grouped kinds, a Dataset otherwise; the kind N hash check and the `input_hashes` warning for every derived kind), `sha256`, `input_hashes`, `git_commit`, `package_version` |
| `lib/control.py` | `load_section`, `read_reduction_manifest`, `load_reduction_inputs`, `ANCHOR_QUANTITY_FOR_RULE`, `read_run_namelist`, `load_run_inputs` (closure comparison, refusals) |
| `lib/gravity.py` | `g_newton`, `G_phi_newton`, `omega_abs`, `g_eff_radial`, `G_phi_eff`, `g_eff_vector` (returns `psi`) |
| `lib/geoid.py` | `U_rigid`, `reference_geoid`, `wind_geoid` (rules `mean_polar_radius`, `north_pole`, `south_pole`, `latitude`, `equatorial_radius`), `radius_at`; SPEC_04 adds `through_anchor` |
| `lib/latitude.py` | `planetocentric_fixed_point`, `planetographic_from_planetocentric` |
| `lib/reduction.py` | `number_density`, `absolute_radius`, `refractivity`, `temperature_from_refractivity`, `quadrature` |
| `lib/geopotential.py` | the field-line integral of the effective gravity along a column (SPEC_03 Step 1) |
| `lib/hydrostatic.py` | the exact log-linear layer integral from the top pressure down (SPEC_03 Step 2) |
| `refrac/anchor.py`, `refrac/reduce.py`, `refrac/product.py` | the reduction to kind N; `casspian-refrac` |
| `forward/production.py` | `produce` on a `Profile` and the loaded inputs; the closure run, `closure_statistics`, `production_record`; `casspian-forward` |
| `tools/run/run_inputs.py` | `casspian-run-inputs`: a run's four inputs from its build control file |
| `tools/wind/`, `tools/composition/`, `tools/gravity/`, `tools/lindal/` | the SPEC_01 tools (kinds W, C, G, R, raw, T, D) |
| `tools/plots/` | `render`, `style`, `figures_product` (F1 to F6 for kind N), `figures_inputs`, `figures_profile` (F5, F6 for kind profile) |

SPEC_04 adds, by its steps: `lib.geoid.through_anchor` and `lib.windfield` (Step 1),
`lib.mesh` (Step 2), `lib.kernel` (Step 3), `forward.transfer`, `forward.estimate` and the
`forward.propagate` hook (Step 4; the hook's object is defined at Step 0), transfer mode in
`forward.production` and `casspian-forward` with F5's across-latitude panel, F6 in transfer
mode and F7 (Step 5); and at Step 0 the three-part kind W, kind C on a latitude grid, the
transfer namelist and the loader's anchor object.
