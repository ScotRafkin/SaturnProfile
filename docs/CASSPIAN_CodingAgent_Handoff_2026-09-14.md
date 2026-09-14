# CASSPIAN coding agent handoff, 14 September 2026

The operational companion to `docs/CASSPIAN_Handoff_2026-09-13.md`. **Read that document first.**
It covers what CASSPIAN is, the roles, the standing rules, the physics as implemented, the
Lindal inputs and decisions, the kind N product, and what comes next. This one covers what a
fresh coding agent needs in order to work in the repository without relearning it: the exact
state, the procedures, the regression suites, and the pitfalls that cost time.

Where this document and the repository disagree, the repository is right.

---

## 1. Exact state at handoff

- **Repository:** `C:\Users\srafkin\SaturnProfile`, branch `main`, remote
  `https://github.com/ScotRafkin/SaturnProfile.git`.
- **Commits:** `HEAD` = `origin/main` = `36d5e96`, "Record the Step 6 sweep". No stash. The
  working tree is clean except the untracked `docs/CASSPIAN_Handoff_2026-09-13.md` (the author's
  handoff) and this file.
- **Specifications:**
  - SPEC_00 v0.15, accepted architecture; section 7.3 is open for SPEC_03.
  - SPEC_01 v0.20, closed.
  - SPEC_02 v0.9, closed at `2cf4457` and swept at `36d5e96`.
  - **SPEC_03 does not exist yet.** The author is writing it.
- **`docs/specs/STATE.md`:** every SPEC_01 and SPEC_02 step is `accepted`.
- **Nothing is in progress. Do not build anything** until SPEC_03 exists, its status line says a
  step proceeds, and `STATE.md` agrees.

**Products under `occul_data/lindal/`**, all clean (no `-dirty`):

| File | Kind | Commit | SHA-256 |
|---|---|---|---|
| `raw/lindal_raw.nc` | raw | `ece58d2` | `3054cc6cf0b54905ee3d7d67bbe28ad12bac18288e1fcd53893957d54d9171eb` |
| `lindal_thermo.nc` | thermo | `ece58d2` | `ff9c12991399f94c08caaae80b1a6a91f45a45feaebc4d059ad563b45a720e85` |
| `lindal_geodesy.nc` | geodesy | `ece58d2` | `16e0c7d551c3026ce5623141b1a2f86187c013a2dbf929f3e330f79f410f787e` |
| `lindal_gravity.nc` | gravity | `ece58d2` | `6db0129cefd01bc9c151020e604cf4563a82381e61c97381d0a7ee383826911e` |
| `lindal_rotation.nc` | rotation | `ece58d2` | `a5c72017a423c50f2a952ddcbc4ccf777c1e94683fa6965ee96f7f43d949a424` |
| `lindal_wind.nc` | wind | `ece58d2` | `18d8cdae4dd8faa5cdcea40d25e6a9dfce4391421b3af4205d2b377f1150e7b9` |
| `lindal_composition.nc` | composition | `6ee8690` | `a1250b60c4b68c42109bb1fbe22244b218af9866aa1e3007f442a9d2ba26654f` |
| `lindal_reduction.toml` | manifest | tracked | `97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee` |
| `lindal_refractivity.nc` | refractivity | `2cf4457` | `888517ccdea0688ea9059f434f07b1384ee117eb57d3c58a5f45ea2094f82a12` |

The manifest anchors on `equatorial_radius` / `radius_equatorial_m` and carries
`[diagnostics] figures = true, format = "png", dpi = 150`. Each hash is recorded in the report
that last produced the file:
- REPORT_01_step9 section 6 for most;
- REPORT_01_step8 section 7 for the composition;
- REPORT_02_step6 section 7 for the manifest and the product.

**Git holds none of the netCDF products.** `.gitignore` excludes `occul_data/*/*.nc`,
`occul_data/*/raw/*.nc` and `occul_data/*/figures/`. Products live on disk only; the reports
record their hashes. A sweep therefore commits a hash record, never a product. The ignore file's
own comment says those rules go "once a reduction is accepted, so that the registered product is
committed". **That has not happened.** It is the author's call; do not change it unasked.

---

## 2. Environment

- **Windows 11.** Python 3.14.3 (`pythoncore-3.14-64`), matplotlib 3.10.8, and numpy, scipy,
  xarray, netCDF4.
- **Install:** the package is installed **editable** from the repository. After adding an entry
  point to `pyproject.toml`, reinstall so the console script exists:

  ```
  pip install -q --no-deps --no-build-isolation -e .
  ```

- **Entry points:**
  - `casspian-lindal-raw`, `casspian-gravity-file`, `casspian-rotation-file`,
    `casspian-wind-from-curve`, `casspian-composition-lindal`, `casspian-lindal-inputs`: the
    SPEC_01 tools.
  - `casspian-refrac <manifest>`: the reduction. It also renders the figures when the manifest
    asks.
  - `casspian-plots <file.nc> [--out DIR] [--format png|pdf] [--dpi N]`: the diagnostics.
- **Shell:** Git Bash for POSIX commands, PowerShell also available. Run Python from the
  repository root; the acceptance scripts use relative paths such as `occul_data/lindal`.
- **Unicode output:** set `PYTHONIOENCODING=utf-8` when a script prints a DataTree or any
  non-ASCII text. The Windows console codec (cp1252) otherwise raises `UnicodeEncodeError`.

---

## 3. The procedure for every step

These rules are recorded as the author's rulings and have been applied throughout.

1. **Gate.** Read the status line at the top of the specification and `docs/specs/STATE.md`.
   A draft specification is not started, and a chat message saying "ready" does not accept a
   draft.
2. **Build only the step.** Acceptance script under `reports/step0<spec>_<N>/`, for example
   `reports/step03_1/accept_step03_1.py`. The SPEC_01 scripts sit under `reports/step<N>/`.
   These directories are ignored by git.
3. **Every check prints its measured value**, and a check the specification did not ask for is
   labeled "beyond the specification". Output goes to `output.txt` beside the script.
4. **Write `reports/REPORT_0<spec>_step<N>.md`:**
   - what was built;
   - decisions: every choice the specification did not make;
   - the acceptance table;
   - findings: wherever a specification number or expectation is wrong, with the measured value
     and the reason;
   - regression;
   - next step.

   State failures plainly. Do not loosen a check to pass it; report the physics (REPORT_02_step6
   finding 1 is the model).
5. **Set the step to `reported` in `STATE.md`. Do not commit. Stop.**
6. **After the review** (`reports/REVIEW_0<spec>_step<N>.md`), follow its order of work
   exactly:
   - make the changes;
   - rerun;
   - commit the author's spec edits and the review, then the step: the acceptance commit;
   - push;
   - run the sweep (section 5);
   - set `STATE.md` to accepted with the commit.
7. **Keep the author's documents separate.** Commit the author's uncommitted documents (specs,
   reviews, `STATE.md` edits) in their own commit before or with the step, never silently mixed
   in. Never overwrite or delete an author file; look first.

**Standing rules for everything written:**
- no em dashes or en dashes anywhere (check every file you write);
- American spelling;
- **no AI attribution in any commit message, document or code**. The author's rule overrides any
  tool reminder to add a `Co-Authored-By` line; no commit in this repository carries one;
- no silent choices: a choice the specification did not make is a decision in the report, or a
  question, never a quiet default.

---

## 4. The regression

Run every earlier suite after any change to `lib`, `refrac` or `tools`, and record the result
lines in the step's `regression.txt` and the report. A loop that works:

```
out=reports/stepXX/regression.txt; : > $out
for s in step1/accept_step1 step1/verify_review_changes step2/accept_step2 step3/accept_step3 \
         step4/accept_step4 step5/accept_step5 step6/accept_step6 step7/accept_step7 \
         step8/accept_step8 step9/accept_step9 step02_1/accept_step02_1 step02_2/accept_step02_2 \
         step02_3/accept_step02_3 step02_4/accept_step02_4 step02_5/accept_step02_5; do
  line=$(PYTHONIOENCODING=utf-8 python reports/$s.py 2>&1 | tail -1); echo "$line   <- $s" >> $out
done
```

**Last full result (at `2cf4457`):** 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01; 6, 7, 9, 9, 7
for SPEC_02 Steps 1 to 5; all pass. The Step 6 acceptance (`reports/step02_6`) passes 7 of 7.

**Timing.** Most suites take seconds. These are slow:

| Suite | Time | Why |
|---|---|---|
| `step7` | 3 to 4 min | wind geoid convergence |
| `step02_2` | about 1 min | |
| `step02_3` | about 3.5 min | reruns the fixed point four times |
| `step02_4` | about 4 min | builds the product through the CLI |
| `step02_5` | about 3 min | |
| `step02_6` | about 12 min | eight extra fixed points |

The whole loop takes about 20 minutes. Run it in the background. **Never run two suites that
write the product at once:** `step02_4`, `step02_5` and `step02_6` all rewrite
`lindal_refractivity.nc` and `figures/`.

**What each suite depends on**, so a later change does not break it silently:
- **`step02_1`** checks input and manifest hashes. It reads REPORT_01_step9 section 6, then any
  rows in REPORT_01_step8 and REPORT_02_step6, later rows overriding earlier. A new manifest or
  input hash must be recorded in a report table in that row format, or add the report to the
  list in the script:

  ```
  | `file` | kind | `64 hex digits` |
  ```

- **`step02_2` and `step02_3`** test values accepted under `mean_polar_radius`. They pin that
  rule in memory with `dataclasses.replace(manifest, anchor_rule="mean_polar_radius",
  anchor_quantity="radius_polar_m")`.
- **`step02_4` and `step02_5`** run on the committed manifest, now `equatorial_radius`. Their
  expectations hold under either rule.
- **`step9`** checks the directory listing and allows `lindal_refractivity.nc` once it exists.

**Side effects of rerunning suites. Clean up before any commit:**
1. `step02_4`, `step02_5` and `step02_6` rewrite `occul_data/lindal/lindal_refractivity.nc`
   from the working tree. It then carries `<HEAD>-dirty`. The sweep fixes that.
2. **`step02_5` overwrites the committed figures attached to REPORT_02_step5** in
   `reports/figures/step02_5_lindal_diag_*.png` with renders of whatever product is on disk.
   Restore them before committing:

   ```
   git checkout -- reports/figures/step02_5_lindal_diag_*.png
   ```

   Better still, change that script to stop copying once its step is accepted, as a reported
   decision.
3. `step02_6` check 1 reruns `casspian-lindal-inputs`. That rewrites `lindal_reduction.toml`
   with identical content, so the hash is unchanged, and keeps T and D.

---

## 5. The sweep, exactly

`lib.io.git_commit` marks a file `-dirty` when `git status --porcelain` shows **anything**,
untracked files included. A clean build therefore needs a completely clean tree.

1. Commit and push everything for the step.
2. Make the tree clean. Move or delete untracked scratch files. If uncommitted work for a later
   stage must survive, `git stash push -u -- <paths>`, then pop it after the build. Confirm the
   porcelain output is empty.
3. Rebuild each product that carries `-dirty` with its own tool. For the product:

   ```
   casspian-refrac occul_data/lindal/lindal_reduction.toml
   ```

   Under the current manifest this also renders `occul_data/lindal/figures/`.
4. Read back `casspian_git_commit` and refuse if it ends in `-dirty`.
5. Compute the SHA-256 with `casspian.lib.io.sha256`.
6. Record the hash in the step report's hash table, in the row format of section 4.
7. Commit that record ("Record the Step N sweep: ...") and push.

**Check what is `-dirty` across the directory:**

```python
from casspian.lib import io as cio
for n, k in [("lindal_thermo.nc","thermo"), ("lindal_composition.nc","composition"),
             ("lindal_geodesy.nc","geodesy"), ("lindal_gravity.nc","gravity"),
             ("lindal_rotation.nc","rotation"), ("lindal_wind.nc","wind"),
             ("lindal_refractivity.nc","refractivity"), ("raw/lindal_raw.nc","raw")]:
    h = cio.read("occul_data/lindal/" + n, k); c = str(h.attrs["casspian_git_commit"]); h.close()
    print(n, c[:12], "DIRTY" if c.endswith("-dirty") else "")
```

---

## 6. Pitfalls that cost time in this repository

1. **Bash heredocs that embed Python with backslashes or triple quotes break.** In several
   different ways:
   - `\n` inside a replacement string became a real newline and split an f-string;
   - bash failed to parse a heredoc containing both `'` and `"""`, and **nothing ran**, while
     the chain reported nothing;
   - a `sed` range deletion ran to the end of a file.

   What works:
   - Write Python patch scripts to a scratch file with the Write tool, then run them.
   - Or use the Edit tool on files already read.
   - Build awkward characters from `chr()`.
   - Assert each replacement matched exactly once.
   - Compile (`python -m py_compile`) after every patch.
2. **The Edit tool refuses a file that has not been read with the Read tool in the session.**
   Reading it with `cat` does not count.
3. **`set -e` is not reliable in a long chained command here.** A failing heredoc Python patch
   did not stop the chain once. Use explicit `|| exit 1` after every step that matters, and
   verify results by reading them back.
4. **Anything shown in a panel or claimed in a report must be measured.** Twice a report sentence
   carried an unmeasured number: a 13 km equator shift that was really 6.5, and a no wind radius
   that was never computed. Check derived figures against the output before writing them.
5. **Matching `STATE.md` rows by `| 6 |` is ambiguous.** SPEC_01 and SPEC_02 both have a Step 6.
   Match on the deliverable text.
6. **Figures need looking at.** Acceptance checks do not see overlapping labels, legends over
   data, or a curve lying on zero. Render, read the PNG, fix, then run the final acceptance.
   Attach figures from a clean product, since the footer shows the commit.
7. **xarray DataTree:** set a root dataset with `tree.dataset = ...`; `tree["/"] = ...` raises.
   Use `node.to_dataset(inherit=False)` to get a group's own variables.
8. **The kind N reader validates hashes.** `lib.io.read(path, "refractivity")` refuses a file
   whose `input_hashes` no longer lists a hash recorded in `reduction_record`. It warns when an
   input on disk differs. Copies made for refusal tests must be edited with that in mind.
9. **Radians everywhere in `lib`.** Degrees are converted once at the `refrac` boundary (the
   manifest tolerance, the kind T label and its uncertainty) and back to degrees once, at the
   product write.
10. **Anchor rules and quantities are paired.** `casspian.lib.control.ANCHOR_QUANTITY_FOR_RULE`
    is the one table; `lib.control` and `casspian-lindal-inputs` both refuse a mismatch. The no
    wind reference geoid is anchored where the rule anchors (`anchor_latitude`), not always at
    the pole.
11. **Long jobs.** The anchored fixed point takes about 30 s under polar anchoring and about 65 s
    under equatorial, where the secant starts farther away. Anything running several of them
    belongs in the background, and the report waits for its result.

---

## 7. Code map, for orientation

| Module | What it holds |
|---|---|
| `lib/constants.py` | CODATA 2018 |
| `lib/schema.py` | kind registry, `validate`, `uncertainty_companion` (`_uncertainty` inserted before the unit suffix; `UNIT_SUFFIXES` includes `_kg_mol`), `WRITER_FILLED_GLOBALS` |
| `lib/io.py` | `write`, `read` (DataTree for C, N and raw; Dataset otherwise; kind N hash check), `sha256`, `input_hashes`, `git_commit`, `package_version` |
| `lib/control.py` | `load_section` for tool control files, `read_reduction_manifest`, `load_reduction_inputs`, `ANCHOR_QUANTITY_FOR_RULE` |
| `lib/gravity.py` | `g_newton`, `G_phi_newton`, `omega_abs`, `g_eff_radial`, `G_phi_eff`, `g_eff_vector` (returns `psi`) |
| `lib/geoid.py` | `U_rigid`, `reference_geoid(..., anchor_latitude)`, `wind_geoid` (rules `mean_polar_radius`, `north_pole`, `south_pole`, `latitude`, `equatorial_radius`; `WindGeoidResult` with polar radii, asymmetry, residual, `equator_radius_m`), `radius_at` |
| `lib/latitude.py` | `planetocentric_fixed_point`, `planetographic_from_planetocentric` |
| `lib/reduction.py` | `number_density`, `absolute_radius`, `mean_over_species`, `refractivity`, `temperature_from_refractivity`, `share_of_remainder_partials`, `quadrature` (NaN terms left out and named), `to_one_sigma` |
| `refrac/anchor.py` | `geoid_setup`, `freeze_anchor` returning `FrozenAnchor` |
| `refrac/reduce.py` | `composition_closure`, `anchor_partials` (total `dr0/dr_anchor`), `reduce_profile` returning `Reduction` with companions, partials and terms |
| `refrac/product.py` | `build_product`; `main` is `casspian-refrac`, rendering figures when asked |
| `tools/plots/` | `render`, `style`, `figures_product` (F1 to F6), `figures_inputs` (W, T, C views) |
| `tools/lindal/build_inputs.py` | kinds T and D, kept when their content is unchanged; the manifest with `[geoid]` and `[diagnostics]` |
| `forward/` | empty; SPEC_03 |

---

## 8. Where the author will pick up, and what to expect

The author is writing SPEC_03, the forward model round trip. The existing handoff's section 8
lists what it must settle:
- the geopotential B5, with a level and layer staggering table and stated numerics;
- the gauge isobar;
- hydrostatic production B7 with `p_b` at the top tabulated level;
- recovery of Table I T and p to a stated tolerance;
- F5 and F6 rendering;
- a forward product kind with `modeled` provenance;
- the `input_hashes` warning for all derived kinds.

It will likely arrive as a draft. **Do not code against a draft.** When a step proceeds:
- expect the F5 and F6 placeholders in `tools/plots/figures_product.py` to be adopted or renamed
  (`geopotential_m2s2(level)`, `pressure_hydrostatic_Pa(level)`, attribute
  `boundary_pressure_Pa`);
- expect `lib.io.read` to gain the general `input_hashes` warning;
- expect a new manifest or namelist, and a `forward/` package with its own entry point.
