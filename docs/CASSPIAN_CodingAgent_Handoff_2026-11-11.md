# CASSPIAN coding agent handoff, 11 November 2026

From the outgoing coding agent, at SPEC_21's closure. Author of record: S. Rafkin (SwRI), who approves everything.
Supersedes `CASSPIAN_CodingAgent_Handoff_2026-10-22.md`; read with `CASSPIAN_CodingAgent_Brief_2026-11-09.md`, whose
rules (materiality, figures, no new keys, Lindal untouched, one composition) stand.

## 1. Before you start

- Read `docs/specs/STATE.md` (the gate: a step's status line there and in its specification decides what may be
  done), then the specification you are given, then its reviews in `reports/`.
- `main` is at `7d61f86`, pushed. SPEC_21 closed at v0.7: closing commit `0226110`, then `STATE.md` with the
  Voyager anchors and the hypothesis 4 wind rebuilt clean at it.
- Uncommitted in the working tree, and not yours to commit unless the author asks:
  - the author's edited documents: `docs/CASSPIAN_CodingAgent_Brief_2026-11-09.md`,
    `docs/CASSPIAN_ReviewingAgent_Handoff_2026-11-09.md`, `docs/casspian_status_2026-10-22.md`;
  - the untracked `docs/CASSPIAN_Manuscript_Update_Brief_2026-11-09.md`;
  - `fits/voyager_fit_2/voyager_fit_2_reference_31s_wind.nc`, the swap fit's wind, a check 5 product left out of
    the closing commit because the review did not list it.
- Next: SPEC_22, the final single answer blending all anchors, with its own profile matrix and sections
  (REVIEW_20_step1 item 5).

## 2. Where things stand

### 2.1 The wind fit, `casspian-wind-fit`

`src/casspian/tools/wind/estimation/`. One control file per fit under `fits/<name>/`, committed with its wind file;
`fits/*/work*/` is ignored.

| Module | What it holds |
|---|---|
| `driver.py` | `read_control` into `Problem`; `extend_top` (wind files reach one decade above 1 Pa, the shear held, provenance 3); `linearize` (forward runs in a process pool, then the radio, IRIS and stability blocks); `iterate` (plain pass, then Anderson(1); ends only on a plain step under `tolerance_K`; prints each solve's cost in its parts); `size_systematic`; `run` (fit, size IRIS's systematic parts, fit again); `write_wind` |
| `state.py` | `Basis`: the B-spline shear state in latitude and ln p; `project`, `x_integral` (linear above the top), `row`, `x_matrix`, `x_matrix_forward` |
| `prior.py` | correlation forms (exponential, Matérn 3/2 and 5/2, Gaussian); `DensePrior`; the temperature-stated prior with the deformation-radius floor |
| `radio.py` | `RadioTargets`: profiles read with their uncertainty; `delivered` (rows `int w Delta s dphi` along the path from the reference); `block` (row per delivered level) and `layer_block` (with `average_layers`: every profile, the reference's shared error included, averaged into the forward wind grid's layers by `layer_mean`, the floor once per layer) |
| `gradient.py` | IRIS gradients at 110, 290 and 730 mbar, the gap rule, systematic parts |
| `stability.py` | the softened hinge on the wind grid's node layers, relative to the anchor's own excess |
| `solver.py` | observation-space solve; `robust` (Huber on whitened residuals per group); `consistency`; `leave_out`; `systematic_term`; `cost_parts` (prior, shared anchor term, each target, each IRIS level, stability; before and after each solve) |

Everything general: the reference is `[radio.reference]`, never assumed to be Lindal's. The control keys added in
SPEC_21 carry their bases in comments in `fits/voyager_fit_2/voyager_fit_2.toml`. `floor_correlation_ln_p` and
`floor_correlation_form` are in code, unused. A control without `average_layers` (as `voyager_fit_1`) behaves as
before.

### 2.2 Hypothesis 4 and the acceptance script

- `fits/voyager_fit_2/voyager_fit_2.toml`, reference 36.3 N (30.8 N), Lindal; targets 31.2 S (26.1 S), 74.2 S
  (71.1 S), 3.0 S (2.4 S). `voyager_fit_2_reference_31s.toml` is the swap with the Voyager 2 egress as reference.
  Its prior stays Lindal's on purpose.
- `tests/step21_1/accept_step21_1.py`. Run from the repository root.
  - The full run takes about six hours, then check 5 a few more.
  - `--resume` starts after the matrix runs from `reports/step21_1/resume.pkl`. **That resume point predates the
    layer averaging and is stale.** Use it only after a full run has rewritten it.
  - `--refit` refits hypothesis 4 under the current control, reports its change in delivered T from the resume
    point's fit (bound 1 K), writes the wind file, and runs check 5 (`check_reference`).
  - Check 5 deletes its runs after comparing them (`fits/voyager_fit_2/work_reference`).
- Timings:
  - A Lindal-referenced linearization takes about 5 minutes; a Voyager-referenced one about 15, since its
    products sit on about 300 levels.
  - Hypothesis 4's fit converges in 4 and 3 passes, the swap's in 4 and 2.

### 2.3 Results the reviews carried forward (for SPEC_22, noted, not chased)

- Reference-independence checks are judged on layer means with the materiality bound (about 1 K), the resolution
  the fit uses (REVIEW_21_step1_check5).
- The Lindal-referenced fit's radio terms are 4 to 7 times the swap fit's at the same targets, though both sit
  well inside their errors (misfit over expected 0.04 to 0.33). Not explained.
- The equatorial shear stays largely the prior's; cloud winds are to constrain it (SPEC_23).

### 2.4 The check 5 history, in one paragraph

The swap fit diverged: its solve asked for about 15 K every pass. The cause was not the operator, whose one-row
steps matched the forward runs to 2 to 6 percent. The Voyager profiles' 300 levels against Lindal's 58 gave
structure below the forward model's 0.115 resolution four times the weight. A correlated floor made it worse
(42 K). Averaging every profile into the wind grid's layers fixed it: the request fell 9 K to 0.4 K in one pass.
A Levenberg-Marquardt trust region was written and removed unused; its design (damping as the prior scaled by
`1 + lambda` about the current state, a step accepted only when the forward-evaluated cost falls) is in the record
of this session only. Full account: `reports/REPORT_21_step1_check5_diagnosis.md`.

## 3. The regression

- `tests/run_regression.sh` lists suites through the older steps. It does not list `step15_1` to `step21_1`; the
  wind fit's suites take hours and are run by their own scripts.
- No regression was run at SPEC_21's closure: no computed value of an accepted product outside the fit could
  change.
- Lindal's registered files were checked unchanged by md5 before and after the rebuild.

## 4. Rebuilding clean, as done at `0226110`

1. Save the committed bytes: `git show HEAD:<file> > <scratch>/committed/<file>` for each product.
2. Anchors, in `occul_data/lindal/` for each of `voyager1_egress`, `voyager1_ingress`, `voyager2_egress`:
   `casspian-voyager-raw <p>_scalars.toml` (run in `raw/`), then `casspian-composition-lindal <p>_build.toml`,
   `casspian-voyager-inputs <p>_build.toml`, `casspian-refrac <p>_reduction.toml`. Mid-chain warnings that an
   `input_hashes` entry does not match the raw file are expected; the finished set is consistent.
3. The fit: `casspian-wind-fit fits/<name>/<name>.toml`.
4. Compare every variable with the saved bytes (xarray, NaN-aware). Values were equal in every variable; bytes
   differ by history and hash attributes. Commit the rebuilt bytes with `STATE.md` (SPEC_20 precedent).

## 5. Things the next agent should know

- **Windows text writes.** `Path.write_text` on Windows writes CRLF. Git normalizes on commit, but a wind file
  records the hash of its control file, so write with `newline="\n"` or as bytes. Convert before committing.
- **Process pools.** Any script that uses the driver's pool needs an `if __name__ == "__main__":` guard (spawn
  start method), or every worker re-runs the script.
- **Long jobs.** Use the tool's background mode, never `&` in Bash. Find processes with `Win32_Process`
  (PowerShell `Get-CimInstance Win32_Process`), and stop a run's workers as well as the parent.
- **Edit scripts.** Bash heredocs fail on some quoting and mangle backslashes; write edit scripts to files in the
  scratchpad. Do not edit a module while a run that spawns workers from it is going, unless the change is inert.
- **Git identity.** No `user.name` or `user.email` is configured in this repository. Commits have used the
  author's identity as recorded in `git log`, passed with `git -c user.name=... -c user.email=...`.
- **Latitudes.** Every latitude planetographic with planetocentric in parentheses.
- **The solve's diagnostics.** Each pass prints the cost in its parts; a solve that raises a target's misfit while
  lowering the total is the signature of weight on unresolved structure (§2.4).

## 6. How the work is judged, and practice

As `CASSPIAN_CodingAgent_Handoff_2026-10-22.md` §6, with the brief of 9 November, in short:

- **General code.** It works on any profile; no rule for one dataset. Explain an unexplained behavior with a
  measurement before acceptance.
- **Bounds.** Physical, general, none set for one case. A missed bound: stop and report with every other check
  measured; do not tune to pass.
- **Diagnosis needs approval.** Get the author's approval before trying a diagnostic approach. Do not proceed
  to a next step without explicit permission.
- **Materiality.** About 1 K in delivered temperature, about 5 percent in wind at probe levels.
- **Figures.** Absolute profiles first, with both hydrogen adiabats. Sections shade the zonal wind under labeled
  temperature contours; no hatching, no IRIS dots.
- **Regression.** Full only when a computed value can change or at a named milestone. Lindal's files untouched.
- **Commits.**
  - Commit only after a step is reviewed and accepted, then push `main` without asking.
  - The acceptance commit holds what the review lists.
  - A second commit adds the hash to `STATE.md`'s step rows, with the products rebuilt clean. The step rows are
    the coding agent's; the header lines are the reviewing agent's.
  - Data files only where precedent tracks them: registered products, fits' wind files, run and section control
    files.
- **Style:**
  - no em or en dashes in prose;
  - American spelling;
  - every measured number with its bound;
  - no attribution of any kind to an AI in commits or code;
  - the work is SwRI funded and proprietary, and nothing leaves the repository.
