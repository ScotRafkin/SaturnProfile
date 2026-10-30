# CASSPIAN coding agent handoff, 3 October 2026

For the incoming coding agent. Author of record: S. Rafkin (SwRI), who approves everything. This
follows `docs/CASSPIAN_CodingAgent_Handoff_2026-09-29.md`; its sections 4 and 5 still hold where
this note does not replace them.

## 1. Before you start

Do not begin until the author says so, and do not start work without a specification. The next
specification, SPEC_08 (a realistic shear case and the first comparison with CIRS), is being
written by a new reviewing agent. Your first deliverable is your reading of it (section 3).

Start from a clean `main`: the documents commit that follows `7652852`, pushed to `origin`.

## 2. Where things stand

- **Specifications.** SPEC_00 to SPEC_07 are closed.
  - SPEC_05: the shear tool `casspian-wind-shear`, the required `[shear]` build step,
    `casspian-new-run`, and eleven named experiments.
  - SPEC_06: the shear read from the wind file's interpolant at the curve's point and label, and
    the tracing integrated per column with the wind file's pressure nodes as breakpoints.
  - SPEC_07, closed and not adopted: a PCHIP wind interpolant. Decision L (linear) stands in the
    code. The work is on branch `spec07-l2-parked` (`89e3e6d`); leave it alone unless a
    specification asks for it.
- **Regression of record:** 30 of 30 suites at their reference counts at `7652852`, about
  115 minutes (6923 s), driven by `tests/run_regression.sh`.
- **Registered products.** They are committed. `python tests/step04_0/sweep.py` rebuilds the whole
  chain, and `--runs` rebuilds only the two runs' files.
- **The outer loop.**
  - The eleven experiment namelists run at a tolerance of 1e-3.
  - `forward/lindal_transfer/lindal_transfer.toml`, which `casspian-new-run` copies, is at 1e-8.
  - `step04_4` and `step04_5` now pin 1e-8 themselves; `step04_5` run 1 asserts it, since that
    run goes through the namelist.
- **A redesign of the vertical is tabled.** It would stagger the variables, after SPEC_00 §7.3
  rule (ii) (`docs/CASSPIAN_Staggering_Design_Note_2026-10-03.md`). Do not act on it.

## 3. Your first deliverable

Read SPEC_08 against the repository and file `reports/REPORT_08_preexecution.md`. Cover:
- every place the specification's description of the code does not match what is on disk;
- every acceptance you cannot write as stated;
- every choice it leaves to you that the author should make.

Measure where a measurement answers a question. Then wait for the rulings, which the reviewing
agent writes into the specification.

## 4. How the work is judged (author's rules, 2 and 3 October)

- **The code must work on any profile, real or invented.** An unexplained behavior is explained in
  the report, with a measurement, before the step is accepted. It is not deferred.
- **Tests are general.** No tolerance is set for one case. Experiments are samples of a general
  test, and bounds are physical: a fraction of a kelvin is noise. Machine epsilon is asked for only
  where the mathematics leaves a value unchanged. If a specification seems to ask for a bound per
  case, say so in your reading.
- **The pressure identity** (`pressure_identity_residual`, level by level in every transfer
  product) is the model's own gauge of whether its numerics limit a run. Report it on every run.
  Its mechanism at a sharp change in shear is the check 7 law of SPEC_06: the anchor layer's mass,
  times the jump in `ln T`, times where in the layer the change falls.
- **When a check misses its bound, stop.** File an interim report with every other check
  measured, and wait. Do not tune numerics, spacings or tolerances to get through.
- **One fix per specification**, each revertible by git. No switch between old and new code paths.
- **Replace, never add.** No new metadata, titles or history lines. No idiot checks (decision N):
  refuse only typos, extrapolation and named numerical failures.

## 5. Practice

- **Regression:** rerun what a change can reach (SPEC_04 v0.22 §0). A change to `lib`, `refrac`,
  `forward` or the registered products means the full set. Name the suites rerun and why the
  others were not.
- **Accepted suites.**
  - A suite whose expectation moves is updated in the same step and named in the report with its
    old and new value.
  - None is edited silently or left failing.
  - Suites are pinned to the spacings and tolerances their bounds were set for.
- **Products from a dirty tree** end `-dirty` and are refused as inputs. Rebuild on a clean tree
  after the acceptance commit. The acceptance scripts relax the refusal for their own process only.
- **Close every netCDF read before opening another.** Do not edit a script while a shell is
  running it.
- **Order after an accepted step:** the author's documents commit, the acceptance commit, the
  sweep if products moved, the push, then `STATE.md`. The step rows of `STATE.md` are yours; the
  header lines are the reviewing agent's.
- **Style:**
  - no em or en dashes in prose;
  - American spelling;
  - every measured number with its bound;
  - no attribution of any kind to an AI in commits or code;
  - the work is SwRI funded and proprietary, and nothing leaves the repository.
