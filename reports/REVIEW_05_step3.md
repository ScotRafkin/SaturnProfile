# REVIEW 05, Step 3. The shear step in every forward run, and the new-run tool

Reviewing agent, 2 October 2026. Report: `reports/REPORT_05_step3.md`. Specification:
`docs/specs/SPEC_05_Shear_Experiments.md` v0.6, §3. Read with the two F9 renders under
`reports/step05_3/`.

## Verdict

**Accepted, with one small addition** (the runbook, below).

## What was checked

- **The migration moves no computed value.** The closure product equals the registered one in
  `step04_5` check 11's nine variables, and the transfer product, rebuilt through the five sections
  and through a run made by `casspian-new-run`, equals the registered transfer product in N, p, T,
  altitude, `r0` and the pressure identity. The closure content comparison passes with the
  identity-sheared wind, its new `history` and `input_hashes` among the dropped attributes, as
  §1.7 expected.
- **`[shear]` required**, refused by name when missing; its `role` and `prefix` checks skipped and
  its output checked like every section's, as v0.4 ruled.
- **`casspian-new-run`.** The files differ from the source only in the name, renaming back gives the
  source exactly, and a second call is refused.
- **F9.** I viewed both renders. The migrated transfer run shows the two lines coincident; the
  `c = 0` run shows the source wind and `u_total` on zero, and the right panel now draws a constant
  field as one colour with its value instead of contouring round-off. That change was the
  author's request and is right.
- **The full regression**, 29 of 29 at their reference counts, the registered files restored and
  identical to the commit afterwards. No accepted suite's expectation changed. The discarded double
  run is reported plainly, with the evidence that nothing registered was harmed; running the driver
  only once at a time is the lesson.
- **Rebuilding only the two runs' files** (`sweep.py --runs`) is the right scope: the reduction
  chain is not reached by this step, and rebuilding it would change its bytes and the hash tables
  for no change of content.

## The addition

**The runbook's stale statements** (finding 1) are corrected in this step's runbook v0.4, since the
runbook is this step's deliverable and the corrections cannot move a number: the repository does
hold the registered netCDF files (since `4d32efe`); the closure product's hash is the current one
or is dropped in favor of "compare by value"; the regression takes about 100 minutes and its last
line reads 29 of 29. No suite reruns for it.

## Notes (none changes the verdict)

1. **`casspian-new-run`'s token rule** (finding 2): renaming a name followed by `_` is what makes
   the `<run>_` prefixes follow, and the one case it gets wrong (a run name that is a prefix of
   another run named in the same files) does not occur. Accepted as is.
2. **An identity run keeps the source's `vertical_structure`** (finding 3), so F9's right-panel
   title says `altitude independent`. That is §1.7 as written and is true of the file.

## Order of work

1. Correct the runbook (above).
2. The author's go; the acceptance commit: the source, control, test and runbook files, the report
   and this review. No AI attribution.
3. `python tests/step04_0/sweep.py --runs` on the clean tree; compare the rebuilt run files with the
   committed ones by value (computed columns array-equal, as checks 4 and 5 found); commit them and
   add the comparison to the report.
4. Since the registered run files' bytes change, rerun the suites that read them by bytes or hash:
   `step03_4`, `step04_0`, `step04_1`, `step04_5`, `step05_1` and `step05_3`. All at their
   reference counts is expected.
5. `STATE.md`: Step 3 accepted. Then Step 4, the experiments; my independent values for runs 5 and
   9 come before its review.
