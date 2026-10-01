# REVIEW 05, Step 1. The shear tool, and the registered transfer product

Reviewing agent, 30 September 2026. Report: `reports/REPORT_05_step1.md`. Specification:
`docs/specs/SPEC_05_Shear_Experiments.md` v0.4. Read from the Windows working tree: the report,
`reports/step05_1/output.txt`, `src/casspian/tools/wind/shear.py`,
`tests/step05_1/accept_step05_1.py`, `.gitignore` and `pyproject.toml`.

## Verdict

**Accepted.** 11 of 11, every number in the report is a line of the printed output, and the code
does what §1.3 to §1.7 say and nothing more.

## What was checked

- **The three layers of §1.5.** The parts (`u_at_pressure`, `position_ln_p`, `position_p`, `ramp`)
  know nothing about any case; `CASES` maps names to functions with their required and optional
  parameters; `assemble` is the one path every case but `identity` takes, and forms `u_shear` as the
  difference, so the sum identity holds by construction. `construct` opens no file; `build` reads,
  constructs, sets `input_hashes` (source and control file) and writes. The Monte Carlo split is
  intact.
- **Replace, never add.** `assemble` overwrites exactly `u_total`, `u_shear`, `value_provenance`
  (2 in every cell), `vertical_structure` (the case name) and, only with `uncertainty_ms`, the
  uncertainty's values and its two descriptive attributes. Check 8 confirms no other global or
  variable attribute differs from the input except the four `lib.io.write` stamps.
- **`u_s` through `WindField`** (ruling 9). Check 2 settles the question left open: at a node the
  interpolant returns the node value exactly.
- **The refusals.** Exactly the seven kinds of check 10, each with a message naming the key and the
  rule, and no file written. `role` and any other stray key are refused by the existing unknown-key
  rule of `load_section`, which is the typo guard and not a role rule.
- **The committed input is not written over.** Its SHA-256 is the same before and after the run.
- **§1a.** The transfer product was rebuilt on a clean tree at `7ef249b`, its 88 computed variables
  array-equal to the previous copy, `r0(10 N)` 60,128,612.966 m and the pressure identity 5.777e-7,
  both the accepted Step 5 values, and committed at `7276f68` with its `.gitignore` exception
  (confirmed on disk). The variable for Step 3 check 3 is
  `reference_surface/reference_surface_radius_m` at the target node.

## Notes, none of which changes the verdict

1. **Checks 4, 5 and 7 measure zero because the script uses the tool's arithmetic.** The report says
   so. The formulas are stated independently from §1.3 and §1.4, and check 6 (a different
   arithmetic form) agrees to 1.5e-16, so there is no concern; the independent values for runs 5
   and 9 at Step 3 remain the real external test of the ramp.
2. **The regression driver fix is right, and its scope is right.** Restoring the committed transfer
   product after every suite is needed now that it is registered, and it can only change what
   `step04_5` leaves on disk; that suite was rerun (11 of 11). No other suite imports changed code.
3. **Step 0's report is not on disk.** `reports/REPORT_05_step0.md` does not exist in the working
   tree. Step 1 does not depend on it, but SPEC_05 closes Step 0 with its review. The author should
   say whether the Linux regression of record has finished; the Step 0 report and its review should
   be filed before Step 2, whose migration reruns most of the regression and would otherwise mix
   with it.

## Order of work

1. The author's go.
2. The acceptance commit: `shear.py`, `pyproject.toml`, `tests/step05_1/accept_step05_1.py`,
   `tests/run_regression.sh`, `REPORT_05_step1.md` and this review. No AI attribution.
3. `STATE.md`: Step 1 accepted, with the commit.
4. Step 2 after Step 0's report is reviewed (note 3).
