# REVIEW 05, Step 2. The column march, faster (decision R, first part)

Reviewing agent, 1 October 2026. Report: `reports/REPORT_05_step2.md`. Specification:
`docs/specs/SPEC_05_Shear_Experiments.md` v0.6, §2a. Read with `reports/step05_2/output_fix2.txt`
and `STATE.md` on the Windows machine.

## Verdict

**Accepted.** The step did what §2a asked and nothing else: no value moved, and the time fell by
an order of magnitude.

## What was checked

- **No value moved.** After fix 1 and again after fix 2, all 518 variables of the 18 registered
  files, rebuilt as candidates in a copy with the `-dirty` refusal relaxed for that script only
  (v0.6), are bit-identical to the committed files; the per-file counts in `output_fix2.txt` add up
  to 518, the closure product's 136 and the transfer product's 177 included. The 1e-14 fallback was
  not needed. No registered file is recommitted.
- **The full regression**, 28 of 28 at their reference counts, `step04_4` included, with the
  registered files restored exactly. This is the evidence that the kept `column` and
  `build_columns` signatures serve every accepted suite unchanged.
- **The refusal inside the array march** (check 3) names the latitude, the geopotential and the
  pressure, as deliverable 4 asks.
- **The three places a last bit could have moved** (the order of the sum over degrees, the
  array and single-value math routines, `np.interp`) were measured before the march was written,
  and the report says plainly that bit-identity on another platform is not measured; the 1e-14
  rule covers it. That is the right statement.
- **Scope.** The changes are in `lib.gravity`, `lib.mesh` and the isobar map and wind table of
  `forward/transfer.py`, plus the one-time validation of the degree list where the gravity file is
  read (`production.py`, `refrac/anchor.py`). No numerical setting, scheme or spacing changed.
- **Step 0's closure** on the Windows result, with the clone's rerun waived and no routine Linux
  runs, is the author's direction (REPORT_05_step0 §10) and is recorded in `STATE.md`.

## The numbers

A production transfer run: 94 s to 12.7 s; the column march inside it: 78 s to 0.05 s. The full
regression: about 8 hours to 95 minutes; `step04_4` from 4.3 hours to 15 minutes.

## Notes for later decision R work (none changes the verdict)

1. **Reading the inputs is now the largest single cost of a run** (7.3 of 12.7 s), with the
   reference-surface march second (4.2 s). For the Monte Carlo both matter differently: the inputs
   are read once per ensemble when the driver passes objects in memory (SPEC_05 §1.5), while the
   reference surface depends on the wind and is marched once per draw.
2. **The kind N build takes about 90 s** (`refrac`, through `lib.geoid` and `refrac.anchor`). It
   is paid in every suite that rebuilds kind N and would be the next target if regression time
   matters again.
3. **The test-side helpers** that still build columns one latitude at a time (`step04_3`'s line
   integrals, the cylinder wind of `step04_2` and `tests/step04_5/fields.py`) now dominate the long
   suites. Moving them to `build_columns` with `wind_all` is a test-side change for whenever it is
   wanted.

## Order of work

1. The author's go.
2. The acceptance commit: the five source files, `tests/step05_2/accept_step05_2.py`,
   `tests/run_regression.sh`, the report and this review. No AI attribution.
3. `STATE.md`: Step 2 accepted. Then Step 3, the shear step in every build.
