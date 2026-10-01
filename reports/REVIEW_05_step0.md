# REVIEW 05, Step 0. The scripts under tests/, the Linux clone, and the regression of record

Reviewing agent, 1 October 2026. Report: `reports/REPORT_05_step0.md`. Read with the Linux logs in
`reports/regression/linux_2026-10-01/` and the files on the Windows machine.

## Verdict

**Accepted with one change** (the ruling on section 5). Step 0 closes when the clone's rerun of the
six suites of section 8 reads 6 of 6 under that ruling.

## What was checked

- **The Linux run of record** (`regression.txt` and the four failing suites' logs, clone at
  `4d32efe`): 21 of 25, the four failures exactly the two causes the report names. `step04_0`
  check 14 fails only because it reruns `step03_3`.
- **Cause 1, line endings.** The diagnosis matches what I found on disk before the fix (the
  Windows working copies of the build control files were CRLF while `.gitattributes` pins LF). Now:
  `occul_data/lindal/lindal_build.toml` and `forward/lindal_transfer/lindal_transfer_build.toml`
  carry no CR; `lindal_refractivity.nc` hashes to `b79c60cb...` and the closure product to
  `f9029197...`, both as section 7 states. The evidence that nothing moved (518 computed
  variables array-equal across all 18 files, only stamps, hashes and resolved paths differing) is
  the right test and is sufficient. The `skip-worktree` handling of the sweep is a sound answer to
  committed products meeting the dirty stamp, and Step 2 should use it.
- **Cause 2, p and T on Linux.** 1 and 2 units in the last place, on 3 and 2 of 66 levels, with
  N, geopotential, refractivity and molar mass identical. That is the platforms' `exp` and `log`
  rounding, not a difference in the calculation.
- **The Windows regression for the change**: 6 of 6, `step03_3` check 7 now silent on the unedited
  kind T, and the new departure lines at 0 ulp on the registering platform.

## Ruling on section 5

One rule on every platform, with no platform detection: `step04_1` check 10 and `step04_5`
check 11 pass when every compared variable agrees with the registered product to **1e-14
relative**. Each check keeps printing whether the variable is bit-identical and its largest
difference in units in the last place, so a change of even one bit on the registering platform is
still visible in the output, and the reviewing agent reads it there. The bound is about 37 times the
largest measured cross-platform difference and many orders below anything physical. Neither suite's
reference count changes.

## Order of work

1. Apply the ruling to the two checks; rerun `step04_1` and `step04_5` on Windows (both at their
   reference counts, 0 ulp printed).
2. The author's go; commit the report, this review, the suite changes and the ruling. No AI
   attribution.
3. The clone pulls and runs the six suites of section 8. Expected 6 of 6, with `step04_1` and
   `step04_5` printing their 1 and 2 ulp departures under the 1e-14 bound.
4. `STATE.md`: Step 0 accepted with that run's result. Then Step 2.
