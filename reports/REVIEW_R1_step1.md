# REVIEW R1, Step 1. The lean pass

Reviewing agent, 30 October 2026. Report: `reports/REPORT_R1_step1.md`. Specification:
`docs/specs/SPEC_R1_Regression_Fixes.md` v0.4.

**Accepted. SPEC_R1 closes.**

- **The lists** follow the rule and §2's rulings item for item. Nothing was kept without a wrong
  number or an obscure failure behind it, and nothing was left unplaced.
- **The regression.** 38 suites in six chunks, every count as `tests/run_regression.sh` now
  lists them, after the two expectation changes. Kind N's hash is `b79c60cb...` before and after
  every chunk, and the closure and transfer files are restored equal. Checked against the chunk
  logs. No computed value changed.
- **`step05_3` checks 4 to 7** ran for the first time since SPEC_11 Step 2: the closure and
  transfer products equal the registered ones. Nothing new.
- **§4, the T and D comparison in `casspian-lindal-inputs`:** accepted. Comparing only the globals
  the writer still writes is what "old files keep theirs and still read" asks, and `step02_6`
  caught the case.
- **The two moved expectations:** accepted.
  - `step03_3` check 7. The registered kind T now warns once, on reading, because its
    `input_hashes` names `lindal_build.toml`, which this step edited. The warning is correct and
    stays until the reduction is next redone. No rebuild now.
  - `step05_3` check 6: `FileExistsError` from `mkdir` is the refusal §2 ruling 3.5 asked for.
- **`tests/step03_0/`** is not in the driver and cannot run on today's reader. Under the rule it
  goes, in the acceptance commit. Its record is `reports/REPORT_03_step0.md`.
- **The `value_source` warning.** SPEC_00 §8 used to list a warning for a missing `value_source`,
  and v0.25 drops it. If any code still issues it, that code goes in the acceptance commit, with
  no rerun.
- **SPEC_00 v0.25** is written: §8 is now the whole of what readers enforce, from the report's
  kept list, and it governs older wording elsewhere. §5, §6.1, §6.2, §6.6, §7.2 and §7.3 are
  amended to match.

**Order of work:** commit and push, with no further regression.
1. The acceptance commit:
   - the code, suites, control files and driver of Step 1, with `tests/step03_0/` removed;
   - `reports/REPORT_R1_preexecution.md`, `REPORT_R1_step1.md` and this review;
   - SPEC_00 v0.25, SPEC_R1 v0.4, SPEC_14 v0.2 and `docs/specs/STATE.md`;
   - `docs/CASSPIAN_CodingAgent_Brief_2026-10-22.md` and `docs/casspian_status_2026-10-22.md`;
   - the untracked reference documents in `docs/` (§2 ruling 4).
2. A second commit with the hash in SPEC_R1's row of `STATE.md`. Push.
3. Then SPEC_14, Step 1.
