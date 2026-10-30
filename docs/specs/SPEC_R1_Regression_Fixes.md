# SPEC R1. The lean pass: checks and keys that are not needed go

Version 0.4, 22 October 2026. Author of record: S. Rafkin. **Status: accepted by the author,
22 October 2026, with the `-dirty` refusal removed; §2 rules on `reports/REPORT_R1_preexecution.md`.**
Depends on SPEC_13 (`cf33cd6`). Precedes SPEC_14.

## 0. Purpose

The regression of 8 October found three failures (`docs/CASSPIAN_CodingAgent_Handoff_2026-10-22.md`
§4). Two are checks that refuse correct data: a closure refused over descriptive attributes, and
a directory listing that names files. Much more data is coming, from many sources. Nothing may
accept only one kind of source or one file name.

**The rule (author, 22 October 2026).** A refusal, a test check or a required key stays only if
without it the code would compute a wrong number without saying so, or fail with no clear
message. Everything else goes. Removing it is the fix, not loosening it.

No computed value changes. The registered products and the committed inputs are not rebuilt.

Protocol as SPEC_04 §0.

## 1. Step 1

**Deliverables.**

1. **Go, known now:**
   - the closure content comparison: `check_closure_inputs`, `CLOSURE_DROPPED_ATTRIBUTES` and
     what only they use;
   - the closure season equality;
   - the refusals on names and labels: `profile_or_run` against the run name or the anchor's
     slug, `role`, `composition_role`, the run prefix on input files, `p_b_location`;
   - `step9` check 6, the directory listing.
2. **Go, by the rule, wherever found:** in `lib.control`, `lib.schema`, `lib.io`, the tools'
   control readers and the suites. The report lists every refusal, required key and check
   removed, one line each with its reason, and every one kept, one line each with the wrong
   number or obscure failure it prevents. What the coding agent cannot place goes in a third
   list and is kept until ruled.
   - **The `-dirty` refusal goes.** The commit and its `-dirty` mark stay written in every
     product, so the record is unchanged. Untracked files in `docs/` stop blocking runs.
   - **Writers** stop writing keys and attributes that no code reads. Old files keep theirs and
     still read.
   - A suite check that only exercised a removed refusal is retired. Its suite's count changes
     in `tests/run_regression.sh` in the same step, old and new count in the report.
3. **`tests/run_regression.sh`** lists the seven newer suites after `step08_1`, in numeric order:
   `step09_1` 3, `step10_1` 2, `step11_1` 5, `step11_2` 2, `step12_1` 4, `step13_1` 3, `step13_2`
   2, as amended by item 2.

**Acceptance.** The full regression, in chunks, since many suites change: every suite at its
reference count as amended. `step05_3` runs checks 4 to 7 for the first time since SPEC_11 Step 2.

The reviewing agent then cuts SPEC_00 §5, the §6 tables, §7.2 and §8 to the kept list (v0.25).

## 2. Ruled on the pre-execution report (22 October 2026)

The classification of its §1 and §2 stands. Then:

1. **Writers (its §4 question 1): the proposal.** Writers stop writing only what existed for a
   removed check: `role`, `composition_role`, `profile_or_run` on the input kinds,
   `thermo_instance` and the source profile set it governed, `p_b_location`, the closure comparison
   keys, and `observation_level_justification`. The records (`production_record`,
   `transfer_record`, `reduction_record`, the embedded namelist) are results, read by people and
   suites, and stay. So do the descriptive globals for now, required by nothing.
2. **Control keys (question 2): the proposal.** `role` and `p_b_location` leave the readers' tables
   and every committed build file, namelist and template, in this step.
3. **Not placed, its §3:**
   1. Path containment: one rule stays, that a run writes only inside its own run directory, so a
      run cannot overwrite registered data. The rest goes: manifest inputs inside the manifest's
      directory, run inputs under `inputs/`. Data may come from anywhere.
   2. The season on every file, and its range on files: goes. A step that needs a file's season
      refuses the file then, where the arithmetic needs it. The namelist's season and its range
      stay.
   3. `package_version`: goes. The installed version is recorded.
   4. Kind N's hash cross-check: goes.
   5. `casspian-new-run` over an existing directory: goes, provided `mkdir` is called so that it
      fails on an existing directory. A kept hypothesis must never be overwritten.
   6. `reject_physical_values`: goes. The kept type checks cover a number given for a path.
   7. The stated latitude global: goes. `refrac`'s `KeyError` names it.
4. **The untracked reference files in `docs/`** are committed in the acceptance commit, beside the
   papers already there, so products are stamped clean.

## 3. Revision history

- v0.1, 22 October 2026: first draft, three narrow fixes.
- v0.2, 22 October 2026: the author's rule, a lean pass over every refusal, check and required
  key.
- v0.3, 22 October 2026: accepted, the `-dirty` refusal's removal included.
- v0.4, 22 October 2026: §2 rules on the pre-execution report.
