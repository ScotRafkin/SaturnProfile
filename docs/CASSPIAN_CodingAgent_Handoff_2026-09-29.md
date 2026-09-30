# CASSPIAN coding agent handoff, 29 September 2026: SPEC_05, the shear tool and experiments

For the incoming coding agent. Author of record: S. Rafkin (SwRI). This follows
`docs/CASSPIAN_CodingAgent_Handoff_2026-09-28.md`, which remains the record of the SPEC_05
Step 0 work (the scripts moved to `tests/`, the regression on the author's Linux clone,
`REPORT_05_step0`).

## 1. Before you start

Do not begin until the author says so. The precondition is that SPEC_05 Step 0 is closed:
`reports/REPORT_05_step0.md` filed and reviewed, the regression of record resolved (every suite at
its reference count, or each failure diagnosed, fixed on `main` and rerun), and `main` committed
and pushed. Start from a clean `main` at that commit.

## 2. What to read, in order

1. This handoff.
2. `docs/specs/SPEC_00_Architecture_and_Data_Files.md` §0 to §2 and `docs/specs/SPEC_04_Transfer.md`
   §0: the protocol, the regression rule (rerun what a change can reach), and decision N (what the
   code checks and what it only records).
3. `docs/specs/SPEC_05_Shear_Experiments.md` v0.3, all of it. It is the only instruction. The
   design note and the two reviews it was built from are background and are not in the repository;
   where anything else you read disagrees with SPEC_05, SPEC_05 governs.
4. `docs/specs/STATE.md`. Add a SPEC_05 table with rows for Steps 0 to 3 (Step 0 as reported in
   `REPORT_05_step0`); the header line above the table is the reviewing agent's.

## 3. Your first deliverable: a pre-execution review

Before writing code, read SPEC_05 against the repository and file
`reports/REPORT_05_preexecution.md`: every place where the specification's description of the code
(file names, functions, attributes, the closure comparison, the run build driver, the schema's kind
W checks, F9) does not match what is on disk; every step whose acceptance you cannot write as
stated; every choice the specification leaves to you that you think the author should make. Short
is fine. SPEC_04's pre-execution review found real errors in the specification and saved a step's
worth of rework. Then wait for the rulings.

## 4. Things SPEC_05 depends on that are easy to get wrong

- **The model does not change.** Nothing under `src/casspian/lib`, `src/casspian/refrac` or
  `src/casspian/forward` is edited in this specification. The code that changes is new
  (`tools/wind/shear.py`, `tools/run/new_run.py`), the run build driver (`tools/run/run_inputs.py`),
  F9's panel (`tools/plots/figures_profile.py`), `pyproject.toml`, the two existing runs' build
  files, and the runbook. If you find you need a model change, stop and say why.
- **Step 0 is done; the work starts at Step 1.** Do not number any shear work as Step 0.
- **`construct` and `build` are separate** (SPEC_05 §1.5). `construct` touches no file; the Monte
  Carlo driver will call it later without writing. Keep it that way even where writing inside it
  would be shorter.
- **General code.** Nothing in the tool assumes Lindal, a 1 bar reference, a grid node at `p_s`,
  or a particular wind grid. The Lindal values in SPEC_05 are the test case's.
- **The Step 2 migration must not move a number.** The closure product stays equal to the
  registered one and the transfer product to the accepted Step 5 one. An accepted suite whose
  expectation changes only because of the new file layout is updated in the same step and named in
  the report with its old and new expectation (SPEC_05 §3); none is edited silently or left
  failing.
- **Run 7 may not complete.** A named failure (isobars crossing, the outer loop not converging)
  is a result. Report it; do not tune numerics, spacings or tolerances to get it through.
- **Bounds come from measurements.** Where SPEC_05 says a bound is set from the measured value,
  report the measurement and let the review set it. Where your number and the specification's
  disagree, give your measurement and ask; do not fit to the stated value.
- **The reviewing agent's independent values** for runs 5 and 9 arrive before Step 3 is reviewed.
  Steps 1 and 2 do not wait for them.

## 5. Carried from earlier agents

- Products from a dirty tree end `-dirty` and are refused as inputs; rebuild on a clean tree.
- Close every netCDF read before opening another (an unclosed HDF5 handle segfaulted at SPEC_04
  Step 4).
- Do not edit a script while a shell is running it.
- Accepted suites are pinned to the spacings their bounds were set for.
- Style everywhere: no em or en dashes, American spellings, every measured number with its bound,
  no attribution of any kind to an AI in commits or code. The work is SwRI funded and proprietary;
  nothing leaves the repository.
- Order after each accepted step: the author's documents commit, the acceptance commit, the push,
  `STATE.md`.
