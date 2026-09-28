# CASSPIAN coding agent handoff, 28 September 2026

For the incoming coding agent. Author of record: S. Rafkin (SwRI). This replaces
`docs/CASSPIAN_CodingAgent_Handoff_2026-09-27.md` as the starting document; that one and
`docs/CASSPIAN_CodingAgent_Handoff_2026-09-17.md` remain as the record of how the previous
agents were started.

## 1. Where the build stands

SPEC_04, the transfer, is closed at the Step 5 acceptance commit `e71d59e` (28 September 2026).
Every specification through SPEC_04 is closed: SPEC_00 v0.22 (architecture and data files),
SPEC_01 v0.28 (the Lindal tool chain), SPEC_02 v0.11 (the reduction and diagnostics), SPEC_03
v0.15 (the forward production and the hydrostatic closure), SPEC_04 v0.22 (the transfer).
`docs/specs/STATE.md` has one row per step with the acceptance commits. The forward model runs
end to end: `casspian-forward forward/lindal_transfer/lindal_transfer.toml` transfers the Lindal
refractivity to 10° N under the closure inputs in about 70 to 100 s and writes a kind `profile`
product with figures F5 to F9. The registered anchor is `occul_data/lindal/lindal_refractivity.nc`
(SHA-256 `920304440ee4...`), the registered closure product `forward/lindal_closure/output/
lindal_closure_profile.nc` (`fcc2e2c4ae71...`); neither is in git, both are rebuilt from control
files.

What comes next, in the author's order: first the regression on the author's Linux clone
(your first task, below; the clone itself is done and reproduces the test case); then SPEC_05, the forward-inputs tools (a wind tool with parameterized vertical
shear, a composition tool with structure in pressure and latitude, a diagnostic decomposition of
a given wind) and a set of named transfer experiments under other assumed states; then the
performance step of SPEC_04 decision R; then the combination specification (A35, the posterior
wind drawn on F9 beside the hypothesis, the reference-surface constant); then the end-to-end
tests and the Monte Carlo wrapper. The author writes each specification before the work starts;
you do not begin SPEC_05 until it exists.

## 2. How the work is done

Read `docs/specs/SPEC_00_Architecture_and_Data_Files.md` §0 to §2 and SPEC_04 §0 first; they
hold the protocol and it has not changed: one step at a time; the specification's status line
and `STATE.md` read before any step; an acceptance script per step under `reports/step0N_M/`
with its printed output; `reports/REPORT_0N_stepM.md` filed by you, `reports/REVIEW_0N_stepM.md`
filed by the reviewing agent; nothing committed until the review says accepted and the author
says go; then the author's documents commit, the acceptance commit, the push. Commits carry no
attribution of any kind to an AI, in the message or the code.

The regression rule (SPEC_04 §0 at v0.22, the author's direction): a change reruns what it can
reach. A change that cannot move a number reruns no calculation; a figure change reruns the
figure checks of the steps it draws for; a change that can move a result reruns every accepted
suite downstream of it, and a change to `lib` is the full set. The report names the suites rerun
and, in one sentence, why the others were not.

What is checked in code (SPEC_04 §0, decision N): the schema, namelist typo guards,
extrapolation, and named numerical failures. Everything else is recorded and warned, never
refused. Do not add consistency checks beyond those; the author has ruled on this twice.

Style, everywhere: no em dashes or en dashes in any document, message, comment or commit (a
minus sign in a number is fine); American spellings; no paragraph or line numbers when
referring to text; every measured number quoted with its bound. Verify, do not assert: a claim
in a report is a number in the printed output. The work is SwRI funded and proprietary; nothing
leaves the repository.

The manuscript `CASSPIAN_AtmosphericModel_Draft9_2.docx` is the draft of record for the equation
labels (A1 to A40, B1 to B8) the specifications cite; it is not in the repository and is not to
be added to it.

## 3. Your first task: the acceptance scripts and the regression on the clone

`docs/RUNBOOK.md` v0.2 is the sequence from `git clone` to the transfer product. The author has
already run its sections 2 to 5 himself on a clean clone of `e71d59e` on his Linux machine
(CentOS 7, glibc 2.17, no C++ compiler, Python 3.12.14 from `uv`, numpy 2.2.6, scipy 1.16.3):
the anchor, the closure product and the transfer product reproduce the registered values to the
last printed digit, and the transfer ran in 91 s. What remains is yours:

1. **Commit the documents and move the acceptance scripts into `tests/`**, on `main` on the
   Windows machine. The author has approved this and the step-by-step instructions are in
   `Claude outputs/INSTRUCTIONS_git_2026-09-29.md` (outside git by design): the runbook and this
   handoff committed first, then the acceptance scripts moved from `reports/step*/` to a
   committed `tests/` directory with everything they write left in the ignored
   `reports/step*/` (SPEC_00 v0.22 rows [Y] and [Z]), `matplotlib` declared, `*.sh` given LF
   endings, and SPEC_00 v0.22 committed with them.
2. **The regression on the clone.** The author pulls on the Linux machine and runs section 6 of
   the runbook (about seven hours; you cannot reach that machine, so the instructions to him are
   the runbook's). He passes you `reports/regression/regression.txt` and the per-suite logs; a
   suite that fails for a platform reason (a path, a fixture, a tolerance at the last digit) is
   diagnosed from its log and fixed on `main`, and the failing suite alone is rerun.
3. **File `reports/REPORT_05_step0.md`** as the runbook's §7 says, with the author's measured
   values from the runbook as the clone's sections 3 to 5, and bring the runbook to v0.3 if the
   regression taught it anything. Do not commit the report until reviewed.

The author will be developing tools for further tests under SPEC_05 after this; the runbook and
the committed scripts are what make that possible.

## 4. Things the last agents learned that you will need

- Products record `casspian_git_commit`; a product written from a dirty tree ends `-dirty` and
  `refrac` and `forward` refuse it as an input. Rebuild on a clean tree; the regression script
  sets registered products aside and restores them after every suite for this reason.
- Close every netCDF read before opening another; an unclosed HDF5 handle segfaulted at Step 4.
- Do not edit a script while a shell is running it; bash reads by byte offset.
- Acceptance suites are pinned to the spacings their bounds were set for (SPEC_04 §0);
  `step04_4` runs at 5,000 m²/s² and takes four and a half hours, `step04_5` at 50,000.
- The reviewing agent measures expected values independently of the repository; when your
  number and the specification's disagree, say so with your measurement and ask for a ruling
  rather than fitting to the stated value. The last three specifications' errors were mostly
  the specification's, and every one was settled by a measurement.
