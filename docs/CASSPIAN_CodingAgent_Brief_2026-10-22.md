# CASSPIAN coding agent brief, 22 October 2026

For the incoming coding agent, from the reviewing agent. Author of record: S. Rafkin (SwRI), who
approves everything. Read `docs/CASSPIAN_CodingAgent_Handoff_2026-10-22.md` first; this brief
updates it.

## 1. What changed since that handoff

- **Its §4 is superseded by SPEC_R1** (`docs/specs/SPEC_R1_Regression_Fixes.md`, v0.3, accepted).
  The author's rule: a refusal, test check or required key stays only if without it the code would
  compute a wrong number without saying so, or fail with no clear message. Everything else goes.
  Removing it is the fix, not loosening it. Much more data is coming, and nothing may accept only
  one kind of source or one file name.
- **The closure item.** The handoff's recommendation (add names to `CLOSURE_DROPPED_ATTRIBUTES`)
  is withdrawn. The whole content comparison goes.
- **The directory listing.** `step9` check 6 holds its own expected list; it does not read
  SPEC_00. The check goes. SPEC_00 v0.25 §2.2 no longer names files.
- **The `-dirty` refusal goes** (author, 22 October). Products still record the commit with its
  `-dirty` mark. Until SPEC_R1 lands, the untracked documents in `docs/` still need
  `.git/info/exclude` for any run, as before.

## 2. The work, in order

1. **SPEC_R1, Step 1.** One report, `reports/REPORT_R1_step1.md`, with three lists, one line per
   item: removed (with the reason), kept (with the wrong number or obscure failure it prevents),
   and not placed (kept until ruled). Retired suite checks with old and new counts. The full
   regression in chunks; logs in `reports/regression/`. Then the review, the author's go, commit
   and push.
2. **SPEC_14, Step 1** (`docs/specs/SPEC_14_Reversed_Shear.md`, v0.2, accepted), after SPEC_R1 is
   pushed. Case (iii) and the run `lindal_iris_v_sl_iii_110`. No regression. Report
   `reports/REPORT_14_step1.md`.

Both specifications are accepted, and their sections in `docs/specs/STATE.md` are open. The step
rows are yours to fill.

## 3. Unchanged

The handoff's §5 and §6: practice, style, commits, latitudes given planetographic with
planetocentric in parentheses. No dashes as parentheticals in prose.
