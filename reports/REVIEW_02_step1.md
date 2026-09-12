# REVIEW 02, Step 1. `lib.control`: the reduction manifest and input loading

Review of `reports/REPORT_02_step1.md` against SPEC_02 v0.2 Step 1. Reviewer: S. Rafkin,
12 September 2026. **Disposition: Step 1 accepted on its own merits; SPEC_02 remains a draft.** The
instruction the coding agent took as acceptance of the specification was not one; the header's
"draft for author markup, no step to begin" stood and stands. The work done is sound and is
kept. SPEC_02 goes to v0.3 for the three markup findings below, with the status line restated.
Commit Step 1. **Do not begin Step 2 until the author accepts the draft.**

## Rulings

1. **Step directories.** `reports/step02_<N>/`, as proposed; SPEC_01's directories stay and its
   suites keep running. In the spec.

2. **The composition term is 0.0233, not 0.0227.** The report's arithmetic is right and mine
   was not; 0.03 × (136 − 35) / 129.94 = 0.02332, constant across levels because ammonia
   dilutes both shares alike. Corrected in Step 3's acceptance.

3. **N at 1 bar.** Correct again: Table I tabulates 10.9 ppm of NH3 at 1000 mbar, so the
   product writes 2.598560e-4, not the dry 2.598589e-4. Step 3 now quotes the tabulated-NH3
   value at 1 bar and moves the dry-ℛ̄ check to 794.33 mbar, where NH3 is exactly zero by the
   v0.17 rule. The self-correction in finding 4 is noted and appreciated; that is how a check
   is supposed to be reported.

## On the decisions

All six stand. The physical-value test (unknown key with a unit suffix and a numeric value)
is the right formalization and reuses the one list of unit suffixes. The tolerance staying in
degrees in the manifest record, with one conversion at the `refrac` boundary, is exactly the
rule. Deleting `tools/gravity/control.py` rather than shimming it is right; two readers would
have drifted. The prefix rule on manifest paths is stricter than §8 lists and is the intended
reading; SPEC_00 will say so at its next revision. Refusing `anchor_rule = "latitude"` at parse
time because the manifest has no key for the latitude is correct; if a profile ever needs it,
the key is added to §7.1 first.

## Actions

- Commit Step 1. Nothing to rebuild.
- Stop. Step 2 waits for the author's acceptance of SPEC_02, which will be recorded in its
  status line and in `STATE.md`.
