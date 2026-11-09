# CASSPIAN coding agent brief, 9 November 2026

From the reviewing agent. Author of record: S. Rafkin (SwRI), who approves everything. Supersedes
`CASSPIAN_CodingAgent_Brief_2026-10-22.md`. Read with `CASSPIAN_CodingAgent_Handoff_2026-10-22.md`
§5 and §6 (practice and style), which stand except where this brief says otherwise.

## 1. What changed since 22 October

- **SPEC_R1 landed.** The `-dirty` refusal is gone; products still record the commit with its mark
  and are rebuilt clean at each acceptance commit. A refusal, check or required key stays only if
  without it the code would compute a wrong number silently or fail obscurely.
- **Materiality (author, 9 November).** A change worth less than about 1 K in delivered temperature
  or about 5 percent in wind at probe levels gets a comment, not a specification. Keep testing that
  the algorithms close and converge; do not chase answers below that size. Recommendations in
  reports should say where they stand against it.
- **Figures start with absolute profiles:** observed with its error band, and delivered, before any
  difference or residual plot.
- **No new keys** unless code reads the value; bases and sources go in comments.
- **Lindal is untouched** by other sources' work. When a shared tool changes, rebuild Lindal's files
  in a copy; equal values mean no full regression.
- **Composition.** Each source is reduced with its own composition and true molecular masses; each
  run and each fit puts all its inputs on one composition.

## 2. The work, in order

1. **SPEC_19's acceptance commit** (`reports/REVIEW_19_step1.md`, rerun section; SPEC_19 v0.5 §5).
   Add the 0.3 K radio floor; refit once; rebuild the three runs; add the absolute profile figure
   (SPEC_19 v0.5 §5 item 6); run `survey_truths.py` and
   `truth_and_noise.py` once as committed. Report the passes, the group ratios and the 110 mbar
   value at 31.2 S. If that value moves by more than 1 K from +3.34 K, stop and report. Otherwise
   commit as the review lists, rebuild clean, put the hash and the three numbers in `STATE.md`,
   push.
2. **SPEC_20, cloud-tracked winds,** when it is accepted. Not yet drafted.

## 3. Known

- **The wind fit** is `src/casspian/tools/wind/estimation/` with `casspian-wind-fit`; SPEC_01's point
  fit stays in `tools/wind/fit.py`. Fits live under `fits/`; their control and wind files are
  committed.
- **The equatorial band** of `voyager_fit_1` is the prior's, not the data's (SPEC_19 v0.5 §5 item 2).
  Do not tune it; cloud winds will constrain it.
