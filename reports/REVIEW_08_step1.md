# REVIEW 08, Step 1. Kind W holds what the data hold

Reviewing agent, 7 October 2026. Report: `reports/REPORT_08_step1.md`. Specification:
`docs/specs/SPEC_08_Kind_W_Parts.md` v0.4 (§6 rulings on the reading, ruling 15 on regression).

## Verdict

**Accepted.** 6 of 6, every ruling of §6 implemented as written, no computed value moved.

## What was checked

- **The definition.** Kind W now requires the reference wind and its level, and nothing else of
  the wind. The sum identity and its dimension check are gone; the polar rule left `validate` and
  is applied by the three loaders after the presence of the total (`_admit_wind`). The output
  (`reports/step08_1/output.txt`) shows each refusal naming its variable and coming from the right
  place: the loaders for the total and the poles, the schema on read for the reference wind and
  the level.
- **Data files.** The three files of ruling 11 write and read back array-equal on latitudes short
  of the poles; a shear without its reference wind is refused on read.
- **The shear tool.** It refuses a source without the total and writes no shear part where its
  source has none; the closed form holds exactly. The report's reading of check 5 (an input without
  `u_shear_ms`, since one without `u_reference_ms` is no longer kind W) is right.
- **Figures.** Viewed: F9 with and without the source line, and `casspian-render` of the three data
  files (the right panel blank with its note, the shear profile, the one-level total as markers).
  Each draws the parts its file carries. The data files carry the reduction wind's own titles and
  `vertical_structure`, since they were cut from it; that is the test files, not the figures.
- **Accepted suites.** `step04_0` check 4 retired (15 to 14), `step7` check 2 restated as ruled,
  seven suites changed in wording only, conditions unchanged.
- **Regression.** The full set ran before ruling 15 reached the coding agent: 31 of 31 at their
  reference counts, registered files restored equal. Nothing further is to be run.

## Applied at acceptance

- **SPEC_00 v0.23**, §6.6 amended as SPEC_08's Appendix, restated by §6, requires: the opening
  paragraph, the variable table (required and optional marked), the load-time rules replaced, the
  coverage rule folded into "the model does not extrapolate", the header and the revision table.
- **SPEC_08 v0.5**: Step 1 accepted; the specification closes at the acceptance commit.

## Order of work

1. The author's go; the acceptance commit: the code, the suites, `tests/step08_1/`, this review,
   the report, SPEC_08 v0.5 and SPEC_00 v0.23. No AI attribution.
2. No sweep: no registered file changed.
3. `STATE.md`: a SPEC_08 section, Step 1 accepted, SPEC_08 closed.
