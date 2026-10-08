# REVIEW 13, Step 2. The runs

Reviewing agent, 8 October 2026. Report: `reports/REPORT_13_step2.md`. Specification:
`docs/specs/SPEC_13_SL_Wind.md` v0.4.

**Accepted once Step 3 (SPEC_13 §2a) is done and the author has viewed its figures. SPEC_13 then
closes.**

- **The gradient from the fitted temperatures closes the round trip:** -0.16 K at 110 mbar, the
  linearization, against -1.50 K with the slope. The model's remaining part is the equatorial
  bridge, 0.6 to 0.9 K at 290 and 730 mbar.
- **Against the egress the derivative does worse at 110 mbar,** +5.3 K against +3.9 K. It carries
  the IRIS change in full, and the radio does not show it.
- **The Sanchez-Lavega cloud wind changes the delivered profile by at most 0.14 K** at the three
  levels. The shear changes with it near the equator, as ruling 3 expects.
- **The correction to SPEC_11's diagnosis is accepted** and recorded in SPEC_11 v0.5 and
  REVIEW_11_step2.
- **SPEC_14 starts from `lindal_iris_v_sl_ii_110`** once it is rerun with the smooth window.

Order of work:
1. Step 3, then the author's view.
2. Commit and push, with no regression:
   - the tools;
   - the data property files and the note;
   - the runs' control files;
   - `tests/step13_1/` and `tests/step13_2/`;
   - the reports, these reviews, SPEC_11 v0.5 and REVIEW_11_step2;
   - the SPEC_13 row in `STATE.md`.
