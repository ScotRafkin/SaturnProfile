# REVIEW 21, Step 1, check 5. The reference swap, and closing SPEC_21

Reviewing agent, 11 November 2026. Report: `reports/REPORT_21_step1_check5_diagnosis.md` §9.
Specification: `docs/specs/SPEC_21_General_Error_Budget.md`, closed at v0.7.

**Check 5 accepted on the materiality rule. SPEC_21 closes; check 5 is not carried.**

- **Matching scales fixed the divergence.** With every profile averaged into the forward grid's layers,
  the swap fit's request fell from 9 K to 0.4 K in one pass and the full swap fit converged in 4 and 2
  passes; every solve lowers the cost at every target. No trust region was needed.
- **Hypothesis 4 refit moves by at most 0.75 K** (3.0 S at 1 mbar), under the 1 K bound; elsewhere
  under 0.25 K. The refit wind is hypothesis 4.
- **The reference does not matter, to within what matters.** Swapping the reference changes the
  delivered profiles by 0.2 to 0.9 K RMS; at the six standard levels every difference is within 0.5 K
  except 3.0 S at 10 mbar, 0.68 K. The largest single-level differences (up to 3.45 K) are at full
  resolution, which the fit no longer constrains; the 0.5 K bound I set was tighter than the
  materiality rule (about 1 K), and compared at a finer scale than the fit uses. Judged as the scales
  rule and materiality now require, check 5 passes. The mirror-image sign between the 36.3 N and
  31.2 S directions says the residual is a small difference in the wind between them, not in a profile.
- **For SPEC_22, noted, not chased:** reference-independence checks are judged on layer means with the
  materiality bound; the Lindal-referenced fit's radio terms are 4 to 7 times the swap fit's, though
  both sit well inside their errors (0.04 to 0.33 of expected).

**Closing commit,** without further review: the layer averaging and cost breakdown in the code; the
refit `fits/voyager_fit_2/voyager_fit_2_wind.nc`; `--refit` in `tests/step21_1/`; a correction note at
the head of REPORT_21_step1 §7 pointing to the diagnosis report (the top edge was a symptom of the
swap's divergence); the diagnosis report and this review with the other SPEC_21 reports and reviews;
SPEC_21 v0.7, the design note v0.9, `STATE.md` with the hash. Rebuild clean; push. No regression;
Lindal's files untouched.
