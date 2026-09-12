# CASSPIAN build state

One row per step of `SPEC_01_Lindal_Tool_Chain.md`. Status is one of `not started`,
`in progress`, `reported`, `accepted`. A step becomes `accepted` only after its report has been
reviewed against the specification.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 0 | Repository reset and skeleton | accepted | `reports/REPORT_01_step0.md`, `reports/REVIEW_01_step0.md` | 2026-09-10 |
| 1 | `.gitattributes`; `lib.constants`, `lib.schema`, `lib.io` | accepted | `reports/REPORT_01_step1.md`, `reports/REVIEW_01_step1.md` | 2026-09-10 |
| 2 | `tools/lindal/build_raw.py`, the raw bundle | accepted | `reports/REPORT_01_step2.md`, `reports/REVIEW_01_step2.md` | 2026-09-10 |
| 3 | `tools/gravity/`, kinds G and R | accepted | `reports/REPORT_01_step3.md`, `reports/REVIEW_01_step3.md` | 2026-09-10 |
| 4 | `lib.gravity` | accepted | `reports/REPORT_01_step4.md`, `reports/REVIEW_01_step4.md` | 2026-09-10 |
| 5 | `lib.geoid` | accepted with changes | `reports/REPORT_01_step5.md`, `reports/REVIEW_01_step5.md` | 2026-09-10 |
| 6 | `lib.latitude` | accepted | `reports/REPORT_01_step6.md`, `reports/REVIEW_01_step6.md` | 2026-09-11 |
| 7 | `tools/wind/`, kind W | accepted | `reports/REPORT_01_step7.md`, `reports/REVIEW_01_step7.md` | 2026-09-11 |
| 8 | `tools/composition/`, kind C | accepted | `reports/REPORT_01_step8.md`, `reports/REVIEW_01_step8.md` | 2026-09-11 |
| 9 | `tools/lindal/build_inputs.py`, kinds T and D, the manifest | accepted | `reports/REPORT_01_step9.md`, `reports/REVIEW_01_step9.md` | 2026-09-11 |

## SPEC_02, `refrac` and the standard diagnostics

SPEC_02 is a draft awaiting the author's markup. Step 1 was built ahead of acceptance and is
accepted on its own; no further step begins until the draft is accepted.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 1 | `lib.control`: reduction manifest and input loading | accepted | `reports/REPORT_02_step1.md`, `reports/REVIEW_02_step1.md` | 2026-09-12 |
| 2 | `refrac/anchor.py`: frozen `phi_c` and `r0` | not started | | |
| 3 | `lib.reduction`, `refrac/reduce.py`: B3.1, B1, B3.3 | not started | | |
| 4 | `refrac/product.py`, `casspian-refrac`: kind N | not started | | |
| 5 | `tools/plots/`, `casspian-plots`: standard diagnostics | not started | | |
| 6 | `refrac/sensitivity.py`: sensitivities and wind Monte Carlo | not started | | |
