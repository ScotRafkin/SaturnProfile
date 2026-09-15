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

SPEC_02 accepted by the author at v0.4 on 12 September 2026. Steps 2 to 5 proceed in order. Current: SPEC_00 v0.15, SPEC_01 v0.20, SPEC_02 v0.9. SPEC_02 closes at the Step 6 acceptance commit and sweep. Next specification: SPEC_03, the forward model round trip.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 1 | `lib.control`: reduction manifest and input loading | accepted | `reports/REPORT_02_step1.md`, `reports/REVIEW_02_step1.md` | 2026-09-12 |
| 2 | `refrac/anchor.py`: frozen `phi_c` and `r0` | accepted | `reports/REPORT_02_step2.md`, `reports/REVIEW_02_step2.md` | 2026-09-12 |
| 3 | `lib.reduction`, `refrac/reduce.py`: B3.1, B1, B3.3 | accepted with changes, applied (SPEC_02 v0.6; SPEC_01 Step 8 amended first) | `reports/REPORT_02_step3.md`, `reports/REVIEW_02_step3.md` | 2026-09-12 |
| 4 | `refrac/product.py`, `casspian-refrac`: kind N | accepted with one change, applied (SPEC_02 v0.7: anchor-rule spread by full rerun) | `reports/REPORT_02_step4.md`, `reports/REVIEW_02_step4.md` | 2026-09-12 |
| 5 | `tools/plots/`, `casspian-plots`: standard diagnostics | accepted with two figure changes, applied (SPEC_02 v0.8) | `reports/REPORT_02_step5.md`, `reports/REVIEW_02_step5.md` | 2026-09-12 |
| 6 | equatorial anchoring: `equatorial_radius` rule, manifest `[diagnostics]`, comparison report | accepted (`2cf4457`); Lindal default anchor rule `equatorial_radius` (SPEC_02 v0.9 decision 9); SPEC_02 closed | `reports/REPORT_02_step6.md`, `reports/REVIEW_02_step6.md` | 2026-09-13 |

## SPEC_03, the forward production and its hydrostatic closure

SPEC_03 accepted by the author at v0.2 on 14 September 2026 (now v0.8). Step 0 accepted and swept; Steps 1 to 4 proceed in order, each after the review of the one before. Current: SPEC_00 v0.16, SPEC_01 v0.24, SPEC_02 v0.10, SPEC_03 v0.8 (REVIEW_03_step2 applied). Next specification: SPEC_04, the transfer (build); SPEC_05, the end-to-end tests.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 0 | pressure grid and Eq. B1 tilt amendments; whole reduction chain rebuilt and swept | accepted (`2149b64`; 13 of 13 as restated at SPEC_03 v0.6; swept, hashes in REPORT_03_step0 section 7) | `reports/REPORT_03_step0.md`, `reports/REVIEW_03_step0.md` | 2026-09-14 |
| 1 | `lib.geopotential` | accepted (`b3efc3d`; 9 of 9; finding 1 restated the Step 4 negative control at SPEC_03 v0.7; no product changed) | `reports/REPORT_03_step1.md`, `reports/REVIEW_03_step1.md` | 2026-09-14 |
| 2 | `lib.hydrostatic` | reported (8 of 8; finding 1 on the linear-T orientation) | `reports/REPORT_03_step2.md` | 2026-09-14 |
| 3 | run directory and inputs (`casspian-run-inputs`), `read_run_namelist`, kind `profile`, `input_hashes` warning for every derived kind | not started | | |
| 4 | `forward/production.py`, `casspian-forward`, the closure run, F5 and F6 | not started | | |
