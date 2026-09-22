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

SPEC_03 accepted by the author at v0.2 on 14 September 2026 (now v0.14). Steps 0 to 4 accepted and swept; **SPEC_03 closed** at the Step 4 record commit `545b385` on 16 September 2026. Current: SPEC_00 v0.20, SPEC_01 v0.28, SPEC_02 v0.11, SPEC_03 v0.14 (Step 3 reviewed and accepted, rulings in SPEC_03 section 8; the whole chain dated 1981-08-26, the occultation, SPEC_00 v0.19 and SPEC_01 v0.28; the seasonal design note is docs/CASSPIAN_Seasonal_Design_Note.md). Next specification: SPEC_04, the transfer (build); SPEC_05, the end-to-end tests.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 0 | pressure grid and Eq. B1 tilt amendments; whole reduction chain rebuilt and swept | accepted (`2149b64`; 13 of 13 as restated at SPEC_03 v0.6; swept, hashes in REPORT_03_step0 section 7) | `reports/REPORT_03_step0.md`, `reports/REVIEW_03_step0.md` | 2026-09-14 |
| 1 | `lib.geopotential` | accepted (`b3efc3d`; 9 of 9; finding 1 restated the Step 4 negative control at SPEC_03 v0.7; no product changed) | `reports/REPORT_03_step1.md`, `reports/REVIEW_03_step1.md` | 2026-09-14 |
| 2 | `lib.hydrostatic` | accepted (`994c787`; 8 of 8; finding 1 restated the linear-T orientation at SPEC_03 v0.8; no product changed) | `reports/REPORT_03_step2.md`, `reports/REVIEW_03_step2.md` | 2026-09-15 |
| 3 | season identifier across the chain (rebuild); run directory and inputs (`casspian-run-inputs`), `read_run_namelist`, kind `profile`, `input_hashes` warning for every derived kind | accepted (`6ae113e`; 16 of 16 after the author's decision on finding 1, every file of the chain dated 1981-08-26; no value changed in any product; swept, hashes in REPORT_03_step3 section 7; the suites of section 5 pass on the swept products) | `reports/REPORT_03_step3.md`, `reports/REVIEW_03_step3.md` | 2026-09-15 |
| 4 | `forward/production.py`, `casspian-forward`, the closure run, F5 and F6 | accepted (`1504056`; 13 of 13 on the SPEC_03 v0.13 bins, the closure namelist dated 1981-08-26; closure product clean at `1504056`, SHA-256 `e6693173...`, F5 and F6 attached; record `545b385`, hash in REPORT_03_step4 section 8; the suites of section 5 pass on the clean tree); SPEC_03 closed | `reports/REPORT_03_step4.md`, `reports/REVIEW_03_step4.md` | 2026-09-16 |

## SPEC_04, the transfer

SPEC_04 accepted by the author at v0.6 on 21 September 2026 (now v0.10). Steps 0 to 5 proceed in order, each after the review of the one before; the coding agent's pre-execution review of v0.5 is ruled on in SPEC_04 section 10, its record is `reports/REPORT_04_preexecution.md`, and the rulings on the step reports are in sections 11 (Step 0) and 12 (Step 1). Current: SPEC_00 v0.20, SPEC_01 v0.28, SPEC_02 v0.11, SPEC_03 v0.14, SPEC_04 v0.10. The manuscript draft of record for the equation labels is `CASSPIAN_AtmosphericModel_Draft9_2.docx`, outside the repository. Next specification after SPEC_04: the combination (A35, the posterior wind, the reference-surface constant), then SPEC_05, the end-to-end tests.

| Step | Deliverable | Status | Report | Date |
|---|---|---|---|---|
| 0 | kind W in three parts, kind C as one structure on (level, latitude), `forward/lindal_transfer/`, the transfer namelist for M anchors, the loader with the anchor object and the propagation hook; the chain rebuilt and swept | accepted (`1c9b310`; 15 of 15 under REVIEW_04_step0 and SPEC_04 v0.9, decisions N and O and the §11 rulings applied and the whole step rebuilt under them; eight findings, finding 1 on the uncertainty column of deliverable 6 and findings 6 to 8 on three accepted acceptances reading retired forms; no value changed in any product; full regression passes; swept, all seventeen products clean at `1c9b310`, hashes in REPORT_04_step0 section 6, record `f38c58b`; the seven deferred suites pass on the swept products at their reference counts) | `reports/REPORT_04_step0.md`, `reports/REVIEW_04_step0.md` | 2026-09-21 |
| 1 | `lib.geoid` `through_anchor`, `lib.windfield`, the anchors' geopotential under the run's wind (`produce` with the wind along the column) | accepted (`780a2de`; 13 of 13 under REVIEW_04_step1 and SPEC_04 v0.10, the §12 rulings applied; deliverable 1 reproduces every stated radius and is 6.4 times cheaper than the `equatorial_radius` rule, the closure production with the wind array bit-identical to the registered product; five findings, finding 2 on the cylinder-extended wind, whose v0.9 values the review withdrew and whose construction decision P places at Step 2 with checks 9 and 14 moved there unloosened, and findings 4 and 5 on `step04_0` checks 0 and 6, both recorded for ratification; full regression passes; swept, the closure product alone rebuilt clean at `780a2de` with `u_column_ms`, hash in REPORT_04_step1 section 6, sweep record `3dcd274`; the four suites that read it pass on it at their reference counts, recorded in REPORT_04_step1 section 5c) | `reports/REPORT_04_step1.md`, `reports/REVIEW_04_step1.md` | 2026-09-22 |
| 2 | `lib.mesh`, the columns, `z_lv` | not started | | |
| 3 | `lib.kernel` | not started | | |
| 4 | `forward.transfer`, `forward.estimate`, the outer loop, the M = 2 identity test | not started | | |
| 5 | production at the target, altitude and datum, kind `profile` in transfer mode, `casspian-forward` transfer mode, F5, F6, F7 | not started | | |
