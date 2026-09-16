# REVIEW 03, Step 4. `forward.production`, `casspian-forward`, the closure run, the product, F5 and F6

Review of `reports/REPORT_03_step4.md` against SPEC_03 v0.12 Step 4. Reviewer: S. Rafkin,
16 September 2026. **Disposition: accepted.** The one failing check fails on a wording defect
in the specification, not on the closure (finding 1); the bins are restated in pressure at
SPEC_03 v0.13 and the check passes on the numbers as measured. Two small changes are applied
before the acceptance commit (findings 3 and 4). SPEC_03 closes at the Step 4 record.

## What was verified independently

The closure was recomputed from the product's own embedded anchor with the reviewing agent's
code of 14 September, written before the specification and sharing nothing with the
repository: the effective gravity from the embedded harmonics, rotation and `u(phi_c)`, `Φ` by
the trapezoid along the field line with the gauge at level 29, the log-linear layer integral
from the anchor's top pressure, `T` by Eq. B6. Against the product: `pressure_Pa` and
`temperature_K` agree to 5.6e-16 relative at every level and `geopotential_m2s2` to 1.2e-9
m²/s². The residual statistics computed from the product are the report's to every printed
figure: 9.740e-4 above 2 mbar, 2.119e-3 at 6.31 mbar, 2.607e-3 between 10 and 100 mbar,
1.714e-3 below 100 mbar excluding the bottom row, −3.281e-3 at the bottom row, mean below 10
mbar −4.895e-4; the temperature residual equals the pressure residual to 4.4e-16. The two
negative controls, computed independently, give −5.072e-3 and −9.632e-3 for the mean below 10
mbar and −7.842e-3 and −1.240e-2 at the bottom row, the report's values. The envelope band
recomputed from the embedded thermo table (50 m over `R T / (M g)`, plus the bottom row's
pressure term) puts 21 of 66 levels outside it, as check 9 counts. `p[0]` equals `p_b` equals
the anchor's top tabulated pressure, 19.952623149688797 Pa. The product's SHA-256 is the
report's. The product carries `mode`, the season and the anchor's season, six NaN companions
with the terms unstated, relative `input_hashes`, the point-latitude marker, `psi_deg`
5.494432°, and the four groups deliverable 3 asks for. F5 and F6 were viewed.

## Rulings on the findings

1. **"The top decade."** Correct, and the fault is the specification's. The expected value
   "within ±9.8e-4 in the top decade" was measured on the ten top levels, above 2 mbar; the
   acceptance's "from 10 mbar down" then left 2 to 10 mbar under no bound, and the partition
   reading put a 1.5e-3 bound on levels the expected values never covered. At v0.13 the
   statistics and the bounds are stated in pressure: above 2 mbar (1.5e-3, measured 9.74e-4),
   from 2 mbar down excluding the bottom row (3e-3, measured 2.607e-3), the bottom row (4e-3,
   measured 3.281e-3), the mean below 10 mbar (1.5e-3, measured 4.9e-4). Every level is under
   one bound, no bound moved, and check 5 passes on the numbers already measured. The choice
   to run the check as decided before the run and report the failure rather than re-read the
   words afterward was the right one.
2. **F5.** Left for the author's figure pass, already an open item; nothing restyled here.
3. **SPEC_00 §7.2 and the resolved namelist.** Correct; the v0.16 words "absolute paths" and
   "into `output/`" contradicted v0.18 §5 and are withdrawn at SPEC_00 v0.20. Decision 6 (the
   resolved form as an attribute of the `namelist` group, relative paths and hashes, no
   separate file) is the rule.
4. **The product's `epoch`.** Adopted: the closure namelist declares `date = "1981-08-26"`, so
   the closure product carries the occultation date like every other file of the run. Applied
   before the acceptance commit; the namelist's hash and the product change, nothing else.
5. **21 of 66 levels outside the band.** Recorded. One reading of F6 for the manuscript, not
   for the code: above 300 mbar the residual is scatter about zero of the size the print
   rounding predicts; below it the residual is a trend, negative and growing with depth to
   −3.3e-3 at the bottom row. That is a statement about the source's own integration (its
   harmonic set, its molar mass, its numerics), sits with the systematic terms of §5, and is
   worth one sentence in the manuscript where the closure is presented. Noted in the expected
   values at v0.13.

## On the decisions

Decisions 1 to 9 are accepted as reported. `produce` on a `Profile` and the loaded inputs
(decision 1) is the interface SPEC_04 needs. The bins of decision 5 are superseded by v0.13's,
which differ only in naming the 2 to 10 mbar bin; both readings need no longer be recorded.
The record of decision 7 carries what the manuscript's closure paragraph will cite.

## Order of work

1. Commit the author's documents (SPEC_00 v0.20, SPEC_03 v0.13, this review, `STATE.md`) in
   their own commit.
2. Apply findings 1 and 4: the acceptance and `closure_statistics` on the v0.13 bins; `date`
   in the closure namelist. Rerun the acceptance; refresh check 5, check 12 (`epoch` now
   1981-08-26) and section 6 in the report, noted at its head.
3. The acceptance commit for Step 4. Push.
4. The closure rerun on the clean tree; F5 and F6 from the clean product attached to the
   report; the product's hash recorded; the regression on the clean tree as section 5 says.
   The record commit.
5. Set `STATE.md` to accepted with the commits and SPEC_03 closed. SPEC_04 is drafted next.

## Record verified and figures accepted (16 September 2026)

The clean-tree closure product on disk carries commit `1504056`, no `-dirty` in any attribute of
any group, `epoch` 1981-08-26, and the SHA-256 the report's section 8 records (`e6693173...`); its
`production_record` holds the reviewed statistics to every digit on the v0.13 bins, and the
regression passes on the clean tree. The "2 mbar edge" of decision 10 is a bookkeeping choice
forced by the specification's double description of one bin ("above 2 mbar", "the ten top
levels"); the grid level is the edge, the row printed 2.00 mbar belongs to the deeper bin, and no
bounded value depends on it.

The author viewed `reports/figures/step03_4_lindal_closure_diag_F5_geopotential.png` and
`step03_4_lindal_closure_diag_F6_hydrostatic.png` and accepted both. On F5: the geopotential and
the height lie on each other and the figure carries little information on a single profile; it is
left as is for the figure pass, and its across-latitude panel, dropped at SPEC_03 until the
transfer builds the surface, is where the figure earns its place (SPEC_04).

Step 4 is accepted at `1504056`, the record at `545b385`; SPEC_03 is closed (v0.14).
