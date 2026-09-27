# REVIEW 04, Step 5. The production at the target, altitude and datum, the product, transfer mode

Review of `reports/REPORT_04_step5.md` against SPEC_04 v0.17 Step 5 (rulings at v0.18, F9 at v0.19). Reviewer: S. Rafkin,
28 September 2026. **Disposition: accepted with changes.** The one failing check fails on the
specification's stated values and bound, both the reviewing agent's, and the delivered
altitudes are right; finding 3 is a floor the specification carried for the fourth run and
omitted for the first. The rulings are in SPEC_04 §18 (v0.18). The author viewed F5 and F7 and
accepted them, and asked for a figure the transfer already has the data for; it is F8.

## What was verified

The delivered product on disk is the report's (SHA-256 `620f2dee...`), 66 levels, the groups
of deliverable 3, one record group, `pressure_identity_residual` derived, the transfer record
carrying the mesh (421 by 82 at 0.05° by 50,000 m²/s²), two passes with the second residual
exactly zero, the largest `|S/g|` 0.09761 per radian against the reviewing agent's 0.0977, the
largest isobar shift −31,240.9 m²/s² at the top level against −31,258.1, `identity_floor`
6.5e-8, the pressure identity 5.777e-7 largest at the bottom row. The delivered altitudes
411,134.9, 98,186.7 and −15,290.3 m and `r0` 60,128,612.97 m are the reviewing agent's
411,132, 98,186, −15,290 and 60,128,613.0. The temperature departures at 10° N and 60° N are
inside 8e-6 of the independent line integrals at every stated level. The regression's
twenty-five rows are at their reference counts and `step04_4` ran at the 5,000 m²/s² it is
pinned to.

On the cylinder-wind altitudes the reviewing agent rebuilt the construction with the reference
surface at decision P's fixed point, marched under the file's own gauge wind, which the v0.10
measurement had not done (it marched under the closure wind). The fixed-point surface at 10° N
is 60,129,725.8 m, six meters from the report's 60,129,732.3 and 1,113 m above the closure
wind's; under the exact hypothesis, the isobars unmoved, the altitudes are 416,319.8, 99,294.8
and −15,456.9 m, twenty meters above the v0.10 values at the top. Given the run's own arrival
geopotential and datum from `diagnose_cylinder_altitude.txt`, the reviewing agent's column
returns 416,353.9, 99,308.7 and −15,450.3 m at 50,000 m²/s², 416,344.1, 99,303.9, −15,451.2 at
25,000 and 416,341.7, 99,305.0, −15,452.2 at 5,000, against the report's 416,353.8, 99,308.7,
−15,450.3; 416,344.0, 99,303.9, −15,451.3; 416,341.7, 99,305.0, −15,452.2. The altitude path
is right to 0.1 m. What remains between the delivered and the exact values is the gridded
cylinder file's isobar shift, decision Q's floor: +180 m²/s² at the top level and −127 at the
datum at the production spacing, and (180 + 127) over `g` is the 33 m measured.

## Rulings

1. **Finding 3, the pressure identity's floor.** The report's diagnosis is right and the
   convergence condition of §7 applies to the tracing's part, measured with the composition
   held on the tabulated axis; the delivered identity is reported at every level and bounded
   at 1e-6 for this state. The closure path is not changed. The specification carried the
   mechanism for the fourth run and not the first; corrected.
2. **Finding 4, check 4.** The stated values are restated at the fixed point (416,320, 99,295,
   −15,457 m) and the bound is decision Q's floor in altitude, 60 m; the 3.494 m/s the report
   inferred is the shifted isobars read as a wind error, and no wind error exists. Nothing
   about the delivered altitudes changes. The 12 m that moves with the geopotential spacing
   when the wind varies along a column is recorded for the performance step.
3. **Findings 1 and 2, decisions 1 to 4.** Accepted as reported. The renamed vertical
   dimensions of the copied subtrees are the right fix and the right scope.
4. **The regression script edited while running.** Recorded; the separate verification stands.
5. **The author's viewing.** F5 and F7 accepted. F8 added at the author's request (SPEC_04 §7
   deliverable 4 at v0.18): `S/g` and the accumulated temperature change along five isobars
   between the gauge and the target, and the delivered `T` against altitude above the 1 bar
   datum and against radius beside the anchor's own; and F9 (v0.19), the wind the run
   assumed, `u_reference(φ)` with the anchors, gauge and target marked and `u_total(φ, p)` as
   filled contours, from the product's `inputs/wind` group. Check 10 covers both.

## Order of work

1. Restate check 4 against the v0.18 values and bound; measure check 2's convergence on the
   tracing's part; render F8 and F9 and add them to check 10; rerun those checks and, since
   the figure module changes, the regression; refresh the report's head, section 3 and
   findings 3 and 4. The author views F8 and F9.
2. Commit the author's documents (SPEC_04 v0.18, this review, `STATE.md`) on `step04_5-wip`
   in their own commit.
3. The acceptance commit for Step 5 on the branch; `main` fast-forwarded; push. No product
   registered, so no sweep. `STATE.md` to accepted; SPEC_04 closed; the Appendix amendments
   to SPEC_00 and SPEC_03 applied at that commit. Then the performance step (decision R).

## Refreshed report and the author's viewing of F8 and F9 (28 September 2026)

Eleven of eleven under v0.18 and v0.19, the regression rerun for the figure module. Check 2's
tracing part falls under halving (1.293e-7 to 1.200e-7) and the delivered identity holds at
5.777e-7 against 1e-6; check 4 passes at 33.79 m of 60 m with the cylinder `r0` 6.3 m from the
reviewing agent's fixed point; check 10 covers F8 and F9. Decision 5, the kernel carried in the
`isobars` group as `shear_kernel_<slug>_per_rad`, is accepted and is deliverable 3's text at
v0.21. The author viewed F8 and F9: F9 accepted; on F8 the staircase in `S/g` is decision L's
derivative of a piecewise-linear wind on the file's 0.5° nodes and not the mesh, the kinks in
the accumulated temperature are its integral, and the altitude panel's offset was the
specification's wording, each profile on its own reference level 90 km apart, restated at
v0.21 on a common datum (SPEC_04 §18 ruling 5). One rerun remains: check 10 with the panel
restated and the author's viewing of it; the author amended the regression rule (SPEC_04 §0 at
v0.22: a change reruns what it can reach, and a figure change reruns the figure checks of the
steps it draws for and no suite that computes), so this rerun is check 10, `step02_5` and
`step03_4`. Then the commits as the order of work states.

