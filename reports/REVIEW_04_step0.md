# REVIEW 04, Step 0. Kind W in three parts, kind C on a latitude grid, the transfer run, the namelist, the loader, the anchor object and the propagation hook

Review of `reports/REPORT_04_step0.md` against SPEC_04 v0.7 Step 0. Reviewer: S. Rafkin,
21 September 2026. **Disposition: accepted with changes.** The one failing check fails on a
baseline written before this step was in prospect (finding 3), not on anything built here.
The changes come from the author's direction of 21 September on consistency checks, applied
at SPEC_04 v0.8 (decision N, §0), from the author's decision O of the same day (kind C one
structure, v0.9), and from the rulings in SPEC_04 §11. Also reviewed:
`reports/REPORT_04_preexecution.md`, filed as the record; its two corrections are adopted
(section 3 below).

## What was verified independently

From the rebuilt files on disk, with the reviewing agent's own code: the SHA-256 of
`lindal_wind.nc`, `lindal_refractivity.nc`, the closure product and the transfer run's wind and
composition are the report's section 6 values. The rebuilt wind carries `u_total_ms`,
`u_reference_ms` and `u_shear_ms` on (361, 61), no `u_cylindrical_ms`, no `decomposition`
attribute, `reference_level_pressure_Pa` 1.0e5, and the sum identity holds to 0.0 m/s. The
rebuilt kind N gives `σ_N / N` = 2.331845e-2 at every one of the 66 levels while `σ_N` runs
from 1.18e-9 to 7.25e-6, which is finding 1. The rebuilt closure product's produced pressure at
the gauge level is 9998.46545058 Pa, the value of the accepted product, and its
`production_record` residuals are the accepted ones to every printed digit. The reviewing
agent's independent transfer, rerun on the rebuilt anchor and wind, returns every v0.6
expected value unchanged (to 10° N: `Δ ln N` 1.1090e-2, 1.0950e-2, 1.0899e-2 at the top,
gauge and bottom; shift −31,258.1 m²/s²; `r0(10°)` 60,128,613.0 m; altitude 411,132 m), which
is the same statement as check 6: the cascade changed no value. The transfer composition is
on (66, 181). The transfer namelist on disk is the specification's with the 0.05° spacing.

## Rulings on the findings

1. **`sigma_ln_N_measurement`.** Correct, and the fault was the specification's. Deliverable 6
   now reads `σ_N / N` times the scale. The implementation is the ruled one.
2. **The `[grid]` refusal.** The check is removed, not relocated. The author's direction is
   that the pipeline assumes its inputs came from the prior steps and the code does not police
   what a user might declare; a range that must cover the anchors is a number the code can
   derive, so at v0.9 the namelist declares `geopotential_spacing_m2s2` and
   `latitude_spacing_deg` only, Step 2 builds the mesh from the anchors' levels and extends it
   when a traced curve needs more, and no refusal remains. Decision 2 is superseded.
3. **Check 15 of `accept_step03_3`.** Refresh the `before/` baseline of `step03_3` at the
   acceptance commit; the check then means "no value changed since SPEC_04 Step 0" and check
   13 passes on the numbers already measured. Not loosened: check 6 of this step is the same
   comparison against the swept copies and passes.
4. **The retrieval refusal.** Removed from the loader; dead code by the schema's prior claim.
   The retrieval leg adds its own when the instance exists.
5. **`accept_step7`.** Accepted: the Appendix's amendment applied to a closed step's acceptance
   record, and recorded.

## The author's direction on checks (decision N)

The author read this report and ruled that there are far too many consistency checks, that the
code is not meant to be idiot-proof, and that the pipeline is to assume it ingests the
appropriate data from the prior steps. SPEC_04 v0.8 §0 states what is checked from here on:
the schema, namelist typo guards (an unknown, missing or malformed key; a scheme not
implemented), extrapolation (an interpolant refuses a point outside its data), and named
numerical failures. Everything else is recorded, warned where the record alone might be
missed, never refused. Applied to this step:

- The loader's coverage checks (wind latitude and pressure against the anchors and the target,
  composition latitude) are removed; the interpolants of Step 1 refuse extrapolation and that
  is where coverage lives.
- A wind or composition at another season is recorded and warned, not refused. This amends
  decision I's "else refused" at the author's direction.
- Kind C becomes one structure, a field on (level, latitude), for every use (decision O):
  the tool's `latitude_grid_deg` is required, the one-latitude form and its absence marker
  are retired, the reduction reads the column at the anchor's latitude by decision L (the
  identical column for a uniform field; no value changes), and the point-composition refusal
  and decision 5's planetographic pairing go with the form. Done inside this step's rebuild,
  so the chain is rebuilt once.
- The `[grid]` range check goes (ruling 2); `[target]` with two values and a closure namelist
  with `[target]` are the generic unknown-or-malformed-key rule, not cases of their own.
- What stays at Step 0: a missing `[target]`, an unknown `[numerics]` table, a scheme not
  implemented, a `weight` other than 0 or 1, and the thirteen closure refusals of SPEC_03,
  which are that specification's. Decision 3 (a nonzero `model_error_correlation_length_deg`
  refused, naming A35) stays, since a recorded-and-ignored value would be the silent kind.
- The kind N reader's refusal (SPEC_02) is of a file whose two hash records disagree, a
  file edited after it was written; a changed input beside it only warns, as for every derived
  kind. It stays. The cascade was the closure's own definition, that its inputs are the
  reduction's, not the reader's; a first draft of this review said otherwise and is corrected.

## On the decisions

Decisions 1, 3, 4 and 6 to 13 are accepted as reported. Decision 2 is superseded by ruling 2
and decision 5 by the direction above. Decision 9 (the swept copies rebuilt in a worktree at
`5f23a8b`) is the right recovery and is recorded so that the sweep's own copies replace them.

## Corrections adopted from REPORT_04_preexecution

The interpolant that gave 1.93991 m/s is PCHIP, and the reviewing agent's v0.5 expected values
were measured with `PchipInterpolator`, not with a cubic spline; SPEC_04 §10 and decision L
are corrected at v0.8. The composition bound is the profile-wide 1.5e-6, as §10 already
carried. The carried-forward item on the derivative at a node the mesh and the file share is
settled in Step 1 deliverable 2 as written at v0.6: the node takes the slope of the interval
to its north, or to higher pressure; the reviewing agent's expected values were measured under
that rule (the largest `|S/g|` lies on the interval south of the 31.0° node, where the slope
is −12.07 m/s per degree).

## Order of work

1. Apply decisions N and O as v0.9 Step 0 states: kind C on (level, latitude) for every use
   (`casspian-composition-lindal` with `latitude_grid_deg` required, the reduction's control
   file and the closure's build file declaring `[-90, 90, 1.0]`, the schema and the readers
   changed, the reduction reading the column at `φ_c` by decision L); the loader's coverage,
   season, point-composition and retrieval refusals removed (season recorded and warned);
   `geopotential_range_m2s2` removed from the namelist and `lib.control` with no replacement,
   the range check gone (the mesh of Step 2 builds and extends itself); `sigma_ln_N_measurement`
   as ruled (already so). Rerun the cascade through the same candidate procedure, no value
   changing. Refresh the `step03_3/before/` baseline. Rerun the acceptance; check 11 becomes
   the four refusals and the one recorded case of v0.9 (the season); check 13 passes; refresh
   the report's sections 1 to 4 and 6 and note it at the head.
2. Commit the author's documents (SPEC_04 v0.8, this review, `STATE.md`, the modified
   `docs/CASSPIAN_Seasonal_Design_Note.md`, which is the author's v0.5 and is to be committed
   with them) in their own commit.
3. The acceptance commit for Step 0. Push.
4. The sweep on the clean tree over `occul_data/` and `forward/`: the chain, the closure inputs
   and product, the transfer inputs rebuilt clean; the deferred suites (`step02_1`, `step02_4`
   to `step02_6`, `step03_1`, `step03_2`, `step03_4`) run on the swept products and recorded;
   the hashes recorded in the report's section 6 in the `step02_1` row format, replacing the
   provisional table; the swept copies replacing decision 9's worktree copies. The record
   commit.
5. Set `STATE.md` Step 0 to accepted with the commits. Step 1 proceeds.

## Refreshed report verified (21 September 2026)

The report refreshed under v0.9 was read and its products verified from disk: the SHA-256 of
the rebuilt kind N, the reduction's composition and wind, the closure product and the transfer
run's composition and wind are the section 6 values; the reduction's composition is a field on
(66, 181) with no absence marker and no variation across latitude; kind N's `σ_N / N` is
2.331845e-2 at every level; the closure product's produced pressure and temperature at the
gauge level are 9998.465450581718 Pa and 83.38720185785154 K, and its top pressure
19.952623149688797 Pa, the accepted values; the reviewing agent's independent transfer, rerun on
the refreshed anchor and wind, returns every v0.6 expected value unchanged. Fifteen of fifteen,
check 14 at 16 of 16 on the refreshed baseline.

Finding 6 (`accept_step8` and `accept_step02_3` moved to the field form of kind C, each first
asserting that nothing varies with latitude) is accepted on the same footing as finding 5.
Decisions 14 to 18 are accepted as reported; `lib/composition.py` as the one home of the
reading rule is the right size, and a flag taken from the nearer node rather than averaged is
correct. **Step 0 is accepted.** Proceed from item 2 of the order of work: the author's
documents commit (the design note included), the acceptance commit, push, the sweep with the
deferred suites and the swept hashes, `STATE.md` to accepted, Step 1.
