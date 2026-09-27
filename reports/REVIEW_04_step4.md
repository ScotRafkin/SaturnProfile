# REVIEW 04, Step 4. The tracing, the transfer, the outer loop, the estimate

Review of `reports/REPORT_04_step4.md` against SPEC_04 v0.13 Step 4. Reviewer: S. Rafkin,
25 September 2026. **Disposition: accepted with changes.** The two failing checks fail on the
test's construction, which the report diagnosed correctly and measured both ways: the M = 2
identity as posed compared two lattices, and on one lattice the transfer is an identity to
3.6e-15. The rulings are in SPEC_04 §15 (v0.15).

## What was verified

The transfer under the closure wind reproduces the reviewing agent's independent values: to
10° N `Δ ln N` within 6.8e-6 at every stated level (bound 1e-4) and the isobar shift within
0.055 percent (bound 1 percent); to 60° N within 9.1e-7 and 0.008 percent. The gauge isobar does
not move, exactly, because `I` is exactly zero on it. Out to 10° N and back under a sheared wind
returns `Φ_k` to 6.6e-8 m²/s² and `ln N_k` to 4.1e-13, falling by 4.9 and 3.6 under halving. The
cylinder-extended wind's traced identity is 1.06e-3 in `ln N` and 163 m²/s² in the isobars,
inside decision Q's bounds and inside the class the two Step 3 kernels gave (4.7e-4 to 1.1e-3;
130 to 224). The loop converges in two passes for the closure wind with the second residual
exactly zero, four for the cylinder wind, six for the sheared wind, every residual below 1e-8.
The two named refusals fire on constructed inputs and a curve reaching the mesh's edge grows
the mesh and repeats. At M = 1 the estimate is the identity to one unit in the last place, and
at M = 2 `φ_r` is the midpoint to 7e-15 degrees, the synthetic anchor's radius sits on the
marched surface to 3e-8 m, and a validation anchor leaves `C` and `φ_r` at their M = 1 values.

## Rulings

1. **The M = 2 identity (finding 1).** The synthetic anchor is written from the Lindal anchor's
   trace on the M = 2 run's own mesh, built first for that run's latitudes; the check is then
   the identity its name claims, bounded at 1e-10 (measured 3.6e-15). The construction as posed,
   the anchor from a separate 60° N run on a lattice that starts at the Lindal latitude, measured
   the tracing's discretization between two lattices (1.2e-4 in `D_12`, 6.0e-5 at the target,
   3.3e-5 in the anchor's `ln N` and 94 m²/s² in its `Φ`); that is a useful number and it is run
   beside the identity and reported, not bounded.
2. **The mesh extension doubling (finding 4).** Ratified: an amount for a rule the specification
   left open, recorded each time, no measured value moved.
3. **The unclosed read handle (finding 5).** A rule for Step 5's driver: reads closed before
   another is opened.
4. **`produce` at the target (finding 6).** At Step 5 `produce` takes the composition on the
   isobar labels as an argument, one path for closure and transfer.
5. **The map's seam (finding 7).** Recorded for the combination specification.
6. **Findings 2, 3 and 8.** Accepted as reported.

Decisions 1 to 12 are accepted; decision 11 is now in the test's text. The production namelist
takes `geopotential_spacing_m2s2 = 5.0e4` from Step 5 (decision R; the reviewing agent's
measurement of the change is 9e-7 in `Δ ln N`), with Steps 0 to 4 recorded as run at 5,000.

## Order of work

1. Restate checks 8 and 10 as v0.15 Step 4 states: the synthetic anchor from the trace on the
   M = 2 run's own mesh, bounded at 1e-10; the lattice variant reported beside it. Rerun those
   runs only (the others are unchanged and their values stand); refresh the report's head,
   section 3 and finding 1.
2. Commit the author's documents (SPEC_04 v0.15, this review, `STATE.md`) in their own commit.
3. The acceptance commit for Step 4. Push. No product changes; `STATE.md` to accepted. Step 5
   proceeds, with the production namelist at 50,000 m²/s², `produce` taking the composition on
   the labels, and every read closed.

## Refreshed report (26 September 2026)

Ruling 1 did what it was for: on the run's own mesh `D_12` fell from 1.2e-4 to 8.2e-8, the
target difference from 6.0e-5 to 4.1e-8, and the synthetic anchor's surface residual to exactly
zero. Checks 8 and 10 still fail against the 1e-10 bound, and the report is right that the bound
was set from a number that never passed through a written anchor; that was the reviewing
agent's reading of the first filing and the report's own wording, and the report's correction is
adopted. The floor the diagnostic isolates is decision G's: an anchor read back from a kind N
file is placed by the field-line integral over its own 66 tabulated levels, which no mesh
refines, and the 0.041 m²/s² gap against the traced characteristic, the same at two resolutions
a factor of ten apart, is that. SPEC_04 v0.16 (§16) bounds the identity through a written anchor
at 1e-7 and the tracing's own identity at 1e-12, separately, removes the script's interpolation
of the anchor's radius from the mesh nodes in the same pass, and has the product's `estimate`
group carry the floor as `identity_floor`, since `D_12` below it says nothing about the
atmosphere. Nothing else in the report changed. **Step 4 is accepted at v0.16 once checks 8 and
10 are rerun under it** (the M = 2 runs only); then the author's documents commit, the
acceptance commit, `STATE.md`, Step 5.

## Verified under v0.16 (27 September 2026)

Eleven of eleven. `D_12` 6.484e-8 and the target 4.1e-8 against 1e-7 through the written anchor,
the tracing's own identity 3.553e-15 against 1e-12 with `Φ_k` returned exactly, the script's
column interpolation removed as ruled (the value is the diagnostic's prediction to three
figures), the lattice variant reported beside it at 1.209e-4. Checks 1, 2, 6, 8, 9 and 10 rerun
and the rest carried from the 25 September run with the rows marked; the regression at every
reference count three times over. **Step 4 is accepted.** Order of work from item 2: the
author's documents commit (SPEC_04 v0.16, this review, `STATE.md`), the acceptance commit,
push, `STATE.md` to accepted, Step 5.
