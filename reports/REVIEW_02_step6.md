# REVIEW 02, Step 6. Anchoring the geoid on the equatorial radius

Review of `reports/REPORT_02_step6.md` against SPEC_02 v0.8 Step 6, SPEC_01 v0.20 and
SPEC_00 v0.15. Reviewer: S. Rafkin, 13 September 2026. **Disposition: accepted. The default
anchor rule for Lindal is `equatorial_radius`** (SPEC_02 v0.9 decision 9, with the comparison
table). SPEC_02 closes at this step's acceptance commit and sweep.

## The decision

The comparison did what it was built to do. Under `equatorial_radius` the declared
uncertainty on `r0` is 11.35 km against 16.38, the shift of `r0` for a 5 percent wind error is
2.60 km against 3.38, and the march predicts Lindal's mean polar radius 3.0 km inside his
±10 km, where the polar anchor predicts his equator at the edge of its ±4 km. The intuition
that the pole is wind-free and therefore the safer anchor does not survive the measurement:
`r0` at 30.8° is the integral of the slope from the anchor to the profile under either rule,
and the mean-polar rule needs the whole pole-to-pole march to fix its mean. Both radii come
from one ellipse fit; the switch chooses the better-constrained end of it. The polar rules
stay in the code and the spread under them stays in `reduction_record`. The reasons and the
table are in SPEC_02 §8 decision 9.

## Rulings on the findings

1. **The asymmetry scales as `r³` with the size of the surface.** Correct, and the derivation
   is the right one: the Eq. B3 slope carries accelerations proportional to `r` over a gravity
   proportional to `1/r²`. The v0.8 expectation of "unchanged to 1 m" was the reviewer's
   error. Check 5 is restated in SPEC_02 v0.9 as the scaled value to 1 m (28,740.31 m
   expected, 28,740.37 measured); rerun it so the report shows seven of seven before the
   commit.
2. **Polar anchoring is the more wind-sensitive.** Agreed, with the geometric reason the
   report gives. The v0.8 expectation of "under one kilometer" was wrong for the reason the
   report identifies. Recorded in decision 9. The 5 percent scaling is a stand-in; the Monte
   Carlo wind draws are to be run under both rules when that wrapper exists.
3. **The three numbers side by side.** This is the table the decision was made from; it is
   copied into decision 9.

## On the decisions

1. Leaving the check failing and reporting the physics was the right call; a loosened check
   would have hidden a real scaling.
2. **Scaling the first secant step.** Accepted. The 1.1e-6 m change in the mean-polar `r0` is
   below every quoted figure and the suites pass; the report says so, which is what matters.
3. **The no-wind reference anchored where the rule anchors.** Accepted, and it was necessary;
   passing an equatorial radius as a polar one would have been a 5,900 km error in the
   reference surface. The no-wind pair under the equatorial anchor (30.799334°, 58,567.166 km)
   is a different reference from the Step 2 no-wind pair, and the record says which rule it
   was built under.
4. **Kinds T and D kept by content.** Accepted. The stale `control_file` hash on those two
   files is the honest consequence and is recorded here; a rewrite for the sake of that
   attribute would have changed hashes that nothing downstream needed changed.
5. **Earlier suites pin the rule they were accepted under.** Accepted, and SPEC_02 v0.9 says
   so beside the Steps 2 to 4 acceptance values.

## Order of work

1. Restate check 5 per finding 1 and rerun; fill in the regression; commit (the acceptance
   commit for Step 6 and for SPEC_02).
2. Sweep: `lindal_refractivity.nc` rebuilt clean under `equatorial_radius` with figures;
   record its SHA-256 in REPORT_02_step6 §7 (the first hash of the product under the new
   default) and the manifest hash as already recorded.
3. Update STATE.md to `accepted` with the commit. SPEC_02 is closed. Nothing further starts
   until SPEC_03 is written.
