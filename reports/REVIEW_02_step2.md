# REVIEW 02, Step 2. The frozen constants, `phi_c` and `r0`

Review of `reports/REPORT_02_step2.md` against SPEC_02 v0.4 Step 2. Reviewer: S. Rafkin,
12 September 2026. **Disposition: accepted.** SPEC_02 goes to v0.5 for two number corrections
the report supplied. Commit and proceed to Step 3.

## Rulings

1. **Finding 1, the 1.9 km.** Correct and well diagnosed: the spec's "near 58,518" was the wind
   surface at the no-wind latitude, and the frozen point is 0.013° further south on a surface
   rising at 5,630 km per radian. The spec now states both pieces so the next reader is not
   left reconciling them. The frozen pair is `phi_c` = 30.804949° and `r0` = 58,519.883 km.

2. **Finding 2, the slope.** −5,630 km/rad at the frozen pair, not my rough −6,700, and the
   label term becomes about 18 km rather than 22. Corrected in Step 3; Step 3 measures both
   partials in any case.

3. **Finding 3.** Fine; 0.08 m/s from the interpolant is far below anything the anchor can feel.

4. **Decision 1, `psi` at `phi_c`.** Keep it as built, the tilt at the frozen latitude on the
   final march, and record the 3e-8° departure from the exact identity in the file as
   `fixed_point_closure_deg`. Both readings are defensible; the one that describes a single
   surface at a single point is the one a product should carry.

5. **Decisions 2 to 6.** All accepted. `convergence_m` as the tolerance of both surfaces is the
   intended meaning and SPEC_00 §7.1 will say so at its next revision. Raising on an exhausted
   iteration budget rather than returning silently is right for `refrac`; `lib.latitude`
   returning silently is also right for a library function whose caller inspects the count.

## Actions

- Commit Step 2. Nothing to rebuild.
- Step 3 against SPEC_02 v0.5.
