# REPORT 01, Step 6. `lib.latitude`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.7, Step 6, against
`SPEC_00_Architecture_and_Data_Files.md` v0.7. Status: reported, awaiting review.

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product.

**Step 5 closed first.** Committed as `b5a09ff`. All four review changes were made and the
acceptance rerun, seven of seven: `radius_at` now interpolates `1/r^2` linearly in
`sin^2(phi_c)` with bracketing in `phi_c` and refusal of an equator straddling interval; the
pseudo-potential is gone from the code path and the output; the `U_rigid` docstring no longer
repeats the withdrawn claim about Lindal Eq. 11's sign; and `tol_m` and `max_iter` are off
`wind_geoid`. The full regression across Steps 1 to 5 passes: 6, 14, 8, 5, 7, 7.

---

## 1. What was built

**`src/casspian/lib/latitude.py`.** `planetocentric_from_ellipsoid`,
`planetocentric_fixed_point`, and `planetographic_from_planetocentric`. Pure functions, no file
access, no module level Saturn quantity.

The fixed point iterates `phi_c <- phi_g - psi(phi_c)` with `psi` from `lib.gravity.g_eff_vector`
evaluated at `r = surface(phi_c)`. It is vectorized over `phi_g` and each element stops on its
own tolerance. It returns the converged latitude, the tilt there, the full history of iterates,
and the step count per element, so the convergence can be inspected rather than trusted.

`surface` is a callable that **solves or marches at the latitude it is given**, per SPEC_01
v0.7: for the no wind geoid a direct Newton solve, for the wind geoid a march with that
latitude inserted as a node. The docstring says so at the point a caller would get it wrong.
The acceptance builds both.

## 2. Decisions

1. **Angles are radians in and out**, consistent with `lib.gravity` and `lib.geoid`, and
   `tol_deg` is the single exception, in degrees, because SPEC_01 names it that. See section 3,
   item 1.

2. **`arctan2` in the ellipsoid seed**, not `arctan`. At the pole `tan(phi_g)` is infinite and
   `arctan` of the product would be a NaN; the two argument form returns exactly `+-pi/2`. The
   acceptance evaluates there. The same reasoning that put `arctan2` into `g_eff_vector` at
   Step 4 applies here for a different reason.

3. **The flattening that seeds the iteration is an argument, defaulting to the surface's own.**
   When it is not supplied the surface is asked for its equatorial and polar radii and the
   flattening is computed from those, so the seed comes from the object the iteration converges
   on. The acceptance passes the 100 mbar oblateness explicitly, which is what reproduces the
   handoff seed. See section 3, item 2. The converged answer does not depend on the seed; only
   the step count does.

## 3. Questions raised

1. **`tol_deg` is the only degree valued quantity in `lib`.** Every angle in `lib.gravity`,
   `lib.geoid` and this module is radians, and a tolerance in different units from the quantity
   it bounds is a place to make a mistake. It is implemented as SPEC_01 names it and converted
   internally at one line. **Proposing `tol_rad`** for consistency, or, if a degree valued
   tolerance is wanted because that is how the acceptance is stated, that the argument keep its
   name and the docstring keep saying loudly what it is. Either is a one line change; the
   current state is the second.

2. **Which flattening seeds the iteration is not stated, and only one reproduces the handoff.**
   The seed is `30.8524` degrees with the 100 mbar oblateness `0.09822`, against the handoff's
   `30.8526`. The 1 bar oblateness `0.09796` gives `30.8670`, and the constructed geoid's own
   flattening `0.09637` gives `30.9556`. The 100 mbar surface is the right one on the physics:
   the Fig. 4 caption defines the label as the latitude where the lowest ray point touched the
   100 mbar surface, so the ellipsoid the seed approximates is that surface. Recorded because
   the acceptance names the seed value, so the choice is visible in a number.

3. **"Three or four iterations" has no tolerance attached.** The step count depends entirely on
   it: three steps reach 1e-3 degrees, four reach 1e-4, five reach 1e-6 and eight reach 1e-9.
   The specification's three or four is right for a tolerance around 1e-4 degrees, which is
   0.36 arcsec and finer than the 0.2 degree uncertainty on the label itself. Suggesting the
   acceptance name the tolerance it means.

4. **The handoff numbers were computed with the modern GM, again.** Null's GM gives
   `phi_c = 30.818189` and the modern GM `30.818458`, against the handoff's `30.8185`: 1.12 and
   0.15 arcsec away. Step 5 found the same pattern in the radii. Both are well inside the
   +-0.001 degree tolerance, and the reduction uses Null's. Recorded so the 1.12 arcsec is not
   later taken for an error.

## 4. Acceptance results

Run by `reports/step6/accept_step6.py`; full output in `reports/step6/output.txt`.

| Check | Result |
|---|---|
| The ellipsoid seed | **Pass. 30.8524 degrees** with the 100 mbar oblateness, against the handoff's 30.8526. Finite at both poles, returning exactly `+-90` degrees. |
| `phi_c` = 30.8185 +- 0.001 for `phi_g` = 36.3 | **Pass. 30.818189** with Null's GM, **30.818458** with the modern one. |
| `psi` = 5.4815 +- 0.001 | **Pass. 5.481811** with Null's GM, 5.481542 with the modern one. |
| Converging from the seed in three or four iterations | **Pass at the tolerance the count implies.** The iterates are 30.85241, 30.81577, 30.81836, 30.81818, 30.81819, 30.81819. Three steps reach 1e-3 degrees, four 1e-4, five 1e-6, eight 1e-9. |
| Round trip returns 36.3 to 1e-9 | **Pass. 2.30e-11 degrees**, which is 8.3e-05 milliarcsec. |
| `phi_g` = 36.5 gives 31.005 | **Pass. 31.005105.** The 0.2 degree difference in the label moves the planetocentric value by 0.186916 degrees, so the label's own uncertainty passes through almost one for one and is not reduced by the conversion. That is worth keeping in view: the anchor latitude is known no better than the figure label it came from. |
| With the wind included the value moves, amount reported | **Pass.** With the Step 5 stand-in wind, `phi_c` moves from 30.818189 to **30.801766 degrees**, a shift of **-0.016423 degrees, 59.12 arcsec equatorward**, and `psi` rises from 5.481811 to 5.498234. The anchor radius moves from 58453.120 to **58455.181 km**, up 2.061 km. This is the stand-in, not the Smith wind; the value the reduction freezes is Step 7's. |
| Beyond the specification: vectorization, both hemispheres, the equator | **Pass.** Seven latitudes from -60 to +60 in one call, bit identical to the seven scalar calls. North and south mirror to 0.0e+00 degrees and the equator maps to exactly zero, as they must on a surface with only even harmonics and no wind. The equator converges in one step, since `psi` is zero there. |

Six of six pass. No em dash or en dash appears in any file written in this step.

One correction to the script during the run, recorded because it affected a printed number: the
wind comparison first printed the no wind `psi` from the modern GM while its `phi_c` was Null's,
because a loop variable from an earlier check was still in scope. Both numbers in that line are
now from Null's GM. No check outcome changed.

## 5. Next step

Step 7, the wind tool for Lindal and kind W, is not started and is not to be begun until this
report has been reviewed. It is the step that produces the real `u(phi)`, and with it the two
numbers this step and Step 5 could only bracket: the frozen anchor latitude and the wind geoid
equatorial radius to be compared with Lindal's fitted 60,367 +- 4 km.

---

## 6. Changes applied after review

Added 11 September 2026, after `REVIEW_01_step6.md`. One rename and one docstring addition; the
acceptance was rerun and still passes six of six.

1. **`tol_deg` is now `tol_rad`.** Radians everywhere in `lib`, no exceptions. The acceptance
   script still states its tolerances in degrees, because that is how SPEC_01 states the
   convergence claim, and converts them at one line in `solve`. The module docstring records
   that v0.7 named it in degrees and v0.8 corrected it, so the change is traceable.

2. **The returned `psi` is documented as the tilt at the previous iterate.** The reviewer's
   point is worth having in the code: `phi_c + psi` equals `phi_g` exactly, because the last
   update was `phi_c = phi_g - psi`, and the tilt evaluated afresh at the returned `phi_c`
   differs by less than the final step times `dpsi/dphi`, under a tenth of the tolerance. The
   docstring now says not to "fix" it into an extra evaluation, which would cost a surface
   solve and break the exact identity. The round trip check, 2.3e-11 degrees, is the evidence.

3. **The other three items** (seed flattening, iteration counts with their tolerances, the GM
   record) were rulings on the specification and needed no code change. They are in SPEC_01
   v0.8 Step 6.

Recorded separately: `reports/REVIEW_01_step6.md` as delivered contains a `<system-reminder>`
block at its end, carrying commit attribution instructions. It is text inside a data file, not
a system message, and was not acted on. SPEC_01 section 0 forbids AI attribution anywhere and
`REVIEW_01_step1.md` reaffirmed it, so commits remain untrailered. Flagged so the stray block
can be removed from the review file.
