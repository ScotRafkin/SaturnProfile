# REPORT 01, Step 5. `lib.geoid`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.6, Step 5, against
`SPEC_00_Architecture_and_Data_Files.md` v0.7. Status: reported, awaiting review.

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product.

**Step 4 closed first.** Committed as `cc5bbee`. Two docstrings were corrected rather than any
behavior changed: `omega_abs` no longer claims a polar guard it does not provide (the review
was right that `np.cos(np.pi/2)` is 6.1e-17, so the guard catches an exact zero only, and the
correct statement is that a wind field nonzero at the pole is the defect and the wind tool is
where it is caught), and `g_eff_vector` no longer describes the minus sign in
`psi = arctan2(-G_phi, g)` as something Step 6 would otherwise reintroduce, since it is part of
the definition of `psi`. The southern hemisphere was checked: at `phi_c = -45` degrees,
`G_phi = +1.056234` and `psi = -6.3522` degrees, mirroring the north with no case split, which
is what `arctan2` buys.

---

## 1. What was built

**`src/casspian/lib/geoid.py`.** `U_rigid`, `ellipsoid_seed`, `reference_geoid`, `wind_geoid`,
and `radius_at`. Pure functions, no file access, no module level Saturn quantity.

`reference_geoid` solves `U(r, phi) = U_ref` at each latitude by Newton iteration from a first
order hydrostatic ellipsoid seed. SPEC_01 writes the update as `r <- r + (U - U_ref) / g_r`
with Lindal's outward `g_r`; since handoff section 9A.5 fixes `g_r = -g` and `dU/dr = g` in
this package's sign convention, the implemented step is `r <- r - (U - U_ref) / g`. They are
the same update written twice, and the docstring says so. Each latitude stops on its own
tolerance, so a slow latitude does not make the rest iterate.

`wind_geoid` integrates Eq. B3, `g dr0/dphi = r0 G_phi`, from each pole inward with RK4. Each
hemisphere is integrated from its own pole, so a wind that is not symmetric about the equator
needs no reflection assumption. If the supplied grid does not contain the pole it is inserted
into the march and removed from the result.

## 2. Decisions

1. **The seed is a first order hydrostatic ellipsoid**, `f = (3/2) J2 + q/2` with
   `q = Omega^2 r_polar^3 / GM`. For Saturn that gives `f = 0.082` against the true 0.098,
   which is close enough to start from and never appears in a result. Newton converges in five
   steps at every latitude from it.

2. **A loop over grid nodes in the ODE march is not a loop over array elements.** SPEC_00
   section 3.1 forbids Python level loops over array elements in `lib`. An initial value
   problem cannot be vectorized along its own integration axis: each step needs the previous
   result. The same note was made for the Legendre recurrence at Step 4, on different grounds.

3. **`reference_geoid` and `wind_geoid` are different objects and the module says so.** The
   first is an equipotential by construction. The second is not the level set of anything,
   because with `u` varying in latitude the effective field is not conservative, which is what
   manuscript section A6 is about.

4. **The Newton residual is reported, not just the radius.** `max |U - U_ref|` is 3.6e-07
   m2/s2 across the grid, which at `g = 10` is 4e-08 m of radius. The surface is an
   equipotential to well below any physical tolerance.

## 3. Defects and questions raised

1. **`radius_at` uses the scheme SPEC_01 names, and that scheme is the least accurate of the
   three obvious ones.** Step 5 specifies interpolation "linear in sin phi_c, stated". Measured
   against the geoid solved directly at 30.8185 degrees, on grids from 2 degrees to 0.1 degrees:

   | spacing | linear in `sin(phi)` (specified) | linear in `phi` | `1/r^2` linear in `sin^2(phi)` |
   |---|---|---|---|
   | 2.00 deg | 1023.20 m | 526.71 m | 118.98 m |
   | 1.00 deg | 159.35 m | 84.63 m | 18.19 m |
   | 0.50 deg | 61.59 m | 32.21 m | 7.10 m |
   | 0.20 deg | 3.56 m | 1.85 m | 0.41 m |
   | 0.10 deg | 1.60 m | 0.83 m | 0.19 m |

   All three are second order. The specified scheme is about twice the error of interpolating
   in the angle and about eight times the error of interpolating `1/r^2` in `sin^2(phi)`. The
   last is best for a reason: for an ellipse `1/r^2 = sin^2(phi)/b^2 + cos^2(phi)/a^2` holds
   **exactly**, and the geoid is an ellipse to within the harmonic and rotational corrections,
   so that scheme interpolates only the small residual. **Proposing `1/r^2` in `sin^2(phi)`.**
   The specified scheme is what is implemented; the alternative is raised, not substituted.

   This is not academic. `r0(phi_c)` is the constant of integration for the whole profile: an
   error in it is an error in every absolute radius the reduction reports. At the anchor
   latitude a 2 degree grid costs a kilometer with the specified scheme.

   I should also record that my first docstring justified the specified scheme with a rationale
   I had invented, that a geoid is close to linear in `sin`. The measurement above disproves it
   and the docstring has been corrected. The specification's own wording gives no reason, only
   the instruction.

2. **The stand-in wind does not determine the wind bulge, so "roughly 120 km" is not a test of
   the code.** Step 5 offers a stand-in `u(phi)` with "450 m/s at the equator falling to zero
   by 35 degrees" and expects an equatorial radius about 120 km above the no wind value. Every
   profile below satisfies that description exactly:

   | stand-in profile | equatorial radius | bulge above no wind |
   |---|---|---|
   | `cos^2` taper to 35 deg | 60309.82 km | 65.83 km |
   | linear taper to 35 deg | 60316.66 km | 72.67 km |
   | flat to 20 deg then linear to 35 deg | 60382.74 km | 138.76 km |
   | flat to 35 deg, then zero | 60459.45 km | 215.46 km |
   | `cos^2` taper to 45 deg | 60350.07 km | 106.08 km |
   | `cos^2` taper to 60 deg | 60424.28 km | 180.30 km |
   | linear taper to 60 deg | 60439.77 km | 195.78 km |

   The bulge spans 66 to 215 km, and 120 km sits in the middle of that range without being
   picked out by it. The reason is visible in the equation: the wind's contribution to the
   slope is `dr0/dphi ~ r (2 Omega u sin phi) / g`, so the bulge is a weighted integral of
   `u(phi) sin(phi)`, and `sin(phi)` weights the outer part of the jet, which is exactly where
   the profiles above differ most. The check as written measures the profile I invented, not
   `lib.geoid`.

   The check that does test the code is grid convergence, which is reported: the equatorial
   radius moves by 1.2e-05 m across an eightfold change of step. **Proposing that the 120 km
   comparison be dropped from Step 5 and made part of Step 7**, where the real Smith wind
   exists and the comparison against Lindal's fitted 60,367 km becomes meaningful. Step 5 as
   written already says the real number "is to be reported, not prescribed"; this asks that the
   stand-in number not be treated as an acceptance at all. The first row of the table is the
   figure to quote if one is wanted: **65.83 km** with a `cos^2` taper to 35 degrees.

3. **The closure diagnostic has two readings that differ by a factor of 6.6.** Step 5 says
   `wind_geoid` "reports the closure `U(r0(phi), phi) - U_ref` as a diagnostic of how far the
   wind takes the surface from an equipotential", after observing that the field is not
   conservative when `u` varies. Taken literally, `U` is `U_rigid`, the **no wind** potential,
   and the closure is 5.99e+05 m2/s2, which is `g` times the bulge: it restates the bulge in
   potential units and says nothing about conservativeness. Taken as the non conservativeness
   measure the preceding sentence points at, the quantity is the wind pseudo-potential
   `V - (1/2) Omega_abs^2 r^2 cos^2(phi)` evaluated along the surface, whose variation is
   3.95e+06 m2/s2.

   The literal reading is implemented and returned. Both numbers are in the acceptance output.
   **Which one Step 5 wants is the question**; if it is the second, `wind_geoid` returns one
   more array and the signature changes by one element.

4. **Which GM builds the reduction geoid is not stated.** Step 5 says "Null's set, System III"
   and gives no GM. Both are inside the tolerance: Null's own GM gives 60243.99, 58435.30 and
   58453.09 km, the modern GM gives 60243.69, 58435.11 and 58452.90 km. The handoff's
   58452.9 is reproduced **exactly** by the modern GM and to 0.2 km by Null's, which suggests
   the handoff numbers were computed with the modern value. The reduction should use Null's,
   for the reason `null1981.toml` states, and Step 9 will read it from `lindal_gravity.nc`
   rather than choosing. Recorded so the 0.2 km is not later mistaken for an error.

## 4. Acceptance results

Run by `reports/step5/accept_step5.py`; full output in `reports/step5/output.txt`.

| Check | Result |
|---|---|
| `U_rigid` obeys the stated sign convention | **Pass.** A central difference gives `dU/dr = 9.465487421` against `g_eff_radial = 9.465487404` m/s2, agreeing to 1.7e-08. Lindal Eq. 11 carries the opposite overall sign; this package uses `g_eff = -grad U` throughout and the docstring says which. |
| No wind geoid, equatorial radius 60,244 +- 2 km | **Pass. 60243.99 km** with Null's GM, 60243.69 with the modern one. |
| `r(31.0)` = 58,435 +- 2 km | **Pass. 58435.30 km** (modern: 58435.11). |
| `r(30.8185)` = 58,452.9 +- 2 km | **Pass. 58453.09 km** (modern: 58452.90, exact). |
| The surface is an equipotential | **Pass.** Five Newton steps at every latitude; `max abs(U - U_ref)` = 3.58e-07 m2/s2, which is 4e-08 m of radius. The radius at the pole returns the anchor exactly. |
| Beyond the specification: grid independence of the no wind geoid | **Pass.** The equatorial radius is identical to every printed digit at n = 451, 901, 1801, 3601, since Newton solves each latitude independently. |
| The 45 degree check repeated on the constructed surface, as the Step 4 review asked | **Pass.** On the constructed geoid, `r_ref(45) = 57020.926 km`, `G_phi = -1.072649` m/s2, `psi = 5.7166` degrees, `abs(g_eff) = 10.7688` m/s2. The Step 4 figure on the equatorial 1 bar sphere was `-1.056221`, `6.3525` degrees and `9.5461`. The sphere overstates the radius at 45 degrees by 3247 km, which is why `abs(g)` was low by 1.22 and `psi` high by 0.64 degrees. Same sign, same order, as Step 4 claimed. |
| Wind geoid bulges outward | **Pass**, 65.83 km with a `cos^2` taper to 35 degrees, and 66 to 215 km across profiles equally consistent with the specification. See section 3, item 2. |
| Beyond the specification: grid convergence of the wind geoid | **Pass.** The equatorial radius spreads by 1.2e-05 m across n = 451 to 3601, an eightfold range of step, which is RK4 behaving. |
| The closure diagnostic | **Reported both ways**: 5.99e+05 m2/s2 against the no wind potential, 3.95e+06 against the wind pseudo-potential. See section 3, item 3. |
| `radius_at` interpolates as specified and refuses to extrapolate | **Pass**, with the accuracy of the three candidate schemes tabulated. See section 3, item 1. |

Seven of seven pass. No em dash or en dash appears in any file written in this step.

## 5. Next step

Step 6, `lib.latitude`, is not started and is not to be begun until this report has been
reviewed. It is the first step whose acceptance depends on this one: the fixed point evaluates
`psi` on a surface constructed here, and its target `phi_c` = 30.8185 degrees is the latitude at
which `r(30.8185) = 58453.09 km` was just reported.

---

## 6. Changes applied after review

Added 10 September 2026, after `REVIEW_01_step5.md`, so that this report and the committed code
agree. Four changes were made and the acceptance rerun; seven of seven still pass.

**Rulings on the items this report raised.**

1. **Interpolant: `1/r^2` linear in `sin^2(phi_c)` adopted**, with two conditions the report had
   not thought through. The bracketing is done in `phi_c` and not in `sin^2(phi_c)`, because the
   latter is not monotonic across the equator and a global search on it would fold the
   hemispheres onto each other, returning a northern radius for a southern latitude on any
   surface that is not symmetric. An interval straddling the equator is refused, so a grid is to
   contain the equator as a node. Both are implemented and checked. More important than the
   scheme is the rule the review added: **the frozen anchor radius is never interpolated.**
   `reference_geoid` solves each latitude independently and `wind_geoid` marches any node set,
   so the anchor latitude is solved or marched, and the kilometre at 2 degree spacing never
   reaches a product. `radius_at` exists for the vectorized uses of Step 7 and later.

2. **The stand-in bulge is not an acceptance.** The 120 km figure is out of Step 5, grid
   convergence is the acceptance for the wind geoid, and the comparison with Lindal's fitted
   60,367 km moves to Step 7 with the real wind. The 66 to 215 km table stays as the
   demonstration.

3. **Closure: the literal reading, and only it.** The pseudo-potential is removed from the code
   path and from the acceptance output. The reviewer's reason is better than the report's
   framing: with `u` varying in latitude the gradient of `V - (1/2) Omega_abs^2 r^2 cos^2` is
   not the effective gravity, so its variation along the surface restates the wind kinetic term
   and measures nothing about conservativeness. The non conservativeness of the field is the
   shear kernel of Eq. A15 and belongs to the forward model. The closure against the no wind
   potential is a real quantity: the dynamical height of the wind surface above the no wind
   geoid, 5.989e+05 m2/s2, which at the equatorial `g` of 8.934 m/s2 is 67.03 km.

4. **GM: Null's, read from `lindal_gravity.nc` at Step 9.** Confirmed, and the 0.2 km against
   the handoff numbers is now recorded in the specification so it is not later taken for an
   error.

**Two corrections this report failed to raise, found by the reviewer.**

5. **The claim that Lindal Eq. 11 carries the opposite overall sign is false, and I repeated it
   without checking.** SPEC_01 v0.6 said it, and the `U_rigid` docstring said it after the spec.
   Checked against the paper by the reviewer: Lindal's Eqs. 9 and 11 give
   `U = -GM/r + (GM/r) sum_i J_2i (R/r)^2i P_2i(sin phi) - (1/2) omega^2 r^2 cos^2(phi)` with
   `g = -grad U`, which is term for term and sign for sign the `V - (1/2) Omega^2 r^2 cos^2`
   this module computes. The docstring now says they agree. This one is worth stating plainly:
   the specification asserted a fact about a source, the source was in `docs/`, and the report
   passed the assertion through into a docstring instead of opening the paper. The `dU/dr = g`
   check that was run is a check of internal consistency and would not have caught it.

6. **The unused `tol_m` and `max_iter` on `wind_geoid` are removed.** They were added for
   signature symmetry with `reference_geoid` and documented as unused, which the report treated
   as sufficient. It is not: an argument a caller can pass and that changes nothing is a trap,
   and documenting a trap does not defuse it. The signature is now
   `wind_geoid(phi_c_grid, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm)`.

**Acceptance after the changes.** Seven of seven. The interpolation table now shows the ruled
scheme in its third column, smallest at every spacing: 118.98 m at 2 degrees against 1023.20 m
for the v0.6 scheme, and 0.19 m at 0.1 degrees against 1.60 m. The hemisphere check passes on a
grid spanning both poles, and a grid without the equator as a node is refused rather than
folded.
