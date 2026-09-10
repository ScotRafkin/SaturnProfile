# REVIEW 01, Step 5. `lib.geoid`

Review of `reports/REPORT_01_step5.md` against SPEC_01 v0.6 Step 5. Reviewer: S. Rafkin,
10 September 2026. **Disposition: accepted with changes.** SPEC_01 is revised to v0.7. The
changes below are to be made and the acceptance rerun before the commit; no product to
rebuild.

## Rulings on section 3

1. **Interpolant: adopt `1/r²` linear in `sin² φ_c`.** The measurement is convincing and the
   reason is right (exact for an ellipse, so only the harmonic and wind residual is
   interpolated). Two conditions on the implementation. First, `sin² φ_c` is not monotonic
   across the equator, so a global `np.interp` on it would fold the hemispheres onto each
   other and would be wrong for any wind geoid that is not symmetric; bracket the target
   between its neighboring nodes in `φ_c` and interpolate between those two, and refuse an
   interval that straddles the equator (a grid is to contain the equator as a node). Second,
   and more important than the scheme: **the frozen anchor radius is never interpolated.**
   `reference_geoid` solves every latitude independently, so the anchor latitude goes into the
   grid and is solved; `wind_geoid` marches any node set, so the anchor latitude is inserted as
   a node. `radius_at` exists for the vectorized uses of Step 7 and later. With that rule the
   kilometer at 2° spacing never reaches a product. Both are now in the spec.

2. **The stand-in bulge is not a test. Agreed.** The 66 to 215 km table is the right way to
   show it. The 120 km figure is removed from Step 5; grid convergence is the Step 5 acceptance
   for the wind geoid; the comparison with 60,367 ± 4 km moves to Step 7 with the real wind,
   reported not prescribed.

3. **Closure diagnostic: the literal reading, and only it.** `U(r₀(φ), φ) − U_ref` against
   the no-wind potential is the dynamical height of the wind surface above the no-wind geoid
   in potential units, which is a meaningful and useful number (it is Lindal's `h(φ)` times
   `g`). The pseudo-potential `V − ½ Ω_abs² r² cos² φ` is not to be evaluated: with `u`
   varying in latitude its gradient is not the effective gravity, so its variation along the
   surface (3.95e6 m²/s²) restates the wind kinetic term `Ω r cos φ u` and measures nothing
   about conservativeness. The non-conservativeness of the field is the shear kernel 𝒮 of
   Eq. A15 and belongs to the forward model, not to the geoid. The spec sentence that invited
   the second reading ("how far the wind takes the surface from an equipotential") is
   rewritten. Remove the pseudo-potential from the acceptance output. Return signature
   unchanged: `(r0, closure)`.

4. **GM for the reduction geoid: Null's, read from `lindal_gravity.nc`.** Agreed with the
   report's reading and with its inference that the handoff numbers were computed with the
   modern GM; the 0.2 km is recorded in the spec so it is not later taken for an error.

## Corrections the report did not raise

5. **`U_rigid` does not carry the opposite sign from Lindal Eq. 11; the spec was wrong and
   the docstring repeats it.** Checked against the paper: Lindal writes
   `U = −GM/r + (GM/r) Σ J_2i (R/r)^2i P_2i(sin φ) − ½ ω_III² r² cos² φ` with `g = −grad U`
   (his Eqs. 9 and 11), which is exactly `V − ½ Ω² r² cos² φ` with this package's `V`. Same
   function, same sign, same convention. The v0.6 sentence claiming otherwise was mine; it is
   withdrawn in v0.7 and the `U_rigid` docstring is to say that the two agree. The acceptance
   check that `dU/dr = g` (1.7e-8 agreement) is the right check and stands.

6. **Drop `tol_m` and `max_iter` from `wind_geoid`.** Arguments accepted "for signature
   symmetry" and unused are a trap for the next reader. The spec signature is corrected.

## Notes on section 2

The Newton update, the per-latitude stopping, the hemisphere-by-hemisphere march from each
pole, the pole insertion, and the reading of the array rule for an initial value problem are
all correct. The 45° repeat on the constructed surface (`G_φ` = −1.0726 m/s², ψ = 5.7166°,
`|g_eff|` = 10.769 m/s²) is accepted as the reference figure for that check from now on; the
3247 km explanation of the Step 4 sphere values is right.

## Actions

- Change `radius_at` to the ruled scheme with `φ_c` bracketing and the equator refusal.
- Remove the pseudo-potential from `accept_step5.py` and its output.
- Correct the `U_rigid` docstring (item 5); remove the unused arguments (item 6).
- Rerun the acceptance, then commit Step 5 as accepted.
- Read SPEC_01 v0.7 Step 5 and Step 6 before starting Step 6. Step 6 is to evaluate ψ on the
  surface solved at the current iterate's latitude, not on an interpolated one, which the
  spec's Step 6 wording already implies and item 1 now requires.
