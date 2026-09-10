# REVIEW 01, Step 4. `lib.gravity`

Review of `reports/REPORT_01_step4.md` against SPEC_01 v0.5 Step 4. Reviewer: S. Rafkin,
10 September 2026. **Disposition: accepted.** SPEC_01 is revised to v0.6 to record the
rulings below. Commit, then proceed to Step 5.

## Independent check

`G_phi_eff` at 45° N on a sphere of radius 60,268 km with u = 0, computed separately from
the closed form Legendre polynomials with the same harmonics and GM: −1.0562344 m/s², ψ =
6.35218°, magnitude 9.54664 m/s². Matches the report to every digit shown. The equatorial and
polar gravities and the Jupiter cross check were verified earlier in the project by the same
route and agree.

## Rulings on section 3

1. **Sign of `G_φ`: the implementation stands.** `G_φ` is the component along increasing
   planetocentric latitude, as the formula defines it, negative in the north. It is not to be
   redefined as equatorward positive. The deciding argument is not Step 6 but Eq. B3:
   `g dr₀/dφ = r₀ G_φ` gives dr₀/dφ < 0 north of the equator with this convention and no
   hemisphere dependent sign anywhere, whereas "equatorward positive" would flip meaning at the
   equator and would break B3 in the south. ψ = arctan2(−G_φ, g) is then simply the definition
   of the tilt of the local vertical toward the pole, not a correction; the docstring's phrase
   "the minus sign is what makes ψ ..." is fine, but the word "reintroduce" in the report
   overstates it. Keep `arctan2`, which the code already uses, so that the southern hemisphere
   is handled without a case split. SPEC_01 v0.6 states the convention in Step 4. The manuscript
   is to be checked for consistency between Eq. A5 and its statement of Lindal Eq. 5; that goes
   on the manuscript notes list, not on the coding agent.

2. **Tolerance on Null's GM: 1e-3 m/s².** Agreed; the shift is g times 5.6e-5 and the polar
   value was not allowed for. Spec revised.

3. **9.29 pairs with 9.102, not 8.96.** Agreed; the spec now states the wrong values with
   their correct counterparts and wind state (9.29 vs 9.102 at u = 0; 9.14 vs 8.951 at u = 450;
   11.68 vs 12.137).

4. **GM from `iess2019.toml`.** Correct reading of SPEC_00 §2.1; the spec now says so.

5. **45° on the equatorial 1 bar sphere.** Accepted as an order of magnitude and sign check.
   Step 5 is to repeat it on the constructed no-wind geoid at r = r_ref(45°) and report the
   value; the spec says so.

## Notes on section 2

The Legendre recurrence for the derivative, the polar limit of `omega_abs`, the trigonometric
form of the centrifugal term, and the reading of the array rule are all correct and need no
change. One remark on the polar limit: `np.where(cos_phi != 0.0, ...)` catches an exact zero
only, and `np.cos(np.pi/2)` is 6.1e-17, not zero, so a caller passing `np.radians(90.0)` takes
the division path and gets `u / (r × 6e-17)`. With u = 0 at the pole that is 0/6e-17 = 0 and
nothing is harmed, which is why the acceptance passes, but a nonzero u at a numerically
near-polar latitude would give an enormous `omega_abs` rather than an error. Not a defect for
any input the Lindal chain produces (the wind tool sends u = 0 at the pole); noted so that
Step 7 makes sure it does, and so that no one relies on the guard for a nonzero u.

## Actions

- Commit Step 4 as accepted. No product to rebuild.
- Read SPEC_01 v0.6 Step 4 (the convention paragraph) and Step 5 before starting Step 5.
- Step 5 adds the 45° repeat on the constructed surface to its acceptance.

## Addendum: manuscript check of the `G_φ` convention

Checked against `CASSPIAN_AtmosphericModel.docx` (current draft) after the review above was
written. Eq. A1 writes `g_eff = −g r̂ + G_φ φ̂`, and Eq. A5 writes
`G_φ = −(1/r) ∂V/∂φ − Ω_abs² r cos φ sin φ`, so the manuscript defines `G_φ` as the component
along increasing planetocentric latitude, which is the convention ruled above and the one the
code implements; Eq. B3 uses the same `G_φ`. The manuscript nowhere states ψ or the
planetographic conversion, so there is nothing in it to correct. The handoff, §9A.2, states
Lindal Eq. 5 as `ψ = arctan(g_φ / g_r)` with `g_r` outward, which with §9A.5's `g_r = −g` is
`ψ = arctan(−G_φ / g)`, again the implemented form; its shorthand iteration line
`φ_c ← 36.3° − arctan(G_φ/g)` drops the sign of `g_r` and is the only inconsistent statement.
That line is a handoff defect, not a manuscript one, and the spec's Step 6 wording
(`φ_c ← φ_g − ψ(φ_c)` with ψ from Step 4) is unaffected. Lindal's own conventions were then
checked against the paper (`Lindal_1985AJ.pdf`, Appendix, Eqs. 2, 3, 5 and 8). His Eq. 3 is
`g_φ = −(GM/r²) Σ J_2i (R/r)^2i dP_2i(sin φ)/dφ − (1/3) ω² r dP_2(sin φ)/dφ` along `u_φ`, the
unit vector in the latitudinal direction; since `dP_2(sin φ)/dφ = 3 sin φ cos φ`, this is
term for term the `G_φ` of Eq. A5 and of the code. His Eq. 2 is `g_r = −GM/r² + (GM/r²)
Σ (2i+1) J_2i (R/r)^2i P_2i + (2/3) ω² r [1 − P_2]`, outward positive, which is `−g` of
Eq. A3, and it confirms the `(degree + 1)` coefficient directly, `2i` being the degree. His
Eq. 5, `ψ = arctan(g_φ / g_r)`, is therefore `arctan(−G_φ / g)`, and his Eq. 8 for the geoid,
`dr = −tan ψ r dφ`, is `g dr/dφ = r G_φ`, which is Eq. B3. Paper, manuscript, and code agree;
the handoff's shorthand line is the single exception.
