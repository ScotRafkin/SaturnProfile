# REPORT 01, Step 4. `lib.gravity`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.5, Step 4, against
`SPEC_00_Architecture_and_Data_Files.md` v0.7. Status: reported, awaiting review.

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product, so there is nothing to rebuild afterward.

**Step 3 closed first.** The three review changes were applied: `index` added to the provenance
vocabulary and carried by `degree`; `J`, `J_uncertainty`, `GM_m3s2` and its companion changed
from `inferred` to `derived`; and the author's `nan` in `null1981.toml` picked up. Step 3 was
committed as `0cb4c38`, its review resolutions recorded in `0e03e2b`, and both tools rerun, so
`lindal_gravity.nc` and `lindal_rotation.nc` carry the clean commit `0e03e2b` and
`J_uncertainty` reads `[1.8e-05, 3.8e-05, nan]`. All four earlier acceptance suites were rerun
afterward and still pass: 6 of 6, 14 of 14, 8 of 8, 5 of 5.

---

## 1. What was built

**`src/casspian/lib/gravity.py`.** Pure functions, NumPy in and out, no file access, no module
level Saturn quantity. `legendre_even` and `legendre_even_derivative`; `potential_V`;
`g_newton`; `G_phi_newton`; `omega_abs`; `g_eff_radial`; `G_phi_eff`; and `g_eff_vector`
returning both components, the magnitude, and the angle `psi`. All accept scalars or
broadcastable arrays.

The Legendre polynomials come from the joint recurrence

    (n + 1) P_{n+1} = (2n + 1) x P_n - n P_{n-1}
    P'_{n+1} = x P'_n + (n + 1) P_n

rather than from the closed form `P'_n = n (x P_n - P_{n-1}) / (x^2 - 1)`, which is singular at
`x = +-1`. That is exactly the pole, where the acceptance evaluates, so the closed form would
have produced a division by zero at one of the two points the step is checked on.

## 2. Decisions

1. **The centrifugal term is written as `Omega_abs^2 r cos^2(phi_c)`** rather than in Lindal's
   form `(2/3) w^2 r (1 - P_2(sin phi_c))`. They are the same quantity, since
   `1 - P_2(sin phi) = (3/2) cos^2(phi)`, and the identity is recorded in the docstring. The
   trigonometric form needs no Legendre evaluation and states the physics directly.

2. **`omega_abs` returns `Omega` at the pole rather than an infinity.** `u / (r cos phi_c)`
   divides by zero there. A zonal wind is zero at the pole by definition, so the limit is
   `Omega`, and the function returns it instead of a NaN that would poison the polar gravity.
   The acceptance evaluates at the pole, so this path is exercised.

3. **A loop over degrees is not a loop over array elements.** SPEC_00 section 3.1 forbids
   Python level loops over array elements in `lib`. The Legendre recurrence loops at most
   twelve times regardless of grid size, and every operation inside it is whole array. Noted
   because it is the first place in the package where the rule could be misread.

4. **Every Saturn number in the acceptance script is read from the static transcriptions**, per
   SPEC_00 section 2.1: the harmonics and `check_values` from `null1981.toml`, the modern GM
   from `iess2019.toml`, the period from `system_iii.toml`, the 1 bar radii from
   `lindal_scalars.toml`. The Jupiter cross check is the one exception, since those values
   appear in no static file; they are quoted from SPEC_01 Step 4 and labeled as such in the
   script.

## 3. Defects and questions raised

1. **The sign convention for `G_phi` is stated two ways, and Step 6 depends on which one
   holds.** SPEC_01 Step 4 defines `G_phi_newton` as `-(1/r) dV/dphi`, which is the component
   along **increasing** planetocentric latitude. On an oblate planet the bulge pulls a mid
   latitude point equatorward, so in the northern hemisphere that component is **negative**;
   the implementation returns `-1.056234 m/s2` at 45 degrees. The acceptance text for the same
   quantity says it is "positive (toward the equator)". Both describe the same physical
   direction. Only the sign convention differs, and there is no arithmetic error on either
   reading.

   This matters at Step 6, where `psi` must satisfy `phi_c = phi_g - psi` with `psi = +5.4815`
   degrees at the Lindal anchor. Taking `psi = arctan(G_phi / g)` with `G_phi` as the
   specification writes it gives a negative `psi` and Step 6 would have to reintroduce the sign.
   **`g_eff_vector` therefore returns `psi = arctan(-G_phi / g)`**, which is positive in the
   north, and the minus sign is documented at the point of use. `G_phi_eff` itself returns the
   quantity the written formula defines, unaltered.

   **The ruling needed:** either `G_phi_eff` keeps the written convention and `psi` carries the
   minus sign, as implemented, or `G_phi_eff` is redefined as equatorward positive and `psi`
   loses it. Kind N stores `psi_deg` as a scalar, so whichever is chosen becomes visible in a
   product file at Step 9.

   *Resolved by `REVIEW_01_step4.md`, and the implementation stands. `G_phi` is the component
   along increasing planetocentric latitude in both hemispheres, never equatorward positive.
   The deciding argument is Eq. B3, not Step 6: `g dr0/dphi = r0 G_phi` holds as written with
   this convention and would need a hemisphere dependent sign with the other. The word
   "reintroduce" above overstates the case, and the docstring has been corrected: the minus
   sign in `psi = arctan2(-G_phi, g)` is part of the definition of `psi` as the tilt of the
   local vertical toward the pole, not a correction applied afterward. `arctan2` is kept so
   the southern hemisphere needs no case split. The reviewer independently reproduced
   -1.0562344 m/s2, psi = 6.35218 degrees and |g| = 9.54664 m/s2 from the closed form Legendre
   polynomials, and checked the convention against the manuscript Eqs. A1, A5 and B3, against
   handoff 9A.2 and 9A.5, and against Lindal's own Eqs. 2, 3, 5 and 8. Paper, manuscript and
   code agree; the single inconsistent statement is a shorthand line in the handoff, which is
   a handoff defect and not this chain's to fix. SPEC_01 v0.6 Step 4 now states the
   convention.*

2. **The acceptance tolerance of 5e-4 m/s2 on Null's own GM is exceeded at the pole.** The two
   GM values differ by 5.592e-05 relative, and the shift in `g` is `g` times that: 5.99e-04 at
   the equator and 6.79e-04 at the pole. All three exceed the stated 5e-4. Nothing is wrong
   with the computation; the tolerance was set without allowing for `g` being 12.14 at the
   pole. **Proposing the bound become 1e-3 m/s2**, which the results satisfy with room. Reported
   rather than adjusted, and the check is marked passing against 1e-3 with the discrepancy
   printed.

   *Resolved by `REVIEW_01_step4.md`: the bound is 1e-3 m/s2, and SPEC_01 v0.6 Step 4 now
   states it with the reason.*

3. **"9.29 and 11.68" pairs unlike numbers.** The wrong coefficient reproduces both, so the
   check discriminates exactly as intended. But 9.29 is the **no wind** equatorial value, while
   the correct number it is set against, 8.96 from Lindal Table II, is the **with wind** value.
   The with wind wrong value is 9.14. The comparison the reader is meant to make is 9.29 against
   9.102, or 9.14 against 8.951. Worth a word in the specification so the next reader is not
   left comparing 9.29 with 8.96 and finding a discrepancy that is only the wind.

   *Resolved by `REVIEW_01_step4.md`: SPEC_01 v0.6 Step 4 now pairs 9.29 with 9.102 at u = 0,
   9.14 with 8.951 at u = 450, and 11.68 with 12.137.*

4. **The GM the specification quotes is in no static file.** Step 4 says `GM = 3.7931206e16`;
   `iess2019.toml` carries `3.7931206234e16`. SPEC_00 section 2.1 forbids hard-coding a value
   that appears in a static file, so the script reads the file. The relative difference is
   6.2e-09 and the effect on `g` is 6e-08 m/s2, far below any tolerance here.

   *Resolved by `REVIEW_01_step4.md`: reading it from the file is the correct reading of
   SPEC_00 section 2.1, and SPEC_01 v0.6 Step 4 now quotes 3.7931206234e16 and says so.*

5. **"At 45 degrees at 1 bar radius" does not say which 1 bar radius.** The geodesy file gives
   an equatorial and a polar radius for that surface, and the true geoid radius at 45 degrees
   lies between them; that radius is Step 5's product and does not exist yet. The check uses the
   equatorial radius and says so. The number reported for `G_phi_eff` at 45 degrees is therefore
   a value on a sphere of that radius, not on the geoid, and should be read as an order of
   magnitude and a sign rather than as a physical constant. Once Step 5 exists this check could
   be repeated on the constructed surface.

   *Resolved by `REVIEW_01_step4.md`: the equatorial 1 bar sphere is accepted for this check,
   and Step 5 repeats it on the constructed no-wind geoid at r_ref(45 degrees), which its
   acceptance now includes.*

## 4. Acceptance results

Run by `reports/step4/accept_step4.py`; full output in `reports/step4/output.txt`.

| Check | Result |
|---|---|
| Legendre polynomials and derivatives to degree 12, vectorized, poles included | **Pass.** Worst absolute departure from the closed forms over `x` in [-1, 1] including both endpoints: **1.110e-16**. `P2(1) = 1.0`, `dP2(1) = 3.0`, evaluated where the closed form for the derivative is singular. |
| `g_eff_radial` at the equator with u = 450 m/s is 8.951 +- 0.005 | **Pass. 8.9512.** Lindal Table II gives 8.96. |
| The same with u = 0 is 9.102 | **Pass. 9.1020.** |
| At the pole, 12.137 | **Pass. 12.1371.** Lindal Table II gives 12.14. |
| With Null's own GM, the same to within tolerance | **Pass against 1e-3**, not against the stated 5e-4: 5.99e-04, 5.99e-04, 6.79e-04. See section 3, item 2. |
| A deliberately wrong coefficient (5, 9, 13) must give 9.29 and 11.68 | **Pass. 9.29** (no wind) **and 11.68.** The departures from the correct values are 0.19 at the equator and 0.46 at the pole, two to four percent, which is the size of error that looks plausible and is why the check exists. |
| `G_phi_eff` is zero at the equator and the pole to round-off | **Pass.** `-0.000e+00` and `-1.282e-16` m/s2. |
| `G_phi_eff` at 45 degrees, magnitude reported | **Pass, with the sign question of section 3, item 1.** `-1.056234 m/s2`, magnitude **1.056234 m/s2**, equatorward. `psi` there is **6.3522 degrees** and the effective gravity magnitude is 9.5466 m/s2. |
| Jupiter cross check: 23.116 and 27.015 | **Pass. 23.1156 and 27.0148.** Lindal gives 23.12 and 27.01. A second planet with different harmonics, radii and rotation through the same code path is independent evidence that the `(l+1)` coefficient is right. |
| Beyond the specification: broadcasting | **Pass.** An array call over five latitudes is bit identical to the five scalar calls, and a scalar call returns a scalar. |

Seven of seven pass. No em dash or en dash appears in any file written in this step.

## 5. Next step

Step 5, `lib.geoid`, is not started and is not to be begun until this report has been reviewed.
Its no-wind acceptance numbers (equatorial radius 60,244 km, `r(31.0)` = 58,435 km) are the
first that depend on `lib.gravity` being right, and its wind geoid needs the `G_phi` sign
question of section 3, item 1 settled, since Eq. B3 integrates `g dr0/dphi = r0 G_phi`.
