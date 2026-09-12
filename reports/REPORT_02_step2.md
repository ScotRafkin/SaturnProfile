# REPORT 02, Step 2. The frozen constants, `phi_c` and `r0`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.4 (accepted by
the author), Step 2, against `SPEC_00_Architecture_and_Data_Files.md` v0.12 and the closed
`SPEC_01_Lindal_Tool_Chain.md` v0.18 (commit `ece58d2`). Status: reported, awaiting review.
**Seven of seven acceptance checks pass.**

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product, so the directory sweep after the acceptance commit will find nothing to rebuild.

---

## 1. What was built

**`src/casspian/refrac/anchor.py`**, built from `lib` only.

- `freeze_anchor(inputs, manifest)` returns a frozen `FrozenAnchor` record: `phi_c`, `psi`, `r0`,
  `u(phi_c)`, the fixed point iterates, iteration count and last step, the tolerance, the anchor
  rule, surface, quantity and radius, the north polar start, both polar radii, the asymmetry,
  the anchoring residual, and the no wind `phi_c`, `psi`, `r0` and iterates. Angles in radians.
- `wind_of_latitude(wind)` is the kind W reference level column as a callable in radians,
  linear in planetocentric latitude on the file's 0.5 degree grid.
- `AnchorConvergenceError` is raised when either fixed point exhausts `max_iterations`.

The sequence is SPEC_02 Step 2, items 1 to 5:
1. `u_of_phi` from kind W.
2. The wind geoid by `lib.geoid.wind_geoid`, anchored on `radius_polar_m` at 10000 Pa by
   `mean_polar_radius`.
3. The fixed point by `lib.latitude.planetocentric_fixed_point` from the kind T label, 36.3
   degrees. Every surface evaluation is a full anchored march with the iterate inserted as a
   node. The seed is the 10000 Pa oblateness.
4. `r0` from one final march with `phi_c` a node.
5. The no wind pair on `lib.geoid.reference_geoid`, same seed and tolerance.

**`src/casspian/lib/control.py`, conformed to SPEC_00 v0.12.** The `[sensitivity]` section is
removed from the vocabulary and from `ReductionManifest`, so a manifest that carries one is now
refused as an unknown section. The one line of `reports/step02_1/accept_step02_1.py` that printed
it was updated; that suite still passes (section 4).

## 2. Decisions

1. **`psi` is evaluated at `phi_c` on the final march, not taken from the fixed point.**
   `lib.latitude.planetocentric_fixed_point` returns the tilt at the previous iterate, which
   makes `phi_c + psi = phi_g` exact. The specification asks for "`psi` at `phi_c`", and the
   final march supplies the radius there at no extra cost, so the record carries the tilt at the
   frozen latitude. The two differ by 3.15e-8 degrees, which is the measured
   `phi_c + psi - phi_g`. Flagged in case the author prefers the exact identity.
2. **`convergence_m` is the geoid tolerance, for both surfaces.** SPEC_00 section 7.1 gives the
   key without saying what it controls. It is passed as `tol_m` to `wind_geoid`, where it stops
   the secant on the north polar start, and to `reference_geoid`, where it stops the Newton
   step. At 1 m it still gives an anchoring residual of 2.6e-6 m, because the map from start to
   outcome is nearly affine and the secant overshoots the tolerance by six orders.
3. **The record's geoid numbers all come from the final march**, the one `r0` was read from, so
   they describe a single surface. With `phi_c` added as a node they differ from a march
   without it only at round-off.
4. **The wind callable adds no rule at the poles.** Kind W carries both poles as nodes with
   exactly zero wind, and `read` checks that, so linear interpolation is zero there by the
   file's own values.
5. **The only degree to radian conversions are in `refrac.anchor`**: the manifest tolerance and
   the kind T label. The record stays in radians. Degrees are for the product (Step 4).
6. **Both fixed points are checked for convergence** and refused if the budget is exhausted.
   `lib.latitude` itself returns silently after `max_iter`.

## 3. Acceptance results

Run by `reports/step02_2/accept_step02_2.py`; full output in `reports/step02_2/output.txt`.

| Check | Measured |
|---|---|
| No wind reference values, Null's GM | **Pass.** GM 3.7929085e16 from the file (Null 1981). `phi_c` = 30.818189 deg (30.8182 +- 0.001), `psi` 5.481811 deg, `r0` = 58,453.120 km (58,453.1 +- 0.5), all as SPEC_01 Step 6. |
| Wind included values, reported | **Pass** on the one tested figure. `phi_c` = **30.804949 deg**, 0.013239 deg below the no wind value (within 0.02). `psi` = 5.495051 deg. `u(phi_c)` = 2.175 m/s. **`r0` = 58,519.883 km**, 66.763 km above the no wind value. See finding 1 for the 1.9 km against "near 58,518". |
| Convergence in at most eight iterations | **Pass.** Wind: 5 iterations, last step 5.48e-7 deg against the 1e-6 deg manifest tolerance. Iterates 30.8524080, 30.8022287, 30.8051059, 30.8049403, 30.8049498, 30.8049493. No wind: 5 iterations, last step 9.16e-7 deg, iterates 30.8524080, 30.8157725, 30.8183598, 30.8181767, 30.8181897, 30.8181888. |
| Anchoring residual below 1e-3 m | **Pass.** -2.62e-6 m. |
| Polar radii and asymmetry to 0.1 km | **Pass.** North 54,423.627 km (54,423.6), south 54,452.373 km (54,452.4), asymmetry 28.745 km (28.7). |
| Beyond the specification: stability and inversion | **Pass.** Halving the march step to 0.025 deg changes `r0` by less than 1e-4 m. The direct relation at `phi_c` on the wind surface returns 36.2999999685 deg, error -3.2e-8 deg. |
| Beyond the specification: refusals | **Pass.** `max_iterations = 2` is refused with the step size named (2.9e-3 deg). A manifest with `[sensitivity]` appended is refused as an unknown section. |

## 4. Findings

1. **The expected 58,518 km is the wind surface read at the no wind latitude.** The acceptance
   splits the 66.763 km between the two frozen radii. On the same anchored march, read at the
   no wind `phi_c` of 30.818189 deg, the radius is **58,518.582 km**: a dynamical height of
   65.462 km at fixed latitude, which matches the specification's "near 58,518" and "about
   65 km". The frozen latitude is 0.013239 deg further south. The surface is rising there at
   5,629.7 km per radian (Eq. B3 at `phi_c`), which adds the remaining 1.301 km. Both terms are
   measured, and their sum closes to the reported `r0`. So nothing disagrees: the "near 58,518"
   omits the latitude shift the wind itself causes.

2. **For the Step 3 markup: the Eq. B3 slope at the anchor is -5,630 km per radian, not about
   -6,700.** SPEC_02 v0.4 Step 3 quotes `dr0/dphi_c = r0 G_phi / g` as "about -6,700 km per
   radian, so roughly 22 km for Lindal". Evaluated at the frozen pair it is -5,629.7 km/rad. An
   ellipse check agrees: through the fitted 100 mbar radii (60,367 and 54,438 km) the slope at
   `phi_c` is -5,587 km/rad, and the ellipse `dphi_c/dphi_g` is 0.923 (the specification's
   0.93). With those and the 0.2 deg label uncertainty, the label term is about 18 km
   (0.923 x 3.49e-3 rad x 5,630 km/rad = 18.1 km), not 22. It would still be larger than the
   10 km anchor term.
   Step 3 will measure both partials. This is noted now so the expected value can be marked
   before that step's acceptance is run. Nothing in Step 2 depends on it.

3. **The wind at the anchor is 2.175 m/s by the file, against 2.092 m/s on the Step 7 curve.**
   The specification's u_of_phi is linear interpolation on the file's 0.5 deg grid. The Step 7
   figure was the PCHIP curve evaluated at 36.3 deg planetographic. The 0.08 m/s difference
   comes from the interpolant. It enters `r0` through `psi` and the slope only, at a level far
   below the 0.1 km the polar radii agree to.

## 5. Regression

The ten SPEC_01 suites and the SPEC_02 Step 1 suite were rerun after the `[sensitivity]` removal
(`reports/step02_2/regression.txt`): all pass, 6, 14, 8, 5, 7, 7, 6, 8, 7 and 6 for SPEC_01 and
6 of 6 for Step 1.

No product under `occul_data/` was rebuilt, so all still carry `ece58d2`. No em dash or en dash
appears in any file written in this step.

## 6. Next step

Step 3, number density, absolute radius and refractivity with the uncertainty companions, is not
started. It waits for the review of this report.
