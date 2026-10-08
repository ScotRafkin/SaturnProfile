# SPEC 11. The Lindal-only wind

Version 0.5, 8 October 2026. Author of record: S. Rafkin. **Status: closed. v0.5 records a correction
to Step 2's diagnosis (§2 item 8).**
Depends on SPEC_08, SPEC_09 and SPEC_10.

## 0. Purpose

Build the wind for the epoch of the Lindal anchor from Voyager data alone: the shear from the IRIS
temperature gradients (SPEC_10's fit), the reference wind from the Voyager cloud tracking (the
Ingersoll and Pollard curve the runs use now). Step 1 builds the winds and the author reviews
them; Step 2 runs the transfer under them.

Protocol as SPEC_04 §0. No regression: a new tool, and no change to the model.

## 1. Step 1

**The tool**, in `src/casspian/tools/lindal/`. Inputs: the cloud wind (a kind W file, as
`casspian-wind-from-curve` writes it), `occul_data/lindal/iris_temperatures.csv`, and the run's
gravity and rotation. Output: a kind W file with all three parts.

1. **The reference wind** is the cloud wind, assigned to 39810.7 Pa (398 mbar, a node of the output
   grid), where the shear is zero.
2. **The gradient** at 110, 290 and 730 mbar is `fit` at its defaults, evaluated at the planetographic
   latitudes of the output's planetocentric grid by the wind tool's rule, and converted to per
   planetocentric radian.
3. **The shear** `(du/dln p)` at each level follows from the gradient by the model's own balance (the
   relation the transfer uses), given the wind and its latitude derivative at that level. Since the
   shear changes the wind, steps 3 to 6 are repeated until the wind changes by less than 0.1 m/s.
4. **In the vertical:**
   - between the three levels, linear in ln p;
   - above 110 mbar, two cases: (i) held at its 110 mbar value; (ii) relaxed,
     `s_110 exp(-x / 2)`, with x the height above 110 mbar in scale heights, `ln(110 mbar / p)`;
   - below 730 mbar, relaxed, `s_730 exp(-y / 1)`, with y the depth below 730 mbar in scale heights.
5. **In latitude:** beyond the data's ends the shear tapers linearly to zero at each pole; across the
   equatorial band, poleward of 5 deg it is computed, and between it is bridged by a shape-preserving
   cubic (PCHIP), which stays between the edge values.
6. **The wind** is the reference wind plus the shear integrated in ln p from 398 mbar, on the cloud
   wind's latitude grid and 1 Pa to 1 MPa at 20 nodes per decade. The uncertainty is the cloud
   wind's, carried unchanged.

**The winds.** The tool builds case (i) or (ii), with the top level at 110 or 150 mbar. Step 2 runs
case (ii) at 110 mbar only; the others are not run (§2 item 6).

**Acceptance (`tests/step11_1/accept_step11_1.py`, under `reports/step11_1/`).** For the author:

1. the shear at the three levels against latitude;
2. the shear at 150 mbar drawn as Conrath and Pirraglia's Fig. 3 is, `-du/dz` in m/s per scale
   height against latitude, to compare with their figure;
3. each wind's `u_total(phi, p)`, and `u(p)` at a few latitudes.

## 1a. Step 2, after Step 1 is accepted

**The run.** `lindal_iris_ii_110`, made from `lindal_transfer`, targeted at 26.4 S planetocentric
(the Voyager 2 egress, 31.35 S planetographic). The wind enters the build as a case `lindal_iris` of
`casspian-wind-shear`, whose function is the tool's `construct`; `casspian-run-inputs` builds
`[composition]` before `[shear]`. The tool's own `build` and `[lindal_wind]` section are removed. No
existing suite is rerun for the reorder.

**Acceptance (`tests/step11_2/accept_step11_2.py`, under `reports/step11_2/`).**

1. **The round trip.** Under the wind, the delivered temperature change from the anchor to the
   target at 110, 290 and 730 mbar against the IRIS fit's change between the same two latitudes.
   Reported; a difference over 0.5 K is explained.
2. **For the author:** the delivered temperature profile at the target with the anchor's; its lapse
   rate with the dry adiabat; the pressure identity, in kelvins, level by level.

## 2. Ruled (author, 7 October 2026)

1. **Below 730 mbar**, the shear relaxes over one scale height rather than stopping, so the wind has
   no corner at 730 mbar. Deeper anchors will need cylinders weighed then.
2. **The equatorial band:** at the equator the axis of the model's balance is horizontal, so the
   gradient sets no vertical shear there and the shear's error grows as `1 / sin(phi)`. The shear is
   computed poleward of 5 deg and bridged between.
3. **The target:** the Voyager 2 egress latitude, so the runs are ready when that profile is
   digitized.
4. **The bridge (8 October):** the cubic matched in slope overshot (about 100 m/s per scale height at
   the top level, above any computed value) and sagged at 290 mbar; it is replaced by PCHIP. The
   Step 1 figures are rerun for the author before Step 1 is accepted.
5. **The band stays at 5 deg.** The shear computed between 5 and 12 deg is 6 to 10 times its
   estimated error there.
6. **No sensitivity runs now** (author): variations come later, possibly in the Monte Carlo. Case
   (ii) is the working wind; case (i) is a bound only.
7. **Composition is fixed.** The tool reads the run's composition only for the mean molar mass in
   `dln p / dr`.

8. **The Step 2 diagnosis, corrected (REPORT_13_step2 §2).** Its path integrals ran over the grid
   nodes inside the path and left out the partial cells at the two ends. Integrated from the
   anchor's latitude to the target's exactly, part (A) at 110 mbar is -1.21 K: -1.15 K the slope
   against the fitted value (A1) and -0.06 K the linearization (A2). The conclusions stand.

## 3. Revision history

- v0.1, 7 October 2026: first draft.
- v0.2, 7 October 2026: split at the author's direction: Step 1 builds the winds, Step 2 runs them.
- v0.3, 7 October 2026: accepted; §2 ruled (relaxation below 730 mbar, the equatorial band).
- v0.4, 8 October 2026: rulings on the Step 1 report: PCHIP bridge, one run in Step 2, the build as
  the coding agent proposed.
- v0.5, 8 October 2026: closed; §2 item 8 records the corrected diagnosis of Step 2.
