# REVIEW 17, Step 1. What the Voyager-epoch observations require of the wind

Reviewing agent, 5 November 2026. Report: `reports/REPORT_17_step1.md`. Specification:
`docs/specs/SPEC_17_Voyager_Constraints.md` v0.3.

**Accepted. SPEC_17 closes.**

**Part A, the 100 mbar surface.**
- **The digitization** is sound: altitude to 0.24 km and latitude to 0.27 deg RMS on the ticks.
  Three dots sit 0.1 to 0.8 deg from their occultations' latitudes. They are recorded as Lindal
  drew them, and the slope correction (0.5 km at most) changes no conclusion.
- **The anchor is reproduced to 0.8 km,** inside Lindal's 3 to 4 km fit, so the reference geoid
  and the anchor isobar agree with Lindal's.
- **South of the equator the delivered surface is too high:** +8.4 km at 31.2 S and +7.1 km at
  74.2 S, outside Lindal's fit. At 3.0 S it is +3.3 km; north of the equator it follows Lindal's
  traced geoid. Both hypotheses agree to 0.1 km, since their shear is the same below 110 mbar.
- **Reading:** the Sanchez-Lavega cloud wind, with the IRIS shear between 398 and 100 mbar, bulges
  the southern 100 mbar surface more than Lindal's altitude-independent Voyager wind. This is a
  constraint on the reference wind and the shear below 100 mbar together, which the fitted wind
  weighs. Pioneer 11's point is drawn and enters nothing, as ruled.

**Part B, the operator.**
- **`balance_weight`** is construct's balance linearized in the shear, on construct's own geometry.
  It is the observation operator. Accepted, at about 9 minutes per run's 66 levels.
- **The linearization holds:** delivered over predicted 0.989 at zero amplitude and 0.969 for a
  4 K correction. Of the 1.1 percent at zero amplitude, 0.5 percent is the wind's own change,
  explained. The 0.6 percent not explained is 0.01 K on a 2 K signal: immaterial, and it is not
  pursued. The second order, 1 percent per 2 K, is w's temperature held fixed; one
  relinearization removes it for large corrections. Nothing reaches north of the perturbed
  segment (1.6e-6 in ln T).
- **The tapers,** inside the segment in latitude and outside the band in pressure: accepted as
  the coding agent read them.

**What the segments ask.**
- **Segment 3 (26.1 S to 71.1 S)** asks little: +1 to +2 m/s per unit ln p at most, steady with
  height. Its -0.74 at 730 mbar is the Fig. 5 depth bias alone.
- **Segment 2 (2.4 S to 26.1 S) at 110 mbar asks -7 m/s per unit ln p.** This is the 31.2 S
  misfit, the radio and IRIS disagreement at 110 mbar, now located in latitude and expressed as
  shear.
- **Segment 1 (30.8 N to 2.4 S)** asks +1 to +3 below 110 mbar, but it is the weakest constraint:
  w passes through zero in it, so its weight integral is 3 to 5 times smaller than segment 3's.
- **Above 30 mbar, segments 1 and 2 swing between -13 and +21 with height.** That is the
  equatorial layering read as shear. A smooth fitted wind should not chase it; those levels enter
  the fit with a large error.

**Order of work:** commit and push, with no regression.
1. The acceptance commit: `balance_weight`; `occul_data/lindal/isobar_100mbar.csv`;
   `tests/step17_1/`; this report and review; SPEC_17 v0.3; SPEC_18 v0.2; `STATE.md`; the status
   note.
2. The hash in SPEC_17's row of `STATE.md`. Push.
3. Then SPEC_18, Step 0.
