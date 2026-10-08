# REVIEW 13, Step 3. Before acceptance

Reviewing agent, 8 October 2026. Report: `reports/REPORT_13_step3.md`. Specification:
`docs/specs/SPEC_13_SL_Wind.md` v0.5, §2a.

**Accepted once the author has viewed the figures. With it Steps 1 and 2 are accepted, and SPEC_13
closes.**

- **The smooth window.** The running maximum over ±4 deg, then a Hann average, is never narrower
  than the adaptive width and has a continuous slope, by construction. The derivative's spikes are
  gone. The fit's RMS about the points is unchanged to 0.001 K. The construction is accepted in
  place of the suggested one, for the reason the report gives.
- **What remains narrow is in the data.** The features near -13 deg are in the slope too, and
  they are amplified by `1 / sin(phi)`. The wind columns near -48, -32 and -13 deg at the top of
  the grid are the derivative's difference from the slope over a few degrees. The relaxation
  carries that difference upward and roughly doubles it.
- **The fitted temperatures move in the IRIS gap.** At 31.2 S, which lies in the gap from -27.7 to
  -38.5 deg, the smooth window changes the fit by up to 1.3 K at 730 mbar. Spanning the gap, the
  IRIS value at the target is uncertain by about that much from the choice of window alone.
- **The round trip** closes to the linearization at 110 mbar (-0.11 and -0.06 K). The model's
  remaining part is the equatorial bridge: 0.6 to 0.7 K at 290 mbar and 0.7 to 0.9 K at 730 mbar.
- **Against the egress:** +5.0, +2.3 and +0.4 K at 110, 290 and 730 mbar. At 110 mbar the
  difference is still the IRIS and radio disagreement.
- **The data property files and SPEC_00 v0.24** are as ruled.

**SPEC_14 starts from `lindal_iris_v_sl_ii_110`**: the derivative, the smooth window, and the
Sanchez-Lavega cloud wind.

Order of work: the author's view; then commit and push, with no regression.
- SPEC_13's code, data property files, runs' control files and tests;
- the reports and the reviews of Steps 1 to 3;
- SPEC_00 v0.24, SPEC_11 v0.5 and REVIEW_11_step2;
- the SPEC_13 rows in `STATE.md`, with the commit hash.
