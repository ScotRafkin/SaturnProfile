# SPEC 14. Hypothesis 2: the shear reversed aloft

Version 0.3, 30 October 2026. Author of record: S. Rafkin. **Status: closed, Step 1 accepted
(REVIEW_14_step1).**
Depends on SPEC_13 and SPEC_R1.

## 0. Purpose

Hypothesis 1 relaxes the shear above 110 mbar to zero, so above about 3 mbar the delivered profile
equals the anchor, up to 15 K warmer than the egress. The radio's egress minus ingress changes sign
near 4 to 5 mbar: warmer to the south below, colder above. Hypothesis 2 reverses the shear there.
**It is tuned to the egress**, so the egress then checks consistency, not independence.

The base is `lindal_iris_v_sl_ii_110` (REVIEW_13_step3). Hypothesis 1 and its runs are kept.

Protocol as SPEC_04 §0. No regression: a new case in the wind tool, and no change to the model.

## 1. Step 1

**Deliverables.**

1. **Case (iii), `above_top = "reversed"`,** in `tools/lindal/lindal_wind.py` beside `"held"` and
   `"relaxed"`. Above the top level the top level's shear is multiplied by f(p), with
   `x = ln(p_hold / p)` and `x_c = ln(p_hold / p_cross)`:
   - f = 1 from the top level up to `p_hold`;
   - f = 1 - x / x_c above `p_hold`, through zero at `p_cross`, up to 0.2 mbar (20 Pa, the
     anchor's top);
   - above 0.2 mbar, f held at its 0.2 mbar value.

   The integral stays closed form. Below the top level nothing changes. Two control values, read
   for this case only: `hold_Pa`, first guess 1500.0, and `cross_Pa`, first guess 450.0. The 0.2
   mbar level is a constant of the case, beside `ABOVE_SCALE_HEIGHTS`.
   `vertical_structure` names the case, its two values and that it is tuned to the egress.
2. **The run `lindal_iris_v_sl_iii_110`**, made by `casspian-new-run --from lindal_iris_v_sl_ii_110`,
   differing in `[shear]` only. Its build file's comments say the case is tuned to the Voyager 2
   egress.

**Acceptance (`tests/step14_1/accept_step14_1.py`, under `reports/step14_1/`).**

1. **Below the top level the wind is the base's,** and so is the shear at the three levels, to
   within the construction's tolerance (v0.3). The largest change in the delivered temperature at
   and below 110 mbar is reported; over 0.1 K is explained.
2. **The pressure identity** of the run, in kelvins, level by level.
3. **For the author:**
   - f(p) and the shear at the anchor's and the target's latitudes against pressure, both cases;
   - the delivered temperature and density against pressure at 31.2 S (26.06 S planetocentric),
     the base and the new run, with the anchor, the egress and its sigma_T, and the IRIS values,
     as SPEC_13's `profiles.png`; a second panel, delivered minus egress with the sigma_T band;
   - a table of delivered minus egress, both runs, at 110, 290 and 730 mbar and at 10, 4.5, 1 and
     0.3 mbar;
   - `u_total(phi, p)` of the new wind, and its difference from the base. The largest |u| at
     each level above 110 mbar and the equatorial wind at 10, 1 and 0.2 mbar are reported;
   - the delivered lapse rate above 110 mbar, both runs.

## 2. Ruled (author, 22 October 2026)

1. **Above 0.2 mbar** f is held, as planned. The equatorial wind there, about 940 m/s at 1 Pa
   against 270 in the base, is accepted: this is one hypothesis of many, its purpose is to work
   the pipeline with an input wind, and later hypotheses bring more physical boundary conditions.
2. **The first guesses stand.** Scaled from the +4.7 K the base carries at 110 mbar, delivered
   minus egress is expected near -0.6, -1.0 and +4.1 K at 10, 1 and 0.3 mbar. A line in ln p
   cannot fit both 1 and 0.3 mbar; later hypotheses and data address that.

## 3. Revision history

- v0.1, 22 October 2026: first draft.
- v0.2, 22 October 2026: accepted; §2 ruled (f held above 0.2 mbar, the first guesses).
- v0.3, 30 October 2026: closed; check 1 to within the construction's tolerance (author, on
  REPORT_14_step1 §6).
