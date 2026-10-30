# REVIEW 14, Step 1. Hypothesis 2

Reviewing agent, 30 October 2026. Report: `reports/REPORT_14_step1.md`. Specification:
`docs/specs/SPEC_14_Reversed_Shear.md` v0.3.

**Accepted. SPEC_14 closes.**

- **Check 1, ruled (author, 30 October): the report's option 1.** The wind below 110 mbar differs
  from the base's by 4.5e-3 m/s and the delivered temperature by 0.03 K, a remainder of the
  construction's stopping rule. The diagnosis shows it falling to 1.8e-7 m/s as the tolerance
  falls, so it belongs to the iteration and not to the case. Neither number matters. The suite asks
  for the difference to be below `TOLERANCE_MS`; with that wording it passes, 3 of 3. No code
  change, and the existing runs are untouched.
- **The case.** The closed form agrees with a numerical integral to 1.3e-10 m/s, and below the top
  level the code is unchanged. `shear.png` shows f(p) and the two columns as specified.
- **Against the egress at 31.2 S (26.06 S planetocentric):** +0.20, -0.77, -2.17 and +2.88 K at
  10, 4.5, 1 and 0.3 mbar, against -1.90, +0.10, +5.24 and +15.10 K for the base. Every value is
  inside the egress's sigma_T. Density aloft moves from -10.0 to -1.8 percent at 0.3 mbar. As
  `profiles.png` shows, case (iii) follows the egress from 10 mbar to the top, where the base
  stayed on the anchor. At 1 mbar the result is 1.2 K colder than the scaled expectation, and at
  0.3 mbar 1.2 K warmer. As expected, a line in ln p does not fit both levels.
- **Unchanged below 110 mbar:** +5.09, +2.30 and +0.40 K at 110, 290 and 730 mbar. The 110 mbar
  difference is still the radio and IRIS disagreement.
- **The winds.** The equatorial wind reaches 472 m/s at 0.2 mbar and 1052 m/s at 0.01 mbar, the
  accepted consequence of holding f (SPEC_14 §2 item 1). Below 1 mbar the change from the base is
  under about 60 m/s.
- **Lapse rate** at the target: -0.55 to +0.44 K/km above 110 mbar, stable throughout and well
  short of the dry adiabat.
- **Pressure identity:** RMS 0.064 K over 66 levels, as the base; largest 0.168 K at 3.99 mbar.

**Order of work:** commit and push, with no regression.
1. The acceptance commit: the case in `lindal_wind.py` and `shear.py`; the run's two control
   files; `tests/step14_1/` with check 1's new wording; this report and review; SPEC_14 v0.3;
   `STATE.md`; the status note.
2. The hash in SPEC_14's row of `STATE.md`. Push.
