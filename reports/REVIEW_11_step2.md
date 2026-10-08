# REVIEW 11, Step 2. The Lindal anchor under the Lindal-only wind

Reviewing agent, 8 October 2026. Report: `reports/REPORT_11_step2.md`. Specification:
`docs/specs/SPEC_11_Lindal_Wind.md` v0.4 §1a.

**Accepted. SPEC_11 closes with this step.**

- **The run** converges in 5 passes at the 1e-8 tolerance. Its pressure identity is at most 0.074 K.
- **The round trip** misses the IRIS fit's change by -1.54, +0.77 and -1.13 K at 110, 290 and
  730 mbar. The report splits each miss into four parts that add to it within 0.01 K. The model's
  own parts are 0.06 to 0.23 K each: the wind as the model reads it (C) and the transfer (D). So the
  balance the tool inverts and the balance the transfer applies agree.
- **The build** is as v0.4 states: the case `lindal_iris`, `[composition]` before `[shear]`, and the
  tool's `build` removed. Two suites were edited and not rerun, `step05_3` check 1 and `step11_1`.
  They run at the next full batch, under the regression rule.

**Carried forward.** These are for the author and do not block closure.

1. **The slope and the value (part A).** The wind is built from fit (a)'s slope. At each latitude
   that slope is the local regression slope, not the derivative of the fitted value. Integrated
   along this path, the two differ by 1.26 K at 110 mbar. Most of the difference is near the
   equator and on the approach to the 110 mbar data gap at the target. It is in the wind, so it is
   in the delivered temperature, not only in the benchmark.
   - Recommended: build the wind from the derivative of the fit's value, so the wind and the
     temperatures it stands for agree. Under SPEC_13, both runs would be rebuilt with it.
2. **The band (part B).** It is up to 0.85 K at 730 mbar, and it is the one part the data cannot
   set. The band stays as ruled. The egress comparison will show whether it matters.
3. **The adiabat.** The report's `c_p` of 7/2 R for H2 is the high-temperature limit, with rotation
   fully excited. At 100 to 150 K hydrogen's rotation is only partly excited.
   - Normal hydrogen's `c_p` there is about 3 R, which puts the dry adiabat near 0.85 to 0.9 K/km
     rather than 0.73 to 0.78.
   - Recommended: normal hydrogen, the adiabat Lindal compares against. He finds the radio lapse
     rate between 0.7 and 1.3 bar about 6 percent below the dry adiabat for 94 percent normal
     hydrogen and 6 percent helium.
   - With it, the anchor lies below the adiabat. The delivered profile exceeds it by about 0.1 K/km
     between 600 mbar and 1.3 bar. That is the stability question deferred to the convective
     adjustment.
4. **The target.** By the model's rule, -26.40 deg planetocentric is -31.49 planetographic, the far
   end of the egress swath. SPEC_12 §2 item 1 proposes Lindal's label, 31.2 S.

Order of work: commit and push the following. No regression.
- the changes to `tools/wind/shear.py`, `tools/run/run_inputs.py` and `tools/lindal/lindal_wind.py`;
- the run's two control files;
- the two edited suites;
- `tests/step11_2/`;
- the report and this review;
- the SPEC_11 rows in `STATE.md`.
