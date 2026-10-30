# REVIEW 16, Step 1. Pioneer 11 refractivity against the transferred Lindal anchor

Reviewing agent, 30 October 2026. Report: `reports/REPORT_16_step1.md`. Specification:
`docs/specs/SPEC_16_Kliore_Refractivity.md` v0.3.

**Accepted. SPEC_16 closes.**

- **The self-check agrees.** Kliore's temperatures and pressures through `reduction.refractivity`
  return his refractivity column to RMS 0.034 n-units. The implied mean refractivity is 125.87
  n-units at STP against our 125.90: the same species refractivities to 0.02 percent. Our
  refractivity arithmetic and Kliore's are the same.
- **The 17.1 mbar row is a misprint for about 17.5 mbar,** shown three ways: the spacing, the gas
  law at the row's own N and T, and interpolation. **Ruled:** the data file is corrected to 17.5
  with an inline comment, as for SPEC_15's two misprints, in the acceptance commit. No rerun: it
  moves nothing beyond its last digit.
- **The radius reading: planetocentric (ruled).** Kliore's latitudes are planetocentric
  (SPEC_15 §2), and this reading needs an offset of only +9.2 km, against +96 km for the other. At
  100 mbar it puts Pioneer 11 at 60,150 km local radius against 60,157 km delivered. The model's
  reference surface at 9.51 S, built from the Lindal anchor and its geoid, agrees with an
  independent spacecraft's radii to 7 km. The next reduction of Kliore's table uses this reading.
- **By radius, from 10 to 182 mbar,** the delivered refractivity matches Kliore's to 2.3 percent
  RMS after the offset. The two hypotheses do not differ there.
- **Above 10 mbar the departure is structure,** the layering SPEC_15 found in temperature, now seen
  in the measured quantity. At equal radius, the delivered refractivity is low by up to 10 percent
  near 9 mbar and high by 10 to 30 percent above 4 mbar. At equal pressure, the delivered
  refractivity is 16 to 18 percent low at 10 mbar. Both say the same thing: from 3 to 15 mbar the
  observed column is colder and more compact than the delivered one.
- **The pressure comparison restates SPEC_15's temperature comparison exactly.** Lindal's rescaling
  preserves Kliore's N, so at equal pressure the refractivity ratio is the inverse temperature
  ratio. The report shows this, and it is right.

**Order of work:** commit and push, with no regression.
1. The acceptance commit: `tests/step16_1/`; the 17.1 mbar correction in
   `occul_data/kliore1980/pioneer11_exit.csv`; this report and review; SPEC_16 v0.3; `STATE.md`.
2. The hash in SPEC_16's row of `STATE.md`. Push.
