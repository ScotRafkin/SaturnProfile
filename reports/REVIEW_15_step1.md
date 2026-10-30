# REVIEW 15, Step 1. Three more radio profiles

Reviewing agent, 30 October 2026. Report: `reports/REPORT_15_step1.md`. Specification:
`docs/specs/SPEC_15_Radio_Targets.md` v0.3.

**Accepted. SPEC_15 closes.**

- **The tracing.** The overlays sit on the ink through the curls and crossings of both figures,
  including Fig. 5's upper branch and Fig. 6's stippled curve. Guided tracing with every hand
  decision coded is the right method for these scans. Calibration is within 0.0085 dex and 0.09 K.
- **The calibration traces.** RMS 0.95 K for 36.3 N against Table I, and 1.51 K for 31.2 S against
  Fig. 4. The departures over 1 K are explained.
  - **A warm bias in Fig. 5 at depth.** From 600 to 1200 mbar, Lindal's drawn line runs 1.5 to 3 K
    warmer than his own Table I. So the 74.2 S profile probably carries the same bias there, and
    its 730 mbar comparison (+2.5 K) understates the delivered warmth by about 2 K. Fig. 6's
    calibration at depth is unknown. The files are not corrected: the bias is the figure's, it
    belongs in the observation error, and the step records it here.
- **Pioneer 11.** Lindal's stippled curve lies 3.31 K below Kliore's table as published and 0.16 K
  above it rescaled. Lindal rescaled the table, and SPEC_15's rule reproduces his. The target
  conversion gives 9.51 S planetocentric, Kliore's own caption value, which confirms §2 ruling 1
  independently.
- **Pressure identity:** 0.25 K at most, at 11.9 S under hypothesis 2.

**What the comparison shows.**
1. **Below about 30 mbar the transfer works at all four targets.** At 730 mbar it removes the
   anchor's warm bias of 4 to 6 K, and at 110 mbar it is within about 2.5 K everywhere except
   31.2 S, where the radio and IRIS disagree.
2. **Above about 30 mbar, near the equator, the observations carry structure no balanced zonal
   wind of this form can make.** At 3.0 S and 11.9 S, layers 10 to 30 K warm sit between 3 and 15
   mbar. Lindal reads them as vertically propagating waves. The untransferred anchor misses them
   by the same amounts, and Pioneer 11 is two years earlier, across the equinox. In the model's
   terms, this is representativeness error plus an epoch difference, not a wind the transfer
   could carry from 36.3 N.
3. **Hypothesis 2 helps where it was tuned and at 74.2 S aloft, and costs a little near the
   equator.** At 74.2 S it improves 1 and 0.3 mbar, from +4.2 and +10.8 K to -3.7 and -2.3 K. At
   3.0 S and 11.9 S it is 1 to 3 K warmer than hypothesis 1 from 20 to 80 mbar. There it holds
   the full 110 mbar shear, which near the equator is the PCHIP bridge, up to 15 mbar instead of
   relaxing it.
4. **At 74.2 S, both hypotheses are 6 to 8 K cold from 3 to 15 mbar,** where that profile is the
   warmest of Fig. 5. A latitude pattern times f(p) cannot reorder the profiles with height, as
   SPEC_14 expected.

**For the next wind fit.** The equatorial profiles above 30 mbar should enter with a large
representativeness error, or not at all, until waves and time are modeled. Otherwise a fitted
shear would bend itself to reproduce waves. The Fig. 5 depth bias goes into those profiles'
observation error.

**Order of work:** commit and push, with no regression.
1. The acceptance commit: `tools/lindal/radio_profiles.py`; the three data files; the six runs'
   control files; `tests/step15_1/`; this report and review; SPEC_15 v0.3; both Kliore papers;
   `STATE.md`; the status note.
2. The hash in SPEC_15's row of `STATE.md`. Push.
