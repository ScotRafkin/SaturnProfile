# REVIEW 13, Step 1. The gradient from the fitted temperatures, and the Sanchez-Lavega wind

Reviewing agent, 8 October 2026. Report: `reports/REPORT_13_step1.md`. Specification:
`docs/specs/SPEC_13_SL_Wind.md` v0.4.

**Accepted once Step 3 (SPEC_13 §2a) is done and the author has viewed its figures.**

- Nothing in a run reads the wind's uncertainty, so the NaNs of the ring gap are harmless (ruling 2).
- The existing paths are unchanged: the Ingersoll and Pollard wind and `lindal_iris_ii_110`'s wind
  are reproduced.
- The Sanchez-Lavega wind passes through all 260 rows. Its ring gap is bridged without reflection.
- The derivative's spikes near -48 and -13 deg and across the IRIS gap come from the fit's window
  widening abruptly. They are not physical, and Step 3 removes them at their source with a smooth
  window.
- `sanchezlavega2000.toml`, and `smith1982_fig4.toml` with it, are cut to the values that matter,
  as Step 3 item 3 states.
