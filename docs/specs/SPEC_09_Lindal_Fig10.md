# SPEC 09. IRIS temperatures at three levels, digitized

Version 0.3, 7 October 2026. Author of record: S. Rafkin. **Status: accepted by the author, §2 settled.**
Depends on SPEC_00 v0.23.

## 0. Purpose

Digitize the Voyager IRIS temperatures against latitude at three pressure levels, and characterize
their spread. SPEC_10 builds the shear from them.

Source: Conrath and Pirraglia (1983), Icarus 53, 286, Fig. 1 (`docs/conrath_etal_IRIS_1983.pdf`),
the same data as Lindal et al. (1985) Fig. 10 in a cleaner scan. The top level is labeled 110 mbar,
as in Lindal and in Conrath et al. (1984); Conrath and Pirraglia label it 150 mbar. Their Fig. 3 (thermal wind shear
at 150 mbar) is SPEC_10's sanity check.

Protocol as SPEC_04 §0. No regression: nothing else runs this code.

## 1. Step 1

**Deliverables.**

1. `src/casspian/tools/lindal/iris_temperatures.py`. Render the figure at 600 dpi, calibrate both
   axes by their ticks, take the center of each dot, and assign it to its level. The broken curve
   (northern top-level data folded onto the south) and the labels are excluded. Any decision by hand
   is coded in the script, so the output is reproduced by running it.
2. `occul_data/lindal/iris_temperatures.csv`, committed: `latitude_planetographic_deg`,
   `pressure_mbar` (110, 290, 730), `temperature_K`, one row per dot.

**Acceptance (`tests/step09_1/accept_step09_1.py`, writing under `reports/step09_1/`).**

1. **The spread.** Per level, a fit in latitude by a Gaussian-weighted running mean of width 4 deg
   (the retrievals' resolution, Conrath and Pirraglia), and the RMS of the points about it,
   reported. The fit is drawn over the points. SPEC_10 differentiates the same fit.
2. **The overlay.** The digitized points over the scan, viewed by the author. A reading error is
   estimated only if the overlay suggests it is not small against the RMS.

## 2. Settled (author, 7 October 2026)

1. **The top level is 110 mbar**, the label of Conrath et al. (1984) and Lindal. SPEC_10 runs the
   shear placed at 150 mbar as a one-sided sensitivity case.
2. **Latitude planetographic:** the points match Lindal's replot, which states it.
3. **The fit's width, 4 deg**, the retrievals' resolution. Neither paper states how its fits were
   made.

The coding agent reads this before starting and reports anything that blocks it.

## 3. Revision history

- v0.1, 7 October 2026: first draft.
- v0.2, 7 October 2026: rewritten at the author's direction: Conrath and Pirraglia (1983) Fig. 1 as
  the source; the spread about a fit added; the radio crosses dropped; the tool in
  `tools/lindal`; nothing kept that the code or the paper already says.
- v0.3, 7 October 2026: accepted; §2 settled (110 mbar, planetographic, 4 deg).
