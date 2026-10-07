# REPORT 09, Step 1. IRIS temperatures at three levels, digitized

Coding agent, 7 October 2026. Specification: `docs/specs/SPEC_09_Lindal_Fig10.md` v0.3. Reading:
`reports/REPORT_09_preexecution.md`, nothing blocking. Working tree on `main` at `2fbcf1d`,
uncommitted. Acceptance: `tests/step09_1/accept_step09_1.py`, **3 of 3 checks pass**; output in
`reports/step09_1/output.txt`. No regression, as §0 directs.

**In short.** The tool finds 697 dots: 243 at 730 mbar, 220 at 290 mbar, and 234 at the top level,
written as 110 mbar. The two axes calibrate to their ticks within 0.067 K and 0.077 deg. About the
Gaussian running mean, the dots scatter by an RMS of 0.49 K at 110 mbar, 0.51 K at 290 mbar and
0.71 K at 730 mbar. The overlay shows every dot found and assigned to its level. One reading is
stated for the author: the 4 deg width is taken as the full width at half maximum (§3).

## 1. Deliverables

- **`src/casspian/tools/lindal/iris_temperatures.py`.** It runs as `python -m
  casspian.tools.lindal.iris_temperatures`, with defaults for the PDF and the output path.
  - **Rendering.** Ghostscript renders page 2 at 600 dpi.
  - **Frame.** Found from the fully inked rows and columns.
  - **Calibration.** A least-squares line through each axis's ticks: 9 temperature ticks, every 5 K
    from 120 to 80 K, at 41.43 px per K; 17 latitude ticks, every 10 deg from 80 to -80, at 16.70 px
    per deg.
  - **Dots.**
    - The labels are blanked, and an opening of radius 4 px removes the broken curve, the ticks and
      the frame.
    - Merged dots are split by area: one dot is 88 px, and k-means splits 254 dots out of
      components holding 2 to 9.
    - A second pass takes back dots printed narrower than the opening's disk: three, measured at 74
      to 80 px and 8 px wide.
    - Two remnants of the broken curve are dropped by the curve region (south of -18 deg, colder
      than 85 K).
  - **Levels.** Each dot is assigned by temperature: above 105 K is 730 mbar; below a
    piecewise-linear line in latitude is the top level; between is 290 mbar.

  Every decision by hand is a named constant with its reason in the module docstring: the figure
  box, the label boxes, the opening radius, the narrow-dot test, the curve region and the two
  boundaries. `running_mean`, the fit of acceptance 1, is in the same module. It characterizes
  the spread and is not the fit the shear is built from, which SPEC_10 defines
  (REVIEW_09_step1).
- **`occul_data/lindal/iris_temperatures.csv`**, with columns `latitude_planetographic_deg`,
  `pressure_mbar` and `temperature_K`, one row per dot. Rows run by level from 110 mbar down and
  then from north to south. Values are written to 0.01 deg and 0.01 K. One pixel is 0.060 deg and
  0.024 K.

## 2. Acceptance

0. **Reproduced.** Running the tool writes the committed CSV byte for byte (697 rows). This check
   is beyond the specification: it is deliverable 1's "the output is reproduced by running it".
1. **The spread.** The fit is a Gaussian-weighted running mean with a full width at half maximum of
   4 deg (sigma 1.70 deg). The RMS is taken at the dots; `reports/step09_1/spread.png` draws the
   fit over them.

   | Level | Dots | Latitudes | RMS about the fit | Largest residual | RMS with 4 deg as sigma |
   |---|---|---|---|---|---|
   | 110 mbar | 234 | 76.67 to -88.87 | 0.491 K | 1.45 K | 0.835 K |
   | 290 mbar | 220 | 76.77 to -89.34 | 0.506 K | 1.79 K | 0.614 K |
   | 730 mbar | 243 | 76.62 to -88.56 | 0.713 K | 2.77 K | 0.779 K |

   The fit is left undefined where no dot lies within 4 deg, so it does not bridge the gaps in the
   data: about -31.5 to -34.5 deg at all three levels, and -45.75 to -46 deg at 110 mbar. It does
   reach up to 4 deg past the last dot at each end, where the mean is one-sided.
2. **The overlay.** `reports/step09_1/overlay.png` shows the scan with every center as a filled disk
   (red 730, green 290, blue 110 mbar) and the blanked label boxes in orange. Every center lies
   inside the frame. For the author to view.

## 3. For the author

- **The width.** "Width 4 deg" is read as the full width at half maximum, the usual statement of
  a resolution, and it is a parameter of `running_mean` (`fwhm_deg`). With 4 deg read as sigma
  instead, the RMS rises to 0.84, 0.61 and 0.78 K (table above), because the wider fit smooths
  through the 110 mbar structure. Neither paper states how its fits were made (§2 ruling 3).
- **The reading error.** Calibration costs at most 0.067 K and 0.077 deg. A dot's center is found
  to well under a pixel where it stands alone, which is 0.024 K. The error is larger only where a
  dot is split out of a merged component (254 of 697 dots), where a center can be off by up to
  about a dot radius, 0.13 K. These are small against the RMS of 0.5 to 0.7 K. The overlay is the
  check (§1 acceptance 2), and no further estimate is made unless the author sees one is needed.
- **Where the count is uncertain.** About three merged components, where 150 mbar dots touch the
  broken curve between 0 and -12 deg, plausibly hold two or three dots each. A remnant of the
  curve can add one to a component's area count there. This changes the number of dots by one or
  two at most, all within 1 K of the fit.
- **The PDF** is untracked in `docs/`. The CSV is the committed product.
- **No console entry** was added to `pyproject.toml`. The tool runs with `python -m`. One can be
  added if the author wants it beside the others.
