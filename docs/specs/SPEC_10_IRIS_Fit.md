# SPEC 10. The fit of the IRIS temperatures and their gradient

Version 0.3, 7 October 2026. Author of record: S. Rafkin. **Status: accepted by the author; Step 1 accepted, (a) chosen.**
Depends on SPEC_09.

## 0. Purpose

Fit the IRIS temperatures of `occul_data/lindal/iris_temperatures.csv` in latitude at each level and
give their gradient, which SPEC_11 turns into shear. The fit must keep the structure the data hold
without inventing structure where they are sparse, and must bridge the gaps. The author chooses the
setting from a side-by-side comparison.

Protocol as SPEC_04 §0. No regression: nothing else runs this code yet.

## 1. Step 1

**The method.** The fit is made in planetographic latitude, as the data are; SPEC_11 places it on
the model's planetocentric grid by the wind tool's rule. At each latitude a weighted polynomial (local linear or local quadratic) is fitted to the
points with Gaussian weights of 4 deg full width at half maximum. Where the data are too sparse for
that window, it is widened until the effective number of points, `(sum w)^2 / sum w^2`, reaches a
minimum. The value is the fit at that latitude, the gradient its slope there, and the slope's
standard error follows from the weighted residuals. The fit is defined from the northernmost to the
southernmost point; beyond them is SPEC_11's.

**Deliverables.**

1. `fit(latitude, temperature, degree, min_points)` in `src/casspian/tools/lindal/iris_temperatures.py`,
   returning the value, the gradient and its standard error on a given latitude grid.
2. **The comparison**, for the author, at each level on a 0.5 deg planetographic grid:
   - (a) local linear, the window widened only to keep 3 effective points;
   - (b) local linear, at least 8 effective points (the reviewing agent's trial);
   - (c) local quadratic, at least 4 effective points.
   
   Per level, the points with the three fits, and below it the three gradients with the standard
   error of (a) as a band; the RMS of the points about each fit.

**Acceptance (`tests/step10_1/accept_step10_1.py`, writing under `reports/step10_1/`).** The author
views the comparison and chooses the degree and minimum, which become the defaults of `fit`.

## 2. For the author

1. The three candidates, or others to add.
2. The fit defined only between the end points; the polar ends left to SPEC_11.

## 3. Revision history

- v0.1, 7 October 2026: first draft.
- v0.2, 7 October 2026: accepted; each level's gradient panel drawn below its temperature panel.
- v0.3, 7 October 2026: the author chose (a), local linear with at least 3 effective points, as the
  default (`reports/REVIEW_10_step1.md`).
