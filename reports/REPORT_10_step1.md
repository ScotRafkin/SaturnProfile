# REPORT 10, Step 1. The fit of the IRIS temperatures and their gradient

Coding agent, 7 October 2026. Specification: `docs/specs/SPEC_10_IRIS_Fit.md` v0.2. Working tree on
`main` at `3ae8dcc`, uncommitted. Acceptance: `tests/step10_1/accept_step10_1.py`, **2 of 2 checks
pass**; output in `reports/step10_1/output.txt`, figures `reports/step10_1/comparison_{110,290,730}mbar.png`,
numbers in `reports/step10_1/comparison.json`. No regression, as §0 directs. Nothing in the
specification blocked the step, so no pre-execution report was filed.

**In short.** `fit` is written, and it is exact on the fields it must reproduce. The three
candidates are drawn for the author:
- (a) and (b) agree wherever the data are dense, and differ where they are sparse;
- (c) has the lowest RMS but invents structure at the edges of the gaps (§3).

The author chose (a), and `degree = 1, min_points = 3` are now `fit`'s defaults (§4).

## 1. The deliverable

`fit(latitude, temperature, degree=1, min_points=3, *, grid, fwhm_deg=4.0)` in
`src/casspian/tools/lindal/iris_temperatures.py`, returning a `Fit`:

| Field | Meaning |
|---|---|
| `value_K` | the fit at each latitude of `grid` |
| `gradient_K_per_deg` | `dT/dphi`, planetographic, per degree |
| `gradient_error_K_per_deg` | the gradient's standard error |
| `fwhm_deg` | the window used at each latitude |

Each field is NaN outside the data's own latitude range.

- **The window.** Gaussian weights of 4 deg full width at half maximum, as SPEC_09 settled. Where
  the effective number of points `(sum w)^2 / sum w^2` is below `min_points`, the width is widened
  to the smallest value that reaches it. That width is solved continuously (Brent's method,
  tolerance 1e-9 deg), so the fit is continuous in latitude where the window widens. `min_points`
  must exceed the number of coefficients, so that a residual is left to estimate the error from.
- **The standard error** is that of a linear smoother whose points have independent errors of one
  variance: `Cov = s^2 A^-1 (X^T W^2 X) A^-1`, `A = X^T W X`. Here `s^2` is the weighted mean square
  residual scaled by `n_eff / (n_eff - (degree + 1))`.
- **Beyond the specification:** `fwhm_deg` is returned too, so the comparison can show where the
  window widened. It can be dropped if the author prefers the three fields §1 names.

## 2. Acceptance

0. **The fit on fields it must reproduce** (beyond the specification; machine precision is asked
   only where the mathematics gives exactness). On the 110 mbar latitudes, with noise-free fields:
   - every candidate returns a linear field's value and slope;
   - (c) also returns a quadratic's.

   The errors are at most 5.7e-14 K and 6.2e-15 K/deg (bounds 1e-9 and 1e-11).

   **The standard error, reported and not bounded.** Over 300 realizations of 0.5 K noise on the
   same latitudes, the mean predicted error against the scatter of the slopes:

   | Candidate | 60 | 20 | 0 | -40 | -80 deg |
   |---|---|---|---|---|---|
   | (a) | 0.98 | 0.99 | 0.99 | 0.97 | 1.02 |
   | (b) | 1.01 | 0.94 | 0.98 | 1.14 | 1.00 |
   | (c) | 1.04 | 1.01 | 1.06 | 1.54 | 1.31 |

   The formula holds for the linear fits. For the quadratic it overstates the error where the data
   are sparse (-40 and -80 deg), because few points against three coefficients make `s^2` noisy.
1. **The comparison**, on a 0.5 deg planetographic grid. Each level's figure has the points with
   the three fits, and below them the three gradients with (a)'s standard error as a band.

   | Level | Candidate | RMS of the points | Window widened | Widest window | Largest abs. gradient | Median standard error |
   |---|---|---|---|---|---|---|
   | 110 | (a) | 0.449 K | 11 % | 8.5 deg | 0.835 K/deg | 0.074 K/deg |
   | 110 | (b) | 0.474 K | 45 % | 14.5 deg | 0.835 K/deg | 0.060 K/deg |
   | 110 | (c) | 0.372 K | 21 % | 10.8 deg | 1.368 K/deg | 0.071 K/deg |
   | 290 | (a) | 0.496 K | 11 % | 8.6 deg | 0.559 K/deg | 0.074 K/deg |
   | 290 | (b) | 0.504 K | 49 % | 13.4 deg | 0.559 K/deg | 0.060 K/deg |
   | 290 | (c) | 0.444 K | 20 % | 10.7 deg | 0.642 K/deg | 0.077 K/deg |
   | 730 | (a) | 0.706 K | 11 % | 8.4 deg | 0.589 K/deg | 0.102 K/deg |
   | 730 | (b) | 0.718 K | 46 % | 13.2 deg | 0.568 K/deg | 0.087 K/deg |
   | 730 | (c) | 0.676 K | 18 % | 10.6 deg | 1.516 K/deg | 0.108 K/deg |

## 3. What the figures show, for the author's choice

- **(c), local quadratic, at least 4,** has the lowest RMS at every level, but it invents structure
  at the edges of the gaps. At -30 deg its 110 mbar temperature spikes to 88.4 K over a stretch
  with no dots, and its gradient reaches -1.37 K/deg. At 730 mbar its gradient reaches +1.2 and
  -1.5 K/deg near -30 and -49 deg. Both are about three times anything the dense data show. Its
  lower RMS is the extra freedom fitting the noise.
- **(a), local linear, at least 3,** follows (b) wherever the data are dense: the two gradients
  overlay to the line width from 75 N to about -25 deg at every level. In the sparse south it keeps
  more structure. Examples are the 110 mbar gradient's swings of 0.3 K/deg between -60 and -85 deg,
  and its 0.43 K/deg step at -48 deg at 290 mbar. Its band shows these are often within one
  standard error of zero. At the northern end (76 N) its gradient rises steeply at 110 and 730 mbar,
  because the end is one-sided and few dots are there.
- **(b), local linear, at least 8,** widens the window at almost half the grid. It bridges the gaps
  and the sparse south smoothly, with the lowest standard error (0.060 to 0.087 K/deg median).
  The price is damped features where the data are sparse, for example the 730 mbar ends and the
  structure between -60 and -85 deg. It also has the highest RMS, by at most 0.03 K over (a).

## 4. Notes

- **Units.** The gradient is per planetographic degree, as §1 places the fit. SPEC_11 converts it
  to the model's planetocentric grid.
- **The SPEC_09 suite after this step.** `fit` was added to the module SPEC_09's tool lives in,
  so `step09_1` was rerun: 3 of 3, the CSV reproduced byte for byte. A first run could not write
  `reports/step09_1/overlay.png` while a viewer held it open (Windows `Errno 22`). It was rerun with
  the figure closed.
- **The defaults, at acceptance (REVIEW_10_step1).** `degree = 1, min_points = 3`, candidate (a).
  `fwhm_deg` is kept. To give the two leading parameters defaults in the order §1 names them,
  `grid` is keyword-only:
  `fit(latitude, temperature, degree=1, min_points=3, *, grid, fwhm_deg=4.0)`.
  - The suite's calls pass `grid=` and it is 2 of 2 after the change.
  - `fit` with no degree or minimum gives (a)'s gradient array-equal at 110 mbar.
- **The ends.** The fit is defined from the northernmost dot (76.8 N at 290 mbar) to the
  southernmost (-89.3 deg at 290 mbar), level by level. Beyond these is SPEC_11's (§2 item 2).
