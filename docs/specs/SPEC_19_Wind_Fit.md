# SPEC 19. The wind fit, and hypothesis 3 at the Voyager epoch

Version 0.5, 9 November 2026. Author of record: S. Rafkin. **Status: closed. Step 1 accepted after
its rerun (REVIEW_19_step1, §5 here).**
Depends on SPEC_17 and `docs/CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` (v0.5).

## 0. Purpose

Build the wind estimation tool of the design note, and its first product: hypothesis 3, the
best-guess shear at the Voyager epoch from the IRIS temperatures and the Voyager radio profiles,
weighted by their errors.

**The tool is built for every dataset to come.** CIRS nadir and limb, Cassini radio occultations,
cloud-tracked winds and an imposed seasonal cycle will each be added as a term, not as a new tool.
So the tool has three parts, and only the first two grow:
- **Observation types,** each a module with three things: a reader for its standard data file, an
  operator that gives the model's equivalent of the data from a wind (with its linearization), and
  an error covariance. This specification implements two: radio targets and temperature gradients.
  CIRS is the second type again; Cassini radio occultations are the first, or anchor pairs; cloud
  winds will be a third.
- **The prior:** a mean and a covariance. Its mean is where an imposed seasonal cycle will enter.
- **The solver,** written once: weighted least squares, linearized, with the posterior covariance,
  the consistency check and leave-one-out. It does not change as types are added.

The state carries an epoch from the start. Here there is one epoch; the time correlation of the
design note §7 adds epochs without changing the solver.

Lean, as always: the three parts are what is built, and no interface beyond what these two types
use. A key is added only if code reads its value; sources and choices go in comments.

Protocol as SPEC_04 §0. No regression: a new tool and new runs, and no change to the model.

## 1. What is estimated, and how

1. **The state:** the shear `s(φ, p)` at the Voyager epoch, on a smooth basis: cubic B-splines in
   planetocentric latitude, knots every 5 deg from pole to pole, and in `ln p`, knots every half
   scale height from 1.3 bar to 0.2 mbar. Above and below, the shear is held at its edge value.
   The polar basis functions are dropped, so the wind is zero at the poles (REPORT_19_step1).
2. **The reference wind is held:** the Sanchez-Lavega cloud wind at 398 mbar, as the base run. It
   enters the state in the next specification, with the cloud winds and the 100 mbar geometry.
3. **The wind:** `u_total = u_ref + ∫ s d ln p` from 398 mbar, as the tool integrates now.
4. **The prior (§2 item 1):** mean zero shear; covariance with an amplitude and correlation lengths
   in latitude and `ln p`. It carries no observation. It supplies smoothness and settles what the
   data leave open: the vertical shape between and beyond IRIS's levels, and the distribution of a
   path integral along the path.
5. **The solve:** minimize `J` of the design note §2, linearized about the current wind. Then the
   wind is rebuilt, the targets rerun and the operators relinearized, until the change in the
   delivered temperature at every target is under 0.1 K.
6. **The composition (author, 5 November 2026):** one composition for every input. Here it is
   Lindal's reduction composition, 94 percent hydrogen and 6 percent helium with his ammonia, with
   true molecular masses. All three Voyager profiles are already on it (SPEC_15). IRIS's own
   assumption cannot be rescaled this way and is part of its error.

## 2. The error budget (ruled, author, 5 November 2026; revised 6 November 2026)

A first budget, revised when warranted: by the consistency check, by leave-one-out, or by new data.
Each number is given with its basis. Code reads them, so they are values in the tool's control file
or the data property files, each with its basis in a comment.

1. **The prior (revised 6 November).**
   - Amplitude: 40 m/s per unit `ln p`. The prior's mean is zero shear (design note v0.5), so the
     prior must hold the whole shear, not a correction to an IRIS-based mean. The equatorial shear
     of hypotheses 1 and 2 reaches about 65 at 110 mbar; at 40 that is within 2 sd. The first value,
     20, was sized for corrections and made the jet a 3.5 sd event.
   - Latitude correlation: 5 deg. The jets are narrower than 10 deg, which flattened the equatorial
     shear at 110 mbar from +65 to +11 (REPORT_19_step1 §4). 5 deg is the knot spacing and close to
     IRIS's 4 deg resolution.
   - `ln p` correlation: one scale height, unchanged.
2. **The radio targets** (Voyager 2 egress at 31.2 S; Voyager 1 ingress at 74.2 S and egress at
   3.0 S), each against the Lindal anchor.
   - **Size:** Lindal's sigma_T curve of the Voyager 2 egress (Fig. 4), for all three. The Voyager 1
     bars drawn at the ends of Figs. 5 and 6 (about 9.9 K at the top, 4 to 6 K at the bottom) agree
     with that curve, as the same instrument and processing would.
   - **Correlation in `ln p` (form revised 6 November):** length 0.5 scale height, form Matérn 3/2,
     `(1 + √3 d/ℓ) exp(-√3 d/ℓ)` with `d = |Δ ln p|`. The exponential expects rough noise that the
     printed profiles do not have, so the group ratio could not test sigma_T (REPORT_19_step1 §5).
     The form is code, with this basis in a comment.
   - **A floor of 0.3 K, uncorrelated (added 9 November).** The linearized operator departs from
     the model by 0.1 to 0.3 K in rough detail, and the Matérn matrix (condition number about
     6e5) gives its smallest modes an sd near 0.07 K, so that departure counted as much as the whole
     radio noise (REPORT_19_step1_v04 §5). The floor is added to the diagonal of each target's
     covariance. Basis in a comment.
   - **The anchor's own error** is common to all three targets, so it enters as one shared term
     across them, the same sigma_T curve. Composition error is common to anchor and targets and
     cancels.
   - **Representativeness at 3.0 S above 30 mbar:** 15 K, the size of the wave layering SPEC_15
     found there.
   - **The Fig. 5 depth bias at 74.2 S, 600 to 1200 mbar:** 2 K added, from REVIEW_15_step1.
3. **IRIS.**
   - **Random only (design note §11 item 3).** Its size is the standard error of the fit's
     gradient at each latitude (SPEC_10), with the smooth window (SPEC_13).
   - **Sampled every 4 deg,** the window's width, so that neighboring gradients are nearly
     independent and the window's smoothing is not counted as information.
   - **Omitted** within 5 deg of the equator, where the balance sets no shear, and where no IRIS
     point lies within 2 deg, such as the gap from -27.7 to -38.5 deg, where the fit only bridges.
   - **The operator is a layer mean (added 6 November).** Conrath and Pirraglia (1983), Fig. 1:
     the values "represent mean values for atmospheric layers approximately one-half scale height
     thick", the retrievals' vertical resolution. The operator averages the point operator uniformly
     in `ln p` over half a scale height centered on the level. Basis in a comment.
   - **A systematic term only if the checks find one** (§3 acceptance 3). None is assumed. The
     candidate form is an error correlated over a long distance in latitude, since that is the error
     that would spoil the radio's long path integrals; Conrath and Pirraglia also note that below
     300 to 400 mbar the retrievals may carry cloud and haze opacity they did not model.

## 3. Step 1

**Deliverables.**
1. **The tool,** in `src/casspian/tools/wind/estimation/` with the command `casspian-wind-fit`, as §0
   and §1 describe (SPEC_01's point fit stays in `tools/wind/fit.py`). It reads one control file
   naming the epoch, the base wind, the anchor, the observations and the error budget, and writes a
   kind W file for the epoch, with the posterior standard deviation of the wind as
   `u_total_uncertainty_ms`.
2. **The two observation types:**
   - **Radio target:** misfit `ln(T_obs / T_delivered)` on the anchor's levels; operator, the path
     integral of `w Δs` with `w` from `balance_weight` (SPEC_17).
   - **Temperature gradient:** `d ln T / dφ` from a fit of IRIS-like data at a level; operator, the
     model's balance at that latitude, linearized in the shear with the wind held, as
     `balance_weight` is, averaged over the layer of §2 item 3.
3. **Hypothesis 3:** the wind `voyager_fit_1` from the base run's wind, and three runs from
   `lindal_iris_v_sl_ii_110` with only the wind changed: to 31.2 S, 74.2 S and 3.0 S. The first fit,
   under v0.2's budget, is recorded in REPORT_19_step1 only.
4. **Files:** `fits/` is a top-level directory for the fits. Its control file and wind file are
   committed; work directories are ignored; the runs' inputs rebuild from the wind file.

**Acceptance (`tests/step19_1/accept_step19_1.py`, under `reports/step19_1/`).**
1. **The synthetic tests, before real data.**
   - **(a) Resolution, reported without a bound:** observations made from hypothesis 1's shear at
     the real observation points, with noise drawn from §2's errors. Recovered against truth, the
     error over its posterior sd, and each group's ratio. Where the basis rather than the prior
     limits the recovery (the basis limit of REPORT_19 §5), it is said; 2.5 deg knots are the next
     step, by ruling.
   - **(b) Statistics, linear (ruled 6 November):** 20 truths drawn from the prior, with data made by
     the operator at the final linearization, each with its own noise. Pass: pooled over the
     ensemble, at least 90 percent of the shear within 2 posterior sd, and each group's mean ratio
     between 0.8 and 1.2. It tests the solver and the budget's statistics. The nonlinear iteration
     is tested by (a), whose truth is physical, and by the start from the prior mean (item 5).
     Truths drawn from the amplitude 40 prior reach unphysical temperatures at 74.2 S, so a
     nonlinear ensemble of them is not run.
2. **Hypothesis 3, for the author:**
   - the shear against latitude and pressure, beside hypotheses 1 and 2, with its posterior
     standard deviation;
   - **the shear and its posterior sd against pressure from 1.3 bar to 0.2 mbar at 50 N, 20 N, 5 N,
     15 S and 50 S planetocentric,** with hypothesis 1, showing how the fit and its uncertainty
     change through IRIS's top level;
   - delivered minus observed at the three targets, all three hypotheses with the first fit beside
     them, with the observed errors, and at 730, 290, 110, 30, 10, 4.5 and 1 mbar in a table;
   - the IRIS residuals at the three levels;
   - **at 110 mbar between 2.4 and 26 S,** where IRIS and the radio disagree: where the fit lands
     between them, in kelvins at 31.2 S, with the path decomposition of REPORT_19 §4.
3. **The consistency check,** per group: each radio target, and IRIS at each level. Each group's
   weighted misfit against its expected value. **If an IRIS group is above 1.5,** or leave-one-out
   misses a radio target by more than 2 stated sd at an IRIS level, the systematic term of §2 item 3
   is sized, by one linear solve, at a latitude correlation of 20 deg and, as a sensitivity, 40 deg:
   its size per level and what it does at 110 mbar at 31.2 S. It is reported and not used until the
   author rules.
4. **Leave-one-out:** each radio target withheld in turn and predicted by the fit to the rest;
   IRIS at each level withheld in turn. The prediction errors against the stated errors.
5. **The solve:** the change at each relinearization, and the pressure identity of each run. The fit
   started once from the prior mean as well as from hypothesis 1's projection; the delivered
   temperatures agree within 0.1 K.
6. **The wind aloft:** the largest change from the base wind, and where, against its posterior sd.

## 4. Amendment to the design note (v0.5)

§5's prior mean was the current hypothesis, the cloud wind with the IRIS shear. With IRIS also an
observation (§3 there), IRIS would count twice. From v0.5 the prior carries no observation: its
mean is zero shear and the held reference wind, or later an imposed seasonal cycle, and every
dataset, IRIS included, is an observation term. Its amplitude must then cover the whole shear
(§2 item 1, revised 6 November).

## 5. Closing (author, 9 November 2026)

1. **Step 1 accepted; SPEC_19 closes.** The tool is the deliverable. Check 1(b) passes, the fit is
   the same from two starts, and the first fit's failure is gone.
2. **Hypothesis 3 is recorded as the Voyager-epoch fit,** `voyager_fit_1`, with this caveat in its
   control file comments: within about 10 deg of the equator IRIS has no rows and the balance sets
   little shear, so the shear there comes from the prior and its stated uncertainty is too small
   (the noise-free resolution test returns a true +54 as +85 at 110 mbar). For delivered
   temperature this is worth about 1 K along any path; for the equatorial wind aloft it is not a
   result. Cloud-tracked winds at more than one level, the next observation type, constrain it.
3. **No further tuning of the Voyager fit.** Materiality (author, 9 November): a change worth less
   than about 1 K in delivered temperature, or about 5 percent in wind at probe levels, gets a
   comment, not a specification.
4. **The iteration** (REPORT_19_step1_v04 §7 item 1) is accepted: Anderson(1) steps while the
   linear change is at or above 0.1 K, ending only on a plain step under 0.1 K; `max_passes` 20.
5. **The radio floor** of §2 item 2 is added at the acceptance commit. The products are rebuilt
   clean with it. If the delivered 110 mbar value at 31.2 S moves by more than 1 K from +3.34 K,
   the step returns for review; otherwise it is recorded in STATE.md.
6. **The absolute profiles (author, 9 November).** Added at the acceptance commit, plotting only:
   for each target, the observed temperature against pressure with its sigma_T band, and the
   delivered profiles of hypotheses 1, 2 and 3 and the untransferred anchor, 1.3 bar to 0.2 mbar,
   log pressure. Every later comparison shows the absolute profiles first, then the differences.
7. **Not adopted, deferred:**
   - the IRIS systematic term: at the computed sizes it discards IRIS, and the excess sits on a few
     rows at the jet's edge, under 1 K on any path;
   - an amplitude that varies with latitude: it affects the prior's own draws, not the fit; taken
     up in the Monte Carlo if its draws need it;
   - a study of the jet's bias, and 2.5 deg knots;
   - the 15 K wave term at 3.0 S, whose ratio (0.42) suggests it is large.

## 6. Revision history

- v0.1, 5 November 2026: first draft.
- v0.2, 5 November 2026: accepted; §2 ruled as proposed, a first budget open to revision.
- v0.3, 6 November 2026: after REPORT_19_step1. Prior amplitude 40 and latitude correlation 5 deg;
  the radio correlation Matérn 3/2; IRIS as a half scale height layer mean; check 1 split into a
  resolution test and a 20-member statistical test; the systematic term sized only if the checks
  ask; the profiles through 110 mbar; the start from the prior mean; the package path and `fits/`.
- v0.4, 6 November 2026: check 1(b) made linear; the nonlinear ensemble is not run.
- v0.5, 9 November 2026: closed (§5); the radio floor of 0.3 K; hypothesis 3 recorded with its
  equatorial caveat; the absolute profile figure.
