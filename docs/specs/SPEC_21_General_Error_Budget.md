# SPEC 21. The general error budget, and hypothesis 4

Version 0.7, 11 November 2026. Author of record: S. Rafkin. **Status: closed (REVIEW_21_step1,
REVIEW_21_step1_check5).**
Depends on SPEC_19, SPEC_20 and `docs/CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` v0.6.

## 0. Purpose

Rebuild `casspian-wind-fit`'s error budget on the general rules of the design note v0.6, so that no
rule is written for one dataset, and fit the Voyager epoch again with it: hypothesis 4, beside
hypothesis 3.

What SPEC_19 and SPEC_20 showed, and this answers:
- **The prior was stated in shear,** uniform in latitude. It was too stiff at the equator, where the
  balance needs large shear, and too loose toward the poles (unphysical draws at 74.2 S; 63 K at
  110 mbar near the north pole).
- **IRIS outweighed the radio.** IRIS's errors were random only and too small (group ratio 2.3 at
  110 mbar against 0.2 to 0.3 for the radio). The two southern radio profiles still miss each other
  by 6.7 K there, about 3 sigma_T.
- **Three hand rules for particular data:** 15 K at 3.0 S above 30 mbar, 2 K for Lindal's Fig. 5 at
  74.2 S, and IRIS omitted within 5 deg of the equator.

Not here: the blended product of all anchors (SPEC_22), cloud winds and a free reference wind
(SPEC_23), time. The state stays the wind; the reference wind stays held.

Protocol as SPEC_04 §0. No regression: the model does not change. Lindal's files are untouched.
Absolute profiles first in every figure set, with both hydrogen adiabats; sections with zonal wind
as dashed contours, no hatching.

## 1. The budget, in general form

1. **The shear prior, stated in temperature.** The input is a temperature-anomaly covariance on
   isobars: amplitude `σ_T(p)`, correlation lengths in latitude and `ln p`, and a correlation form
   smooth enough to have a curvature (twice differentiable). The model's own linearized balance
   (`balance_weight`, SPEC_17) carries it to the shear's prior covariance. Near the equator the
   balance requires the temperature gradient to vanish and ties the shear to the curvature, so the
   construction must keep the shear finite there. A suggested construction, for the pre-execution
   report to confirm or replace: a stationary prior on `d ln T / dφ`, conditioned on zero at the
   equator, divided by `w`.
1a. **A deformation-radius floor (v0.4):** the latitude correlation length at each latitude is the larger
   of the input length and the local deformation radius (`c / |f|` away from the equator,
   `√(c a / 2Ω)` at it, in a form joining the two), with `c = N ×` the prior's vertical correlation
   length in height and `N` from the anchor's profile. A balanced structure cannot be narrower.
1b. **No superadiabatic profile (v0.4, form v0.5):** a smooth one-sided penalty on layer slopes between
   the forward run's wind-grid nodes, with rows formed from the wind as the forward run sees it. It
   holds the excess over the dry adiabat (normal hydrogen, `c_p` from composition) to no more than the
   source anchor's own excess in the same layer: the transfer may not make a profile less stable than
   its source allows. Excess already in the source is reported, not penalized. Reported where active;
   inherited by every Monte Carlo member. If it still prevents convergence, it becomes a reported check
   on the fit and a penalty only in the Monte Carlo draws, after a report.
2. **The reference wind** becomes its own part of the state with its own prior, held at zero spread
   in this specification.
3. **The basis follows the lengths:** knots no farther apart than half of each correlation length.
4. **Representativeness for point observations:** `σ_rep(p) = σ_0 (p_0 / p)^(1/2)`, capped, added to
   every observation made at one place (radio profiles), none to zonal means (IRIS). The observation
   type knows which it is; no key per dataset.
5. **Source errors at ingest.** Each radio profile's temperature uncertainty is written into its own
   data at ingest, from its source: Lindal's Fig. 4 sigma_T, his reliability statement for Voyager 1,
   his Fig. 5 depth bias at 74.2 S. The fit reads them; it holds none of them. Kind T's uncertainty
   column, NaN since SPEC_18 for these profiles, is filled.
6. **An uncertain level.** IRIS's top level carries an error equal to its operator's change between
   110 and 150 mbar, the two published labels. The level's range is IRIS's property at ingest.
7. **IRIS's systematic part,** correlated over a long distance in latitude, sized by the consistency
   check: fit without it, size it per level to bring each IRIS group to its expected value, fit again
   once.
8. **Robust weighting:** Huber's rule on normalized misfits, with its threshold an input.
9. **No exclusions.** IRIS is not omitted near the equator; where an observation carries no
   information its operator gives it none.
10. **Every amplitude, length, form and threshold is an input** in the fit's control file, with its
    basis in a comment. Nothing in code.

**One reference, every profile (design note §6).** The four Voyager profiles enter as in SPEC_19: one
reference anchor with its error shared, the others against it, which carries the same information as
every pair once. The reference must not matter (acceptance 5).

## 2. First values (ruled, author, 9 November 2026)

Each a first guess from physics, to be varied by the Monte Carlo and estimated from Cassini data.
Ruled as proposed until results are seen; the prior and the representativeness may be too
optimistic (author). The temperature prior serves the radio profiles until CIRS, whose data will
estimate the atmosphere's variability.

| Input | First value | Basis |
|---|---|---|
| `σ_T(p)` | 2 K at 1 bar, growing as `(1 bar / p)^(1/2)`, capped at 15 K | meridional contrasts on isobars of a few K near the clouds, 5 to 10 K near the tropopause, 10 to 15 K in the stratosphere |
| temperature correlation, latitude | 5 deg, or the local deformation radius where larger (v0.4) | the jets' scale; IRIS resolves 4 deg; balance below the deformation radius is impossible |
| temperature correlation, `ln p` | 1 scale height | as SPEC_19 |
| temperature correlation form | Matérn 5/2 | the least smooth form with a curvature |
| reference wind prior | held, zero spread | until SPEC_23 |
| `σ_rep(p)` | 0.5 K at 1 bar, growing as `(1 bar / p)^(1/2)`, capped at 15 K | about 1.5 K at 110 mbar, 5 K at 10 mbar, 15 K by 1 mbar: the size of the layering in the Voyager profiles aloft |
| radio errors | Lindal Fig. 4 sigma_T; Voyager 1 profiles twice it; 2 K added at 74.2 S, 600 to 1200 mbar | Lindal's statements; the factor of 2 is the reviewing agent's judgment |
| radio correlation in `ln p` | Matérn 3/2, 0.5 scale height, 0.3 K floor | as SPEC_19 |
| IRIS systematic, latitude correlation | 20 deg | longer than the jets, so it spoils long integrals, not local gradients |
| IRIS top level range | 110 to 150 mbar | Lindal and Conrath 1984; Conrath and Pirraglia 1983 |
| Huber threshold | 2 normalized units | common choice; reported, not tuned |

## 3. Step 0: pre-execution report

`reports/REPORT_21_preexecution.md`, before building:
1. **The prior's construction** (§1 item 1): how the temperature covariance is carried to the shear,
   how the equatorial limit is taken, and a check that its draws are finite and smooth there.
2. **The cost** of the finer basis (2.5 deg knots in latitude) and of the second fit for IRIS's
   systematic part.
3. **The ingest change** for the profiles' uncertainties, and that Lindal's files stay untouched
   (his anchor's uncertainty is written where the fit reads it, not into his registered files, if
   that is the only way).
4. **Anything that blocks.**

## 4. Step 1

**Deliverables.**
1. `casspian-wind-fit` with the budget of §1, its control file with §2's values as ruled.
2. **Hypothesis 4:** the wind `voyager_fit_2`, its runs from `lindal_iris_v_sl_ii_110` with only the
   wind changed to the three other Voyager profiles, and its sections from all four anchors.

**Acceptance (`tests/step21_1/accept_step21_1.py`, under `reports/step21_1/`).**
1. **Synthetic, before real data.**
   - **(a) Resolution, reported:** hypothesis 1 as the truth.
   - **(b) A balanced equatorial anomaly, reported:** a 20 K warm layer, half-width 10 deg in latitude,
     about one scale height thick near 15 mbar, with its shear from the prior's own construction.
     What the fit recovers from the real observation points.
   - **(c) Statistics, linear:** 20 truths drawn from the new prior. Pass: at least 90 percent of the
     shear within 2 posterior sd, and each group's mean ratio between 0.8 and 1.2.
   - **(d) The prior's draws are physical:** the largest temperature anomaly of the 20 draws at each
     Voyager latitude and at 87.5 N and S, against `σ_T`; and at the equator, the largest shear and
     the wind change from 10 to 1000 mbar (REVIEW_21_step0); the equatorial deformation radius used,
     and where the stability constraint was active (v0.4).
2. **Hypothesis 4, absolute first:**
   - the anchor matrix as SPEC_20 drew it, under hypothesis 4, beside hypothesis 3, with the
     table of delivered minus observed at 730, 290, 110, 30, 10 and 1 mbar;
   - the four sections with zonal wind;
   - the shear and the wind with their posterior sd, beside hypotheses 1 and 3;
   - the IRIS residuals at the three levels.
3. **The checks:** the consistency check per group before and after IRIS's systematic part, with
   the part's size per level; leave-one-out, each radio profile and each IRIS level; which
   observations the robust weighting downweighted, and by how much.
4. **The solve:** passes, the start from the prior mean, the pressure identity of each run.
5. **The reference does not matter:** the fit repeated with the Voyager 2 egress as reference. The
   delivered temperatures at the four profiles agree with the Lindal-referenced fit within 0.5 K.

## 5. Revision history

- v0.1, 9 November 2026: first draft.
- v0.2, 9 November 2026: accepted; §2 ruled as proposed.
- v0.3, 9 November 2026: Step 0 accepted; the wind added to acceptance 1(d).
- v0.7, 11 November 2026: closed. Observations compared at the forward model's resolution, every
  profile averaged into the wind grid's layers, the floor uncorrelated once per layer (design note
  v0.9); hypothesis 4 refit (at most 0.75 K change); check 5 accepted on the materiality rule (0.2 to
  0.9 K RMS, within 0.5 K at the standard levels but 0.68 K at 3.0 S, 10 mbar).
- v0.6, 10 November 2026: Step 1 accepted on check 5; wind files extend a decade above the
  forward grid's top; 1c judged on 500-draw means.
- v0.5, 10 November 2026: the stability constraint at the forward model's resolution, relative to
  the source anchor, smooth.
- v0.4, 9 November 2026: after the first Step 1 run, a deformation-radius floor on the latitude length
  and a no-superadiabatic constraint; no shear clip.
