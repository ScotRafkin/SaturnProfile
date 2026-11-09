# CASSPIAN design note: the best-guess wind

Reviewing agent and author, 5 November 2026. **Version 0.5. Status: accepted by the author, 5 November 2026; v0.5 amends §5 (SPEC_19 §4).**
Companion to `CASSPIAN_Seasonal_Design_Note.md`, whose §5 it builds on and whose §6 it amends
(§6 here).

## 0. The goal

The model delivers one profile for one input. That profile is only as good as the wind it is
transferred with. The goal is a best-guess wind at each epoch: the wind that, used by the model,
reproduces every observation of that epoch to within its error, and changes with time no faster
than the atmosphere can. If it does, the model can be trusted where and when there are no
observations, with an uncertainty that grows with distance from them in space and time. Varying
that wind about its uncertainty is then the Monte Carlo.

The best-guess wind is an input, built outside the model by a tool and written as a kind W file
per epoch, like any other hypothesis (SPEC_00 principle 5). The tool uses the model: its misfits
come from running the transfer, and its operator is the model's own balance (SPEC_17).

## 1. What is estimated

The whole wind, at each epoch `e` that has observations:
- the **reference wind** `u_ref(φ; e)`, the wind at the reference level;
- the **reference level** `p_ref`, itself uncertain, since cloud-tracked features sit somewhere
  between the upper haze and the ammonia cloud;
- the **shear** `s(φ, p; e)`, the change of the wind from the reference level.

`u_total` follows from the three. The estimate is a correction `x(e)` to a prior wind `x_b(e)` (the
current hypothesis), on a smooth basis in latitude and `ln p`. **All epochs are estimated
together** (§5, §7), so each epoch's wind is informed by the others in proportion to how close
they are in time.

## 2. The estimate

One weighted least-squares problem over all epochs:

    J(x) = (x - x_b)^T B^-1 (x - x_b) + Σ_k (y_k - H_k(x))^T R_k^-1 (y_k - H_k(x))

- `B` is the prior's error covariance: how far, and over what distances in latitude, `ln p` and
  time, the true wind may depart from the prior (§5). Its correlation in time is what couples
  the epochs.
- For each observation `k`: the data `y_k`, the operator `H_k` that gives the model's equivalent
  of the data from the wind at the observation's epoch, and the error covariance `R_k` (§4).

**Solved** by linearizing `H_k` about the current wind and iterating if the correction is large.
SPEC_17 shows the transfer is linear to about 1 percent for small corrections and 3 percent for
a 4 K correction, so one or two passes suffice.

**Outputs**, all from the same solve:
- the best-guess wind at each epoch, written as kind W;
- its **posterior covariance** `(B^-1 + H^T R^-1 H)^-1`, which the Monte Carlo samples;
- the **consistency check**: each term of `J` against its expected size. Misfits larger than the
  stated errors allow mean an error is understated, and the check says which (§8);
- **leave-one-out**: each observation withheld in turn and predicted from the rest. Its
  prediction error is the honest measure of trust where there are no observations.

## 3. The observations

Every observation type is three things: a data file in its standard form, an operator, and an
error model. Adding a type adds a term to `J`. Nothing else changes.

| Type | Constrains | Operator `H` | Status |
|---|---|---|---|
| Radio occultation pairs at the same or near epoch | the shear's latitude integral between the two latitudes, at every level | each profile transferred to the other's latitude, linearized by `balance_weight` (SPEC_17) | Voyager: Lindal ingress and egress, Voyager 1 ingress and egress |
| IRIS temperatures | the shear itself along three isobars, at every latitude | the balance, locally; vertically a layer average over the retrieval's weighting function | in hand |
| CIRS cross-sections | the same as IRIS, at many levels and epochs | as IRIS | not yet |
| Cloud-tracked winds | the reference wind at its level | the wind at `p_ref` | Voyager curves in hand as priors; Cassini and HST not yet |
| The 100 mbar isobar's radius | the reference wind and the shear below 100 mbar together | the isobar's radius from the wind geoid (SPEC_17 Part A) | Lindal Fig. 9 in hand |
| Probe in situ wind, Doppler winds | the wind at their levels | the wind there | future |

## 4. Errors

The weights are the errors, so the errors decide the answer.

- **What matters is the error of what the operator compares.** A radio pair constrains a
  difference between two latitudes; an error common to both profiles cancels (the composition
  assumed in the inversion, for example). Only the independent part counts.
- **Correlation length decides which scale an error spoils.** An error independent from one
  latitude to the next makes local gradients noisy but averages out over a long path. An error
  that varies slowly with latitude (a retrieval assumption that changes across hemispheres)
  leaves local gradients nearly right but biases the long integral. So:
  - IRIS is trusted for the shear's latitude structure; its slowly varying systematic part is
    given a long correlation length;
  - the radio, which measures the long integrals directly, then dominates the integrals while
    IRIS shapes everything within them.

  The fit does this by itself once the correlations are stated; no hand-built rule is needed.
- **Representativeness.** What a zonal-mean balanced wind cannot carry (waves, local anomalies)
  enters as error, not as signal. The equatorial profiles above 30 mbar are the present case.
- **Model error.** The model's own approximations (the equatorial bridge, 0.6 to 0.9 K) take their
  share as error before the wind takes the rest.
- **Where the numbers live.** An error that code reads is a value in a data file or a data
  property file, with its basis in a comment. Nothing else is a key.

## 5. The prior

- **The prior carries no observation (v0.5).** Its mean is zero shear and the held reference wind,
  or later an imposed seasonal cycle. Every dataset, IRIS included, is an observation term (§3), so
  nothing counts twice. The prior supplies smoothness and settles what the data leave open.
- **Its error covariance `B`** has an amplitude and correlation lengths in latitude, in `ln p` and
  in time.
  - The vertical length makes IRIS's influence fade away from its levels over about the
    retrieval's own thickness; above 110 mbar and below 730 mbar the radio and the prior carry
    the shear alone.
  - The time correlation (§7) is what lets anchors close in epoch inform each other's wind.
  - A smooth prior also settles which of the many path-equivalent corrections the fit takes, and
    keeps it from fitting noise with unphysical winds.
- **The equatorial band** keeps its own treatment: there the balance sets no shear from a
  temperature gradient, and the prior and the wind observations carry it.

## 6. Anchors: every radio profile is one

**No radio profile is special (author, 5 November 2026).** Any radio occultation can serve as an
anchor; Lindal's is first only because it was reduced first. Once a profile is reduced to kind N,
it is an anchor like any other.

- **Anchors at the same epoch constrain the wind through each other.** Each is transferred to the
  other's latitude, and the misfit there is a measurement of the shear's integral between them.
  To first order the transfer from i to j and from j to i carry the same information, so each pair
  enters `J` once, not twice. With several anchors, the pairs that share a latitude segment locate
  the correction in latitude, as SPEC_17's segments did from one anchor.
- **The wind tool is where all wind information is combined.** The model's combination of Eq. A35
  (seasonal note §6), which corrects the kernel from the anchors of one run, does the same thing
  for one run's anchors only. With the tool in place, a run takes the best-guess wind of its epoch
  and its posterior uncertainty, and the model's combination blends the anchors' constants
  (Eq. A30) without correcting the wind again. Otherwise the same anchors would count twice.
- **A radio profile not yet reduced to kind N** enters the tool as a target of the reduced anchors
  (the present Voyager 2 egress and Voyager 1 profiles). It enters once, and leave-one-out
  withholds it from wherever it is.

## 7. Time

**Each epoch has its own best-guess wind (author, 5 November 2026),** built from that epoch's
anchors and observations, and informed by other epochs in proportion to their closeness.

- **The epochs are estimated together.** `B` correlates the wind at two epochs by a function of
  their separation in season, falling to the interannual floor at large separations. So an anchor
  near in time pulls on the other epoch's wind almost as hard as on its own, and one far in time
  hardly at all. There is no threshold in either direction.
- **The correlation time is a continuum with height:** long in the deep troposphere, where little
  changes, short aloft, where the season acts fast. Recommended: tie it to the radiative time
  constant `τ(φ, p)` of the seasonal propagator (seasonal note §5.2), so the wind is allowed to
  change in time as fast as the temperature that balances it, and no faster.
- **The wind in time is never derived from a propagated thermal field** (seasonal note §5.1, v0.4).
  It comes from observations at their own epochs, coupled through `B`. Temperature in time remains
  the anchors' business: each anchor is propagated by the seasonal propagator to a run's season
  before it is transferred.
- **Interannual variability**, which no seasonal rule predicts, is part of `B`'s amplitude: it is
  the floor below which no observation from another year can pin the wind.

## 8. The checks

- **Within an epoch:** the model, run with the epoch's best-guess wind and its anchors, reproduces
  every observation of that epoch within its stated error.
- **Across epochs (author, 5 November 2026):** two epochs close in time must not have wildly
  different winds. In the joint estimate this is the time-correlation term of `J`: if the data
  pull near epochs apart further than `B` allows, that term is large and the consistency check
  reports it. It means one of three things, and the check does not hide which: an observation's
  error is understated, the time correlation is too long (the atmosphere really changed), or the
  model is missing something.
- **Leave-one-out** gives the predictive error away from observations; the posterior covariance
  gives how it grows with distance in space and time.
- **The Monte Carlo** samples the wind jointly with everything that built it, `p_ref` among it (§9),
  and delivers the envelope.

## 9. The Monte Carlo, sampled jointly

**The wind and the anchors are not independent (author, 5 November 2026).** The anchors build the
wind, and the model is the operator that does it. So the wind's uncertainty already contains the
anchors' errors. Drawing anchor perturbations and wind perturbations independently would count the
anchor errors twice, and would pair anchor and wind states that cannot occur together: a member
whose anchor is warmer, carried by a wind that does not know it. The envelope would be too wide,
in the wrong places.

**Each member rebuilds its own wind.** Per member:
1. **Perturb every observation within its error:** the anchors' reduction parameters, giving that
   member's kind N instances; IRIS; the cloud winds; the isobar geometry; the targets.
2. **Run the wind tool** on the perturbed inputs, giving that member's best-guess wind at each
   epoch.
3. **Run the model** with that member's anchors and that member's wind, and the forward-model
   parameters being varied (`p_ref`, the propagator's rule and its parameters, and the like).

Each member is internally consistent: its wind is the best wind for its anchors. The spread across
members is the uncertainty, with the correlations built in. This is the driver of the manuscript's
plan (reduction parameters varied first, then forward-model parameters for each instance of N),
with the wind tool as the step between them. In the language of data assimilation it is the
perturbed-observations ensemble.

**The model as the operator is not circular, on two conditions.** Each observation counts once
(§6). And the model's own errors are carried as error and sampled, not fitted away: its
approximations (the equatorial bridge) enter `R` as model error, and its declared alternative
rules are drawn by the Monte Carlo. Errors in the model's physics cannot be found by fitting the
model to data.

**Cost.** The wind fit is nearly linear (SPEC_17). Linearized about the best guess, a member's wind
correction is a fixed linear map, the gain, applied to that member's perturbations: one matrix
product per member, not a fresh fit. A member whose perturbations fall outside the linear range
gets a fresh fit, and the record says which members did. The forward runs remain the Monte
Carlo's real cost, as they always were.

**Later, the same thing explicitly.** If anchors are allowed to move within their own errors, their
corrections join the wind in the fit's state, and the posterior covariance holds the anchor and
wind correlations directly. Sampling that joint posterior is equivalent to the per-member rebuild
in the linear case. The per-member rebuild comes first, because it needs nothing new.

## 10. The path

1. **The Voyager-epoch shear blend:** IRIS and the radio targets of the Lindal anchor, the
   reference wind held at the cloud wind. Hypothesis 3. Its first ruling is the error budget.
2. **The reference wind:** the cloud winds and the 100 mbar geometry; `p_ref` sampled.
3. **More anchors at the Voyager epoch:** the other Voyager profiles reduced to kind N, entering as
   anchor pairs; the model's combination reduced to Eq. A30 (§6).
4. **Time:** the joint estimate across epochs with the time correlation of §7; the seasonal
   propagator for anchors; Kliore at 1979; CIRS cross-sections and Cassini cloud winds.
5. **The Monte Carlo,** each member rebuilding its own wind (§9).

Each is one specification or a few.

## 11. Ruled (author, 5 November 2026)

1. **The model does not correct the wind (§6).** The wind tool uses the anchors; the model's
   combination keeps the inverse-variance blend of Eq. A30. This amends the seasonal note §6.1 item
   3 and the manuscript's Eq. A35.
2. **The time correlation (§7)** is set by the radiative time constant `τ(φ, p)`: published values
   for Saturn first, replaced by the propagator's own fitted values once it exists.
3. **The IRIS error model** is random by default. A systematic part is added only where comparison
   with independent observations shows a bias. The fit's consistency check is that comparison; the
   hypothesis 3 specification sets how it sizes the systematic term.
4. **`p_ref` is a Monte Carlo parameter.** Each member draws its own reference level, and its wind
   is rebuilt with it in the wind tool (§9 step 2). The range and distribution, by default uniform in
   `ln p` from the upper haze to the ammonia cloud top, come from the cloud-structure literature
   when the Monte Carlo is specified. If the envelope shows the level matters, the 100 mbar geometry
   is the constraint that can later estimate it.

## 12. Revision history

- v0.1, 5 November 2026: first draft.
- v0.2, 5 November 2026: every radio profile is an anchor, same-epoch pairs constrain the wind
  through each other, once per pair; a best-guess wind per epoch, all epochs estimated together
  with a time correlation; the check across epochs; the model's combination keeps Eq. A30 only.
- v0.3, 5 November 2026: the Monte Carlo sampled jointly, each member rebuilding its own wind (§9);
  the model does not correct the wind (author).
- v0.4, 5 November 2026: accepted; §11 ruled (the time correlation from `τ`, IRIS random by
  default, `p_ref` a Monte Carlo parameter).
- v0.5, 5 November 2026: §5, the prior carries no observation; IRIS enters as an observation only
  (SPEC_19 §4).
