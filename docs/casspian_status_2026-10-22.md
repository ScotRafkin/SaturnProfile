# CASSPIAN status, 22 October 2026 (updated 9 November)

Reviewing agent. Supersedes `casspian_status_2026-10-08.md`. The results of SPEC_13 and the first
comparison are in that note and are unchanged; this note records what has happened since.

## Where things stand

- **Closed:** SPEC_00 to SPEC_19 and SPEC_R1. SPEC_17 and SPEC_18 closed on 5 November, SPEC_19 on
  9 November. SPEC_00 is at v0.25: its §8 is now the whole of what readers enforce.
- **Accepted:** `docs/CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` v0.5, the plan for the
  best-guess wind (below).
- **Next:** cloud-tracked winds as the third observation type (below, Next).
- **The coding agent is new.** It works from `docs/CASSPIAN_CodingAgent_Handoff_2026-10-22.md`
  and `docs/CASSPIAN_CodingAgent_Brief_2026-10-22.md`.

## SPEC_R1, the lean pass

**Why.** The regression of 8 October on `cf33cd6` found three open items. Two were checks
refusing correct data: the closure refused over five descriptive attributes that SPEC_00 v0.24 had
dropped, and `step9` check 6 failed on two data files its own list did not name. The third was
the regression driver missing the seven newest suites. Much more data is coming, from many
sources, and nothing may accept only one kind of source or one file name.

**The rule (author, 22 October).** A refusal, a test check or a required key stays only if without
it the code would compute a wrong number without saying so, or fail with no clear message.
Everything else goes. Removing it is the fix, not loosening it.

**What goes.**
- The closure content comparison, in full, and the closure season equality.
- The refusals on names and labels: `profile_or_run`, `role`, `composition_role`, run and manifest
  prefixes, the namelist named for its run, `p_b_location`.
- The `-dirty` refusal. Products still record the commit with its `-dirty` mark.
- `step9` check 6, the directory listing.
- About 25 descriptive requirements of the schema: required `units`, `long_name`, `provenance`,
  `title`, `source` and `epoch`, the bare-name refusal, the source profile globals, the dimension
  markers, `positive` and `direction`, among others.
- From the not-placed list: path containment except one rule, the season required on every file,
  `package_version`, kind N's hash cross-check, `reject_physical_values`, the stated latitude
  global.
- The suite checks that only exercised removed refusals, with their counts amended in
  `tests/run_regression.sh`.

**What stays.**
- Unknown keys refused (a typo), missing or mistyped values, and the cross-checks that prevent a
  wrong number: rotation of W against R, the anchor rule and its quantity, the geodesy and thermo
  levels, composition levels, the gauge isobar.
- The wind's total present and zero at the poles; coverage spanning the run; no extrapolation.
- Required variables, dimensions and groups; monotonic coordinates; mole fractions summing to one;
  float64; `casspian_kind` and the schema version; `harmonic_convention`.
- The named numerical failures.
- A run writes only inside its own run directory, and `casspian-new-run` never overwrites an
  existing directory.

**Writers** stop writing only what existed for a removed check, plus
`observation_level_justification`. The records (`production_record`, `transfer_record`,
`reduction_record`, the embedded namelist) are results and stay. The descriptive globals stay
written, required by nothing. `role` and `p_b_location` leave every committed build file and
namelist.

**Housekeeping in the same commit.** The untracked reference files in `docs/` (the 10-03 handoff,
the Sanchez-Lavega files, the Conrath et al. IRIS PDF) are committed, so products are stamped clean.
The seven newer suites join the driver.

**Accepted (30 October).** 38 of 38 suites at their amended counts, kind N's hash unchanged, no
computed value changed. `step05_3` checks 4 to 7 passed, their first run since SPEC_11 Step 2.
The registered kind T now warns once on reading, because its `input_hashes` names the edited
`lindal_build.toml`; it stays until the reduction is next redone. `tests/step03_0/` is removed.

## SPEC_14, hypothesis 2: the shear reversed aloft

- **The case.** `above_top = "reversed"` in `tools/lindal/lindal_wind.py`. Above 110 mbar the 110
  mbar shear is multiplied by f(p): f = 1 to `hold_Pa` (1500 Pa); then linear in ln p through
  zero at `cross_Pa` (450 Pa, where the two radio profiles cross); held at its 0.2 mbar value above
  0.2 mbar, the anchor's top.
- **The run.** `lindal_iris_v_sl_iii_110`, made from `lindal_iris_v_sl_ii_110` and differing in
  `[shear]` only. Tuned to the Voyager 2 egress, so the egress checks consistency, not
  independence.
- **Ruled (author, 22 October).** The equatorial wind above 0.2 mbar, about 940 m/s at 1 Pa
  against 270 m/s in the base, is accepted: this hypothesis exists to work the pipeline with an
  input wind, and later hypotheses bring physical boundary conditions. The first guesses stand.
- **Result (REPORT_14_step1, accepted 30 October).** Delivered minus egress at 31.2 S:

  | Level (mbar) | 730 | 290 | 110 | 10 | 4.5 | 1 | 0.3 |
  |---|---|---|---|---|---|---|---|
  | Hypothesis 1 base (K) | +0.40 | +2.32 | +5.00 | -1.90 | +0.10 | +5.24 | +15.10 |
  | Hypothesis 2 (K) | +0.40 | +2.30 | +5.09 | +0.20 | -0.77 | -2.17 | +2.88 |
  | Egress sigma_T (K) | 2.73 | 2.29 | 2.27 | 2.88 | 3.47 | 5.42 | 8.74 |

  Aloft, every value is inside sigma_T; density at 0.3 mbar moves from -10.0 to -1.8 percent. At
  110 mbar the radio and IRIS disagreement remains, about 5 K. The equatorial wind reaches 472 m/s
  at 0.2 mbar and 1052 m/s at 0.01 mbar, as accepted. The lapse rate stays stable.
- **Check 1** was ruled to the construction's tolerance: the wind below 110 mbar differs from the
  base's by 4.5e-3 m/s, from the stopping rule, which is immaterial.

## SPEC_15, three more radio profiles

- **Data.** Lindal Figs. 5 and 6 traced along hand-placed guides: Voyager 1 ingress at 74.2 S and
  egress at 3.0 S (`occul_data/lindal/`). Pioneer 11 from Kliore et al. (1980, JGR) Table 1,
  transcribed as published (`occul_data/kliore1980/`), rescaled from 90/10 to 94/6 by Lindal's
  rule where compared. Lindal's own Fig. 6 curve agrees with the rescaled table to 0.16 K, so he
  rescaled it the same way. Epochs are in file comments: Voyager 1 on 13 November 1980, Pioneer 11
  on 1 September 1979.
- **A figure bias.** From 600 to 1200 mbar, Lindal's Fig. 5 is drawn 1.5 to 3 K warmer than his
  Table I. The 74.2 S profile probably carries it. It is recorded as observation error, not
  corrected.
- **Result.** Delivered minus observed (K):

  | Target | Profile | 730 | 290 | 110 | 10 | 4.5 | 1 | 0.3 mbar |
  |---|---|---|---|---|---|---|---|---|
  | 74.2 S | hypothesis 1 | +2.52 | +0.12 | +0.40 | -8.10 | -6.30 | +4.21 | +10.76 |
  | | hypothesis 2 | +2.52 | +0.11 | +0.50 | -5.88 | -7.26 | -3.70 | -2.29 |
  | | anchor | +5.56 | -0.23 | -4.78 | -10.12 | -7.54 | +3.66 | +10.33 |
  | 3.0 S | hypothesis 1 | +0.08 | +2.21 | +0.76 | -10.30 | +6.84 | +16.48 | |
  | | hypothesis 2 | +0.08 | +2.20 | +0.84 | -8.27 | +6.32 | +10.04 | |
  | | anchor | +5.29 | +3.26 | -2.27 | -10.97 | +6.99 | +17.29 | |
  | 11.9 S | hypothesis 1 | | | +2.38 | +16.86 | +27.36 | | |
  | | hypothesis 2 | | | +2.52 | +20.06 | +26.09 | | |
  | | anchor | | | -4.15 | +15.10 | +26.80 | | |

- **What it shows.**
  - Below about 30 mbar the transfer works at every target: it removes the anchor's 4 to 6 K warm
    bias at 730 mbar and is within about 2.5 K at 110 mbar, except at 31.2 S.
  - Above about 30 mbar near the equator, the observed profiles carry warm layers of 10 to 30 K,
    which Lindal reads as waves. No balanced zonal wind of this form makes them, and the epochs
    differ. They are representativeness error and time, not a wind the transfer could carry.
  - Hypothesis 2 improves 74.2 S aloft and is 1 to 3 K worse near the equator from 20 to 80 mbar,
    where it holds the bridged shear to 15 mbar.
  - At 74.2 S both are 6 to 8 K cold from 3 to 15 mbar; a pattern times f(p) cannot reorder the
    profiles with height.

## SPEC_16, Pioneer 11 refractivity

- **The self-check agrees.** Kliore's temperatures and pressures through our refractivity
  arithmetic return his measured column; the implied mean refractivity matches ours to 0.02
  percent. The 17.1 mbar row is a misprint for 17.5, corrected.
- **The geometry agrees.** On the planetocentric reading of his radius mapping, Pioneer 11's
  100 mbar level is 7 km from the surface the model delivers at 9.51 S from the Lindal anchor.
- **The profile does not, above 10 mbar.** From 10 to 182 mbar the refractivity matches to 2.3
  percent; above, the observed column is colder and more compact, the layering SPEC_15 found.
- **Kliore is Pioneer 11, September 1979,** two years before the anchor and across the equinox.
  It is not used to correct the Voyager-epoch wind, which would mix time into a spatial
  correction (author, 30 October). It becomes a second anchor and a cross-validation once time
  propagation exists.

## SPEC_17 and SPEC_18 (accepted 30 October)

- **SPEC_17** measures what the Voyager-epoch observations ask of the wind, and fits nothing.
  - The Voyager 1 profiles are near-epoch with the Voyager 2 anchor; as distinct times later.
  - Part A: Lindal's Fig. 9 100 mbar altitudes against the delivered 100 mbar surface, the
    reference wind's test.
  - Part B: the linearized balance weight, the observation operator; every path starts at the
    anchor, so the three Voyager targets split the shear constraint into three latitude
    segments; the linearization tested at two amplitudes.
- **SPEC_18** reduces Kliore to a kind N with a closure run. Three layers: source-specific ingest,
  general tools from the standard kinds on, the model unchanged. Kliore's table has no 100 mbar
  level, so his kind T gets one interpolated level (about 0.01 K). Lindal's answer must not
  change: every changed tool in his chain rebuilds his files in a copy, equal in value, and with
  that no full regression is run.

## SPEC_17 and SPEC_18, results (5 November)

- **SPEC_17.** The anchor's 100 mbar surface matches Lindal's Fig. 9 to 0.8 km; south of the equator
  the delivered surface is 7 to 8 km high, a constraint on the reference wind and the shear below
  100 mbar. The operator (`balance_weight`) is linear to about 1 percent. The segments ask little of
  26 to 71 S; segment 1 (30.8 N to 2.4 S) is weakly constrained; above 30 mbar near the equator the
  requests swing with the wave layering; and 2.4 to 26 S at 110 mbar asks -7 m/s per unit ln p,
  the radio and IRIS disagreement located in latitude.
- **SPEC_18.** Kliore is a working kind N with a closure run, and Lindal's answer is unchanged. The
  profile tools are split into source-specific ingest and shared writers. The closure's drift is
  Kliore's own: his molar mass times gravity is 0.6 percent below ours, likely 2.2 amu. When Kliore
  becomes an anchor, re-deriving his T and p from his refractivity with the exact mass is to be
  considered.

## The best-guess wind (design note v0.5, accepted 5 November)

- **One weighted least-squares fit,** outside the model, using the model as the operator: a prior
  wind and its covariance, and one term per observation type (data, operator, error model). Adding
  data adds a term.
- **Every radio profile is an anchor.** Same-epoch anchor pairs constrain the shear integral between
  them, once per pair. The model does not correct the wind; its combination keeps Eq. A30 only.
- **IRIS shapes the shear in latitude; the radio sets the long integrals.** The fit weighs them by
  their errors, including how those errors correlate in latitude. IRIS errors are random by default;
  the consistency check finds any bias.
- **A best-guess wind per epoch,** all epochs estimated together, coupled by a time correlation set
  by the radiative time constant. Near epochs that pull the winds far apart are flagged.
- **The Monte Carlo samples jointly:** each member perturbs every observation and rebuilds its own
  wind. `p_ref` is one of its parameters.
- **Later data:** Cassini radio occultations (1.5 bar to 0.1 mbar) as anchors, with same-day
  ingress and egress pairs; CIRS nadir and limb as gradient observations, limb giving the shear
  aloft; both need the time machinery.

## Direction (author and reviewing agent, 30 October)

- **The goal.** The model gives the best answer for a given input, using every observation that
  constrains its physics. Every observation incorporated pulls the answer toward itself at its
  time and place, within its error. Away from observations, the framework and the observations
  together give the minimum-error answer, with uncertainty growing with distance in space and time.
- **The wind is the whole wind:** the shear, from temperatures, and the reference wind, from cloud
  tracking or anything that measures a wind at a level. The reference level is itself uncertain
  and is a Monte Carlo variable.
- **A misfit is mostly wind error, but not only.** Observation error and bias, the missing time
  step, and the model's own approximations take their share first.
- **Order: improve the wind, then blend.** A fitted wind (adjusted within its uncertainty to fit
  the observations) carries corrections everywhere. Optimal interpolation then handles what a
  balanced wind cannot, with honest weights, on a background already close.
- **Hold profiles out.** Each profile used to build the wind stops being an independent test;
  leave-one-out gives the honest error.
- **Not yet in the model:** cloud tracking beyond Voyager, and transport in time (the seasonal
  propagator, CIRS seasonal fields).

## SPEC_19, the wind fit (closed 9 November)

- **The tool:** `casspian-wind-fit` (`src/casspian/tools/wind/estimation/`). Observation types
  (reader, operator, error covariance), a prior, and one solver: weighted least squares, iterated,
  with the posterior covariance, the consistency check and leave-one-out. Two types so far: radio
  targets and temperature gradients. New data add a type or a term, not a tool.
- **The budget as closed:** prior 40 m/s per unit ln p, 5 deg in latitude, one scale height in ln p;
  radio sigma_T from Lindal Fig. 4, Matérn 3/2 over 0.5 scale height, a 0.3 K floor; IRIS random,
  every 4 deg, as a half scale height layer mean.
- **The first budget failed honestly.** With the prior's mean moved to zero, a 20 amplitude and a
  10 deg correlation could not hold the equatorial jet, and the fit's 110 mbar answer came from that.
- **Hypothesis 3** (`voyager_fit_1`), delivered minus observed at 110 mbar: +3.3 K at 31.2 S (from
  +5.0 under hypothesis 1), +0.6 K at 3.0 S, -3.6 K at 74.2 S. Better than hypotheses 1 and 2 at
  most levels, most of all aloft. The radio and IRIS disagreement at 110 mbar is now 1.5 sigma_T at
  31.2 S.
- **The equator is a physics limit.** IRIS has no rows within about 8 deg, and there the balance
  ties shear to the curvature of temperature, not its gradient. The equatorial shear is the prior's
  and its uncertainty is understated. About 1 K in delivered temperature; for the equatorial wind
  aloft, not a result.

## Materiality (author, 6 and 9 November)

- **The yardstick.** A Monte Carlo mean within a few kelvin of what a probe measures, and winds within
  about 10 percent, would be an excellent result.
- **The rule.** A change worth less than about 1 K in delivered temperature or about 5 percent in wind
  at probe levels gets a comment, not a specification. Algorithm checks (closure, convergence,
  synthetic recovery) remain, since they test the tool, not the answer.
- **What the fit does.** The observations, within their errors, are the answer where and when they
  exist. The model fills between them with balanced, smooth temperature gradients and shears. The
  prior is what keeps the shear physical and smooth; balance alone does not.
- **Where the large errors are,** largest first by the reviewing agent's estimate:
  1. **Helium.** Radio temperatures scale with the assumed mean molecular mass. Lindal used about 6
     percent helium by number; later estimates (Conrath and Gautier 2000) are near 10 to 13 percent.
     About 7 to 10 K in absolute temperature. It cancels in the latitude differences the wind fit
     uses, not in the profile a probe measures. To be checked against current values and made a
     Monte Carlo parameter.
  2. **Epoch.** Voyager is 1980 to 1981; a probe arrives decades later, in another season. Stratospheric
     seasonal swings reach 10 K or more at high latitudes.
  3. **The deep wind** below the clouds, for probe drift. The fit holds the cloud wind and relaxes the
     shear below 730 mbar; Cassini gravity (Iess et al. 2019) has the cloud-level winds extending
     thousands of km deep, which supports that to first order.
  4. **Equatorial winds aloft,** which temperatures cannot set (above).

## Next

In order:
1. **Cloud-tracked winds** as the third observation type, with the reference wind entering the
   state: Voyager and Cassini winds, including those measured at two levels near the equator
   (Cassini imaging in different filters), which constrain the equatorial shear directly; the 100
   mbar geometry; `p_ref` sampled.
2. **Cassini radio occultations and CIRS** (nadir and limb), with the time machinery they need: the
   seasonal propagator and the joint estimate across epochs, Kliore at 1979.
3. **The Monte Carlo,** each member rebuilding its own wind, with helium and `p_ref` as parameters.
4. **More Voyager anchors** (the Voyager 2 egress and Voyager 1 profiles as kind N anchor pairs) only
   if they move an answer by more than the materiality rule.

**Deferred from SPEC_19:** the IRIS systematic term, a latitude-varying prior amplitude (for the
Monte Carlo, if its draws need it), the jet's bias study, 2.5 deg knots, the 3.0 S wave term.

**Later**, unchanged from the 8 October note: the IRIS profile catalog; the dry adiabat for normal
hydrogen in lapse-rate plots; convective adjustment mixing entropy, momentum and composition;
cylinders at depth.

## Parked

Decision L (the linear wind interpolant), the branch `spec07-l2-parked`, and the tabled redesign of
the vertical on a staggered grid.
