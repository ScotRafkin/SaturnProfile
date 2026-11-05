# CASSPIAN status, 22 October 2026 (updated 30 October)

Reviewing agent. Supersedes `casspian_status_2026-10-08.md`. The results of SPEC_13 and the first
comparison are in that note and are unchanged; this note records what has happened since.

## Where things stand

- **Closed:** SPEC_00 to SPEC_16 and SPEC_R1. SPEC_R1, SPEC_14, SPEC_15 and SPEC_16 closed on 30
  October. SPEC_00 is at v0.25: its §8 is now the whole of what readers enforce.
- **Committed:** SPEC_R1 at `15965af`, SPEC_14 at `cfe8030`.
- **In progress:** SPEC_17, the Voyager-epoch constraints on the wind.
- **Accepted, waiting:** SPEC_18, Kliore et al. (1980) as a second registered profile.
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

## Next

1. SPEC_17, then SPEC_18.
2. **A fitted wind at the Voyager epoch** from IRIS, the Voyager radio profiles and the 100 mbar
   geometry, on SPEC_17's operator, with leave-one-out. It needs a ruling first on the 110 mbar
   radio and IRIS disagreement. The equatorial profiles above 30 mbar enter with a large error.
3. **Time:** the seasonal propagator, CIRS seasonal fields and Cassini cloud tracking. The CIRS
   shear hypothesis of the earlier queue comes here. Steps 2 and 3 may swap. Kliore then enters as an anchor at its own time.
4. **The combination step:** Kliore and others as construction anchors, and optimal interpolation.
5. **Later**, unchanged from the 8 October note: the IRIS profile catalog; stratospheric IRIS near
   1 mbar; the dry adiabat for normal hydrogen in lapse-rate plots; how NaN uncertainties are used;
   the Monte Carlo; convective adjustment mixing entropy, momentum and composition; cylinders at
   depth; a seasonal wind prior; the assimilation pulling hypothesis 1 toward the radio egress.

## Parked

Decision L (the linear wind interpolant), the branch `spec07-l2-parked`, the tabled redesign of the
vertical on a staggered grid, and the radio versus IRIS disagreement at 110 mbar (about 5 K).
