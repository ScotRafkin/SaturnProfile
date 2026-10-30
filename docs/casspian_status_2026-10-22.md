# CASSPIAN status, 22 October 2026 (updated 30 October)

Reviewing agent. Supersedes `casspian_status_2026-10-08.md`. The results of SPEC_13 and the first
comparison are in that note and are unchanged; this note records what has happened since.

## Where things stand

- **Closed:** SPEC_00 to SPEC_14 and SPEC_R1 (both 30 October; REVIEW_R1_step1, REVIEW_14_step1). SPEC_00 is at v0.25:
  its §8 is now the whole of what readers enforce.
- **SPEC_R1** is committed at `15965af`. SPEC_14's acceptance commit is next.
- **Next:** more radio targets (below).
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

## Next

1. SPEC_14's acceptance commit and push.
2. **More radio targets, in one digitizing step:** Lindal Fig. 6 (Voyager 1 egress near 3.0 S, in
   the equatorial band; Pioneer 11 near 11.9 S) and the 74.2 S profile of Fig. 5. Then compare the
   hypotheses at several targets. Lindal rates the Voyager 1 profiles less reliable.
3. **Hypothesis 3:** shear from a CIRS cross-section of the same season (Aug 2009 to Nov 2010,
   before the storm), used outside about +/-15 to 20 deg.
4. **Later**, unchanged from the 8 October note: the IRIS profile catalog; stratospheric IRIS near
   1 mbar; the dry adiabat for normal hydrogen in lapse-rate plots; how NaN uncertainties are used;
   the Monte Carlo; convective adjustment mixing entropy, momentum and composition; cylinders at
   depth; a seasonal wind prior; the assimilation pulling hypothesis 1 toward the radio egress.

## Parked

Decision L (the linear wind interpolant), the branch `spec07-l2-parked`, the tabled redesign of the
vertical on a staggered grid, and the radio versus IRIS disagreement at 110 mbar (about 5 K).
