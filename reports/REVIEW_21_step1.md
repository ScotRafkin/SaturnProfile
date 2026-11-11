# REVIEW 21, Step 1. The general error budget, and hypothesis 4

Reviewing agent, 10 November 2026. Report: `reports/REPORT_21_step1.md`. Specification:
`docs/specs/SPEC_21_General_Error_Budget.md` v0.5, closed at v0.6.

**Accepted, on check 5 finishing as ruled below. SPEC_21 closes with its acceptance commit.**

- **Every fit converged under v0.5.** The stability rows now match the forward run (actual over predicted
  0.99 to 1.03); the constraint acts only on opening passes and nowhere at a final fit. The Voyager 2
  egress's own excess at 794 to 891 mbar is reported, not penalized, as ruled.
- **Hypothesis 4 improves the matrix materially.** At 110 mbar the two southern profiles now miss each
  other by 2.7 K (from 6.7 under hypothesis 3); 31.2 S from Lindal is -0.8 K (from +3.3); the 10 to 12 K
  misses at 10 mbar involving 3.0 S fall to 0.1 to 2.1 K. IRIS residuals sit within their stated errors;
  only IRIS 110 mbar needed a systematic part, sized by the check (1.12 after). Huber downweighted nothing.
- **Hypothesis 4 reproduces the 3.0 S inversion with balanced equatorial shear** (the matrix figure): the
  transfers from every other anchor follow it. That is the behavior the author asked of a balanced model
  if the feature is real; whether it is real stays for Cassini, and one equatorial profile is thin
  support. The equatorial wind aloft that does this must be reported (below).
- **The budget is generous, not optimistic.** The overall ratio is 0.14, the radio groups 0.10 to 0.16,
  and leave-one-out errors 0.10 to 0.21: the stated uncertainties are about 2.5 to 3 times the misfits in
  sd. Equal scaling would not move the best guess, only its uncertainty. Not tuned on four profiles, as
  ruled; to be set from CIRS statistics. Recorded as a finding for the Monte Carlo, whose envelope it
  would otherwise widen.
- **1b** does not recover an equatorial layer at 15 mbar and says so (posterior sd 230 to 430): honest,
  and the reason wind observations come next. **1d:** the floor brought the prior's equatorial draws from
  thousands to about ±1900 m/s of wind change; the posterior at the equator is what the Monte Carlo
  samples, and it is far tighter (58 to 223 m/s per unit `ln p` at 110 mbar).
- **1c** fails its bound only on IRIS 290 mbar, accepted as sampling (500-draw mean 1.04, 4.5 effective
  degrees of freedom). From now on the 1c bound is judged on the 500-draw means, with the 20-draw values
  reported.

**Departures:** `radio.temperature_on_wind_grid` read a variable name written for Lindal's anchor; now it
takes the run's single refractivity variable. Accepted: it was exactly a rule for one dataset.

**Check 5, ruled:**
1. **The wind file's top.** Every wind file the fit writes extends at least one decade above 1 Pa (the
   forward grid's top), the shear held at its top value there, which is the state's own rule above its
   highest knot. The forward run must never stop because the column's top moved by 0.3 percent. Not a
   damped first step: that would hide the edge, not remove it.
2. **Then finish check 5** from the resume point: the Voyager 2 egress-referenced fit, its delivered
   temperatures at the four profiles against the Lindal-referenced fit, bound 0.5 K.
3. **Add to the report:** hypothesis 4's equatorial wind at 30, 10, 3 and 1 mbar with its posterior sd,
   beside hypothesis 3's, since it now carries the 3.0 S structure.

**Closing:** if check 5 passes and the Lindal-referenced values are unchanged by item 1, the acceptance
commit follows with no further review: the code, `fits/voyager_fit_2/`, the runs' control files,
`tests/step21_1/`, the reports and reviews of Step 0 and Step 1, SPEC_21 v0.6, the design note v0.8,
`STATE.md` with the hash and check 5's largest difference. Rebuild clean; push. If check 5 fails its
bound, report it.
