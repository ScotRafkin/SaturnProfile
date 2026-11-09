# REVIEW 19, Step 1. The wind fit, and hypothesis 3 at the Voyager epoch

Reviewing agent, 6 November 2026. Report: `reports/REPORT_19_step1.md`. Specification:
`docs/specs/SPEC_19_Wind_Fit.md`, reported against v0.2, revised to v0.3.

**Not accepted. The tool is accepted; the budget is revised (SPEC_19 v0.3) and the step is rerun.**

- **The tool works.** With truths drawn from its own prior it recovers the shear within its
  posterior uncertainty and every group is near its expected value (check 1b). It converges in 3
  passes, the gradient operator returns hypothesis 1's data to a fifth of the stated error, and the
  runs reproduce the fit's last runs exactly.
- **The failure is the prior, and the error is the reviewing agent's.** Design note v0.5 moved the
  prior's mean from the IRIS shear to zero without rescaling the prior. Its 20 m/s per unit `ln p`
  was sized for corrections to an IRIS-based mean, and its 10 deg is wider than the jets. So the
  prior could not hold the equatorial jet, the fit flattened it at 110 mbar (+65 to +11 at 5 N), and
  that one place added about +2.9 K along the path to the radio targets. The fit's 110 mbar answer
  was this, not a weighing of IRIS against the radio. Hypothesis 3 is not adopted.
- **The sharp change near 100 mbar is the top of the IRIS data,** not a fault: the posterior sd of
  the shear rises from 0.7 at 110 mbar to 5.2 at 70 mbar at 50 N. It is drawn as the retrievals'
  half scale height layer in v0.3, which is what IRIS measures, and shown as profiles for the author.
- **The IRIS consistency excess is mostly the prior.** The synthetic case with no systematic error
  gives 8.1, 2.1 and 1.5. A systematic term sized now would absorb the prior's failure. It is sized
  after the refit, and only if the checks still ask.
- **The radio's low ratios test nothing yet.** The exponential correlation expects rough noise the
  printed profiles lack. The form changes; the ruled length stays.

**On the report's points:**
1. **The package beside SPEC_01's `fit.py`:** accepted; SPEC_19 names the path.
2. **The exponential correlation:** replaced by Matérn 3/2 at the same length.
3. **The prior on coefficients (point sd 18.5 for 20):** accepted; the point sd is reported.
4. **Polar functions dropped:** accepted.
5. **First guess:** kept, with one start from the prior mean as a check.
6. **Linear leave-one-out and systematic sizing:** accepted.
7. **Check 1's bound** was the agent's, honestly tightened. v0.3 splits it: hypothesis 1's truth as a
   resolution test without a bound, and 20 prior draws as the statistical test with one.
8. **`fits/`:** accepted. Commit the control file and the wind file.
9. **`pip install -e .`** before the rerun, so the command itself is exercised.

**Order of work:** rerun under SPEC_19 v0.3, report, no commit until acceptance. No regression.

## The rerun under v0.4 (REPORT_19_step1_v04), 9 November 2026

**Accepted. SPEC_19 closes (v0.5 §5).**

- **The tool is sound.** Check 1(b) passes (95.6 percent within 2 sd, group ratios 0.93 to 1.18). The
  fit from the prior mean agrees within 0.04 K. The first fit's failure is gone: IRIS ratios 2.32,
  1.88, 1.04, from about 10, 5 and 3.
- **Hypothesis 3 is within a few kelvin of the radio,** better than hypotheses 1 and 2 at most levels
  and most of all aloft. At 110 mbar: +3.3 K at 31.2 S, +0.6 K at 3.0 S, -3.6 K at 74.2 S.
- **The equatorial jet is a physics limit, not a fault.** IRIS has no rows within about 8 deg, and
  near the equator the balance ties shear to the curvature of temperature, not its gradient. The
  shear there is the prior's, and its stated uncertainty is too small (true +54 returned as +85).
  It moves delivered temperature by about 1 K; for the equatorial wind aloft it is not a result.
  Cloud winds at more than one level, the next observation type, constrain it. Recorded as a
  caveat with `voyager_fit_1`.
- **Materiality (author).** Items below about 1 K in temperature or 5 percent in wind get a comment,
  not a specification. On that basis the agent's items 2, 3, 4 and 6 are deferred (v0.5 §5 item 7).

**On the report's points:**
1. **The iteration** (Anderson(1), plain step to end, `max_passes` 20): accepted. Two starts, one
   answer.
2. **The radio operator floor** (item 5): accepted at 0.3 K, uncorrelated, on each target's diagonal.
   It removes a numerical weight on rough operator detail and may cut the passes.
3. **The Gaussian form of the systematic term:** moot; the term is not adopted.
4. **The test files** run as scratch copies (`survey_truths.py`, `truth_and_noise.py`): run them as
   committed once before the acceptance commit.
5. **Products** rebuilt clean at the acceptance commit.

**Order of work:** no regression.
1. Add the floor, run `casspian-wind-fit` once, and rebuild the three runs. Report three numbers:
   passes, the group ratios, and 110 mbar at 31.2 S. If that value moves by more than 1 K from
   +3.34 K, stop and report.
   Add the absolute profile figure of SPEC_19 v0.5 §5 item 6 (plotting only).
2. Otherwise the acceptance commit: `src/casspian/tools/wind/estimation/`, `pyproject.toml`,
   `fits/voyager_fit_1/` (control file with the caveat comment, wind file), the three run control
   files, `tests/step19_1/`, both reports and this review, SPEC_19 v0.5, `STATE.md`, the status note,
   and `docs/CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` (v0.5).
3. Products rebuilt clean at that commit; the hash and the three numbers in SPEC_19's row of
   `STATE.md`. Push.
