# CASSPIAN reviewing agent handoff, 8 October 2026

For the next reviewing agent. Supersedes `CASSPIAN_ReviewingAgent_Handoff_2026-10-03.md`. Read this,
then `casspian_status_2026-10-08.md`. Physics and history are in the specs; this note is how to work
and where things stand.

## 1. The role

- **You are the reviewing agent and the spec agent.** You write the specifications
  (`docs/specs/SPEC_xx_*.md`), rule on the coding agent's reports, and write the reviews
  (`reports/REVIEW_xx_stepN.md`). You keep `docs/specs/STATE.md`: a section per spec, a row per step.
  You never write code and never run git.
- **The author** is S. Rafkin, a planetary atmospheric scientist. He decides. Ask him when a choice
  is his; recommend, and say plainly when you think he is wrong. He wants pushback.
- **The coding agent** codes, runs the acceptance, writes `reports/REPORT_xx_*.md`, and commits and
  pushes after your review and the author's go. Messages to it go through the author; end each
  review with a short message he can paste.

## 2. Mechanics

- **The repository** is on the author's Windows machine, `C:\Users\srafkin\SaturnProfile`, mirrored
  at GitHub `ScotRafkin/SaturnProfile`. Reach it through the device bridge:
  - list with `device_list_dir`;
  - read by staging (`device_stage_files`); files land under
    `/mnt/user-data/uploads/SaturnProfile/...`;
  - write with `device_commit_files`.
- **The stale-copy fault.** `device_commit_files` sometimes writes an earlier copy of a file that
  was edited after a previous commit from the same path.
  - Write each new version to a fresh path (`/mnt/user-data/outputs/vN/`) and commit from there,
    with `force: true` for an existing file.
  - Then stage the file back and `cmp` it. Every time.
  - Never edit and commit in the same parallel batch.
- **The project** holds the specs as `claude/SPEC_xx_<name>_draft.md`, plus status notes and
  handoffs. `project_write` needs `local_path` inside `/home/claude/`, so copy there first.
- **Reports and figures.** Read every figure the report cites. The author looks at them too and
  often sees things first.

## 3. The author's rules

- **Not enterprise software.** Regression is by judgment of risk. A change unlikely to move a result
  that matters gets no regression. Full runs come after a batch of changes, when a result could
  plausibly move, or at a milestone he names. Push after every accepted commit.
- **Lean documentation.**
  - Keys and attributes are kept only for values code reads. Descriptive text goes in comments
    (SPEC_00 v0.24 §5, kind W so far). No `*_source`, `justification` or `note` keys.
  - Specs are short and say what to do, not what the code or the paper already says.
  - No captions or metadata nobody asked for.
- **Prose without dashes** as parentheticals, in everything you write for him.
- **Comparisons are absolute:** temperature and density at the target. Changes are informative but
  not sufficient. Every latitude is given planetographic, with planetocentric in parentheses.
- **Every wind is a hypothesis.** Each has its own kind W and its own run, made with
  `casspian-new-run --from`, and differs from its parent in one thing. Earlier hypotheses are kept,
  never overwritten. Every kind file that drives a run stays in that run's build. Files that can be
  regenerated stay out of git.
- **No nonphysical states,** even harmless ones. Fix them at their source, not by filtering.
- **Data is data:** files hold data where observed. Plots show data where it exists. NaN means no
  information.
- **New options need his leave.** The rule is "replace, never add".
- **The shear aloft has always been a guess.** Tuning it against the radio is acceptable when it is
  labeled as a new hypothesis.

## 4. Where things stand

**Specs:**
- **Closed:** SPEC_08 to SPEC_13.
- **SPEC_00:** v0.24.
- **SPEC_13:** v0.5, committed at `cf33cd6`, with SPEC_00 v0.24, SPEC_11 v0.5 and
  REVIEW_11_step2.

**The full regression (8 October, on `cf33cd6`).** The script ran 31 suites, step1 to step08_1,
and 29 matched their reference counts. The coding agent ran the seven newer suites by hand, and all
passed. Logs are in `reports/regression/`. The coding agent's handoff,
`docs/CASSPIAN_CodingAgent_Handoff_2026-10-22.md` (§4), has the detail. Three things are open, and
they come before SPEC_14:
1. **The closure run is refused (step05_3, a real break).**
   - **The failure.** The run's loader (`control.load_run_inputs`) compares each input with the copy
     embedded in the registered anchor product, after removing a fixed list of bookkeeping
     attributes. The rebuilt wind lacks five descriptive attributes that SPEC_00 v0.24 dropped:
     `epoch_note`, `observation_level_justification`, `solar_longitude_source`,
     `source_latitude_convention_source`, and `value_source` on `reference_level_pressure_Pa`.
     Every value is equal.
   - **Reach.** Any closure run built from freshly rebuilt inputs. Transfer runs are not affected.
     step04_5 passes only because it runs the closure from copies of the committed inputs.
   - **The fix is the spec's choice:** add the five names to that list; ignore every descriptive
     attribute; or compare values only. Recommended: ignore descriptive attributes, which is what
     v0.24 means. Do not rebuild the registered product.
   - step05_3 checks 4 to 7 never ran, so they run after the fix.
2. **The directory listing (step9 check 6, benign).** It has failed since SPEC_09.
   `iris_temperatures.csv` (SPEC_09) and `voyager2_egress.csv` (SPEC_12) sit in
   `occul_data/lindal/`, and SPEC_00 §2.2's listing does not name them. Add them to §2.2. Confirm
   whether the suite reads its expected list from the spec or holds its own.
3. **The regression script lacks the newer suites.** step09_1, step10_1, step11_1, step11_2,
   step12_1, step13_1 and step13_2 all passed when run by hand, step11_1 included (5 of 5). What is
   missing is their place in `tests/run_regression.sh`. The coding agent's handoff gives their
   reference counts and the order among them. The full set then takes about 2 h 56 min, more than
   one background job allows on this machine, so the script runs in chunks.

A small fix step, with no new physics, settles all three.

**The runs**, all from the Lindal anchor (36.3 N, 30.81 N planetocentric) to Lindal's egress label,
31.2 S (26.06095 S planetocentric):

| Run | Gradient | Fit window | Cloud wind |
|---|---|---|---|
| `lindal_iris_ii_110` (hypothesis 1) | the fit's slope | adaptive | Ingersoll and Pollard |
| `lindal_iris_v_ii_110` | derivative of the fitted temperatures | smooth | Ingersoll and Pollard |
| `lindal_iris_v_sl_ii_110` (SPEC_14's base) | derivative of the fitted temperatures | smooth | Sanchez-Lavega |

All three use case (ii), with the IRIS shear at 110, 290 and 730 mbar and the reference at
398 mbar:
- relaxed above 110 mbar over 2 scale heights;
- relaxed below 730 mbar over 1;
- PCHIP across |phi_c| < 5 deg;
- tapered to the poles.

**Results** are in the status note. In short:
- The transfer works and closes its round trip, except the equatorial bridge, 0.6 to 0.9 K.
- Below about 300 mbar it beats the untransferred anchor.
- From 50 to 300 mbar it is about 5 K too warm: IRIS has the south warmer at 110 mbar and the
  two radio profiles do not.
- Above about 3 mbar it equals the anchor, which is up to 15 K warmer than the egress, because the
  shear has relaxed to zero.

## 5. Key files

- **Data:**
  - `occul_data/lindal/iris_temperatures.csv` (SPEC_09);
  - `occul_data/lindal/voyager2_egress.csv` (SPEC_12; Lindal Fig. 4, with sigma_T);
  - `data_static/winds/vasavada_saturn_winds.txt` (Sanchez-Lavega 2000 Table II, checked row for
    row);
  - `sanchezlavega2000.toml` and `smith1982_fig4.toml`.
- **Tools:**
  - `tools/lindal/iris_temperatures.py`: `fit`, with `window` adaptive or smooth;
  - `tools/lindal/lindal_wind.py`: `construct`, with `gradient` and `window`;
  - `tools/lindal/egress_profile.py`;
  - `tools/wind/shear.py` (the case `lindal_iris`), `tools/wind/curve.py` and
    `tools/wind/build_wind.py` (the table path).
- **Papers in `docs/`:**
  - Lindal et al. 1985;
  - Conrath and Pirraglia 1983;
  - Sanchez-Lavega et al. 2000, `Sanchez-Lavega_etal_cleanPDF.pdf`.

## 6. Next: the fix step of §4, then SPEC_14, hypothesis 2

Case (iii), "reversed". The base is `lindal_iris_v_sl_ii_110`. Above 110 mbar the 110 mbar shear is
multiplied by f(p):
- f = 1 up to p_hold, first guess 15 mbar;
- linear in ln p from there, through zero at p_cross, first guess 4.5 mbar, where the two radio
  profiles cross;
- continuing at the same slope above, and held constant above 0.2 mbar, the anchor's top.

**Rough expectation:** about +3 K at 10 mbar, about -7 K at 1 mbar and about -14 K at 0.3 mbar,
against the radio's +3.7, -4.9 and -14.7.

**The spec:**
- the new value of the case's above-110 behavior, with its two parameters;
- one run made from the base;
- the comparison figure as SPEC_12 and SPEC_13 draw it, with the IRIS values;
- the u_total figure, watching the equatorial winds aloft.

It is labeled as tuned to the egress.

## 7. The queue after SPEC_14

1. **More radio targets, in one digitizing step:**
   - Lindal Fig. 6: Voyager 1 egress near 3.0 S (inside the equatorial band), and Pioneer 11 near
     11.9 S (1979);
   - Fig. 5's 74.2 S profile.

   Then compare the hypotheses at several targets. Lindal rates the Voyager 1 profiles less
   reliable.
2. **Hypothesis 3:** shear from a CIRS cross-section of the same season, Aug 2009 to Nov 2010,
   before the storm. Use it outside about +/-15 to 20 deg. Interannual change is large in the
   tropics.
3. **Later:**
   - the IRIS profile catalog;
   - stratospheric IRIS near 1 mbar;
   - the dry adiabat for normal hydrogen (about 0.85 to 0.9 K/km) in lapse-rate plots;
   - how NaN uncertainties are used;
   - the SPEC_00 metadata audit of the other kinds;
   - the Monte Carlo;
   - convective adjustment that mixes entropy, momentum and composition together;
   - cylinders at depth;
   - a seasonal wind prior;
   - the assimilation pulling hypothesis 1 toward the radio egress, as kind T.
