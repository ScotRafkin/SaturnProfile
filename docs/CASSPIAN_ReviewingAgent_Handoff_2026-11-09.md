# CASSPIAN reviewing agent handoff, 9 November 2026

For the next reviewing agent. Supersedes `CASSPIAN_ReviewingAgent_Handoff_2026-10-08.md`. Read this,
then `casspian_status_2026-10-22.md` (updated 9 November) and
`CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` (v0.5). Physics and history are in the specs.

## 1. The role

- **You are the reviewing agent and the spec agent.** You write the specifications
  (`docs/specs/SPEC_xx_*.md`), rule on the coding agent's reports, and write the reviews
  (`reports/REVIEW_xx_stepN.md`). You keep the header lines of `docs/specs/STATE.md`; the step rows
  are the coding agent's. You never write code and never run git.
- **The author** is S. Rafkin, a planetary atmospheric scientist. He decides. Recommend, and say
  plainly when you think he is wrong. He wants pushback, including on your own earlier rulings.
- **The coding agent** codes, runs the acceptance, writes `reports/REPORT_xx_*.md`, and commits and
  pushes after your review and the author's go. Messages to it go through the author; end each
  review or ruling with a short message he can paste.

## 2. Mechanics

- **The repository** is on the author's Windows machine, `C:\Users\srafkin\SaturnProfile`, mirrored
  at GitHub `ScotRafkin/SaturnProfile`. Reach it through the device bridge: `device_list_dir`;
  read by `device_stage_files` (files land under `/mnt/user-data/uploads/SaturnProfile/...`); write
  by `device_commit_files`.
- **The stale-copy fault.** Write each new version to a fresh `/mnt/user-data/outputs/vN/`, commit
  with `force: true`, then stage the file back and `cmp` it. Every time. Restage any file you are
  about to edit first: the coding agent edits `STATE.md` too.
- **The project** holds the specs as `claude/SPEC_xx_<name>_draft.md`, plus status notes and
  handoffs. `project_write` needs `local_path` inside `/home/claude/`.
- **Reports.** The coding agent may name a rerun's report with a suffix (`REPORT_19_step1_v04.md`);
  list the directory rather than restaging an old name. Read every figure the report cites.

## 3. The author's rules

- **Materiality (9 November).** A Monte Carlo mean within a few kelvin of what a probe measures, and
  winds within about 10 percent, would be an excellent result. A change worth less than about 1 K in
  delivered temperature or about 5 percent in wind at probe levels gets a comment, not a
  specification. Algorithm checks (closure, convergence, synthetic recovery) stay. Test every
  recommendation, the coding agent's and your own, against this before writing it into a spec.
- **The lean rule (SPEC_R1).** A refusal, check or required key stays only if without it the code
  would compute a wrong number silently or fail obscurely. No new keys unless code reads the value;
  sources, bases and choices go in comments. Nothing may accept only one kind of source or one file
  name.
- **Regression by risk.** When tools change, Lindal's original results must come out the same; if
  they do, no full regression. Push after every accepted commit.
- **Lindal is untouched** by other sources' work. `refrac` changes only to fix an error. A new
  source gets its own ingest so it can use the same tools.
- **Composition.** Temperature always comes from refractivity with the appropriate composition, with
  true molecular masses. Reproducing a source uses that source's composition; each run and the wind
  fit choose one composition for all their inputs.
- **Time is not space.** Spatial wind corrections are not mixed with time corrections. Pioneer 11
  (1979) waits for the time machinery.
- **Every wind is a hypothesis,** with its own kind W and run; earlier ones are kept.
- **Comparisons are absolute** (temperature and density at the target). Every comparison figure set
  starts with the absolute profiles, observed with its error band and delivered, before any
  difference or residual plot. Every temperature profile panel carries the dry adiabats for normal
  and equilibrium hydrogen, `c_p` from composition, started from the deepest observed point.
  Sections show zonal wind as dashed contours; no hatching. Latitudes planetographic, with planetocentric in parentheses.
- **Prose without dashes** as parentheticals, in everything you write. Lean documentation.

## 4. Where things stand

- **Closed:** SPEC_00 to SPEC_19 and SPEC_R1. SPEC_19 closed on 9 November at v0.5, committed at
  `5493fee` with the 0.3 K radio floor (5 passes; 110 mbar at 31.2 S +3.33 K).
- **Closed: SPEC_20 v0.5** (REVIEW_20_step1). Three new anchors; `casspian-section`; the anchor matrix.
  At its acceptance commit: heat capacity from composition with both hydrogen adiabats; SPEC_15's
  egress gravity fixed; zonal wind as dashed contours on the sections, hatching removed. A single
  anchor carries everything it holds; the product blends all anchors (author). At 110 mbar the
  southern profiles miss each other by 6.5 to 6.8 K: the IRIS weighting, to SPEC_21.
- **Error budget rules for later fits (author, 9 November):** a source's own reliability statement
  inflates its errors (Lindal's Voyager 1 profiles, a factor of 2 as a start); structure a balanced
  zonal model cannot represent is representativeness error; robust weighting in the solver when
  Cassini data arrive. The Voyager fit stays as is.
- **The 3.0 S inversion** (28 K between 20 and 11 mbar) is in Lindal's Fig. 6 as drawn. Pioneer 11 has
  a similar single step near 4 mbar. Artifact of the method or a real equatorial feature: Cassini
  will tell. The author is skeptical of it.
- **The wind fit** (`casspian-wind-fit`, `src/casspian/tools/wind/estimation/`): observation types
  (reader, operator, error covariance), a prior, one solver. Two types: radio targets and
  temperature gradients. The budget as closed is in SPEC_19 §2 and the status note.
- **Hypothesis 3** (`fits/voyager_fit_1/`): at 110 mbar, delivered minus observed +3.3 K at 31.2 S,
  +0.6 K at 3.0 S, -3.6 K at 74.2 S. Within about 10 deg of the equator the shear is the prior's and
  its stated uncertainty too small: IRIS has no rows there, and the balance ties shear to the
  curvature of temperature, not its gradient.

## 5. Next

The plan in `casspian_status_2026-10-22.md` (Next and "The general error budget"), in order:
1. **Design note v0.6:** the general error budget, no rule for one dataset: state the wind; shear
   prior from a temperature-anomaly covariance through the balance; every amplitude, length and
   form an input; one representativeness error by pressure for every observation; robust
   weighting; an IRIS systematic term with the 110 versus 150 mbar ambiguity. Draft it for the
   author first.
2. **SPEC_21:** that budget built, every Voyager radio profile an anchor in the fit, hypothesis 4.
3. **SPEC_22:** the blended product, every anchor combined (Eq. A30): the final single answer, the
   profile matrix with every anchor and leaving each out, and the sections.
4. **SPEC_23:** cloud winds (kind W, one new operator, the reference wind estimated).
5. **Cassini radio occultations and CIRS** with time; the equatorial oscillation as climatology.
6. **The Monte Carlo,** with helium and `p_ref`. **Check helium first:** Lindal used about 6 percent
   by number; later estimates are higher (Conrath and Gautier 2000); confirm before ruling.

**Deferred from SPEC_19** (v0.5 §5 item 7): the IRIS systematic term, a latitude-varying prior
amplitude, the jet's bias study, 2.5 deg knots, the 3.0 S wave term.

**Parked:** decision L (`spec07-l2-parked`); the vertical on a staggered grid.

## 6. Lessons from SPEC_19

- **A prior's size follows its mean.** Moving the mean to zero without rescaling the amplitude broke
  the first fit. When a spec changes what a term represents, recheck its numbers.
- **Synthetic tests with a physical truth** (hypothesis 1) find what draws from the prior cannot.
- **The coding agent's next-step lists are thorough and often below materiality.** Rule on each
  against the yardstick; defer, don't spec, what falls under it.
