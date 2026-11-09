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
  difference or residual plot. Latitudes planetographic, with planetocentric in parentheses.
- **Prose without dashes** as parentheticals, in everything you write. Lean documentation.

## 4. Where things stand

- **Closed:** SPEC_00 to SPEC_19 and SPEC_R1. SPEC_19 closed on 9 November at v0.5.
- **Pending:** SPEC_19's acceptance commit. The coding agent adds a 0.3 K floor to the radio
  covariance, refits once, and reports the passes, the group ratios and the delivered 110 mbar value
  at 31.2 S. If that value moves by more than 1 K from +3.34 K, it comes back for review; otherwise
  it is recorded in STATE.md and committed (REVIEW_19_step1, rerun section).
- **The wind fit** (`casspian-wind-fit`, `src/casspian/tools/wind/estimation/`): observation types
  (reader, operator, error covariance), a prior, one solver. Two types: radio targets and
  temperature gradients. The budget as closed is in SPEC_19 §2 and the status note.
- **Hypothesis 3** (`fits/voyager_fit_1/`): at 110 mbar, delivered minus observed +3.3 K at 31.2 S,
  +0.6 K at 3.0 S, -3.6 K at 74.2 S. Within about 10 deg of the equator the shear is the prior's and
  its stated uncertainty too small: IRIS has no rows there, and the balance ties shear to the
  curvature of temperature, not its gradient.

## 5. Next

1. **SPEC_20, cloud-tracked winds** as the third observation type, with the reference wind entering
   the state. Voyager first (Sanchez-Lavega et al. 2000 is in `data_static/winds/`), then Cassini,
   including the equatorial winds measured at two levels in different imaging filters (García-Melendo
   et al. 2011), which constrain the equatorial shear directly. The 100 mbar geometry of
   SPEC_17 Part A; `p_ref` sampled later. Find the sources and their error statements before
   drafting; ask the author for papers the tools cannot reach.
2. **Cassini radio occultations and CIRS** with the time machinery.
3. **The Monte Carlo,** with helium and `p_ref` as parameters. **Check helium first:** Lindal used
   about 6 percent by number; later estimates are higher (Conrath and Gautier 2000). The reviewing
   agent estimated 7 to 10 K in absolute temperature from that difference; confirm against current
   values before ruling.
4. **More Voyager anchors** only if they pass the materiality rule.

**Deferred from SPEC_19** (v0.5 §5 item 7): the IRIS systematic term, a latitude-varying prior
amplitude, the jet's bias study, 2.5 deg knots, the 3.0 S wave term.

**Parked:** decision L (`spec07-l2-parked`); the vertical on a staggered grid.

## 6. Lessons from SPEC_19

- **A prior's size follows its mean.** Moving the mean to zero without rescaling the amplitude broke
  the first fit. When a spec changes what a term represents, recheck its numbers.
- **Synthetic tests with a physical truth** (hypothesis 1) find what draws from the prior cannot.
- **The coding agent's next-step lists are thorough and often below materiality.** Rule on each
  against the yardstick; defer, don't spec, what falls under it.
