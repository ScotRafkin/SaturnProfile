# CASSPIAN handoff: the staggered vertical

3 October 2026, by the outgoing reviewing agent. For the agent that takes up the redesign of the
model's vertical numerics, and through it the coding agent that builds it.

**Tabled (author, 3 October 2026).** The redesign waits until the realistic shear case and the
first comparison with CIRS are done on the current numerics. Use this handoff when the comparison
or the Monte Carlo reopens it, and refresh "The state you start from" first.

## The project in one paragraph

CASSPIAN is a one dimensional Saturn atmosphere reference model for probe mission design, written
for S. Rafkin (SwRI), the author, who approves everything.

- **Anchors.** Occultation profiles (Lindal's Voyager profile first) are reduced to refractivity on
  their own levels.
- **The wind hypothesis.** Supplied as kind W: a reference-level wind plus shear along the local
  vertical.
- **The transfer.** The model carries each anchor's isobars to a target latitude through the shear
  kernel (manuscript Eqs. A15, A16, A24, A27, A28) and produces the target's profile
  hydrostatically (B4 to B7).

The code is in `C:\Users\srafkin\SaturnProfile` (GitHub `ScotRafkin/SaturnProfile`):

| Location | What |
|---|---|
| `docs/specs` | the specifications, and `STATE.md` |
| `reports/` | reports and reviews |
| `tests/` | acceptance scripts |
| `docs/RUNBOOK.md` | the runbook |

## Roles and protocol

- **The reviewing agent** writes specifications and reviews reports.
- **The coding agent** reads each specification first and files a pre-execution reading. The
  reviewing agent rules on it in the specification.
- **Each step:** the coding agent codes it, runs its acceptance script and files a report. The
  reviewing agent reviews it.
- **Commits** happen only after review and the author's go.

Rules:

- One fix per specification, each revertible by git.
- Rulings go into the specification as numbered sections; they are not written over earlier text.
- The author asked (3 October) that a specification not be rewritten while the coding agent is
  working from it.

## Read first

1. `claude/CASSPIAN_Staggering_Design_Note_2026-10-03.md` (project). This is the problem and the
   proposed arrangement.
2. `claude/casspian_status_2026-10-03.md` (project). Where everything stands, and the author's
   standing rules.
3. `docs/specs/SPEC_06_Kernel_On_The_Isobar.md`, `reports/REPORT_06_step1.md` and
   `reports/REVIEW_06_step1.md`. The principle that the mesh never resamples the wind, and the
   check 7 law.
4. `docs/specs/SPEC_07_Wind_Interpolant.md` (§8 and §9) and `reports/REPORT_07_step1.md` (§4, §5,
   §7). The corner at a node, and the dense-node reference.
5. `SPEC_00` §7.3, rules (i) and (ii), and `SPEC_04` §4, "The staggering". What was asked for and
   what was built.
6. The code: `forward/transfer.py` (`trace`, `transfer`, `produce_at_target`, `outer_loop`),
   `forward/production.py`, `lib/hydrostatic.py`, `lib/kernel.py` (`ColumnIntegral`) and
   `forward/estimate.py`.

## The task

Turn the design note into a specification once the author has decided its §8 questions:

1. the arrangement;
2. the rule for the level temperature;
3. the estimate's form (A29 to A35, which reaches the manuscript);
4. the order;
5. the flag's tolerance.

Before the specification, have the coding agent confirm the note's §3 claim by measurement: that
the pressure identity is the disagreement between the tracing's layer integral of the shear and the
transfer's two level samples. Run 6 and a dense run show it directly.

## What not to do

- **No tolerance set for one case.** Tests are general and every experiment is a sample of them
  (author, 3 October). The bounds per run proposed at the end of SPEC_07 were withdrawn for this
  reason.
- **No smoothing of a hypothesis inside the model.** The model reads what it is given. Resolution is
  the hypothesis's own, and a level the model cannot resolve is flagged, not hidden.
- **No stability diagnostic or convective adjustment yet.** The author deferred both until the
  slope changes are dealt with.

## Practice

- Files written to the author's machine go through a fresh outputs filename each time and are read
  back.
- Project documents are written with the Projects tool.
- No dashes in prose.
- No AI attribution in commits.
- SwRI proprietary: nothing leaves the repository.

## The state you start from

- **`main` at `7652852`.**
  - SPEC_00 to SPEC_07 closed. SPEC_07 was not adopted; its record is at `5f2e05c`.
  - `step04_4` and `step04_5` pin their own outer-loop tolerance, 1e-8. `step04_5` run 1 asserts
    it instead, since that run goes through the namelist.
- **Regression of record** on the clean tree at `7652852`: 30 of 30 suites at their reference
  counts, 6923 s.
- **SPEC_07's L2 work** is on branch `spec07-l2-parked` (`89e3e6d`, pushed).
  - It includes `tests/step07_1/`. The reports on `main` name those scripts, which exist only on
    the branch.
  - The coarse L2 experiment products are no longer on disk; the branch rebuilds them. The dense
    comparison is kept in `reports/step07_1/dense/comparison.json` (ignored).
- **The author's answers to the design note's §8:** not yet given. Settle them with the author
  first.
