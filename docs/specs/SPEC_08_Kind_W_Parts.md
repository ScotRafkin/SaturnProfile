# SPEC 08. Kind W holds what the data hold

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.5, 7 October 2026. Author of record: S. Rafkin. **Status: accepted by the author at v0.2
on 7 October 2026, the decisions of §3 confirmed as drafted; v0.3 rules on the coding agent's
reading (§6), which amends §1: the reference wind and its level are required; v0.4 drops the
full regression (§6 ruling 15); v0.5: Step 1 accepted (`reports/REVIEW_08_step1.md`), the Appendix
applied as SPEC_00 v0.23; closed at the acceptance commit.**
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.22, the closed `SPEC_04_Transfer.md` v0.23,
`SPEC_05_Shear_Experiments.md` v0.7, `SPEC_06_Kernel_On_The_Isobar.md` v0.4 and
`SPEC_07_Wind_Interpolant.md` v0.6 (closed, not adopted), which it does not repeat. It starts from
`main` as pushed at `84c870a`.

---

## 0. Why, and how this specification is to be worked

**Why (author, 6 October 2026).** The next work builds a wind field for the Lindal epoch from
observations: Lindal et al. (1985) Fig. 10, which gives a temperature gradient in latitude at three
pressure levels and so a shear and nothing else, and the Sanchez-Lavega et al. (2000) table, which
gives a wind at one level. Each is to be written as a kind W file holding the data where they are
and nowhere else. Kind W cannot hold that today. It requires all three parts on every node, checks
that they add up, and requires both poles with zero wind on every read. The model needs none of
this from a data file: it reads only the total, and only from a file it runs on.

**The remedy.** Each of the three parts becomes optional, and a file holds whichever the data give,
on its own latitudes and pressures. The sum identity check is removed: it guards against a writer
error, not a physical one (author: an idiot check). The polar rule moves from every read to the
places where the model takes a wind, where it guards a real singularity. The model's loaders refuse
a wind file that lacks what the model reads, and name it.

**Isolated on purpose.** One step, its own commit, revertible by git. No model arithmetic changes,
so no product value can move.

**Protocol.** SPEC_04 §0 and SPEC_05 §0 apply: report and review, nothing committed before the
author's go, no AI attribution, no dashes, the regression rule of SPEC_04 v0.22 with the author's
direction of 1 October. Decision N: the code refuses schema faults, typos, extrapolation and named
numerical failures, and nothing else.

**The coding agent's reading.** Before execution the coding agent reads this draft and reports
(`reports/REPORT_08_preexecution.md`) on the points of §4. The reviewing agent rules on the reading
in a section of its own before the step starts.

---

## 1. Kind W, redefined

**What a kind W file is.** Zonal wind data on the file's own latitude and pressure grid, in up to
three parts along the local vertical:

| Variable | Required | Meaning |
|---|---|---|
| `u_total_ms(latitude, pressure)` | no | the wind; with `u_total_uncertainty_ms` beside it when present (the existing companion rule) |
| `u_reference_ms(latitude)` | no | the wind at `reference_level_pressure_Pa` |
| `u_shear_ms(latitude, pressure)` | no | the change of the wind from `reference_level_pressure_Pa` along the local vertical |
| `reference_level_pressure_Pa` | no | the level the reference wind is given at and the shear is measured from |
| `value_provenance(latitude, pressure)` | yes | unchanged |
| the coordinates | yes | unchanged |

The required global attributes are unchanged. A file holds the parts its source gives. A file with
the reference and the shear but no total has no total: no reader forms one. The grid is the
source's own; nothing requires the poles, a range of pressure or a number of levels. (SPEC_00 §6.6
already lets the reduction's wind have one pressure level.)

**What is no longer checked.**

- **The sum identity** `u_total = u_reference + u_shear`: removed from the schema and from the two
  forward loaders. A tool that writes all three forms them consistently, as every present tool
  does; nothing checks it.
- **No rule ties the parts to each other** (for example a shear without its reference level). Each
  reader that needs a part names it when it is absent (§1, the model; §2 deliverable 5, the shear
  tool).

**What the model reads, and where it refuses.** Three loaders take a wind for the model:
`lib.control.load_reduction_inputs` (the reduction), `lib.control.load_run_inputs` (closure mode)
and `lib.control._load_transfer_inputs` (transfer mode). Each refuses, naming the variable, a wind
file without:

1. `u_total_ms`, the field every path reads (`WindField`, `refrac.anchor.wind_of_latitude`);
2. `reference_level_pressure_Pa`, the level whose wind sets the anchor's geoid in the reduction,
   the closure production's wind at `phi_c` (`forward.production.produce` through
   `wind_of_latitude`), and the reference-level value the transfer production records.

Each loader then applies the polar rule as it stands: both poles nodes of the latitude grid, and
`u_total_ms` exactly zero there at every level. The forward loaders already call it; the reduction's
loader gains the call. The rule is removed from `lib.schema.validate`, so a data file that stops
short of the poles reads and writes. Coverage needs no rule of its own: the model refuses to
extrapolate (`WindField`).

**The name `u_reference_ms` is kept.** The author's "mean wind" is this part; renaming it would
touch every file and suite for no change of meaning.

---

## 2. Step 1

**Deliverables.**

1. **`lib.schema`.** In the kind W table, `u_total_ms`, `u_reference_ms`, `u_shear_ms` and
   `reference_level_pressure_Pa` are not required. `check_wind_components` is deleted with its
   calls. `check_wind_poles` is no longer called by `validate`; it stays in the module for the
   loaders. The kind's `notes` restated to §1.
2. **`lib.control`.** The three loaders of §1 refuse, by name, a wind without `u_total_ms` or
   `reference_level_pressure_Pa`, and apply the polar rule. The two forward loaders drop their
   `check_wind_components` call.
3. **Nothing else in `lib`, `refrac` or `forward` changes.** `WindField`, the reduction and the
   productions read the same variables as now, from files the loaders have admitted.
4. **F9** (`tools/plots/figures_profile.py`). The left panel draws the source line
   (`u_reference_ms`) only when the file carries it; the line of `u_total` at the reference level is
   drawn always, as now. A run's wind always carries the reference level (§1).
5. **The shear tool** (`tools/wind/shear.py`). It forms its output from what its input has:
   - `u_total_ms` is required, since every case reads `u_s` from it; an input without it is refused
     by name;
   - `u_reference_ms` and `reference_level_pressure_Pa` are carried when present;
   - `u_shear_ms = u_total - u_reference` is written only when the input carries `u_reference_ms`;
     otherwise the output has no shear part.
6. **`casspian-wind-from-curve`** is unchanged: it writes all three parts, consistently.
7. **SPEC_00 §6.6** amended by the reviewing agent at acceptance (Appendix).

**Acceptance (`tests/step08_1/accept_step08_1.py`).** Files are made in memory or under
`reports/step08_1/`, never over a committed file.

1. **Nothing moves.** The registered wind files (the reduction's `lindal_wind.nc`, the closure and
   transfer runs' source and run winds) read as before. The full regression passes, with the suites
   named below restated. Since no arithmetic changes, no registered file is rebuilt and there is no
   sweep.
2. **A file holds what the data hold.** Three kind W files are written and read back through
   `lib.io`, each on a grid that stops short of both poles:
   - shear only: `u_shear_ms` and `reference_level_pressure_Pa` on three pressure levels;
   - total only, on one pressure level;
   - reference only, with its level.
   Each reads back array-equal.
3. **The model takes only what it can run on.** For each of the three loaders, a copy of the
   registered wind with `u_total_ms` removed, and one with `reference_level_pressure_Pa` removed,
   is refused, the message naming the variable.
4. **The polar rule at the model.** A copy of the registered wind with the north pole set to
   1e-6 m/s reads through `lib.io`, and is refused by each of the three loaders.
5. **The shear tool.** On an input without `u_reference_ms`, the case of SPEC_05 Step 1 check 5
   writes a file with `u_total_ms` equal to the closed form to 1e-12 relative and no `u_shear_ms`;
   on an input without `u_total_ms` it is refused by name. The `step05_1` suite passes unchanged in
   its values.
6. **F9.** Drawn for the transfer run (both lines) and for a run whose wind has no `u_reference_ms`
   (the total line only); the author views both.

**Accepted suites whose checks change.**

- `step04_0` check 4 (a copy with the sum identity broken is refused): **retired**; the reference
  count goes from 15 to 14.
- `step7` check 2 (a copy with a nonzero pole wind is refused by `lib.io.read`): **restated**, the
  copy refused by `load_reduction_inputs`; its comment about the sum identity firing first removed.
- **Wording only, conditions unchanged:** `step7` check 6, `step05_1` checks 3 and 8, `step04_2`
  check 8 and the report lines of `step04_4` and `step04_5` that say the reader checks the sum
  identity. Each computes the identity itself from the file, which remains true of every file the
  present tools write.
- Any other suite the reading finds (§4) is named in the report with its old and new expectation.
  No accepted suite is edited silently or left failing.

**Regression.** `lib` changes: the full set, about 100 minutes, on the registered products as
committed.

---

## 3. Decisions (confirmed by the author, 7 October 2026)

1. **What runs the model:** `u_total_ms` and `reference_level_pressure_Pa` (§1). The reference level
   is read by the reduction, the closure production and the transfer production's record, so a file
   without it cannot run any present path. Making the transfer independent of it would be a change
   to `forward`, not wanted here.
2. **No rule ties the parts together** (§1). The alternative, refusing a shear without a reference
   level, is a check on the writer of the kind removed with the sum identity.
3. **`step04_0` check 4 retired**, rather than restated to show that a broken sum now reads.
4. **The shear tool requires the total** and writes a shear only when its input has a reference
   (§2 deliverable 5), the author's rule that what is undefined stays undefined.
5. **`u_reference_ms` keeps its name** (§1).

---

## 4. For the coding agent's reading

1. Every reader of `u_total_ms`, `u_reference_ms`, `u_shear_ms` and `reference_level_pressure_Pa`
   in `src`, confirming §1's list or adding to it.
2. Every accepted suite, diagnose script or fixture comparison that asserts the sum identity or the
   polar rule through `lib.io.read` or `validate`, beyond those named in §2.
3. Whether anything in `lib.io` (writing, or reading a kind `profile`'s embedded `inputs/wind`
   group) relies on `validate` running the two checks.
4. Whether any registered file or fixture could change under this step.

---

## 5. Revision history

- v0.1, 6 October 2026: first draft, from the author's direction of 6 October: kind W holds any of
  its three parts, on the data's own grid; the sum identity is not checked; a file without the
  total cannot run the model.
- v0.2, 7 October 2026: accepted by the author as drafted; the decisions of §3 confirmed. Next in
  order: SPEC_09 (digitizing Lindal et al. 1985 Fig. 10), SPEC_10 (the Lindal-only wind and a run
  with it), SPEC_11 (Sanchez-Lavega et al. 2000, the consistency test and the combined wind).
- v0.3, 7 October 2026: the coding agent's reading (`reports/REPORT_08_preexecution.md`) ruled on
  (§6). The author's ruling on what kind W is: the reference wind and its level are required, the
  total and the shear optional; data that are not a wind at a stated level are not kind W. The
  dimension check goes with the sum identity; figures draw the parts a file has.
- v0.4, 7 October 2026: no full regression for this step, and the standing rule for when the
  full set is run (§6 ruling 15, author).
- v0.5, 7 October 2026: Step 1 accepted; the Appendix, as §6 amends it, applied to SPEC_00 (v0.23);
  SPEC_08 closes at the acceptance commit.

---

## 6. Rulings on the coding agent's reading (7 October 2026)

`reports/REPORT_08_preexecution.md`, measured by `tests/step08_1/preexecution_measure.py` on a
simulated schema. A careful reading; every finding was checked against the code and is correct.
The author's ruling on what kind W is comes first, because it settles R1 and R2.

**What kind W is (author, 7 October 2026).** A kind W file is a wind at a stated level: it requires
`reference_level_pressure_Pa` and `u_reference_ms`. `u_total_ms` and `u_shear_ms` are optional.
Data that are not a wind at a stated level (a shear alone, a temperature gradient) are not kind W;
they stay in a format of the tool's choosing, such as a CSV with its note, as the digitized sources
in `data_static/` do now. This amends §1's table: `u_reference_ms` and
`reference_level_pressure_Pa` are required, enforced by the schema's required variables and by
nothing else. §1's "no rule ties the parts to each other" stands for the two optional parts.

An example, for the record. The Sanchez-Lavega et al. (2000) table is a wind at one level: as kind
W, its reported wind is `u_reference_ms` at its assigned level; where `u_shear_ms` is written it is
zero there, and where `u_total_ms` is written it equals the reference. It runs the model only once a
declared vertical structure fills the pressure grid, as `casspian-wind-from-curve` does now.

1. **R1** is settled by the ruling: every kind W file carries the level, so `WindField` and the
   shear tool need no change for it. The shear tool refuses an input without `u_total_ms`, by name
   (§2 deliverable 5).
2. **R2.** The shear tool writes `u_shear_ms = u_total - u_reference` only where its input carries
   `u_shear_ms` (values replaced, the variable never created); otherwise the output has no shear
   part. Every input carries the reference, so a shear without a reference cannot arise.
3. **R3. Dropped** (author): the dimension half of `check_wind_components` goes with the sum
   identity. Once the pipeline is automated and called in order, there is no one for it to catch.
4. **R4, accepted.** In all three loaders, directly after the rotation check: the presence of
   `u_total_ms`, then the polar rule. The level needs no loader check, since the schema now requires
   it on read. Check 3 restated: a copy without `u_total_ms` is refused by each loader, the message
   naming it; a copy without `reference_level_pressure_Pa` or `u_reference_ms` is refused on read
   by the schema, naming the variable.
5. **R5, accepted:** the loaders pointed at the copies through the loaded manifest and namelists
   (`dataclasses.replace`); the copies written with `xarray.to_netcdf`, attributes kept.
6. **R6, accepted:** F9 drawn from the registered transfer product held in memory with
   `u_reference_ms` dropped from its `inputs/wind` group.
7. **R7, accepted:** the restated `step7` check 2 perturbs `u_total_ms` only and hands the copy to
   `load_reduction_inputs`.
8. **R8, accepted.** Docstrings that become false are restated, text only. Deliverable 3 means code.
9. **R9, accepted:** the two `step04_3` docstrings restated, the `step04_1` and `step02_6` texts
   left; each named in the report.
10. **Figures draw what a file has (author).** Every figure that draws a kind W file draws the parts
    it carries, and omits a line or leaves a panel blank where a part is absent: F9, and
    `casspian-render` of a kind W file (`_figure_wind`), brought into this step from the reading's
    note. F9's own case is already covered by R6.
11. **Check 2 restated** to the ruling. Written and read back array-equal, on latitudes that stop
    short of both poles:
    - the reference wind at one level (the Sanchez-Lavega form);
    - the reference wind with the shear on three levels, no total;
    - the reference wind with the total on one level.
    A file with the shear and the level but no reference wind is refused on read, naming
    `u_reference_ms`.
12. **Check 6 extended:** `casspian-render` draws check 2's three files; the author views them with
    F9's two cases.
13. **Corrections to §1 and §2, from the reading.** The transfer production does not record the
    reference level; only `WindField` reads it, in its constructor. Decision 1 stands for that
    reason. The regression of record takes about 115 minutes, not 100.
14. **Noted, no change.** The pre-SPEC_04 winds under `tests/step04_0/fixtures/swept/` now read;
    no suite reads them through `lib.io`. Uncertainty companions for the reference and shear parts
    belong to SPEC_10 or SPEC_11.

15. **Regression (author, 7 October 2026).** No full regression for this step: it changes no
    computed value, so the full set would show nothing the step's own checks do not. Run
    `step08_1` and the two suites whose checks it edits, `step04_0` and `step7`, and nothing else;
    the wording-only edits of §2 need no run. This replaces §2's regression paragraph. **The
    standing rule from here:** the full set is run only when a change can alter a computed value
    (the model's arithmetic or numerics, or the registered products), or at a milestone the author
    names. Otherwise a step runs its own checks and the suites whose checks it edits.

**The Appendix is applied as amended here.** In its opening paragraph, "in up to three parts" reads
"a wind at a stated level, with the reference wind and its level required and the total and the
shear optional"; in its table, `u_reference_ms` and `reference_level_pressure_Pa` are required; in
its load-time rules, the second reads "The file carries `u_total_ms`; otherwise refuse, naming it."

---

## Appendix. Amendments to SPEC_00, to be applied at acceptance

**Section 6.6 (kind W).** The opening paragraph: "Zonal wind data on the file's own latitude and
pressure grid, in up to three parts along the local vertical: the wind, the reference level wind,
and the change of the wind from the reference level. A file holds the parts its source gives, where
the source gives them. The model runs only on a file that carries the total and its reference
level."

The variable table: `u_total_ms`, `u_reference_ms`, `u_shear_ms` and `reference_level_pressure_Pa`
marked optional; `u_shear_ms` defined as "the change of the wind from `reference_level_pressure_Pa`
along the local vertical"; the clause "zero everywhere when the hypothesis is altitude independent"
kept.

"Rules the model enforces at load time" replaced by:

- The file's `rotation_system_name` and rate match the kind R file; otherwise refuse.
- The file carries `u_total_ms` and `reference_level_pressure_Pa`; otherwise refuse, naming the
  missing variable.
- Both poles are nodes of the latitude coordinate and `u_total_ms` is exactly zero there at every
  level; otherwise refuse. The rule is applied where the model takes a wind, not on every read.
- The model does not extrapolate the wind; a point outside the file's grid is refused.
- The prose provenance is required and is not read by the model (unchanged).

The sum identity rule is deleted. Section 6.6's paragraph on the reduction instance is kept.

**Revision table:** v0.23, the amendments above, at SPEC_08's acceptance.
