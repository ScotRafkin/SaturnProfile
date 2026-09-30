# SPEC 05. The shear tool and the first shear experiments

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.3, 29 September 2026. Author of record: S. Rafkin. **Status: accepted by the author at
v0.3 on 29 September 2026. Step 0 closes with the review of `REPORT_05_step0`; Step 1 begins after the coding
agent's pre-execution review is ruled on.** Built from the design
note `claude/SPEC_05_shear_design_2026-09-29.md` and its review
`claude/SPEC_05_shear_design_review_2026-09-29.md` (the review is the later word where they
differ), with the author's rulings of 29 September. Depends on `SPEC_00_Architecture_and_Data_Files.md`
v0.22, the closed `SPEC_04_Transfer.md` v0.23, and `docs/RUNBOOK.md` v0.3, which it does not
repeat. Composition experiments (helium, NH3) are not in this version; the author will add them
after the shear runs are read, because what these runs show may shape that part.

---

## 0. How this specification is to be worked

**Protocol.** SPEC_04 §0 applies unchanged: one step at a time; the status line and `STATE.md`
read before any step; `REPORT_05_step<N>.md` by the coding agent and `REVIEW_05_step<N>.md` by
the reviewing agent; nothing committed until the review says accepted and the author says go;
no attribution to an AI in any commit; no dashes; American spellings; no silent choices; every
measured number quoted with its bound. Acceptance scripts are committed under
`tests/step05_<N>/` and write only under the ignored `reports/step05_<N>/` (SPEC_00 v0.22 rows
[Y] and [Z]); each accepted suite is added to the driver table of `tests/run_regression.sh`
with its reference count.

**Regression rule.** SPEC_04 §0 at v0.22: a change reruns what it can reach. The report names
the suites rerun and, in one sentence, why the others were not.

**What is checked (decision N).** The code refuses only schema faults, namelist or control-file
typos (an unknown case, a missing or unused parameter, a malformed value), extrapolation, and
named numerical failures. Everything else is recorded. This specification adds no other check.

**Step 0 is already done.** In this project Step 0 has always been the groundwork that gives the
new specification a clean base (SPEC_03 and SPEC_04 Step 0 both rebuilt what came before). For
SPEC_05 that groundwork is the clone on the author's Linux machine, the runbook, the move of the
acceptance scripts to `tests/`, and the regression of record, reported in
`reports/REPORT_05_step0.md`. **The work of this specification starts at Step 1.** Do not number
any shear work as Step 0.

**What this specification covers.** A tool that writes a kind W file under a named shear
hypothesis (Step 1); its place in every forward run's build, and a tool that makes a new named
run (Step 2); a first set of named transfer experiments under those hypotheses (Step 3). **The
model does not change.** What is left for later: the constant-on-cylinders case (§1.6), a
latitude-confined case, decay above a height or radius, a modify mode, the in-memory call from
the Monte Carlo driver (§1.5), and the composition tools and experiments.

**Two meanings of "profile" kept apart.** The *wind grid* is the latitude and pressure grid of
whatever kind W file the shear tool is given, however that file was made. The *anchor column*
is the range of levels of whatever anchor a run transfers from. Neither is fixed by this
specification, and nothing in the code assumes a particular source for either. Every pressure
bound in the shear tool is a bound on the wind grid of its input, never on any anchor column.
For the Lindal test case used in the acceptance and the experiments below, the wind grid runs
from 1 Pa to 1 MPa at ten levels per decade (the `pressure_grid_Pa` of the Lindal runs' wind
build section), and the anchor column from about 20 Pa to 1.294 bar (129,422 Pa), with six of Lindal's levels
below 1 bar.

---

## 1. The design, as settled with the author

### 1.1 What the model reads

Checked in the repository. The model reads only `u_total_ms`: the reduction, `lib.windfield`,
the kernel, the tracing and the production all read it and nothing else. `u_reference_ms` and
`u_shear_ms` are read in two places only: the sum-identity check (`lib/schema.py`,
`check_wind_components`, which refuses a file whose parts do not add up to 1e-12 of the field's
size) and F9's left panel. `reference_level_pressure_Pa` is read by the reduction
(`refrac/anchor.py`, the geoid wind); `WindField.reference_wind` exists but no forward code calls
it; the forward run's reference surface uses `u_total` at the gauge isobar. `value_provenance`
is read by nothing in the model. The model interpolates `u_total` linearly in latitude
(decision L) and linearly in ln p, and refuses outside the file's grid.

Every shear hypothesis is therefore a kind W file written by a tool. No model code changes.

### 1.2 The three parts of kind W under a shear case

- **Source data, carried unchanged.** The input file's `reference_level_pressure_Pa` and
  `u_reference_ms(phi)` describe where the input came from (for Lindal, the cloud-tracked wind
  assigned to 1 bar; for an invented input, whatever its maker assigned). The tool copies them to
  the output and never uses them.
- **The case's wind.** The case builds `u_total_ms(phi, p)` from the shear reference pressure
  `p_s` and `u_s(phi)` alone (§1.3). `p_s` is a control-file value, independent of the input's
  reference pressure; it is not written to the kind W file.
- **The shear, formed by the tool.** `u_shear_ms = u_total_ms - u_reference_ms`, formed in the
  one assembly function every case passes through, so the sum identity holds by construction in
  memory and on disk. `u_shear` need not be zero at the reference pressure: a nonzero value there
  is how far the hypothesis departs from the source at the source's own level, and a reader
  should see it.

A consequence stated here so nobody infers otherwise: a wind file used by a reduction is read by
the reduction at its reference pressure as the source's wind, so a shear case there would move
the anchor itself. Whatever tool or script makes a reduction's wind file (for the Lindal anchor
it happens to be `casspian-lindal-inputs`), it is made without a shear case. Shear cases belong
to the wind files of forward runs, where the shear step is part of every build (§3).

### 1.3 Reading u_s

`u_s(phi)` is the input's `u_total_ms` interpolated to `p_s` at every latitude node, linearly in
ln p between the two bracketing pressure nodes (the model's own vertical rule). `p_s` need not be
a node of the wind grid, and the tool does not insert nodes (author, 29 September). `p_s` must
lie within the wind grid, endpoints included, or the tool stops and names the grid's range.

The output is written on the input's latitude and pressure grids, unchanged. A kink of the
hypothesis that falls between nodes is smeared over one grid cell, because that is what the
hypothesis is on its grid; a kink at a node is represented exactly, since the model's rule is
linear in ln p. A shape that is curved in ln p (linear in p, §1.4) is represented to the wind
grid's resolution, which is set by `pressure_grid_Pa` in the wind build section. Run 7 measures
that (§4).

### 1.4 The cases

The control file names a case; a table in the tool maps the name to its case function. Each case
function takes `u_s`, the grids and its own parameters and returns `u_total(phi, p)`. Case
functions are built from small shared parts that know nothing about any case. The author's case
numbers are kept as labels in documents; the control file uses the names.

| No. | Name | Parameters | `u_total(phi, p)` |
|---|---|---|---|
| -1 | `identity` | none | the input file returned unchanged |
| 0 | `uniform` | `shear_reference_pressure_Pa`, `scale` (c, default 1) | `c u_s(phi)` at every p |
| 1 | reserved | | constant on cylinders; not in this version (§1.6) |
| 2 | `decay_above` | `shear_reference_pressure_Pa`, `shape`, `stop_pressure_Pa`, `stop_fraction` | `u_s(phi) F(p)`, change above `p_s` |
| 3 | `increase_below` | the same four | `u_s(phi) F(p)`, change below `p_s` |

**Cases 2 and 3, one definition.** Let `p_stop` be the stop pressure and `f` the stop fraction.
Define the position `x(p)` between `p_s` and `p_stop` by the shape:

- `linear_ln_p`: `x = ln(p / p_s) / ln(p_stop / p_s)`; requires `p_stop > 0`.
- `linear_p`: `x = (p - p_s) / (p_stop - p_s)`; `p_stop = 0` is allowed.

Then `F = 1` for `x <= 0` (the side of `p_s` the case does not change), `F = 1 + (f - 1) x` for
`0 < x < 1`, and `F = f` for `x >= 1` (held beyond `p_stop`). The change is a fraction of
`u_s(phi)`, so the jet's latitude structure scales; a uniform subtraction in m/s would leave the
latitude derivative the kernel reads nearly unchanged.

`decay_above` requires `p_stop < p_s`; `increase_below` requires `p_stop > p_s`. These are the
cases' definitions, checked as malformed parameters. `f` is not restricted: `decay_above` with
`f > 1` is a wind that grows with height, and nothing refuses it. `p_stop` need not lie within
the wind grid; it only sets the shape (in the Lindal test case, case 3's nominal 10 bar happens
to sit exactly on the wind grid's bottom edge).
`linear_ln_p` with `p_stop = 0` is refused with a message naming `linear_p`.

A parameter the named case does not use is refused as a typo (author, 29 September: keys only
where code reads them).

**Poles.** The input is exactly zero at both poles (the schema's `check_wind_poles`), so
`u_s`, and with it every case's `u_total`, is exactly zero there. The acceptance checks it.

### 1.5 Structure of the tool, and the Monte Carlo later

One module, `src/casspian/tools/wind/shear.py`, console entry `casspian-wind-shear`, in three
layers: the parts, the case functions with their table, and the assembly and write. The public
functions are:

- `construct(source, case, parameters) -> xarray.Dataset`: takes an in-memory kind W dataset and
  returns the complete output dataset, three parts, attributes and provenance included. No file
  is read or written inside it.
- `build(control_path, section="shear") -> Path`: reads the control section and the source file,
  calls `construct`, writes with `lib.io.write`. This is what the command line and the run build
  call.

The split is a design constraint for later, not a feature of this version: the Monte Carlo
driver will import `construct`, skip the write, and hand the returned object to the model. The
switch that does so and the model's entry that accepts a kind W object belong to the Monte Carlo
specification. Nothing in this version adds either.

### 1.6 The cylinder case is left out

Constant on cylinders (case 1) is not in this version (review, Revision 1). Knowing where a
cylinder meets the `p_s` surface needs the model's reference surface, therefore an anchor, a
fixed point between that surface and the wind being built, a per-hemisphere construction, an
equatorial gap rule, and at M of 2 or more a choice of anchor. The test the cylinder wind
performs is already in the regression (`tests/step04_2` builds it; `tests/step04_4` and
`tests/step04_5` check the transfer against decision Q's floor). When deep-wind hypotheses call
for it, the specification that brings it back states: the anchor that sets the geometry, the
fixed point, per-hemisphere U(s), the pin at `p_s` (at the input's reference pressure it
reproduces decision P), and the equatorial gap rule (the author's ruling of 29 September: fill
with the wind on the outermost cylinder that still meets the `p_s` surface, and flag the cells).

### 1.7 What the output file carries besides the wind

- **`u_total_uncertainty_ms`**: set by the optional control key `uncertainty_ms` (m/s),
  accepted by every case except `identity`, which returns its input untouched. With no value,
  the input's field passes through unchanged, its `long_name` stating that it is the source's
  uncertainty; it is never scaled with the wind (author, 29 September: a smaller wind is not a
  better known one). With a value, that value replaces the input's at every node, the
  `long_name` saying it was set in the shear control file. The model does not read this field
  until the combination specification.
- **`value_provenance`**: for every case except `identity`, every cell takes a new value 5,
  meaning `constructed_by_shear_case`, appended to the file's `flag_values` and `flag_meanings`
  (they are written per file; the schema does not change). The author accepts this as metadata to
  be revisited in a later audit.
- **Globals**: `vertical_structure` becomes `shear case <name>`; every other global of the input
  is carried unchanged (`observation_level_*`, `method`, `epoch`, the season, the rotation).
  `input_hashes` lists the source wind file and the build control file, so the file traces to the
  exact case and parameters that made it without separate labels. One `history` line records the
  case name and its parameters.
- **`identity`** copies the input dataset and adds only its `history` line and its own
  `input_hashes`. The closure comparison (`lib.control.CLOSURE_DROPPED_ATTRIBUTES`) drops both,
  so a closure run's identity-sheared wind still compares content-identical to the anchor's
  embedded copy.

---

## 2. Step 1: the shear tool

**Deliverables.**

1. `src/casspian/tools/wind/shear.py` as §1.4, §1.5 and §1.7 describe, with parts tested on their
   own: the ln p interpolation to `p_s`, the two shapes of `x(p)`, and the ramp `F`.
2. Console entry `casspian-wind-shear` in `pyproject.toml`, taking the control file and
   `--section` (default `shear`), as `casspian-wind-from-curve` does.
3. The `[shear]` control section's keys: `role`, `prefix`, `source` (path to the input kind W),
   `output`, `case`, `title` (optional), the case's own parameters of §1.4, and the optional
   `uncertainty_ms` of §1.7. Paths resolve as in every other build section.

**Acceptance (`tests/step05_1/accept_step05_1.py`), on the transfer run's own wind file
`forward/lindal_transfer/inputs/lindal_transfer_wind.nc` as input.**

1. `identity`: every variable array-equal to the input; the closure comparison's content check
   passes against the input.
2. `uniform`, `c = 1`, `p_s = 1e5`: `u_total` array-equal to the input's (the input is altitude
   independent); `u_reference` and the reference pressure array-equal to the input's.
3. `uniform`, `c = 0`: `u_total` exactly zero; `u_shear = -u_reference` exactly; the file reads
   back through the schema (sum identity).
4. `uniform`, `c = 1`, `p_s` between two nodes (for example 1.1e5 on a synthetic input sheared in
   pressure made in the script): `u_s` equals the ln p interpolant at every latitude to 1e-12
   relative.
5. `decay_above`, `linear_ln_p`, `p_s = 1e5`, `p_stop = 700`, `f = 0`: at every node, `u_total`
   equals the closed form `u_s F(p)` to 1e-12 relative; `F` is 1 at and below `p_s`, 0 at and
   above `p_stop`.
6. `decay_above`, `linear_p`, `p_stop = 0`, `f = 0`: `u_total = u_s p / p_s` above `p_s` to 1e-12
   relative.
7. `increase_below`, `linear_ln_p`, `p_stop = 1e6`, `f = 1.5`: closed form to 1e-12 relative;
   `u_total = 1.5 u_s` at 1e6.
8. Every output of checks 2 to 7: sum identity through the schema on read; `u_total` exactly zero
   at both poles; `value_provenance` 5 everywhere; source data array-equal to the input's.
9. `construct` on the in-memory input returns a dataset whose every variable is array-equal to the
   file `build` writes and reads back, for the case of check 5.
9a. In every output of checks 2 to 7, `u_total_uncertainty_ms` is array-equal to the input's;
    with `uncertainty_ms = 25` on the case of check 5 it is 25 at every node.
10. The named refusals, and only those: unknown case; missing parameter; a parameter the case
    does not use; `p_s` outside the wind grid; `p_stop` on the wrong side of `p_s`; `linear_ln_p`
    with `p_stop = 0`.

**Regression.** A new module that no existing code imports: its own suite only.

---

## 3. Step 2: the shear step in every forward run, and the new-run tool

**The author's ruling (29 September).** The shear step is part of every forward run's build, so
the pipeline is always the same; a run that wants the source wind as it is sets `case =
"identity"`.

**Deliverables.**

1. **`casspian-run-inputs` runs five sections**, in the order gravity, rotation, wind, shear,
   composition. `[shear]` is required, with the same checks the driver already makes (role
   `forward`, the run's prefix, output under the run's `inputs/` carrying the prefix).
2. **File names.** `[wind]` writes the source wind as `inputs/<run>_wind_source.nc`; `[shear]`
   reads it and writes `inputs/<run>_wind.nc`, which is the path every namelist's `[inputs] wind`
   already names. No namelist changes.
3. **The two existing runs migrated.** `forward/lindal_closure/lindal_closure_build.toml` and
   `forward/lindal_transfer/lindal_transfer_build.toml` gain `[shear]` with `case = "identity"`,
   and their `[wind] output` becomes `inputs/<run>_wind_source.nc`. Their products must not move:
   the closure product equal to the registered one by the bit-identity check `step03_4` already
   makes, and the transfer product's values equal to the accepted Step 5 product's.
4. **`casspian-new-run`** (`src/casspian/tools/run/new_run.py`): `casspian-new-run <name> --from
   <run directory>` makes `forward/<name>/` with `<name>_build.toml` and `<name>.toml` copied from
   the named run, every prefix, output name, `[run] name` and title carrying the new name, and
   nothing else changed. It refuses if `forward/<name>/` exists. It makes no inputs; the person
   edits `[shear]`, `[target]` or whatever the experiment changes, then runs
   `casspian-run-inputs` and `casspian-forward` as the runbook says. Its control files are
   committed; its `inputs/` and `output/` are ignored like every run's.
5. **F9's left panel.** The label `u_reference at <p> mbar` becomes `source wind, assigned to <p>
   mbar`, and a second line, `u_total at <p> mbar`, is always drawn beside it (author, 29
   September: where the two agree they lie on top of each other).
6. **Runbook.** `docs/RUNBOOK.md` names the five sections and `casspian-new-run`, v0.4.

**Acceptance (`tests/step05_2/accept_step05_2.py`).**

1. A build file without `[shear]` is refused, naming the five sections.
2. `casspian-run-inputs` on both migrated build files writes both wind files; the `_wind.nc` file
   is array-equal to `_wind_source.nc` in every variable.
3. The closure run's content comparison against the anchor's embedded copies passes.
4. The closure product equals the registered product under the existing bit-identity check.
5. The transfer product's delivered N, p, T, altitude, `r0` and pressure identity equal the
   accepted Step 5 product's (`array_equal`).
6. `casspian-new-run` from `lindal_transfer` makes a directory whose two files differ from the
   source only in the name fields; a second call with the same name is refused; the new run
   builds and runs to the same values as check 5.
7. F9's left panel carries both lines, for the migrated transfer run (coincident) and for a
   `uniform, c = 0` run (the total on zero).

**Regression.** `tools/run` and the build files change how every forward input is made:
rerun every suite that builds forward inputs or reads a forward product (the coding agent names
them from the driver table); the F9 change reruns the figure checks of SPEC_04 §0. The
reduction chain is not reached.

**Accepted suites whose expectations change.** Some accepted suites may assert the contents or
file names of a run's `inputs/` directory, or compare a rebuilt input set against a fixture made
before this step (`step03_3` does the latter). A suite whose expectation changes only because of
this migration is updated in this step, each one named in the report with its old and new
expectation, and its reference count in the driver table kept unless a check is added or
removed. No accepted suite is edited silently, and none is left failing.

---

## 4. Step 3: the named experiments

**Runs.** Each is a named run made by `casspian-new-run` from `lindal_transfer`, differing only
in its `[shear]` section and, for run 3b, its target. Target 10 N. `p_s = 1e5` Pa for every case
but `identity`. The input to the shear step is the run's own `_wind_source.nc`.

| Run | Directory | Case | Parameters | Purpose |
|---|---|---|---|---|
| 2 | `shear_r2_uniform` | `uniform` | `c = 1` | reproduces the accepted SPEC_04 transfer |
| 3a | `shear_r3a_nowind` | `uniform` | `c = 0` | exact null at 10 N |
| 3b | `shear_r3b_nowind_anchor` | `uniform` | `c = 0`, target at the anchor's `phi_c` | the null's comparison |
| 3c | `shear_r3c_half` | `uniform` | `c = 0.5` | linearity of the latitude term (optional) |
| 4 | `shear_r4_decay12` | `decay_above` | `linear_ln_p`, `p_stop = 20` Pa, `f = 0` | about 12 percent per scale height |
| 5 | `shear_r5_decay20` | `decay_above` | `linear_ln_p`, `p_stop = 700` Pa, `f = 0` | about 20 percent per scale height |
| 6 | `shear_r6_decay40` | `decay_above` | `linear_ln_p`, `p_stop = 8000` Pa, `f = 0` | about 40 percent per scale height |
| 7 | `shear_r7_decay_linp` | `decay_above` | `linear_p`, `p_stop = 0`, `f = 0` | stress test |
| 7f | `shear_r7f_decay_linp_fine` | as run 7 | wind grid 20 levels per decade | grid representation of a shape curved in ln p |
| 8 | `shear_r8_increase25` | `increase_below` | `linear_ln_p`, `p_stop = 1e6`, `f = 1.25` | increase below |
| 9 | `shear_r9_increase50` | `increase_below` | `linear_ln_p`, `p_stop = 1e6`, `f = 1.5` | Galileo-like increase |

(Run 1, `identity`, is Step 1's check 1 and Step 2's migration; it is not repeated here.) Rates
check: `ln(1e5/20) = 8.52`, `ln(1e5/700) = 4.96`, `ln(1e5/8000) = 2.53` scale heights, so 11.7,
20.2 and 39.6 percent of `u_s` per scale height. Run 3b's target is the anchor's own `phi_c` to
full precision, the value the SPEC_04 Step 5 acceptance used. Run 7f differs from run 7 only in
`[wind] pressure_grid_Pa` (1 Pa to 1 MPa, twenty per decade); the shear tool writes on its
input's grid.

**Every run from 2 to 9 is made at geopotential spacings 5e4 and 2.5e4** (the namelist's
`[grid] geopotential_spacing_m2s2`), and reports the pressure identity at both, the number of
outer-loop passes, and the wall time. Bounds on the identity are set from these measurements;
nothing is loosened first. For scale: the Step 5 synthetic sheared run gave 3.3e-4 at 5e4 and
1.5e-4 halved, against 5.8e-7 for the closure run. Decision Q's truncation floor (latitude
roughness of the source curve) applies to every run and is reported, not bounded.

**Expected outcomes.**

- **Run 2** reproduces the accepted transfer product: N, p, T, altitude, `r0` array-equal.
- **Runs 3a and 3b, the exact null.** With zero wind the kernel is exactly zero (no latitude
  term, no vertical term, and no composition term for the uniform composition), isobars are flat
  in geopotential, and N is carried along them unchanged. The delivered N, p and T of 3a and 3b
  are equal to round-off, targeted at 1e-12 relative; if the measurement is larger, report it and
  ask for a ruling. The reference surface is the no-wind surface: `r0(10 N) = 60,092,307.69` m
  (the reviewing agent's independent march, DOP853 at relative tolerance 1e-13, the same code
  that returns the accepted closure-wind value 60,128,612.966 m to the millimeter) against
  60,128,612.97 m under the closure wind.
  Altitudes differ between 3a and 3b because gravity does; that is geometry, not a failure.
- **Run 3c.** The latitude term is `2 Omega_abs cos(phi) du/dphi` with `Omega_abs` containing
  `u`, so halving the wind gives a temperature change about 1.5 percent short of half of run 2's
  at 10 N (`u / (r cos phi)` is about 3 percent of `Omega` there). Measured and reported, not
  bounded.
**The transfer is local in the vertical.** Along each isobar the change in ln T is set by the
kernel at that level, which depends on the wind and its shear at that level. So each case changes
the delivered temperature only at the levels where it changes the wind.

- **Runs 4 to 7.** The change sits between `p_s` and `p_stop` and nowhere else. Above `p_stop`
  the wind is zero, the kernel is zero, and the delivered temperature equals the anchor's;
  relative to run 2 those levels are therefore warmer by run 2's own change (about 1.1 percent).
  Below `p_s` only the latitude term acts, as in run 2. Between the two, the temperature at 10 N
  falls relative to run 2, the more so the faster the decay: the reviewing agent's rough scaling
  for run 5 is about 10 percent between the anchor's latitude and 10 N, ten times the closure
  run's, colder toward the equator. For run 5 the gauge (100 mbar) lies inside the decay zone and
  the column top (0.2 mbar) above it. Every decay case also weakens the wind at the gauge, so the
  reference surface moves: the delivered radius at 10 N changes by up to the wind's 36 km share of
  the flattening, scaled by the decay at the gauge (about half of it for run 5).
- **Run 7 is a stress test.** Linear decay in p to zero is 100 percent per scale height at 1 bar,
  so the largest change is just above 1 bar, several tens of percent by the rough scaling; below
  1 bar the case changes nothing. The transfer is not linearized and may return an answer, or may
  stop with a named failure (isobars crossing, the outer loop not converging). Either is a result;
  the report records which, and nobody tunes the numerics around it. Run 7f then says whether the
  ten per decade grid's representation of the curved shape matters.
- **Runs 8 and 9.** At the levels they touch the local rate is about 11 percent (run 8) and 22
  percent (run 9) of `u_s` per scale height, the wind decreasing upward as in runs 4 and 5, so the
  change there is of the same order as run 5's and colder toward the equator. What is small is the
  number of levels: in the Lindal column only the six levels below 1 bar are touched. Their role
  now is to show the case works; deeper anchors will give them leverage.
- **Independent values.** The reviewing agent computes runs 5 and 9 independently (this needs
  the vertical-shear term added to its own transfer code first) and states them in a later version of this
  specification, before Step 3 is reviewed, at several levels inside and outside each case's
  zone rather than as a single number.

**Deliverables.**

1. The eleven run directories' control files, committed.
2. `tests/step05_3/run_experiments.py`: builds and runs every experiment at both spacings (7f at
   5e4 only), writing under `reports/step05_3/`; and `tests/step05_3/accept_step05_3.py`: the
   checks below.
3. F7, F8 and F9 for every run, as `casspian-forward` already draws them.
4. One comparison figure in the report, made by the Step 3 script under `reports/figures/` and
   not added to the package: temperature at 10 N minus run 2's against pressure for runs 3a to 9,
   one line per run, and a table of `r0(10 N)` per run.
5. `reports/REPORT_05_step3.md` with a results table: per run and spacing, the largest delivered
   T change from run 2 and the pressure where it occurs, and the change at 1 bar, at the gauge and
   at the top of the column; altitudes at the same levels; `r0(10 N)`; the pressure identity; passes; wall time; and for run 7 whether it
   completed or which named failure stopped it.

**Acceptance checks.** 1: run 2 array-equal to the accepted transfer product. 2: runs 3a and 3b
equal in N, p, T to the measured round-off. 3: run 3a's `r0(10 N)` equal to 60,092,307.69 m within
0.1 m. 4: the pressure identity of every
completed run at both spacings, bounded from the measurements. 5: runs 5 and 9 against the
reviewing agent's independent values (stated before Step 3 is reviewed). Every other outcome is reported, not checked.

**Regression.** Step 3 adds run directories and scripts and changes no code: its own suite only.

---

## 5. Decisions of v0.1, ruled by the author (29 September 2026)

1. Case names `identity`, `uniform`, `decay_above`, `increase_below` in the control file, the
   author's numbers kept as labels in documents, number 1 reserved: accepted.
2. Key names `scale`, `shear_reference_pressure_Pa`, `stop_pressure_Pa`, `stop_fraction`,
   `shape`: accepted.
3. `f` unrestricted in cases 2 and 3, only the side of `p_stop` checked: accepted.
4. `inputs/<run>_wind_source.nc` and `inputs/<run>_wind.nc`: accepted.
5. Experiment directories `shear_r<N>_<tag>` under `forward/`: accepted.
6. F9's second line: drawn always, not only when the two differ (author's change).
7. The 1e-12 relative target for the null of runs 3a and 3b, replaced by the measured value with a
   ruling if larger: accepted.

The uncertainty (§1.7, author, 29 September): carried unchanged and never scaled with the wind,
unless the optional control key `uncertainty_ms` gives a value, which then replaces it.

---

## 6. Revision history

- v0.1, 29 September 2026: first draft, from the design note and its review.
- v0.2, 29 September 2026: the author's review. §0 and §1.2 made general (no source of the wind
  or the anchor assumed; the Lindal values stated as the test case's); the uncertainty carried
  unchanged and not scaled, or replaced by the optional key `uncertainty_ms` (§1.7); F9 draws both lines always; §5 ruled.
- v0.3, 29 September 2026: the reviewing agent's review of v0.2
  (`claude/SPEC_05_draft_v0.2_review_2026-09-29.md`). Expected outcomes of runs 4 to 9 restated
  by level (the transfer is local in the vertical; runs 8 and 9 are the size of run 5's change on
  six levels, not a few tenths of a K); the report table by the largest change and where it
  occurs; the null's radius stated (60,092,307.69 m, bound 0.1 m); the anchor column's bottom
  1.294 bar; accepted suites whose expectations change with the Step 2 migration updated in the
  step and named.

---

## Appendix. Amendments to SPEC_00, to be applied at acceptance

- **Section 2.3 (run build files).** A forward run's build file carries five sections: gravity,
  rotation, wind, shear, composition. `[wind]` writes `inputs/<run>_wind_source.nc`; `[shear]`
  writes `inputs/<run>_wind.nc`, the file the namelist reads.
- **Section 6.6 (kind W).** `u_reference_ms` and `reference_level_pressure_Pa` are the source's
  wind and the level the source assigned it to, carried unchanged by the shear tool; `u_total` at
  that pressure may differ from `u_reference`, and `u_shear` is their difference. The sum identity
  is unchanged. A file written by the shear tool may carry `value_provenance` value 5,
  `constructed_by_shear_case`.
- **Console tools.** `casspian-wind-shear` and `casspian-new-run` added to the list of tools.
