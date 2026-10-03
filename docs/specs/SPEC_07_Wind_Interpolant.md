# SPEC 07. A wind interpolant with a continuous slope

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.6, 3 October 2026. Author of record: S. Rafkin. **Status: closed, not adopted (§9).**
The decisions of §5 were confirmed by the author at v0.2; v0.3 ruled on the coding agent's reading
(§7); v0.4 on the interim report of Step 1 (§8); v0.5 closes the specification without adopting
decision L2, the work kept on a branch and the findings carried to the redesign of the vertical. Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.22, the closed
`SPEC_04_Transfer.md` v0.23, `SPEC_05_Shear_Experiments.md` v0.7 and `SPEC_06_Kernel_On_The_Isobar.md`
v0.4, which it does not repeat. It starts after SPEC_06 Step 1 is committed, its run files swept and
`step05_2`, `step05_3` and `step05_4` rerun (REVIEW_06_step1, order of work).

---

## 0. Why, and how this specification is to be worked

**Why (author, 2 and 3 October 2026).** SPEC_06 made the model read the wind only from the wind
file's interpolant and integrate it with the file's nodes as breakpoints. What remains is the
interpolant itself. Decision L reads the wind linearly between nodes, in latitude and in `ln p`, so
its first derivatives are constant in each cell and jump at every node. By the thermal wind a jump
in `(du/dln p)` is a jump in the horizontal temperature gradient, and integrated in latitude it
becomes a step in the delivered temperature at the node's pressure. That is the stepping in run 7f's
lowest 50 km, and through the check 7 law of SPEC_06 (identity at a break equal to the anchor
layer's mass times the jump in `ln T` times `s - s0`) it is also what sets the remaining pressure
identity. Along an isobar the jump in `(du/dphi)` at every 0.5 deg latitude node is the staircase
the author saw in F8 (SPEC_06 §5 ruling 3). Neither is numerical instability: both are the
derivative of a piecewise linear hypothesis, integrated exactly. A real or invented wind can change
slope anywhere (the standing rule, SPEC_06 §6), and the model should not turn every node of the file
into a temperature step.

**The remedy.** Decision L revisited: a shape preserving cubic (PCHIP) between the same nodes, whose
first derivatives are continuous, with no overshoot. The node values do not change and nothing in
any wind file changes. The model and the reduction change together, because the reduction's wind
(`refrac.anchor.wind_of_latitude`) must be the model's (the at-anchor identity of SPEC_04 Step 5).

**What is not here.** Denser anchor levels for the production and the transfer (SPEC_08, the
anchor sub-levels). The superadiabatic profiles of the extreme cases (a stability diagnostic, and
whether to adjust): later, after the slope changes are corrected (author, 3 October). The
composition's interpolation, which stays linear.

**Isolated on purpose.** One step, its own commit, revertible by git. No switch between the old
and new interpolant in the code (as SPEC_06 §2 decision 2).

**Protocol.** SPEC_04 §0 and SPEC_05 §0 apply: report and review, nothing committed before the
author's go, no AI attribution, no dashes, the regression rule of SPEC_04 v0.22 with the author's
direction of 1 October. Tolerances are physical, not numerical: a fraction of a kelvin is noise.
Machine epsilon is asked for only where the mathematics leaves a value unchanged.

**The coding agent's reading.** Before execution the coding agent reads this draft and reports
(`reports/REPORT_07_preexecution.md`) on the points of §4, with measurements where they are asked
for. The reviewing agent rules on the reading in a section of its own before the step starts.

---

## 1. The interpolant (decision L2)

**Construction.** On the file's grid of latitude nodes `phi_j` and pressure nodes `p_k`, with
`x = ln p`:

1. **Latitude first, at every pressure node.** Each row `u(phi_j, p_k)`, `k` fixed, is interpolated
   in latitude by the one dimensional PCHIP of Fritsch and Carlson in the form
   `scipy.interpolate.PchipInterpolator` uses (node slopes the weighted harmonic mean of the
   adjacent secant slopes, `w1 = 2h_k + h_(k-1)`, `w2 = h_k + 2h_(k-1)`, zero where the secants
   differ in sign or either is zero, and scipy's three point rule at the ends). The rows are fixed
   data, so their coefficients are formed once, when the `WindField` is built. This gives
   `v_k(phi)` and its derivative `v_k'(phi)` at the point.
2. **Then `ln p`, at the point.** The values `v_k(phi)` are interpolated in `x` by the same PCHIP,
   with node slopes formed from the `v_k(phi)` of the interval and its two neighbors (four rows at
   most, so the cost per point stays independent of the grid's size).
3. **The derivatives are the interpolant's own** (decision L's rule kept). `(du/dln p)_phi` is the
   derivative of step 2's cubic in `x`. `(du/dphi)_p` is the derivative of the same expression in
   `phi`, carried through step 2: the Hermite form differentiated with `v_k'` in place of `v_k` and
   the slopes' own derivatives, `dm/dphi`, from the partials of the harmonic mean with respect to
   the two secants (a tangent of the same formulas, not a separate fit and not a finite difference).
4. **Evaluated in local form only** (each piece as the node value plus powers of `x - x_k`), so
   that a node returns the file's value bit for bit and a column with no vertical shear returns
   exactly the same value at every pressure. The Hermite basis form is not used: `h00 + h01` is not
   exactly one in floating point (§7 ruling 1).
5. **Unchanged from decision L:** extrapolation refused (decision N, `_refuse_outside`); the index
   rule at a node (`searchsorted(..., side="right") - 1`, last interval at the last node); the API of
   `wind_at`, `wind_derivatives`, `wind_on_mesh` and `reference_wind`.

**Properties, which the step's checks confirm.**

- **Node values exact.** Every node returns the file's `u_total`; `reference_wind` at a latitude
  node returns the reference row.
- **No overshoot.** Each one dimensional PCHIP stays within its interval's two end values, so the
  wind at any point lies within the four corner values of its cell. No new maximum or minimum is
  created anywhere; a jet's peak stays at the node value.
- **`(du/dln p)` is continuous everywhere,** in pressure and in latitude, since for every latitude
  the profile in `ln p` is one C1 PCHIP of continuous data. This is the derivative that sets the
  temperature steps: they go.
- **`(du/dphi)` is continuous across every latitude node and every pressure node.** It can still
  jump on one kind of locus: where a vertical secant `v_(k+1)(phi) - v_k(phi)` changes sign along
  latitude, that is, where the wind's vertical profile turns over exactly at a pressure node, the
  harmonic mean switches to zero and its derivative jumps. The jump is about one half of the
  difference in `(du/dphi)` between the interval's two pressure nodes (measured, §8), and it is a
  slope change in the shear term along that locus, not a step in any delivered quantity. A wind with no vertical shear, or one
  monotonic in pressure at every latitude, has no such locus.
- **Linear data reproduced exactly.** Where three consecutive nodes are collinear in either axis,
  the interval between the inner two is the line. A constant column is constant.
- **A slope change at a node is rounded over the node intervals either side** (the shape
  preserving slope at a corner is the smaller, often zero). The node values are unchanged, so the
  shear integrated across each interval, `u` at its upper node minus `u` at its lower node, is
  unchanged too; what changes is how that shear is spread inside the interval. A sharper corner is
  had by placing nodes closer around it in the hypothesis: structure finer than the file's grid is
  not in the hypothesis (SPEC_06 deliverable 3).

**Where it applies.** Every read of the wind: `lib.windfield.WindField` (the model, through the
callers listed in §3 deliverable 2) and `refrac.anchor.wind_of_latitude` (the reduction),
which becomes `WindField(wind).reference_wind`, so that one interpolant serves both and the
at-anchor identity holds by construction. The composition is not touched.

---

## 2. The behaviours expected

1. **No temperature step at a wind node.** The delivered temperature is continuous across every
   wind pressure node (it may change slope there, since PCHIP's curvature jumps at nodes), and the
   shear term along an isobar has no step at the 0.5 deg latitude nodes.
2. **Nearly the same answer** (author, 3 October). Where the wind is smooth on its nodes, every
   result is what decision L gave, to within small numerical differences of the size the two
   interpolants differ between nodes. Where the mathematics leaves a value unchanged (node values,
   a wind with no vertical shear in `ln p`, the at-anchor identity), it is unchanged to machine
   epsilon, at the accepted bound of the check that already tests it. In the shear experiments, away
   from the intervals adjacent to a slope change, the delivered temperature is also nearly the same,
   because the shear integrated across every interval is unchanged.
3. **Nothing non-physical introduced.** No wind outside the node values of its cell: no new
   maximum, minimum or reversal of the wind anywhere.
4. **Corners rounded where the hypothesis has them,** within the node interval on the sheared side,
   by design (§1). The shear there peaks at 4/3 of the zone's, and no interpolant through the nodes
   can avoid some overshoot at a corner (v0.4, §8 ruling 4). Reported, not bounded.
5. **The pressure identity smaller where a wind node interval spans several anchor layers;** where it
   does not, it can move either way with where the change lands in the layer (v0.4, §8 ruling 7).
   What remains is the anchor's level spacing, for SPEC_08.

---

## 3. Step 1: the interpolant, in the model and the reduction

**Deliverables.**

1. **`lib.windfield.WindField`** reads by decision L2 (§1). Its module docstring and every docstring
   that states the linear rule are reworded (`tools/wind/shear.py`'s `u_at_pressure` included).
2. **No caller changes** beyond `refrac.anchor.wind_of_latitude`, which returns
   `WindField(wind).reference_wind` and keeps its name and its signature. The callers that reach the
   new interpolant, named in the report: `tools/wind/shear.py`; `lib/kernel.py`; `forward/transfer.py`
   (the column, `_wind_table`, the geoid gauge wind, the columns, the production); `forward/production.py`;
   `refrac/anchor.py`.
3. **The registered files rebuilt** after the acceptance commit: the full sweep,
   `python tests/step04_0/sweep.py`, since the reduction's wind moves `phi_c` and `r0` and so every
   file of the chain; then `--runs`. Every value that moves is compared with its predecessor and the
   largest change per variable reported.
4. **The SPEC_05 experiments rerun** through `tests/step05_4/run_experiments.py` at both spacings,
   with the comparison figure redrawn against the SPEC_06 values; run 7f's F7 and F8 viewed by the
   author.
5. **The figures read the wind through `WindField`** (v0.3, §7 ruling 5):
   `tools/plots/figures_product.py` and `tools/plots/figures_inputs.py` stop interpolating the wind
   with `np.interp` and read it from the run's `WindField`, so that what is drawn is what the model
   used. There is then one interpolant of the wind in the code base.
6. **Decision L amended where it is stated** (SPEC_04 Step 1 deliverable 2, and SPEC_00 where it
   names the rule) by the reviewing agent at this specification's closure, as SPEC_04's amendments
   were applied at its closure.

**Acceptance (`tests/step07_1/accept_step07_1.py`).** Runs in copies under `reports/step07_1/`, from
the working tree with the `-dirty` refusal relaxed for the script only, as in SPEC_05 and SPEC_06.
"Every wind file" means the registered source and run wind files and the eleven experiment files.

1. **The rule.** On every row and every column of every wind file, the one dimensional rule agrees
   with `scipy.interpolate.PchipInterpolator`, values and first derivatives, to 1e-14 relative.
2. **Nodes exact.** `wind_at` at every node of every wind file returns `u_total` bit for bit;
   `reference_wind` at every latitude node returns the reference row bit for bit.
3. **No overshoot.** On a sample of at least 20 points per cell along each axis of every wind file,
   every value lies within its cell's four corner values, to 1e-12 of the file's largest `|u|`.
4. **Continuous slopes, and the interpolant's own.**
   - Across every pressure node and every latitude node of every wind file, the one sided limits of
     `(du/dln p)` agree, and those of `(du/dphi)` agree, to 1e-9 of the file's largest value of each.
   - The largest jump of `(du/dphi)` anywhere off the poles (the turnover loci of §1), found on the
     same sample, is reported with its place, and is smaller than decision L's largest jump of
     `(du/dphi)` across a latitude node on the same file. No present file has such a locus (§7
     ruling 1), so a synthetic wind file with one is added (a vertical profile whose `du/dln p`
     changes sign across a pressure node at some latitudes and not at others): there the measured
     jump is compared with §1's prediction, about one half of the difference in
     `(du/dphi)` between the interval's two pressure nodes, and with decision L's on the same file.
   - Away from those loci, both derivatives agree with centered differences of `wind_at` (steps
     1e-6 rad and 1e-6 in `ln p`) to 1e-6 relative.
5. **No vertical shear, no vertical derivative.** Under the closure wind, `(du/dln p)` is exactly
   zero at every sampled point, and the SPEC_04 Step 3 closed form check passes at its accepted
   bound (its closed form taking `(du/dphi)` from the interpolant, if the suite forms it itself).
6. **Nearly the same answer, the chain** (behaviour 2). Against the registered files before the
   sweep:
   - `phi_c` within 1e-4 deg and `r0` within 100 m in every reduction, and the record's polar
     radii and polar asymmetry within 100 m (the reduction's march crosses the equatorial jet, where
     the two interpolants differ most, §7 ruling 2);
   - the closure product and the Lindal transfer product: delivered temperature within 0.1 K at
     every level, `ln N` and the pressure on the labels within 1e-3 relative, radii and altitudes
     within 100 m;
   - the at-anchor identity (`step04_5`, the transfer at the anchor's own latitude against the
     closure product) at its accepted bound.
   The measured changes are reported beside the predictor, the largest `|u_L2 - u_L|` on each file
   and at `phi_c`. They are expected to be far below the bounds.
7. **No steps** (behaviour 1).
   - Run 7f at 5e4, along the 999 mbar isobar: the largest step in `S/g` between consecutive curve
     nodes at least ten times smaller than in SPEC_06 Step 1, and no pattern tied to the wind's
     latitude nodes. Old and new drawn together.
   - Run 7f's lowest 50 km: at every anchor layer straddling a wind pressure node, the layer's
     `dT/dln p` lies within the range of its two neighboring layers' values widened by 10 percent
     (a continuous profile passes even where its slope changes; a step does not). Old and new F7
     drawn together.
   - Every experiment run: the largest pressure identity not above its SPEC_06 value; the two
     spacings within 10 percent of each other (SPEC_06 check 6, rerun); the largest run's remaining
     identity explained by the layer law of SPEC_06 check 7, applied to a profile now continuous,
     with a measurement.
   - In the experiments with no slope change in the hypothesis (runs 2 and 3a to 3c), the delivered
     temperature within 0.1 K of SPEC_06 at every level. In the others, at every level farther than
     one wind node interval from a slope change, within 0.5 K; the levels inside those intervals are
     reported, not bounded (behaviour 4).
8. **The quadrature is adequate for a cubic interpolant.** Run 7f at 5e4 three ways: two Gauss
   points per piece (the model as delivered), four, and two with the column's isobar map knots added
   to the breakpoints (§7 ruling 7). The model as delivered within 0.1 K of each of the other two at
   every level; the three measured differences reported. The map knots are measured here, not
   adopted.
9. **The cost.** The Lindal transfer run (14 s after SPEC_06) and run 8 (16.5 s) each no more than
   twice as long; the times reported. Holding the latitude step's rows per column is allowed if it
   is needed, provided every value is bit for bit what the uncached reading gives.

**Regression.** `lib`, `refrac` and `forward` reach every product: the full set, about 100 minutes,
on the swept products. An accepted suite whose stated value moves only because of decision L2, within
check 6's bounds, is updated in this step and named in the report with its old and new value; a
suite that forms its own expectation by the linear rule is restated to the interpolant's own. A value
that moves by more than check 6's bounds stops the step and is reported. Reference counts stay unless
a check is added or removed. No accepted suite is edited silently or left failing. The suites the
reading named are handled as §7 ruling 4 says: the hash table recorded in this step's report and
`step02_1` and `step02_6` pointed at it; the read-only fixtures kept, not regenerated.

---

## 4. For the coding agent's reading

1. **The construction.** Whether §1 gives what it claims: `(du/dln p)` continuous everywhere,
   `(du/dphi)` continuous across nodes, jumps only on the turnover loci. Where those loci are in
   every wind file, and the predicted jump there.
2. **The predictor.** The largest `|u_L2 - u_L|` on every wind file, and at `phi_c` of every
   reduction, so that check 6's measured changes can be read against it before the run.
3. **The experiments' grids.** For each experiment, the wind node spacing in `ln p` near `p_s` and
   the stop pressure, and the number of anchor layers per node interval there; this predicts how
   far each identity falls (behaviour 5).
4. **The accepted suites** whose expectation is formed by the linear rule or stated to a bound
   tighter than check 6's, by name.
5. **Any other interpolation of the wind** in `src` or `tools` outside `lib.windfield` and
   `refrac.anchor`.
6. **The derivative in latitude** by the tangent of §1 item 3: its cost against a centered
   difference, and anything in it the coding agent would do differently.

---

## 5. Decisions (confirmed by the author, 3 October 2026)

1. **The order of construction** (§1): latitude first, then `ln p`. It makes the vertical
   derivative, the one that sets the temperature steps, continuous everywhere, and leaves the only
   possible slope jumps in `(du/dphi)` on the loci where the wind turns over in the vertical exactly
   at a node. The other order would put them in `(du/dln p)` instead.
2. **The reduction changes with the model** (§1, §3 deliverable 3). The cost is the full sweep;
   the alternative, the reduction left linear, breaks the at-anchor identity.
3. **The shear cases keep their node values** and their corners are rounded within one node
   interval (behaviour 4). If an experiment needs a sharper corner, its hypothesis is written with
   nodes closer around it; no change to the shear tool now.
4. **The bounds:** behaviour 2 at 0.1 K, 1e-3 relative in pressure and 100 m for smooth winds and
   the chain; 0.5 K in the experiments away from corners (SPEC_06 check 5's bound); twice the
   cost.
5. **The composition stays linear.** Its own specification if it is ever needed.

---

## 6. Revision history

- v0.1, 3 October 2026: first draft, from SPEC_06 §5 ruling 3 (the 0.5 deg staircase is decision L),
  REVIEW_06_step1 (the identity's mechanism), and the author's direction of 3 October: PCHIP, with
  the caveat that it must still give nearly the same answer, differences of small numerical size
  and machine epsilon tolerated.
- v0.1 revised the same day (author): the stability diagnostic taken out. The slope changes are
  corrected first; the superadiabatic profiles wait for a later specification. One step.
- v0.2, 3 October 2026: the five decisions of §5 confirmed by the author as written, to be revisited
  only if results call for it.
- v0.3, 3 October 2026: the coding agent's reading ruled on (§7): local form only (§1 item 4); the
  figures read the wind through `WindField` (deliverable 5); check 4 off the poles and with a
  synthetic turnover; check 6 adds the polar radii; check 8 three ways; check 9 on SPEC_06's times;
  the hash table and fixtures in the regression.
- v0.4, 3 October 2026: the interim report of Step 1 ruled on (§8): `phi_c`'s bound restated; the
  corner at a node stated as the mathematics requires; a dense-node reference measurement added
  before the author decides between L2 and L; check 7 restated in part.
- v0.5, 3 October 2026: closed, not adopted (§9), at the author's direction.
- v0.6, 3 October 2026: §9 addendum, the redesign of the vertical tabled (author).

---

## 7. Rulings on the coding agent's reading (3 October 2026)

`reports/REPORT_07_preexecution.md`, measured on the 27 kind W files in the tree at `a238306` with
decision L2 built from scipy beside decision L (`reports/step07_1/preexecution.json`). A careful
reading; every point is accepted, with the rulings below.

1. **The construction.** Confirmed as §1 states it. Two corrections taken into the text: the local
   form only, since a Hermite basis is not bit exact for a constant column (§1 item 4); and the
   loci, which the reading found only at the poles, where kind W holds the wind at exactly zero and
   both sides are zero. The reason no present file has a genuine locus is worth keeping: every shear
   case is `u_s(phi) F(p)`, so adjacent vertical secants agree in sign wherever `u_s` is not zero.
   A Monte Carlo hypothesis need not have that form, so check 4 searches off the poles and adds a
   synthetic file with a turnover, as the reading suggests, so that §1's prediction is measured
   once.
2. **The predictor.** Small where the runs read the wind: 1.04 m/s at most between 9 and 32 deg,
   -0.23 m/s at `phi_c`, an estimated -5 m in `r0(10 N)` against the 100 m bound. The largest
   differences, 5.3 m/s and up to 16 m/s in run 7, sit at the equatorial jet's sharpest bins south
   of the equator, where the jet's curvature is largest. Neither reading is the more correct there:
   both pass through the same nodes, and L2 stays within each cell's corner values (check 3). The
   reduction's march crosses that region on its way to the south pole, so check 6 now bounds the
   record's polar radii and asymmetry too. The estimate of `phi_c`'s move is not needed before the
   run; check 6 measures it.
3. **The experiments' grids.** Accepted, and the prediction stands as behaviour 5 says it: at `p_s`,
   where each adjacent wind interval holds five anchor layers (three and two in run 7f), the
   identity should fall several times; at run 5's stop pressure, one layer, little; run 5 may remain
   the largest. That remainder is the anchor's spacing, for SPEC_08.
4. **The accepted suites.**
   - `step02_5`: the expectation is restated to scipy's `PchipInterpolator` on the reference row,
     formed in the suite, so it stays independent of `lib.windfield`.
   - `step04_1` check 12: restated to what L2 promises, the one sided values at a node agreeing and
     the derivatives agreeing with scipy's; the sheared check's analytic `du/dln p` takes the
     reference row by scipy's PCHIP in latitude.
   - The cylinder-extended synthetic wind in `step04_2`, `step04_3`, `step04_4` and
     `step04_5/fields.py` is the tests' own construction of a file and is kept. Where a suite then
     compares the model's reading of that file between nodes with the construction, the comparison
     stays at decision Q's floor if it passes there, and is otherwise restated to the file's nodes;
     each named in the report.
   - **The hash table.** The new table is recorded in this step's report (`REPORT_07_step1.md`, in
     the sweep's section), and `step02_1` and `step02_6` are pointed at it, as SPEC_05 Step 0 did.
   - **The read-only fixtures** (`tests/step03_3/fixtures/`, `tests/step04_0/fixtures/swept/`) are
     not regenerated: they are the record of what those steps produced. A comparison against them
     of a value that L2 moves is restated to check 6's bounds against the same fixture, and a value
     L2 cannot reach stays exact. Named in the report.
   - `step05_2` and `step04_5` check 11 pass after the sweep's commit, as the reading says.
5. **The figures' own interpolation.** Accepted as deliverable 5: the figures read the wind through
   `WindField`. A figure showing decision L while the model uses L2 would be wrong in the way the
   author would least expect.
6. **The derivative in latitude.** The tangent, as the reading proposes: the latitude step by scipy's
   `PchipInterpolator(phi, u, axis=0)` formed once with its `derivative()`, the `ln p` step by hand
   on the four rows in local form, `dm/dphi` in closed form and zero on the harmonic mean's zero
   branch. The end rule's clamps act only on the grid's top and bottom rows, outside the model's
   range; the tangent follows the branch taken.
7. **The isobar map's knots as breakpoints.** A good catch. Along a column the wind is read at
   `p(Phi)`, and `ln p(Phi)` is the map's piecewise linear interpolation, so `S/g` changes slope at
   every map knot as well as at every wind node. Under decision L the wind nodes' breaks dominated;
   under L2 the knots may be what check 8 sees. Measured three ways in check 8, not adopted here:
   adding them is a second change, and this step must stay one revertible change. If the model as
   delivered misses check 8's 0.1 K against the knotted rule, the step stops there and the knots
   come back to the author as a step of their own. They are also the place SPEC_08's sub-levels act.
8. **The cost's baseline.** SPEC_06's times, 14 s for the Lindal transfer run and 16.5 s for run 8
   (check 9). Caching the rows per column is allowed, bit for bit.

---

## 8. Rulings on the interim report of Step 1 (3 October 2026)

`reports/REPORT_07_step1.md` (interim). The step stopped where the regression paragraph says it
must, at check 6, and was filed with every other check measured. Checks 1 to 5, 8 and 9 pass: the
rule is scipy's to 4e-16, every node exact, no value outside its cell, both slopes continuous, an
exact zero vertical derivative without vertical shear, the cost 1.4 to 1.5 times SPEC_06's. Away
from the corner at `p_s`, run 7f's delivered profile has lost the zigzag of its lowest 40 km (the old
and new F7), which is what this specification was written to remove.

1. **Check 4's measurement** (one sided values read 1e-12 from the node, since PCHIP's curvature
   jumps there) is accepted. The synthetic turnover's jump measured at one half of the rows'
   difference, not my "one third"; §1 corrected to the measurement.
2. **The map's knots** (§7 ruling 7) are settled: 0.03 K, not needed.
3. **`phi_c`.** Explained by the predictor: the wind at the anchor changes by -0.227 m/s, and the
   geoid's slope term `2 Omega Delta u sin(phi_c) / g` gives 2.1e-4 deg against 2.4e-4 measured.
   My bound of 1e-4 deg was set before the predictor at `phi_c` was known; what `phi_c` reaches,
   `r0` (22.9 m) and the products (1e-3 to 5e-3 K), is far inside its own bounds. Restated to 1e-3
   deg, about 1 km along the surface. L2 is not held back at the anchor: that would break the
   at-anchor identity.
4. **A corner at a node, stated as the mathematics requires.** Take a corner at a node with no shear
   on one side and shear `delta` on the other. The node values fix the mean shear over each
   interval: zero on the flat side, `delta` on the sheared side. If the shear is continuous at the
   node, with some value `m` there, then either `m > 0` and the shear must go negative somewhere on
   the flat side (the wind reverses), or `m < delta` and it must exceed `delta` somewhere on the
   sheared side. **No interpolant through the nodes renders a corner with a continuous shear and
   no overshoot.** By the thermal wind, a corner in `u` is a step in temperature; a smooth rendering
   of it must overshoot somewhere. PCHIP takes `m = 0`: none of it on the flat side (the wind never
   reverses) and all of it on the sheared side, where the shear peaks at 4/3 of the zone's two
   thirds of the way across the interval. That is the -85 K at level 57 against the zone's -74, and
   the near-zero change at level 59, just past the rounded corner. Decision L renders the same
   hypothesis as the step itself. Neither is the hypothesis's own temperature, which is a
   discontinuity. Behaviour 4 stands with its size now known.
5. **Which reading is closer to the hypothesis: measured before the author decides.** Away from the
   corners, L and L2 differ by up to 7.2 K in run 7, whose wind is linear in `p` and so curved in
   `ln p` everywhere in the zone. One of the two is that far from the hypothesis. A reference
   answers it:
   - Runs 5, 7 and 7f rebuilt with the same hypothesis on ten times the pressure nodes (the source
     built on the denser grid, the shear case unchanged), at 5e4, in copies under
     `reports/step07_1/dense/`. Run under L2 (the working tree) and under L (a worktree at
     `a238306`); dense L and dense L2 should agree away from the corners, and their difference is
     reported as the reference's own uncertainty.
   - For each run, at every anchor level, coarse L minus dense and coarse L2 minus dense in the
     delivered temperature change from run 2: the largest and the mean magnitude away from the
     corner intervals, and the values inside them, with the pressure identity of each.
   - Nothing in the code changes for this measurement.
6. **Check 7, restated in the parts that are test design** (whatever the author decides):
   - The isobar test is taken on the anchor level nearest 500 mbar in run 7f, inside the zone and
     away from any corner: the largest step in `S/g` at least ten times smaller than SPEC_06's on
     the same isobar, and the mean step at the wind's latitude nodes against elsewhere reported.
   - The continuity test excludes the wind intervals adjacent to a slope change (behaviour 4).
   - The identity is reported with the layer law, not bounded against SPEC_06 (ruling 7).
   - The bound away from the corners waits for ruling 5: SPEC_06's values under L are not the right
     reference if L is the one in error on a curved hypothesis. It is restated against the dense
     reference.
7. **The identity at the stop pressures.** Run 4 rises to 7.1e-3 and run 6 to 2.4e-2, the latter
   above `step05_4` check 4's accepted 1e-2. The layer law explains it: under L, run 6's break sat
   at `s` = 0.495 of its anchor layer, almost on the zero crossing, so its small value was the luck
   of placement; L2 moves the change within the layer (effective `s - s0` about 0.28, 208 Pa of a
   741 Pa scale). It is the anchor's spacing, SPEC_08's subject. Behaviour 5 is corrected: the
   identity falls where a wind interval spans several anchor layers (five to six times at `p_s`) and
   can move either way where it does not. How `step05_4` check 4 is carried is for the author with
   the decision of ruling 8.
8. **For the author, after ruling 5's measurement:** keep L2, with corners rendered as ruling 4
   states and realistic hypotheses written without corners (the shear tool's transitions given a
   width of a few node intervals, when the realistic case is built); or revert to L. If L2 is kept,
   `step05_4` check 4 is either restated for runs 4 and 6 until SPEC_08 restores it, or SPEC_08 is
   taken first. The sweep and the full regression wait for that decision.

---

## 9. Closure (3 October 2026, author)

**Closed, not adopted.** Decision L stands in the code. The work of Step 1 is kept on a branch, not
merged; the specification, the reading and the report are kept on `main` as the record.

**Why.** The dense reference (REPORT_07_step1 §7) settled the question this specification asked
and showed it was the wrong question.

- Wherever the hypothesis is resolved, L2 is the closer reading: within 0.15 K of the reference
  over most of run 7's zone, where L is off by up to 4.4 K, alternating level to level.
- At a change in shear sharper than the anchor's layers, neither reading converges. The identity,
  the level next to the change and the offset carried below it depend on where the change falls in
  an anchor layer, not on the interpolant.

The author's reading (3 October): every change in the wind or the temperature is a corner, and the
numerical trouble grows with its sharpness against the grids it passes through. The cause is where
the variables live. The model transports `Phi` and `ln N` on the same levels and rebuilds the
hydrostatic column from them, which is a collocated arrangement. SPEC_00 §7.3 rule (ii) asked for an
explicit staggered arrangement, and it was never written. The redesign of the vertical on a
staggered arrangement (reviewing agent's design note of 3 October, in the author's project) is
handed to a new agent. The interpolant is judged again inside it, on general tests rather than
bounds set per run.

**Carried forward:**

- the corner at a node (§8 ruling 4: no interpolant through the nodes renders a corner with a
  continuous shear and no overshoot);
- the dense reference and its measurements;
- the cost of L2 (1.4 to 1.5 times SPEC_06's);
- the map's knots as breakpoints (not needed, 0.03 K);
- the `phi_c` sensitivity to the wind at the anchor (2.4e-4 deg for 0.23 m/s).

**Not carried forward:** check 7's bounds, which were set per run. No tolerance is set for one
case.

**Addendum (author, 3 October 2026).** The redesign of the vertical is tabled. A realistic shear
case and the first comparison with CIRS come first, as originally planned, with the current
numerics.
- **Why it is good enough for that.** The trouble measured here grows with how sharply the shear
  changes against the anchor's layers. A hypothesis written smoothly, on wind nodes dense enough
  where it curves, stays in the regime where it is small.
- **How it is watched.** The pressure identity every product carries, level by level, tells
  whether the numerics limit a given run.
- **What reopens the redesign.** The comparison, or the Monte Carlo, which will draw hypotheses
  nobody inspects. The design note and its handoff are kept for then.
