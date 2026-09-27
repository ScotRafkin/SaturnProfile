# REPORT 04, Step 4. The tracing, the transfer, the outer loop, the estimate

Coding agent, 26 September 2026; refreshed under the REVIEW_04_step4 addendum and SPEC_04 v0.16.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.16, section 6, decisions C, D, F, G, H, K, M
and Q, and the section 15 and 16 rulings. v0.14's decision R touches no step; section 7 carries
what it changes for Step 5 and finding 4 answers the one sentence of it that reaches back here.
Acceptance script `reports/step04_4/accept_step04_4.py`, output `reports/step04_4/output.txt`,
console `reports/step04_4/run.txt`. The diagnostic of finding 1 is
`reports/step04_4/diagnose_anchor_placement.py`.
Regression `reports/step04_4/run_regression.sh`, output `reports/step04_4/regression.txt`.

**Eleven of eleven checks pass.**

The two quantities section 16 ruling 1 separates are seven orders of magnitude apart, which is the
point of separating them:

* **The transfer's own identity is round-off.** The Lindal anchor's curves carried to 60 N on the
  M = 2 mesh and taken straight back to the gauge, with no anchor file and no placement in
  between, return its own `C_i` to **3.553e-15** against a bound of 1e-12 and its own `Phi_k`
  **exactly**.
* **Through a written anchor it is 6.484e-08**, against a bound of 1e-7. That is decision G: an
  anchor read from a kind N file is placed by the field-line integral of `|g_eff|` over its own 66
  tabulated levels, while the trace that wrote it is the characteristic in `(phi, Phi)`, and no
  mesh refines the anchor's levels. Finding 1 carries the measurement.

The three earlier fillings of checks 8 and 10 are what this one is built on and they are kept: the
lattice difference of the first (1.209e-04), which section 15 ruling 1 removed by writing the
anchor on the run's own mesh, and the interpolation of the second (8.191e-08), which section 16
ruling 1 removed by integrating the column on the traced `Phi`. The lattice variant is still run
beside the identity and reported, not bounded, at 1.209e-04.

What the step establishes:

* **The transfer reproduces the reviewing agent's independent values.** To 60 N, `d ln N` agrees
  to **9.133e-07** against a bound of 1e-4, three of the four levels to 2e-7, and the isobar shift
  to **0.008 percent** against 1 percent. To 10 N, 6.807e-06 and 0.055 percent.
* **The gauge isobar does not move**, exactly, at every latitude, because `I` is exactly zero on
  it; and at M = 1 the estimate is the identity to one unit in the last place of `ln N`.
* **The scheme is reversible.** Under a sheared wind, out to 10 N and back returns `Phi_k` to
  6.601e-08 m2/s2 and `ln N_k` to 4.103e-13, falling by 4.93 and 3.59 under halving.
* **The cylinder wind holds its identity at decision Q's floor**, 163.3 m2/s2 and 1.064e-03,
  inside the class the two Step 3 kernels gave.
* **At M = 2 the estimate is what its equations say**: `phi_r` the weighted centroid to
  0.0e+00 degrees, the reduced chi-square 3.689e-12, the synthetic anchor on the marched reference
  surface exactly, and a validation anchor leaving `C` and `phi_r` at their M = 1 values.
* **Two numerical failures are named and refused**, curves that cross and a loop that will not
  converge; a curve reaching the mesh edge is not a failure but a report, and the loop grows the
  mesh and repeats.

Five defects were found and fixed across the three fillings, all five mine; findings 1 to 5 carry
them. Three were found by re-reading the code or by chasing a number rather than by a check, and
each now has a check or a measurement that would have caught it.

**How this pass was run.** SPEC_04 v0.16 changes the M = 2 construction and nothing else, and the
review asks for the M = 2 runs only, so checks 1, 2, 6, 8, 9 and 10 were rerun and checks 3, 4, 5
and 11 carry their rows from the run of 25 September, kept as `output_v015.txt` and `run_v015.txt`,
each row marked `CARRIED` in `output.txt` with the file it came from. Checks 1 and 2 were rerun
although nothing in them changed, because checks 8 to 10 are built on their states; both reproduced
their accepted values. Check 7 reports the pass counts of this pass's runs and points at the
carried rows for the rest. The regression was run in full after the change.

## 0. Before this step

HEAD `19ffd57`, Step 3 accepted at `bb75d16`. Every product on disk is the swept one: the
reduction chain clean at `1c9b310` and the closure product at `780a2de`. This step registers
nothing, changes no existing module, and writes nothing inside the repository outside
`reports/step04_4/`, which git ignores. `forward/transfer.py` and `forward/estimate.py` are new.

## 1. What was built

**`forward.transfer.trace(mesh, shear_integral, geopotential_m2s2, phi_from, phi_to)`.** RK4 of
`dPhi/dphi = -I(phi, Phi)` (Eq. A24, B6.2) on the mesh's own latitude nodes, with `I` bilinear
between them. The bilinear evaluation is `lib.kernel.transfer_kernel` with its composition
arguments omitted, which is what that function is with no composition term, so the tracing and the
transfer read the mesh the same way. Levels may be given in any order and come back in it. Two
outcomes the specification names are reported rather than assumed: a curve that leaves the mesh's
geopotential range stops the trace, which returns what it has with the side and the excess in
`reached`, and a pair of curves that cross is refused in the specification's own words.

**`forward.transfer.transfer(mesh, s_over_g, curves, ln_N, composition_latitude_rad,
composition_slope)`.** Eq. A28: the trapezoid of `K = S/g + (d ln(R_bar/m_bar)/dphi)_p` along each
curve on the curve's own nodes, taken outward from the starting latitude in both directions, as
`lib.kernel.shear_integral` takes its vertical integral, so that the value at the start is `ln N`
unchanged and the integral carries its own sign whichever way the curve runs. The composition
slope is one row per level on the composition file's latitude grid, which is what
`lib.kernel.composition_term` returns for that level's label; omitted, `K` is the shear term
alone, which is what a composition uniform in latitude gives and what this run has.

**`forward.transfer.outer_loop`.** The fixed point of SPEC_00 section 7.2. One pass places every
anchor on its own levels under the wind the map gives along its column (`place`, Step 1
deliverable 3), builds the reference surface through the first anchor and the column at every
latitude node, forms `S`, `S/g` and `I` on the nodes, traces every anchor's levels from its own
latitude to both edges of the mesh, and reads the new map off the curves. It repeats until the
largest `|d ln p|` on the mesh, between the map that built a pass and the map that pass's curves
give, is below `relative_tolerance_ln_p`, and refuses at `max_iterations` naming the residual. A
curve that reaches a geopotential edge does not end the run: the mesh is extended and the pass
repeats, both the excess and the amount added recorded.

**`forward.transfer.IsobarMap`.** `p(phi, Phi)`: at each latitude node, `ln p` against `Phi`
through that latitude's knots, linear between them and continued at the slope of the last interval
outside them. The knots are the isobar labels and the geopotential the curves put them at, the
nearest anchor's in latitude; on the first pass they are the anchors' own productions at their own
latitudes, the same at every latitude, which is the barotropic guess.

**`forward.estimate`.** `gauge_latitude(anchors)` is decision K's weighted centroid,
`1 / mean(sigma_i^2)` over each construction anchor's levels, and the one place `phi_r` is
computed; `sigma_ln_N(anchor)` is the measurement and season terms in quadrature, read from the
anchor object and from no file (decision J). `estimate(arrivals, sigma_K, phi_r)` forms the union
of the construction anchors' arrival levels, interpolates each anchor's `C_i`, `sigma_i` and label
onto it within its own span and marks it absent outside, weights by `1 / (sigma_i^2 + P_i)` with
`P_i = (sigma_K (phi_i - phi_r))^2`, and returns `C` (A30), its variance, the reduced chi-square
per level with the inflation where it exceeds one, `D_ij` (A33) for every pair including validation
anchors, and the isobar label of each union level with the anchors' disagreement beside it.

## 2. Decisions

1. **The isobar map is continued beyond an anchor's outermost labels at the slope of its last
   interval in `ln p` against `Phi`.** The mesh reaches one geopotential spacing past the anchor's
   levels at each end, so the map is asked for a pressure there. Step 2 decision 8 continues the
   flat map that way and section 14 ruling 2 named the clamped alternative as the Step 3 fault; one
   rule now serves both, in `IsobarMap` and in the acceptance script's flat map, which is the same
   object.

2. **Between anchors the map takes the nearest anchor's knots, and ties take the first listed.**
   The specification states that rule for the barotropic guess; it is kept for the map the curves
   give, so that the guess and the converged map are one function and the loop is its fixed point.
   The rule puts a seam at the midpoint of two anchors' latitudes, measured in finding 7.

3. **Each anchor's levels are traced from its own latitude to both edges of the mesh.** The map is
   wanted at every latitude node, so the two traces together cover the mesh; they are joined at the
   anchor's node, which both carry.

4. **The convergence measure is the largest `|d ln p|` over the mesh nodes** between the map that
   built a pass and the map that pass's curves give. The state returned is the one the curves were
   traced on, with the map that built it, so that the pressure, the kernels and the curves in it
   are consistent with each other; the map the curves give differs from it by less than the
   tolerance, which is what convergence means here.

5. **The union of arrival levels is the construction anchors'.** A validation anchor is
   propagated, placed, traced and reported like any other, and its `C_i` and every `D_ij` are on
   the union, but it does not decide which levels the estimate is made on. That is what leaves `C`
   at the values it would have without it, which is the check the specification asks for.

6. **The isobar label of a union level is the inverse-variance weighted mean, in `ln p`, of the
   labels the anchors present carry**, with the same weights as the estimate, and the spread
   between them is recorded per level. At M = 1 it is the anchor's own label exactly.

7. **The mesh is built once, on the first pass, and afterwards only extended.** A pass that moves
   an anchor's levels outside it extends it by what they need, as a traced curve does; the latitude
   nodes never change, so the map's knots, which are one array per latitude, stay as they are.

8. **`IsobarMap` holds its knots per latitude rather than as one rectangle**, because two anchors
   need not carry the same number of levels and the map at a latitude is one anchor's.

9. **The numbers `outer_loop` works from are its arguments, not the namelist.**
   `outer_loop(inputs, anchors, *, gauge_latitude_rad, target_latitude_rad, gauge_isobar_Pa, p_b,
   latitude_spacing_rad, geopotential_spacing_m2s2, relative_tolerance_ln_p, max_iterations,
   field=None)`: `inputs` carries the loaded composition, gravity, rotation and wind, and every
   number the loop uses is named. The driver of Step 5 reads them from the namelist; the acceptance
   varies the spacings and the wind file, which is what the specification's checks ask of it.

10. **The production at the target is formed in the acceptance script, not in the library.** It is
    Step 5's deliverable. The M = 2 identity test compares `p` and `T` there and decision H reports
    the pressure identity at every level of every run, so both are formed there through
    `lib.hydrostatic` on the geopotential the tracing gives, with the composition on the isobar
    labels. Finding 6 says why `forward.production.produce` cannot be used as it stands.

11. **The synthetic anchor of the M = 2 test carries the Lindal anchor's relative refractivity
    uncertainty**, not its absolute column, so that `sigma_ln_N` is identical on the two anchors and
    `phi_r` is the midpoint the specification states. Copying the absolute column instead moves
    `phi_r` off it, because the transferred `N` is not the anchor's.

12. **This step's own suite is not inside the regression loop.** It reads only the registered
    products and writes only inside `reports/step04_4/`, which git ignores, so the tree it sees is
    the tree the loop restores; running it twice would cost four hours and buy nothing. Its row in
    the regression table is the acceptance run itself, and the porcelain line at the end of the
    regression is the evidence that it left the tree as it found it.

13. **The mesh is an argument of `outer_loop`, and left out it is built from the anchors as
    before.** Section 15 ruling 1 has the script build the M = 2 mesh first and form the run on
    that same mesh, which the loop could not do while it built its own. Given a mesh, the loop
    still grows it when a curve or an anchor's levels need more, and refuses one that has no node
    at an anchor, the gauge or the target.

14. **The three M = 2 runs share one mesh, and it carries both gauge latitudes as nodes.** The
    identity, the validation case and the season case are compared level for level, so they have
    to be on one mesh; giving the season case its own would put back a smaller copy of the lattice
    difference ruling 1 removes. Its gauge latitude is computable before anything is traced,
    because the season column is the script's own and both anchors share the measurement column,
    and the script asserts both gauge latitudes against what `gauge_latitude` gives from the
    anchors as they arrive.

## 3. Acceptance results

| # | Check | Result |
|---|---|---|
| 1 | closure wind to 10 N at the namelist's spacings | PASS, `d ln N` 6.807e-06 of 1e-4; shift 0.055 of 1 percent |
| 2 | closure wind to 60 N | PASS, `d ln N` 9.133e-07; shift 0.008 percent |
| 3 | both spacings halved, carried | PASS, 4.574e-06 and 4.483e-06 of 2e-5; 0.035 and 0.178 of 0.2 percent |
| 4 | cylinder wind at decision Q's floor, carried | PASS, 163.3 of 300 m2/s2; 1.064e-03 of 2e-3 |
| 5 | curves that cross are refused, carried | PASS, both constructions refused, all runs ordered |
| 6 | a curve reaching the mesh edge, and the extension | PASS, both sides, 4 extensions in the loop |
| 7 | the loop's pass count for every wind file | PASS, 2 to 6 passes, every residual under 1e-8 |
| 8 | M = 2, the identity on the run's own mesh | PASS, `\|D_12\|` 6.484e-08 of 1e-7; the tracing's own 3.553e-15 of 1e-12 |
| 9 | a validation anchor | PASS, `phi_r` exact, `C` to 1.8e-15 |
| 10 | a nonzero season column after the hook | PASS, target 4.076e-08 of 1e-7, the same floor as 8 |
| 11 | a sheared wind, out and back, carried | PASS, 6.601e-08 m2/s2 and 4.103e-13, falling 4.93 and 3.59 |

**The closure wind, checks 1 and 2.** On the mesh of 421 by 783 nodes the loop converges in two
passes with the second residual exactly zero, which is what a wind that does not vary with pressure
must give and is not special-cased. To 10 N: `d ln N` +1.108357e-02 at the top level,
+1.098165e-02 at 10 mbar, +1.094319e-02 at the gauge and +1.089265e-02 at the bottom, each within
6.8e-06 of the stated values; the isobar shift -31,240.9 m2/s2 at the top and +11,333.5 at the
bottom, 0.055 percent from the stated -31,258.1 and +11,339.7, and exactly zero at the gauge. To
60 N, on 588 latitudes: +2.512071e-03, +2.492366e-03, +2.484913e-03, +2.475105e-03, within
9.1e-07, and -7,117.6 and +2,585.3 m2/s2, 0.008 percent.

Two measurements beyond the acceptance list are reported with check 1. At M = 1 the estimate is the
identity: `phi_r` is the anchor's own latitude, the union is its own 66 levels, and `C - ln N` is
1.8e-15, one unit in the last place; no `D` exists and the reduced chi-square is absent at every
level. And the pressure identity of decision H, the produced `p` at the target against the label
each isobar carries, is **5.789e-07**, beside the 4.1e-7 section 1 states.

**The halved spacings, check 3.** To 10 N the departure from the stated values falls from
6.807e-06 to 4.574e-06 and the shift from 0.0549 to 0.0348 percent, ratios 1.49 and 1.58. To 60 N
it rises, from 9.133e-07 to 4.483e-06 and from 0.0080 to 0.1779 percent, ratios 0.20 and 0.04.
Both are inside the 2e-5 and 0.2 percent bounds and no order is claimed: the 60 N run at the
namelist's spacings sits on a cancellation, and what falls or rises here is the departure from a
measurement, not from an exact value.

**The cylinder wind, check 4.** The file is rebuilt here by the Step 2 construction, the inversion
curve sampled at the mesh spacing over 3,962 latitudes, and the fixed point converges in five
passes. On the construction, at the anchor's own levels with no file grid in between, `u` is 9.488,
3.716, 2.167 and 1.925 m/s against decision P's stated 9.490, 3.717, 2.168 and 1.926, largest
departure **1.712e-03 m/s** against Step 2's 1e-2 bound. Read back through the written file at the
same points it is 9.552, 3.890, 2.167, 1.970, up to **0.174 m/s** away: a field constant on
cylinders, stored on 361 latitudes by 61 pressures and read back by decision L's interpolant, is
not constant on cylinders any more. That is decision Q's truncation measured at the wind rather
than at the kernel, and it is reported and not bounded. Under that field the traced isobars hold
the anchor's own to **163.3 m2/s2** of 300 and the transferred `ln N` to **1.064e-03** of 2e-3,
1.25 and 2.24 times the Step 3 line integrals' estimates and inside the 130 to 224 and 4.7e-4 to
1.1e-3 class SPEC_04 v0.13 states.

**The two refusals and the mesh edge, checks 5 and 6.** A constructed kernel that reverses sign
twice in the vertical, `I = 4e6 sin(2 pi Phi / 2e5)` per radian on a coarse mesh, carries the lower
of two curves past the upper one and is refused at 25.000000 degrees; an anchor with two levels at
one geopotential, a profile with a layer of no thickness, is refused at the starting node. The
three runs' own 66 curves stay strictly ordered at every node. A constructed kernel of 3e6 per
radian everywhere stops the trace at the first node with side `above` and excess 11,799.4 m2/s2,
and with its sign reversed at the other edge with the same excess; extending and repeating as the
loop does carries the curve across in four extensions. In the loop, the closure wind with
`u_total` reversed in sign turns the shift outward and the loop records four extensions and
converges in two passes. The closure and cylinder runs need none: their isobars move inward, the
top level down and the bottom up, so the anchor's own levels stay the outermost.

**The estimate at M = 2, checks 8, 9 and 10.** All three runs are on one mesh of 1,007 latitudes by
783 geopotential nodes, built before anything was traced for this run's required latitudes
(section 15 ruling 1) and carrying both gauge latitudes as nodes, and each converges in two passes
with the second residual exactly zero. `phi_r` is 45.402784 degrees, the midpoint of 30.805568 and
60, to **0.0e+00 degrees**, and what the anchors as they arrive give through `gauge_latitude` to
1.1e-16 radians. The union at `phi_r` has 131 levels from the two anchors' 66 each; the reduced
chi-square is 3.689e-12 against a bound of 1, the isobar labels agree between the anchors to
7.457e-08 in `ln p`, and the synthetic anchor's radius at the gauge isobar sits on the reference
surface marched from the Lindal anchor **exactly**.

The two quantities section 16 ruling 1 bounds apart:

* **The transfer's own identity**, with no anchor file and no placement in between: the Lindal
  anchor's curves carried to 60 N on this mesh and taken straight back to the gauge return its own
  `C_i` to **3.553e-15** against 1e-12, and its own `Phi_k` to **0.000e+00 m2/s2**, exactly.
* **The identity through a written anchor**: `|D_12|` **6.484e-08** against 1e-7, and the target
  profile departs from the M = 1 product by **3.242e-08** in `ln N`, 7.583e-08 in `p` and
  4.376e-08 in `T`, each against 1e-7. The difference between the two is decision G's placement,
  which finding 1 measures.

The lattice variant of section 15 ruling 1, the same anchor written from the separate M = 1 run to
60 N and substituted into this run's estimate, gives `|D_12|` 1.209e-04 and a target difference of
6.046e-05, and differs from the anchor this run wrote by 3.307e-05 in `ln N` and 93.6 m2/s2 in
`Phi_k` at 60 N: the tracing's discretization between two lattices at the namelist's spacings,
reported and not bounded.

With the weight of the synthetic anchor set to 0 the gauge returns to the construction anchor's own
latitude **exactly**, `C` equals its own `C_i` to 1.8e-15, and the reduced chi-square is absent at
every level, with the validation anchor still placed, traced and reported, its `C_i` from
-16.803458 to -8.075440 and `D_12` from -6.484e-08 to +1.970e-09. With a season column of 2e-3 set
by the script after the hook, which is 0.09 times the anchor's mean measurement term, `phi_r` moves
from 45.402784 to 45.349290 degrees, toward the anchor with the smaller uncertainty, and the
synthetic anchor's mean weight falls to 0.9927 of its M = 2 value; the target profile is the M = 1
product's to 3.229e-08 in `ln N`, 7.583e-08 in `p` and 4.388e-08 in `T`, the same floor as check 8
and for the same reason. The hook itself leaves the column zero and says so, `season_term` absent
and `propagation` "none: propagator not implemented".

**The sheared state, check 11.** The file is written in its three parts and read back through the
kind W reader, whose sum identity holds exactly and whose poles are zero; on the anchor's column
it runs 4.013 m/s at the top level against 2.167 under the closure wind. Out to 10 N and back
returns `Phi_k` to 3.255e-07 m2/s2 and `ln N_k` to 1.473e-12 at twice the namelist's spacings,
and to **6.601e-08** and **4.103e-13** at the namelist's own, falling by 4.93 and 3.59. The
curves stay ordered on both meshes and the loop converges in six passes on each, against a bound
of ten.

**The loop over every file, check 7.** Two passes for the closure wind at both targets and both
spacings and for all three M = 2 runs, with the second residual exactly zero; four for the cylinder
wind, residuals 2.977e-04, 1.939e-06, 1.085e-08, 3.836e-11; two for the reversed closure wind; six
for the sheared wind at both spacings, the residual falling by about forty per pass. Every run ends
below the namelist's 1e-8.

**What it cost.** The eleven runs took 14,889 s, four hours and eight minutes, on one core. The
column build dominates: every pass integrates `r` and `z` by RK4 at every latitude node, and
`lib.gravity`'s effective gravity is evaluated once per stage on scalars. A pass on the 421 by 783
mesh is about 200 s and on the 1,006 by 783 mesh of the M = 2 runs about 500 s. This is recorded
because the author asked what a transfer run costs (REPORT_04_step1 section 5b): a single
converged transfer at the namelist's spacings is two to six of those passes, so 7 to 20 minutes
for the closure wind and up to an hour for a wind with vertical shear whose first pass grows the
mesh.

## 4. Findings

**1. The M = 2 identity took three fillings to state correctly, and what it is bounded by now is
two things measured apart: the transfer, at round-off, and the placement of an anchor read back
from a file, at 6.5e-08.** Ruled on at section 15 ruling 1 and section 16 ruling 1; closed.

The first filing bounded one quantity at 1e-5 and measured 6.0e-05, and I reported the cause as the
test comparing two latitude lattices, which it was. Section 15 ruling 1 had the anchor written from
the Lindal anchor's trace on the M = 2 run's own mesh, which removed that: `|D_12|` fell from
1.209e-04 to 8.191e-08 and the synthetic anchor's reference-surface residual from 2.980e-08 m to
exactly zero, because the surface is now marched on the mesh the anchor came from. The lattice
variant stays, run beside the identity and reported at 1.209e-04, which is the tracing's
discretization between two lattices and a number worth having.

The bound ruling 1 set, 1e-10, came from a measurement of mine that had not passed through a
written anchor, and the second filing reported that correction and the floor underneath. The floor
was two parts, separated in `reports/step04_4/diagnose_anchor_placement.py`:

| | `place`'s `Phi` against the traced | `D_12` |
|---|---|---|
| the column interpolated from the mesh's nodes, 0.5 deg by 5e4 | 1.0183 m2/s2 | 3.3832e-06 |
| the same, at the namelist's 0.05 deg by 5,000 | 4.4896e-02 | 8.1914e-08 |
| the column integrated on the traced `Phi`, 0.5 deg by 5e4 | 4.0868e-02 | 6.4850e-08 |
| the same, at the namelist's spacings | 4.0855e-02 | 6.4837e-08 |

* **Mine, and it fell with the mesh.** The anchor's `radius_m` and its field-line height were read
  off the state's column by interpolating from the mesh's geopotential nodes onto the traced `Phi`,
  an O(dPhi^2) step of the script's. Section 16 ruling 1 has it removed and it is: the 60 N column
  is integrated on the traced `Phi` itself, which `lib.mesh.column` allows because it takes the
  nodes it is given, so this is a change of nodes and not of method. `|D_12|` is now **6.484e-08**.
* **Not mine, and it does not fall with the mesh.** What is left is decision G. An anchor read from
  a kind N file is placed by the field-line integral of `|g_eff|` along the profile's own tilted
  vertical (SPEC_03 Step 1), taken as a trapezoid over the anchor's **66 tabulated levels**, while
  the trace that wrote it is the characteristic `dPhi/dphi = -I` integrated on the mesh. The gap is
  **0.0409 m2/s2**, the same to four figures at two resolutions a factor of ten apart in both
  spacings, zero at the gauge level where both constructions put `Phi = 0` by definition, and
  growing outward. `D_12` is that times `d ln N / d Phi`, about 2.2e-06 per m2/s2 on this profile.

Section 16 ruling 1 bounds the two apart, 1e-7 through a written anchor and 1e-12 for the tracing's
own, and both are now measured in check 8: **6.484e-08** and **3.553e-15**, the second with the
traced `Phi_k` returning exactly. SPEC_04 section 1 states the mechanism in advance, that the radial
and field-line constructions of `Phi` "are values of one scalar field wherever the effective gravity
is conservative and differ by the shear's circulation around the drift loop otherwise"; at 60 N
`psi` is 4.7394 degrees.

**What was mine in this, plainly.** The 1e-10 of ruling 1 was set from my number and my wording:
the first filing called it "the second anchor written from this run's own trace", which reads as the
file path when the measurement substituted the traced values directly. While chasing the floor I
also said twice, in the course of the work, that the residual was resolution-independent and that
the drift geometry was the whole of it; the first was wrong and the second was three quarters wrong,
and the table above is what settled both. Section 16 ruling 2 carries the consequence downstream:
two anchors whose isobars coincide disagree in `D_12` at 6.5e-08 for no physical reason, invisible
against the anchors' declared uncertainties of 2.3e-2, and the product's `estimate` group states it
once as `identity_floor`.

**2. Check 4 compared a gridded read-back with a value measured on the construction. Mine.** The
first filing of this script read the cylinder wind back through the written kind W file at the
anchor's latitude and levels and compared that with decision P's stated column values, which are
measured on the construction itself. It failed by 0.173 m/s. The construction reproduces the stated
values to 1.712e-03 m/s. Before finding the right answer I chased two wrong ones, and both are
worth recording because they are ruled out by measurement: whether the anchor's flat map must be
rebuilt from the current `Phi_k` on every pass of the fixed point, as Step 2's script does and this
one did not, and whether the inverse map must clamp, as Step 2's script does, or continue, as
section 14 ruling 2 says. Four variants of the construction, held or rebuilt and clamped or
continued, agree to 0.001 m/s on the anchor's column, so neither is the difference. The check now
compares the construction with the stated values and reports the read-back beside it as what the
file's own grid costs, 0.174 m/s, which is decision Q's truncation at the wind.

**3. The lower edge test could not fire. Mine.** `trace` tested the geopotential edges with the
smallest overshoot over the levels instead of the largest, so a curve leaving the bottom of the
mesh was never detected while one leaving the top was. Found by re-reading the module before the
recorded run, not by a check, which is why check 6 now exercises both sides on constructed inputs
and says so in its output.

**4. Extending the mesh by the excess at the stopping node makes the pass repeat once per node.
Mine.** The excess a stopped trace reports is the overshoot at the one node it stopped on and says
nothing about how much further the curve climbs over the latitudes it has not reached. Under the
sheared wind the first pass extended **eleven times**, each repeat rebuilding every column, and the
pass cost 1,084 s on a 213 by 416 mesh. A second stop on the same side within a pass now adds twice
what the one before it added: six extensions and 766 s on the same case, four in check 6's loop.
The specification's own words are "the caller extends the mesh and repeats", which this is; the
doubling is how much it adds each time, and both the excess and the amount added are recorded.

Decision R, which arrived at v0.14 while this step's suite was running, says that no step of the
specification is to be optimized while it is being verified, so this change is reported for
ratification rather than presented as settled. It is not only a speed change: the mesh it leaves
is slightly larger, 418 geopotential nodes against 416 on the sheared case, because the last
extension overshoots what the curve needed. Every value in this report was measured with the
doubling in place. If the reviewing agent would rather the timid rule stood, the cost is that the
sheared run's first pass rebuilds every column eleven times instead of six, and no measured value
in section 3 moves except the extension list in check 6 and the mesh sizes the sheared runs report.

**5. Ruled on at section 15 ruling 3, a rule for Step 5's driver. A `lib.io.read` handle left
open segmentation faults the interpreter when the same file is opened again. Mine, and a hazard
for Step 5.** The first run of this script printed an attribute from `cio.read(path, kind)`
without closing the handle, and the third open of that file killed the process with SIGSEGV
rather than an exception. `lib.control` already loads every file through a read, load and close;
the script now does the same through one helper. Step 5's driver opens the anchors and the four
inputs and writes a product, so it is worth stating once: in this HDF5 build, a read that is not
closed is not merely untidy.

**6. Ruled on at section 15 ruling 4: `produce` takes the composition on the labels at Step 5,
one path for closure and transfer. `forward.production.produce` cannot produce at the target,
and Step 5 must give it the composition on the isobar labels.** `produce` takes the composition
column level by level, which is right when the levels are the anchor's own and wrong when they
are the union's: at M = 2 the union has 131 levels and the composition file 66, and the call
fails on the shape. Section 10 finding 5 already states the rule for transfer mode, that the
composition is read log-linearly in pressure onto each isobar's label. The acceptance forms `p`
and `T` through `lib.hydrostatic` with the composition on the labels, which is four lines of the
same arithmetic; Step 5 should add the argument to `produce` rather than leave two paths.

**7. Recorded for the combination specification at section 15 ruling 5. The map's seam at the
nearest-anchor boundary, measured.** Decision 2 above takes the nearest anchor's knots at each
latitude, so the map changes hands at the midpoint of two anchors' latitudes. The two anchors'
maps differ there by **7.943e-08 in `ln p`** at M = 2, which is a step in the map at one node.
Under the closure wind it reaches nothing: `du/dln p` is exactly zero, and the map's slopes
enter `S` only multiplied by it. Under a wind with vertical shear it would put a spike in `S` at
that latitude. The combination specification, which reconciles the anchors with the wind
hypothesis, is where a map that does not change hands abruptly belongs; this is recorded, not
solved.

**8. Convergence toward the stated values is not monotone, and check 3 reports it.** At 60 N the
departure from the stated `d ln N` rises from 9.133e-07 at the namelist's spacings to 4.483e-06
with both halved, and the shift from 0.0080 to 0.1779 percent. Both are inside the bounds. The
stated values are themselves a measurement, made with a different scheme integrated to 1e-12, so
the departure from them is the difference of two discretizations and need not fall monotonically;
the run at the namelist's spacings sits on a cancellation. No order is claimed anywhere in this
step.

## 5. Regression

The full regression runs every accepted suite on the registered products of the Step 3 acceptance
commit `bb75d16`. This step adds two modules and changes none, so no accepted behavior can move;
the run is the procedure's confirmation of that. This step's own suite is not in the loop, for the
reason decision 12 gives, and its row is the acceptance run recorded in
`reports/step04_4/output.txt`.

**Every accepted suite passes at its reference count.**

| Suite | Result | Suite | Result |
|---|---|---|---|
| `step1/accept_step1` | 6 of 6 | `step02_2` | 7 of 7 |
| `step1/verify_review_changes` | 14 of 14 | `step02_3` | 9 of 9 |
| `step2/accept_step2` | 8 of 8 | `step02_4` | 9 of 9 |
| `step3/accept_step3` | 5 of 5 | `step02_5` | 7 of 7 |
| `step4/accept_step4` | 7 of 7 | `step02_6` | 7 of 7 |
| `step5/accept_step5` | 7 of 7 | `step03_1` | 9 of 9 |
| `step6/accept_step6` | 6 of 6 | `step03_2` | 8 of 8 |
| `step7/accept_step7` | 8 of 8 | `step03_3` | 16 of 16 |
| `step8/accept_step8` | 9 of 9 | `step03_4` | 13 of 13 |
| `step9/accept_step9` | 6 of 6 | `step04_0` | 15 of 15 |
| `step02_1` | 6 of 6 | `step04_1` | 13 of 13 |
| | | `step04_2` | 13 of 13 |
| | | `step04_3` | 9 of 9 |
| | | `step04_4` | 11 of 11, the acceptance run |

The registered kind N is restored bit for bit, SHA-256
`920304440ee4554648c3aa6321b713741ee26eb0eb1ff223f5a81a8f9b20fd37` before and after, and the
closure inputs, the closure product and the transfer inputs are restored equal. The porcelain at
the end carries six paths, none of them written by a suite: this step's two new modules, this
report, and the author's SPEC_04 v0.15, REVIEW_04_step4 and `STATE.md`.

The loop was run three times, once after each filing of checks 8 and 10, and gave the same
twenty-four rows every time. The changes between them are `outer_loop`'s mesh argument
(decision 13) and the acceptance script's own construction; the first is additive and no accepted
suite reaches it, and the second is not library code at all. `regression.txt` holds the last run.

## 6. Hashes

Nothing is registered and no product changed, so there is no sweep. `reports/REPORT_04_step1.md`
section 6 and `reports/REPORT_04_step0.md` section 6 remain the record of the products' hashes.

## 7. Next step

Step 5: the production at the target, the altitude and the datum, kind `profile` in transfer mode,
`casspian-forward` in transfer mode, and F5, F6 and F7.

What this step leaves for it, beyond its own deliverables: the production at the target, which the
acceptance forms through `lib.hydrostatic` and which needs the composition on the isobar labels
(finding 6); the pressure identity of decision H, measured here at 5.789e-07 at the target and to
be reported in the product at every level of every run; and the `transfer_record` group, which
`outer_loop` fills with the pass count, the residual per pass and every extension with its excess
and the amount added.

**What SPEC_04 v0.14 changes for Step 5.** Decision R gives Step 5 two things to carry. The
production namelist takes a coarser geopotential spacing, the author's direction being that the
mesh is oversampled in geopotential by an order of magnitude, with the expected values remeasured
there; on this step's own evidence the geopotential spacing is indeed not what limits the transfer,
since halving it along with the latitude spacing moved the departure from the stated values by a
factor of 1.5 at 10 N and the wrong way at 60 N (check 3, finding 8), while the latitude spacing at
0.5 degrees costs 2.0e-04 in `d ln N` against 6.8e-06 at 0.05. And the reference-surface march step
stops being the 0.05 degree default argument in `lib.geoid`: either a key of
`[numerics.reference_surface]` beside the scheme, or the wind file's latitude nodes united with the
mesh's and the anchor's. The second reaches this step's cost, though not its values:
`through_anchor` marches on its own dense grid, `arange` over both poles at `march_step_deg`,
unioned with the caller's nodes, so every pass of every run here marched about 3,600 latitudes
whatever the mesh's
spacing was, and the two coincide only because the namelist's spacing happens to be the same 0.05
degrees as the default. At the coarse spacings the constructed cases use, the march is the same
3,600 nodes as at the fine ones.

**What v0.15 adds for Step 5**, from the section 15 rulings: `produce` takes the composition on the
isobar labels as an argument, one path for closure and transfer, which finding 6 measured the need
for (ruling 4); the driver closes every file it reads before opening another (ruling 3); and the
production namelist's geopotential spacing is 50,000 m2/s2, with Steps 0 to 4 recorded as run at
5,000.

One thing is open and one is closed. **Open:** nothing in this step. Section 16 settled the last
of it, and what it leaves for Step 5 is one number to carry: the product's `estimate` group
states the identity floor once as `identity_floor`, because `D_12` below 6.5e-08 says nothing
about the atmosphere (section 16 ruling 2).

**Closed:** the transfer itself. Under the closure wind it reproduces an independent implementation
to 9.1e-07 in `d ln N` over 29 degrees of latitude, three of four levels to 2e-07; out and back
under a sheared wind it is reversible to 4.1e-13; the gauge isobar does not move at any latitude;
and at M = 2, with no file between the trace and the estimate, it is an identity to 3.6e-15.
