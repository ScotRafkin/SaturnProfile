# REPORT 07, pre-execution. The coding agent's reading of SPEC_07 v0.2

Coding agent, 3 October 2026. Specification: `docs/specs/SPEC_07_Wind_Interpolant.md` v0.2, §4.
Working tree on `main` at `a238306` (SPEC_06 closed, its run files swept). Measurements by
`tests/step07_1/preexecution_measure.py`, which builds decision L2 from
`scipy.interpolate.PchipInterpolator` exactly as §1 describes (latitude first on every pressure
row, then `ln p` at the point) and sets it beside decision L; nothing in `src` is used for L2.
Output: `reports/step07_1/preexecution.json`. "Every wind file" is the 27 kind W files in the tree:
`occul_data/lindal/lindal_wind.nc`, the two registered runs' source and run winds, and the eleven
experiments' source and run winds.

**In short.** §1 gives what it claims. No wind file in the tree has a turnover locus. The predictor
is small where the runs read the wind (1.04 m/s at most between 9 and 32 deg, -0.23 m/s at `phi_c`,
an estimated -5 m in `r0` at 10 N). Six accepted suites form an expectation by the linear rule and
four more compare against committed hashes or fixtures that the full sweep will move; they are
named in item 4. Three points are added for the ruling: the isobar map's knots as quadrature
breakpoints, the hash tables, and the cost's baseline.

## 1. The construction

**It gives what §1 claims.** For a fixed latitude the profile in `x = ln p` is one PCHIP of the
values `v_k(phi)`, so `(du/dln p)` is continuous in `x` everywhere; the `v_k(phi)` are C1 in
latitude, so `(du/dln p)` is also continuous across latitude nodes. `(du/dphi)` is the derivative
of the Hermite form with `v_k'` and `dm/dphi`; it is continuous wherever the node slopes `m_k` are
C1 in latitude, which fails only where the harmonic mean switches branch: an adjacent pair of
vertical secants changing from one sign to opposite signs, or one of them passing through zero.

Two properties follow that the specification relies on and I confirm:

- **Node values bit for bit.** PCHIP in local form (scipy's piecewise polynomial in powers of
  `x - x_k`, constant term the node value) returns the node value exactly at `t = 0`; scipy takes
  the right-hand interval at an interior node, where `t = 0` too.
- **A column with no vertical shear is exactly constant.** Identical rows give identical `v_k(phi)`
  bit for bit (each row's PCHIP is the same arithmetic on the same numbers), every vertical secant
  is exactly zero, so every slope is zero and every higher coefficient is zero: the value is
  `v_k(phi)` at every pressure and `(du/dln p)` is exactly zero, which check 5 needs. This holds
  only in local form; a Hermite basis written as `h00 v_k + h01 v_(k+1)` would not be bit-exact,
  since `h00 + h01` is not exactly one in floating point.

**The turnover loci, measured.** On each file the latitude-interpolated rows were sampled 20 times
per latitude cell and every interior pressure node tested for a change in the sign agreement of
its two vertical secants:

| Files | Loci found | Where | Largest jump of `(du/dphi)` | Decision L's largest jump across a latitude node |
|---|---|---|---|---|
| the 20 with no vertical shear (the Lindal wind, both registered runs' source and run winds, every experiment's source, runs 2, 3a, 3b and 3c) | none | | | 0 to 4.86e3 m/s per rad |
| runs 4 to 9, 7f (sheared) | 17 to 187 each | every one at a pole, +90 or -90 deg | 0.0107 m/s per rad, finite-difference noise | 4.86e3 to 7.29e3 m/s per rad |

Every locus found is at a pole, where kind W holds the wind at exactly zero, so the secants there
are exactly zero and the rule's zero branch is taken; the wind is zero on both sides, and the
measured jump is the finite difference's noise. **No file has a genuine turnover locus.** The
shear cases cannot have one: their wind is `u_s(phi) F(p)`, so every vertical secant is
`u_s(phi)` times a fixed number, and two adjacent secants agree in sign at every latitude where
`u_s` is not zero (PCHIP is homogeneous, so the rows' interpolants inherit the product form). A
locus needs a wind whose vertical profile turns over at a node somewhere and not elsewhere; none
of the present files does. Check 4's reported jump will be the polar zero, and its comparison
against decision L's 4.86e3 is not informative on these files. I suggest a synthetic file with a
turnover (a sign change of `du/dln p` crossing a pressure node along latitude) in check 4, so that
§1's predicted jump is measured once.

## 2. The predictor

The largest `|u_L2 - u_L|` on 8 samples per cell along each axis:

| Files | Largest `|u_L2 - u_L|` | Where |
|---|---|---|
| sources, closure, transfer, runs 2, 3a, 3b (no shear) | 5.344 m/s (runs 3a and 3b's own winds: 0) | -8.84 deg, every pressure (no vertical structure) |
| run 3c (half) | 2.672 m/s | the same place |
| run 4 | 5.344 m/s | the same place |
| run 5 | 5.634 m/s | -7.22 deg, 930 mbar |
| run 6 | 8.764 m/s | -7.22 deg, 930 mbar |
| run 7 | 16.014 m/s | -7.22 deg, 930 mbar |
| run 7f | 9.793 m/s | -7.22 deg, 965 mbar |
| run 8 | 6.706 m/s | -8.84 deg, 1075 mbar |
| run 9 | 8.069 m/s | -8.84 deg, 1075 mbar |

The largest differences are in latitude, at the equatorial jet's sharpest bins south of the
equator, and in the sheared runs where the latitude difference meets the corner at `p_s`. **Where
the runs read the wind** (9 to 32 deg) the largest is 1.04 m/s (10.35 deg), the mean `|du|` 0.17
m/s and the mean `du` -0.03 m/s, of either sign. **At `phi_c`** (30.806 deg), on the reference level
of every file whose source is the Lindal wind: -0.2272 m/s.

**What it predicts for check 6.** The geoid's rise along the gauge isobar from `phi_c` to 10 N,
`int (2 Omega u r sin(phi) + u^2 tan(phi)) dphi`, changes by -48 m2/s2 with L2, about -5 m in
`r0(10 N)` (bound 100 m). The reduction's `phi_c` moves through the wind at `phi_c` only: -0.23
m/s of the 2.17 m/s there, a tenth of a wind that is itself small (`u / (r cos phi_c)` is about
4e-8 rad/s against `Omega` 1.6e-4); I expect far below 1e-4 deg but have not estimated it. The closure and
transfer temperatures are expected to move by hundredths of a kelvin.

## 3. The experiments' grids

The wind grid is ten nodes per decade (`d ln p` = 0.230) in every run but 7f (twenty, 0.115). The
anchor levels are about 0.046 apart in `ln p` near 1 bar and about 0.19 to 0.23 apart above 10
mbar. The anchor levels inside each wind interval adjacent to the breaks:

| Run | `p_s` (on a node) | stop pressure |
|---|---|---|
| 2, 3a to 3c | 1000 hPa: 5 levels in the interval above, 5 below (no break: uniform) | none |
| 4 | 1000 hPa: 5 and 5 | 20 Pa, the column's top, between nodes: 1 level in 20 to 25 Pa |
| 5 | 1000 hPa: 5 and 5 | 700 Pa, between nodes 631 and 794 Pa: 1 level (632.3 Pa, at the lower end) |
| 6 | 1000 hPa: 5 and 5 | 8000 Pa, between nodes 7943 and 10000 Pa: 2 levels |
| 7 | 1000 hPa: 5 and 5 | 0 Pa (linear in `p` to zero; no node break above) |
| 7f | 1000 hPa: 3 above, 2 below | as 7 |
| 8, 9 | 1000 hPa: 5 and 5 | 1 MPa, beyond the column's bottom (1.294 bar): no level |

**What it predicts for behaviour 5.** Under L2 the shear's change at a node is rounded over the
interval on either side, so the temperature break that SPEC_06 check 7 measured at one layer is
spread over the anchor layers those intervals hold. At `p_s` (5 and 5 layers, 3 and 2 in 7f) the
identity should fall several times; at run 5's stop pressure (1 layer) and run 6's (2) much less;
run 4's stop is at the column's top. The largest identities of SPEC_06 were run 5 at its stop
pressure (1.5e-2) and runs 7, 7f at `p_s` (1.1e-2, 1.6e-2); I expect 7 and 7f to fall well below
1e-2 and run 5 to remain the largest.

## 4. The accepted suites

**Forming an expectation by the linear rule** (to be restated to the interpolant's own, §3
Regression):

- `step02_5` (line 171): `u_direct` is `np.interp` of the reference row at `phi_c`, set against the
  model's value.
- `step04_1` check 12: `wind_derivatives` constant across a cell and its one-sided values at a node;
  and the sheared check near line 385, whose analytic `du/dln p` takes the reference row by
  `np.interp` in latitude.
- `step04_1` near line 220: `reference_wind` against `wind_of_latitude` bit for bit, which holds by
  construction once `wind_of_latitude` returns `reference_wind` (no edit expected).
- `step04_2` (lines 363, 440), `step04_3` (424, 481), `step04_4` (498) and `step04_5/fields.py`
  (111): the cylinder-extended synthetic wind is built from the reference row by `np.interp` in
  latitude. That is the test's own construction of a wind file, not a read, so it need not change;
  but where the suite then compares the model's reading of that file with the construction between
  nodes, the comparison is to the linear rule. To be checked suite by suite at execution, against
  decision Q's floor at which these are bounded.
- `step04_3` (line 177) forms the closed form from `wind_derivatives` itself and needs nothing (check 5).

**Compared against committed hashes or fixtures that the full sweep moves** (the reduction's wind
moves `phi_c` and so every chain file's bytes):

- `step02_1` and `step02_6` read the hash table of `reports/REPORT_05_step0.md`. After the full
  sweep every chain file's hash changes and both fail unless a new table is recorded and the two
  suites pointed at it, as SPEC_05 Step 0 did. The specification should say where the new table
  lives (this step's report).
- `step03_3` compares against `tests/step03_3/fixtures/` and `step04_0` against
  `tests/step04_0/fixtures/swept/` (read-only fixtures from SPEC_03 and SPEC_04). Their comparisons
  of reduction products will see the sweep's changes; their bounds are to be read at execution, and
  a fixture comparison that is exact is a ruling for the review (restate, or regenerate the fixture).
- `step05_2` compares the 18 registered files built as candidates against the committed ones; it
  passes after the sweep's commit, as under SPEC_06.

**Stated tighter than check 6's bounds:** `step04_5` check 11 (the closure production at 1e-14
relative to the registered product) passes after the sweep, since it rebuilds against the swept
product; `step04_5`'s at-anchor identity is check 6's own. No other bound tighter than check 6's
reaches a value L2 moves.

## 5. Other interpolations of the wind

- `tools/plots/figures_product.py` (lines 55 and 63): the product figures interpolate the wind
  themselves, `np.interp` in `ln p` and in latitude. They draw the wind, so they should read it
  through `WindField`; otherwise the figure shows decision L while the model uses L2.
- `tools/plots/figures_inputs.py` (line 159): the inputs figure takes a profile by `np.interp` in
  latitude. Same remark.
- `tools/wind/shear.py`: `u_at_pressure` reads through `WindField`, so it moves with L2 when `p_s`
  is not a wind node; in every present experiment `p_s` is the 1000 hPa node, so no shear file
  changes.
- `tools/wind/curve.py` and `fit.py` build wind files from data (already PCHIP in places); they do
  not read a kind W file. `lib/geoid.py` interpolates radii, not the wind.

## 6. The derivative in latitude by the tangent

**Cost.** The tangent adds, per point, the four rows' latitude derivatives (scipy's derivative of
the rows' PCHIP, formed once) and the partials of two harmonic means: a few tens of operations on
arrays already in hand, about half again the cost of a value. A centered difference costs two more
full evaluations and gives about half the digits, and at a node it straddles the corner the
tangent resolves. The tangent is the better choice.

**What I would do.** The latitude step by `scipy.interpolate.PchipInterpolator(phi, u, axis=0)`
itself, formed once, with its `derivative()`: check 1 then holds by construction. The `ln p` step
by hand on the four rows, in local form, with `dm/dphi` closed form; on the harmonic mean's zero
branch `dm/dphi` is zero, consistent with `m` there. The end rule's clamps (scipy's three point
rule, clipped where it changes sign or exceeds three secants) apply only at the grid's top and
bottom rows (1 Pa and 1 MPa), outside the model's range, and the tangent follows the branch taken.

**The cost of the run.** A transfer run now spends about 3.3 s of 20 s (profiled) reading the wind,
in 71,000 small calls, most from the tracing's quadrature at a mesh latitude or a midpoint between
two. L2 may cost three to five times as much per call. Since the tracing reads the wind at only
the mesh's latitude nodes and midpoints, the latitude step's rows `v_k(phi)` can be held per
column, which removes most of the cost; I would do that if check 9 needs it. The baseline is
better taken as SPEC_06's (the Lindal transfer run 14 s, run 8 16.5 s) than SPEC_05 Step 2's 12.7 s.

## 7. Added for the ruling

1. **The isobar map's knots as quadrature breakpoints.** SPEC_06's column integral breaks at the
   wind file's pressure nodes, where the wind's derivatives change. But the integrand is read
   along geopotential, and `ln p(Phi)` in a column is the map's piecewise linear interpolation
   through the anchor levels, so `S/g` also changes slope at every map knot. A two-point Gauss
   rule across such a knot is not exact. Under decision L this was masked by the larger breaks at
   the wind nodes; under L2 it may be what check 8 measures. Adding the column's map knots to the
   breakpoints costs about one more piece per anchor layer. I suggest check 8 run three ways:
   two points, four points, and two points with the map knots added.
2. **The hash tables** (item 4): where the new table is recorded, and `step02_1` and `step02_6`
   pointed at it.
3. **The cost's baseline** (item 6): SPEC_06's run times rather than 12.7 s.
