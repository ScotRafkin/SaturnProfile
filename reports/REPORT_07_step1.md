# REPORT 07, Step 1 (interim). The wind interpolant with a continuous slope

Coding agent, 3 October 2026. Specification: `docs/specs/SPEC_07_Wind_Interpolant.md` v0.3, §3, with
the rulings of §7. Working tree on `main` at `a238306`; nothing of this step is committed.
**Interim: the step stopped at check 6**, as the regression paragraph directs when a value moves by
more than check 6's bounds: `phi_c` moved by 2.4e-4 deg against 1e-4. The sweep and the full
regression were not run. Filed for a ruling, with what the other checks found. Section 7 adds the
dense reference of v0.4 §8 ruling 5, for the author's decision between L2 and L.

**Result.** `tests/step07_1/accept_step07_1.py`: 6 of 9 as run, 7 of 9 with check 4's measurement
corrected (section 2). Decision L2 is implemented as §1 states it and does what it promises: the
rule is scipy's to 4e-16, every node is exact, no value leaves its cell's corner values, both
derivatives are continuous, a wind with no vertical shear has a vertical derivative of exactly zero,
the cost is 1.4 to 1.5 times SPEC_06's, and the quadrature is adequate to 0.03 K. Three findings
need the ruling: `phi_c` (section 3), the rounding of a corner at a node, which moves shear inside
the interval and so moves the delivered temperature by tens of kelvins at the corner (section 4), and
the identity, which falls five to six times at `p_s` but rises at two stop pressures (section 5).

## 1. What changed

| File | Change |
|---|---|
| `src/casspian/lib/windfield.py` | Decision L2. The latitude step is `scipy.interpolate.PchipInterpolator(phi, u, axis=0)`, formed once with its `derivative()`; the `ln p` step by hand on the four rows about each point's interval, in local form, with scipy's slope rule (`_find_derivatives`, `_edge_case`) and coefficients (`CubicHermiteSpline`); `(du/dphi)` by the tangent, the harmonic mean's partials in closed form and zero on its zero branch. The last node of either axis returns the file's value exactly (the last interval's polynomial reaches it only to rounding). Rows are evaluated once per distinct latitude in a call; no cache across calls (check 9 did not need one). The API and the refusal outside the grid are unchanged. |
| `src/casspian/refrac/anchor.py` | `wind_of_latitude(wind)` returns `WindField(wind).reference_wind`. |
| `src/casspian/tools/plots/figures_product.py`, `figures_inputs.py` | The figures read the wind through `WindField` (deliverable 5). |
| `src/casspian/lib/kernel.py`, `tools/plots/figures_profile.py`, `tools/wind/shear.py` | Docstrings that stated the linear rule. |
| `tests/step02_5/accept_step02_5.py` check 6 | The direct `u` at `phi_c` is scipy's PCHIP of the 1 bar column, formed in the suite (§7 ruling 4). |
| `tests/step04_1/accept_step04_1.py` checks 12 and 13 | Check 12: the one sided slopes at a node agree with each other and with scipy's PCHIP derivative. Check 13: the sheared field's `du/dln p` is the interpolant's own formed independently by scipy (rows in latitude, then `ln p`), since PCHIP is not linear in the data and the shear term's PCHIP alone is not the answer. Neither suite has been run yet. |
| `tests/step07_1/accept_step07_1.py`, `preexecution_measure.py` | The acceptance and the reading's measurements. |

The callers that reach the new interpolant: `tools/wind/shear.py`, `lib/kernel.py`,
`forward/transfer.py` (the column, `_wind_table`, the geoid gauge wind, the columns, the production),
`forward/production.py` (through `wind_of_latitude`), `refrac/anchor.py`, and the figures.

## 2. The checks

| Check | Measured | Result |
|---|---|---|
| 1. the rule | every row and column of the 16 wind files against scipy: values and both derivatives within 4e-16 of the file's largest | pass |
| 2. nodes exact | every node of every file, and `reference_wind` at every latitude node, bit for bit | pass |
| 3. no overshoot | 20 points per cell per axis: no value beyond its cell's corner values (0.0 in every file) | pass |
| 4. continuous slopes | as run, **fail**: across latitude nodes `(du/dphi)` differed by 3.8e-8 of its largest, against 1e-9. That was the measurement: the one sided values were read 1e-10 from the node, and PCHIP's second derivative jumps there, so the difference is that jump times the distance. It falls tenfold with each tenfold smaller step (3.8e-10 at 1e-12, 3.8e-11 at 1e-13; the same in every file and in `ln p`); the script now reads 1e-12 from the node, where every file passes. Centred differences agree to 5e-9. No file has a locus off the poles. The synthetic turnover file has 12; the largest jump of `(du/dphi)` is 13.3 m/s per rad at 11.09 deg and the 2512 Pa node, against decision L's 64.9 on the same file, and its ratio to the rows' difference in `(du/dphi)` across the two intervals is 0.500, not §1's "of order one third" | pass, corrected |
| 5. no vertical shear | under the closure wind `(du/dln p)` exactly 0 on check 3's sample | pass |
| 6. the chain before the sweep | **`phi_c` +2.385e-4 deg (bound 1e-4)**; `r0` -22.9 m; polar radii +2.6 and +2.1 m, mean +2.3 m, asymmetry -0.5 m; closure product T 1.1e-3 K, `ln N` 0, pressure 7.2e-6, radius 22.9 m; transfer product T 5.1e-3 K, `ln N` 3.7e-5, pressure on the labels 7.2e-6, radius 15.4 m, altitude 12.9 m | **fail: stop** |
| 7. no steps | section 4 and 5 | **fail** |
| 8. the quadrature, three ways (run 7f at 5e4) | four Gauss points: 0.020 K from the model as delivered; two with the map's knots: 0.029 K | pass |
| 9. the cost | the Lindal transfer run 20.1 s (14.0), ratio 1.44; run 8 25.2 s (16.5), ratio 1.53 | pass |

Check 8 settles §7 ruling 7: the map's knots move the answer by 0.03 K at most, so they are not
needed. (The two variant runs took 14 and 21 minutes because the test's own quadrature is not
vectorised as the model's is; that is the test, not the model.)

## 3. `phi_c` moved by 2.4e-4 deg

Everything else in the chain moved by far less than its bound; `phi_c` by 2.4 times its bound. The
predictor accounts for it. At `phi_c` the reference wind changes by -0.227 m/s out of 2.17 m/s
(the anchor lies where the wind is weak and changing, so a tenth of it moves). The reduction turns
the occultation's geometry into a planetocentric latitude through the geoid's slope at the anchor,
and that slope carries the term `2 Omega u sin(phi) / g`: with `Omega` 1.638e-4, `Delta u` 0.227
m/s, `sin(phi_c)` 0.512 and `g` about 10.4 m/s2, `2 Omega Delta u sin(phi_c) / g` is 3.7e-6 rad,
2.1e-4 deg, against 2.4e-4 measured. Nothing in the reduction changed; the wind it reads at the
anchor changed by the difference of the two interpolants there.

`r0` moved 22.9 m and the closure product's radii the same; the products' temperatures moved by
1e-3 to 5e-3 K and their pressures by 7e-6. The bound of 1e-4 deg was set before the predictor at
`phi_c` was known; 2.4e-4 deg is 250 m along the surface.

## 4. A corner at a node: the shear moves inside the interval

Where the hypothesis has a corner exactly at a wind node (every experiment's `p_s`, at the 1000 hPa
node), PCHIP's slope at that node is zero: the secant on one side is the zone's and on the other
zero, so the harmonic mean takes its zero branch. Between that node and the next, the cubic leaves
the node with zero slope and must still end at the next node's value, so it is steeper in the middle
of the interval than the secant. The wind stays within its corner values (check 3), but the
vertical shear inside the interval is moved: smaller near the corner, larger in the middle. The
delivered temperature at 10 N integrates that shear along the isobars, and the anchor levels inside
the interval read it. Run 7f at 5e4 (the change from run 2, K):

| Level | Label (hPa) | SPEC_06 (L) | now (L2) |
|---|---|---|---|
| 55 | 831.4 | -69.99 | -69.40 |
| 56 | 870.8 | -69.84 | -71.29 |
| 57 | 911.5 | -74.39 | **-85.04** |
| 58 | 953.5 | -73.10 | **-80.37** |
| 59 | 998.7 | -71.25 | **-4.68** |
| 60 | 1045.4 | +2.22 | -0.09 |

Level 59 lies 0.13 percent above the 1000 hPa node, inside the interval where L2 brings the shear to
zero, so it now reads almost no shear and almost no change; levels 57 and 58 read more than the
zone's. The same in every decay run: the 999 mbar isobar of run 7f now carries `S/g` near zero
along its whole path (`reports/step07_1/F8_isobar_L_against_L2.png`), where under decision L it
carried the zone's shear. The delivered profile's slope at layer 56 to 57 is -261 K per unit `ln p`
against its neighbours' 0.65 and 144.7 (check 7's continuity test, the one layer of 16 that fails;
under SPEC_06, 9 of 16 failed); the delivered temperature at level 57 is 44.3 K.

This is behaviour 4, "corners rounded where the hypothesis has them, within one node interval either
side, by design", but its size was not foreseen. In every decay run the level just above `p_s`
(998.7 mbar) now reads almost none of the zone's change: run 4 -9.90 to -0.27 K, run 5 -16.77 to
-0.48, run 6 -32.08 to -0.80, run 7 -68.65 to -2.41, run 7f -71.25 to -4.68. And the levels in the
interval above it read more than the zone: run 5 -22.4 K at level 56 against the zone's -17.2, run 6
-42.8 against -33.7, run 7 -86.8 against -72.8. The shear integrated across the interval is unchanged, as
§1 says; what the anchor levels inside it read is not.

**Check 7's first part** (run 7f's isobar ten times smoother than SPEC_06's) is therefore measured on
an isobar whose shear L2 has taken away: the largest step is 1.45e-2 against 5.89e-2 (4.1 times), and
the mean step at the wind's latitude nodes against elsewhere is 1.77. On this isobar the comparison
does not measure smoothness. An isobar away from the corner would.

**Away from the slope changes** (check 7's fourth part): runs 2, 3a to 3c within 5e-3 K of SPEC_06
(pass); in the others the levels farther than one wind interval from `p_s` or the stop pressure
differ by 0.10 to 0.18 K (runs 8, 9), 0.7 to 0.8 K (runs 4, 6), 1.8 K (run 5), 2.4 K (run 7f) and 5.4
to 7.2 K (run 7), against 0.5 K. Two causes, both measured. In run 5 the difference decays with
depth as `1/p` (-1.79, -1.33, -1.00, -0.74 K at 10, 12.6, 15.9, 19.9 mbar): it is the production's
absolute pressure offset from the corner's layer, carried down (SPEC_06 check 7), and it changes
because the identity there changes. In run 7 (`u` linear in `p`, so curved in `ln p` everywhere in
the zone) the two interpolants differ between every pair of nodes, up to 16 m/s (the predictor); that
is behaviour 2's "small numerical differences of the size the two interpolants differ", not small
here because the hypothesis is curved on a coarse grid.

## 5. The identity

| Run | SPEC_06, 5e4 / 2.5e4 | now, 5e4 / 2.5e4 (level) |
|---|---|---|
| 2, 3a, 3b, 3c | 5.4e-7 to 5.8e-7 | 5.4e-7 to 5.8e-7 (unchanged) |
| 4 | 1.64e-3 / 1.64e-3 | **7.06e-3 / 7.09e-3 (1)** |
| 5 | 1.51e-2 / 1.51e-2 | 8.29e-3 / 8.31e-3 (17) |
| 6 | 6.41e-3 / 6.40e-3 | **2.38e-2 / 2.53e-2 (28)** |
| 7 | 1.10e-2 / 1.08e-2 | 2.24e-3 / 2.28e-3 (58) |
| 7f | 1.64e-2 | 3.04e-3 (57) |
| 8 | 1.45e-3 / 1.45e-3 | 2.50e-4 / 2.50e-4 (64) |
| 9 | 2.89e-3 / 2.89e-3 | 5.01e-4 / 5.02e-4 (64) |

Every run is mesh independent (the two spacings within 6 percent, run 6; 2 percent or less in the
others). At `p_s` the identity falls five to six times (runs 7, 7f, 8, 9), as behaviour 5 predicted
for intervals holding five anchor layers. At run 5's stop pressure it nearly halves. **At run 4's
and run 6's stop pressures it rises**, four times, against behaviour 5 ("no larger in any run").

The layer law of SPEC_06 check 7 reads it. Run 6's stop pressure, 8000 Pa, lies between the nodes
7943 and 10000 Pa; under decision L the file's break sat at the 7943 Pa node, at `s` = 0.495 in the
anchor layer from 7243 to 8728 Pa, almost at the law's zero crossing, and that layer gained +27 Pa.
Under L2 the change is spread from that node into the interval below it, so the layer's temperature
change lies mostly in its lower half: the layer now gains +208 Pa, against a law scale `M |jump in ln
T|` of 741 Pa (an effective `s - s0` of about 0.28). Run 4's stop pressure (20 Pa) lies in the first
interval of the column's top, and its largest identity is now at level 1 there. The identity at a
stop pressure is the anchor layer's integration of a temperature change that L2 moves within the
layer; whether it rises or falls depends on where the change lands, which is SPEC_08's.

## 6. For the ruling

1. **`phi_c`'s bound.** The move is the predictor's: a tenth of the weak wind at the anchor, 2.1e-4
   deg estimated, 2.4e-4 measured, `r0` and both products far within theirs. Restate the bound (for
   example to what the predictor gives), or hold L2 back at the anchor.
2. **Corners at nodes.** PCHIP's zero slope at a corner moves the shear inside the adjacent interval,
   and anchor levels there read tens of kelvins differently; the level 0.13 percent above `p_s`
   reads almost no shear at all. That is shape preservation doing what it does on a corner, and it
   is what the hypothesis on its nodes allows. Options: accept it as behaviour 4 with its size now
   known; write the shear cases with nodes closer around their corners (§5 decision 3); or a
   different node slope at a corner (not PCHIP's), which is a change of interpolant.
3. **Check 7 as written** cannot pass under 2: its isobar test is taken inside a rounded corner, its
   continuity test meets the rounding at one layer, its "not above SPEC_06" meets the stop pressures,
   and its 0.5 K away from corners meets the identity's offset carried down and the curved hypothesis
   of run 7. Each part is measured above.
4. **The rest stands:** checks 1 to 5, 8 and 9. The map's knots are not needed (0.03 K).

Nothing is committed. The experiment products under `forward/` and `reports/step05_4/` are now L2's
(ignored files); SPEC_06's are kept under `reports/step07_1/spec06/`.

## 7. The dense reference (SPEC_07 v0.4, §8 ruling 5)

`tests/step07_1/dense_reference.py`. Runs 5, 7 and 7f were rebuilt with the same hypothesis on ten
times the pressure nodes: 601 nodes for runs 5 and 7 (100 per decade) and 1201 for run 7f (200 per
decade). Every coarse node is among them, and 1000 hPa is a node. They ran at 5e4 with figures off,
in copies under `reports/step07_1/dense/`. Under L2 they ran from the working tree (109, 166 and
402 s). Under L they ran from a worktree at `a238306` (25, 32 and 39 s), outside the repository,
whose stamps are clean. The `-dirty` refusal was relaxed in those processes only. No code changed.

Every value is the delivered temperature change from run 2 under the same reading, at the 66 anchor
levels. The reference is dense L2. Dense L minus dense L2 is the reference's own uncertainty. The
corner intervals are the coarse wind intervals on either side of each slope change: `p_s` in every
run, and in run 5 also the stop pressure (700 Pa). Results are in `comparison.json` and
`comparison.txt`, with the logs in `run.txt`.

### Away from the corners

Largest (mean) magnitude, K:

| Run | coarse L | coarse L2 | dense L (uncertainty) |
|---|---|---|---|
| 5 | 1.08 (0.11) | 0.71 (0.07) | 0.54 (0.01) |
| 7 | 4.39 (0.83) | 1.21 (0.11) | 3.02 (0.17) |
| 7f | 2.26 (0.55) | 0.31 (0.03) | 2.76 (0.23) |

The uncertainty column has two parts, and they must be read separately:

- **Inside the zone** (above the corner intervals), dense L and dense L2 agree to 0.49 K in run 7
  (level 45; 0.3 K elsewhere), 0.24 K in run 7f and 0.01 K in run 5. There the reference holds.
  - **Coarse L** is off by up to ±4.4 K in run 7 and ±2.2 K in run 7f, alternating in sign below
    about 200 hPa. This is the staircase of a linear reading in `ln p` on a hypothesis that is
    linear in `p`.
  - **Coarse L2** is within 0.15 K down to 630 hPa in both runs. The rest of its error is next to
    the corner: 0.4 to 1.2 K in run 7 (levels 50 to 54, 660 to 793 hPa) and 0.31 K in run 7f.
  - **Run 5** (linear in `ln p`, so L is exact between its corners): both coarse readings are
    within 0.01 K except just below the stop corner. There, at level 17 (10 hPa), coarse L is 1.08 K off and
    coarse L2 0.71 K off.
  - **On run 7's 7.2 K question:** L is the one in error. Coarse L minus coarse L2 there is up to
    5.4 K (level 54), and nearly all of it is L's.
- **Below `p_s`** (levels 60 to 65, a region with no shear), the dense readings themselves differ by
  2.5 to 3.5 K (runs 7 and 7f) and 0.54 K (run 5), so the reference does not hold there.
  - The change there is a uniform offset, carried down from the corner. It tracks the identity.
  - Dense L2 gives -1.06 K in run 7 against -0.18 K in run 7f: the same hypothesis on two dense
    spacings.
  - Coarse L2 (+0.8 K from dense L2 in run 7, +0.09 K in run 7f) lies closer to dense L2 than
    either L reading does. No reading settles that region.

### Inside the corners

Corner levels, dense L2's change; then coarse L, coarse L2 and dense L minus it, K:

| Run | Level (hPa) | dense L2 | coarse L | coarse L2 | dense L |
|---|---|---|---|---|---|
| 5 | 15 (6.3) | +1.30 | -8.24 | -0.27 | -0.00 |
| 5 | 16 (8.0) | -15.69 | +1.48 | +4.09 | -0.01 |
| 5 | 55 to 58 (831 to 954) | -17.0 to -17.2 | +0.01 | -5.17 to +4.84 | -0.00 |
| 5 | 59 (998.7) | -4.12 | -12.64 | +3.65 | -12.66 |
| 7 | 55 to 58 | -69.1 to -73.3 | -3.83 to +2.63 | -15.71 to +16.46 | +0.20 to +0.30 |
| 7 | 59 (998.7) | -22.27 | -46.38 | +19.86 | -51.19 |
| 7f | 57, 58 | -72.4, -73.3 | -1.95, +0.15 | -12.60, -7.12 | +0.12, +0.06 |
| 7f | 59 (998.7) | -38.10 | -33.15 | +33.42 | -35.50 |

Levels 60 to 64 carry the offset below `p_s` discussed above; level 14 agrees to 0.00 K throughout.

- **Inside the corner intervals, above `p_s`** (levels 55 to 58, and 15 and 16 in run 5), the
  reference holds: dense L and dense L2 agree to 0.3 K.
  - **Coarse L** is closer there: within 4 K in run 7 and 2 K in run 7f. The exception is level 15
    in run 5, where L is 8 K off.
  - **Coarse L2** is off by up to 16.5 K (run 7) and 12.6 K (run 7f). This is ruling 4's
    overshoot, now measured against the hypothesis on a finer grid.
- **At level 59** (998.7 hPa, 0.13 percent above `p_s`), the reference does not exist.
  - Even on the dense grid, the level lies inside the last dense interval above `p_s` (977 to
    1000 hPa), so the corner's rendering decides it. Dense L and dense L2 differ there by 12.7 K
    (run 5), 51 K (run 7) and 36 K (run 7f).
  - This is ruling 4: the hypothesis's own temperature is a step, and a level that close to the
    step reads whatever the rendering puts there. It is not resolved by more nodes, only moved.

### The identities

Largest `|p/p_label - 1|` (level):

| Run | coarse L | coarse L2 | dense L | dense L2 |
|---|---|---|---|---|
| 5 | 1.51e-2 (16) | 8.29e-3 (17) | 2.92e-3 (60) | 2.49e-3 (59) |
| 7 | 1.10e-2 (60) | 2.24e-3 (58) | 1.81e-2 (60) | 1.22e-2 (59) |
| 7f | 1.64e-2 (60) | 3.04e-3 (57) | 2.04e-2 (60) | 9.55e-3 (59) |

At `p_s`, the dense grids raise the identity.
- **Runs 7 and 7f:** dense L2 is four and three times coarse L2.
- **Why:** the dense wind interval (0.023 in `ln p`) is half an anchor layer (0.046). The corner's
  temperature step now falls inside a single anchor layer rather than across five, and the layer law
  gives the offset by where in the layer it falls. This is behaviour 5 as ruling 7 corrected it.
- **Run 5:** the largest identity moves from the stop pressure (coarse) to `p_s` (dense). The dense
  wind nodes around 700 Pa spread the stop corner over more anchor layers.

The identity at a corner is the anchor's spacing against the corner's width, SPEC_08's subject. It
is not a property of either reading.

### For the author's decision (ruling 8)

1. **Away from corners and inside the zone, L2 is the closer reading.**
   - Coarse L2 is within 1.2 K of the reference in run 7 and 0.31 K in run 7f, mostly next to the
     corner and 0.15 K or less elsewhere.
   - Coarse L is off by up to 4.4 K and 2.2 K, alternating level to level. The 7.2 K of ruling 5 was
     L's staircase.
2. **In the corner intervals, L is the closer reading.**
   - Coarse L2's overshoot reaches 16.5 K against the hypothesis on a ten times finer grid.
   - L's step reaches 4 K, with run 5's level 15 the exception.
   - Ruling 8's remedy for L2 removes that region from realistic cases: corners are written with a
     width of a few node intervals.
3. **At `p_s` itself and below it, neither reading has a reference.** The level next to the corner
   and the offset carried below it depend on the rendering and on the spacing, dense or coarse.
4. **For check 7 (ruling 6):** dense L and dense L2 agree within the zone to 0.5 K.
   - That supports the far bound proposed below: coarse L2 against dense L2, away from the corner
     intervals, above the first corner interval.
   - Measured: 0.71 K (run 5, at 10 hPa below the stop corner), 1.21 K (run 7) and 0.31 K (run 7f).
   - **Proposed bound:** 1.5 K for runs 5 and 7, 0.5 K for run 7f.
   - The region below `p_s` is not bounded, for the reason in 3.
   - **The check is not yet edited:** this waits for the author's decision.

The L worktree is kept until the decision; `git worktree remove` clears it. Nothing is committed and
nothing is swept. The dense products under `reports/step07_1/dense/` are ignored files.

## 8. Closing (SPEC_07 v0.5, §9)

SPEC_07 is closed, not adopted. The cleanup the author gave the go for:

| Item | Result |
|---|---|
| 1. Park L2 | branch `spec07-l2-parked` at `89e3e6d`, pushed to `origin`: `lib/windfield.py`, `refrac/anchor.py`, the figure and docstring edits (`kernel.py`, `figures_inputs.py`, `figures_product.py`, `figures_profile.py`, `wind/shear.py`), the `step02_5` and `step04_1` restatements, and `tests/step07_1/` |
| 2. Restore `main` | `src/` and `tests/` as at `a238306` (no diff); the decision L worktree removed with `git worktree remove`; the ignored products under `reports/step07_1/` left in place |
| 3. The record | `5f2e05c`: SPEC_07 v0.5, `REPORT_07_preexecution.md`, this report, and a SPEC_07 section in `docs/specs/STATE.md` |
| 4. Outer loop pinned | `7652852`: `step04_4` and `step04_5` set 1e-8 themselves; the namelist is unchanged |
| 5. `reports/New folder` | removed (empty, untracked) |
| 6. Regression of record | on the clean tree at `7652852`: **30 of 30 suites at their reference counts**, no path dirty after any suite, the registered kind N file, the closure and the transfer inputs and products restored equal; 6923 s |

Not as expected, listed rather than worked around:

1. **`step04_5` run 1 cannot take a pinned tolerance.** It is the namelist's own run through the
   driver (`tr.run`), so a test cannot hand it a value without changing the namelist. There the pin
   is an assertion that the namelist carries 1e-8, as the suite already does for its spacings. Runs
   2 to 5, and every outer loop in `step04_4`, take 1e-8 from the suite. Both suites still read the
   iteration cap (50) from the namelist.
2. **The coarse L2 products are gone from disk.** The regression's `step05_4` reran the
   experiments under decision L, overwriting the ignored products under `forward/` and
   `reports/step05_4/` that §7 read as coarse L2. The comparison itself is kept in
   `reports/step07_1/dense/comparison.json`, and the branch rebuilds them.
3. **The reports on `main` name `tests/step07_1/` scripts that are now only on the branch.**
4. **`main` is not pushed.** It is ahead of `origin/main` by the SPEC_06 commits and these; the go
   covered pushing the branch only.
