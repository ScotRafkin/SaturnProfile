# REPORT 06, Step 1. The shear read at the wind file's resolution, in the transfer and the tracing

Coding agent, 2 October 2026. Specification: `docs/specs/SPEC_06_Kernel_On_The_Isobar.md` v0.4,
§1 with the rulings of §4 to §6. Working tree on `main` at `d5f4147`; nothing of this step is
committed. Filed for acceptance; two interim reports preceded it (section 9; their full text is
kept in `reports/step06_1/REPORT_06_step1_interim.md`).

**Result.** `tests/step06_1/accept_step06_1.py` passes 7 of 7. The model no longer lets its mesh
resample the wind: the transfer reads the shear on the isobar at the curve's label, and the tracing
integrates it in each column with the wind file's pressure nodes as breakpoints. The sawtooth
along the isobar is gone (131 times smaller); run 8's leak above `p_s` is 0.0000 K (SPEC_05 -5.87
K); every run's pressure identity and its 1 bar temperature are the same at both mesh spacings;
and the two runs that cycled under SPEC_05 converge at 1e-8. What remains of the identity is the
production's integration across the temperature break a wind file node puts between two anchor
levels, measured and explained (check 7): its size is the layer's mass times the jump in `ln T`,
its sign set by where in the layer the node falls. The full regression: of 30 suites, the three
that compare against the registered transfer product fall short, since this step moves it by
2.3e-6 K and the sweep after acceptance rebuilds it; every other suite is at its reference count
(section 7).

## 1. What changed

| File | Change |
|---|---|
| `src/casspian/lib/kernel.py` | `IsobarKernel` and `isobar_kernel` (the mesh's `r`, `g` and the two isobar slopes, the `WindField`, `Omega`, formed once per state); `shear_at(kernel, phi, Phi, p)`, Eq. A15 with the wind read from the interpolant at `(phi, p)` and the geometry bilinear on the mesh; `shear_on_isobar`, the same at the curve's label; `ColumnIntegral` and `column_integral`, `I` (Eq. A16) cell by cell, each cell piece by piece between the geopotentials where the column crosses a wind file pressure node, a two-point Gauss rule on each piece; between nodes the node below plus the partial cell; between mesh latitudes linear in latitude. `shear_kernel`, `shear_integral` and `transfer_kernel` unchanged. |
| `src/casspian/forward/transfer.py` | `IsobarMap.geopotential_at`, the inverse of `ln_p_at` (the breakpoints). `TransferState.isobar_kernel`; `shear_integral` is now a `ColumnIntegral`. `trace` keeps its signature and reads a node array bilinearly as before (the synthetic traces of `step04_4`). `transfer(kernel, curves, ln_N, label_pressure_Pa, ...)` (deliverable 7); `carry_to`, the anchor arrivals and F8's `shear_kernel_<slug>_per_rad` pass the labels; that variable's description says where the term is read. Run time unchanged (run 8, 16.5 s). |
| `tests/step06_1/` | `accept_step06_1.py` (the seven checks), `stop_position.py` (check 7's measurement), `cycled_runs.py` (the two runs that cycled, at 1e-8). |

**Accepted suites edited**, each mechanical or restated by the specification:

| Suite | Edit | Why |
|---|---|---|
| `step04_4/accept_step04_4.py` (6 calls), `step04_4/diagnose_anchor_placement.py` (2), `step04_5/accept_step04_5.py` (1) | `tr.transfer(state.mesh, state.s_over_g, curves, ...)` becomes `tr.transfer(state.isobar_kernel, curves, ..., labels, ...)`, the labels each call already gave `composition_term` | deliverable 7, the signature. One diagnose script calls `transfer`, not two |
| `step04_4/accept_step04_4.py` check 11 | the report line's ratios name a finer value of exactly zero instead of dividing by it | the sheared out-and-back `ln N` is now exactly 0 at the finer spacing (1.8e-15 at the coarser) and the line crashed; the check's condition is unchanged and passes, the geopotential falling from 2.8e-9 to 4.7e-10 m2/s2 under halving |
| `step05_2/accept_step05_2.py` | `compare` takes a 0-d variable (`np.atleast_1d`) | a scalar of a registered file differed for the first time and the path crashed |
| `step05_4/accept_step05_4.py` check 4 | REVIEW_05_step4's bound of 1e-2 on the identity restated as SPEC_06 v0.4 restates check 6: every run completes and the largest identity agrees between 5e4 and 2.5e4 within 10 percent, the size reported | §6 ruling 1: the identity at a break is the production's integration over the anchor's levels, which no mesh setting reaches; runs 5, 7 and 7f measure 1.1e-2 to 1.6e-2. For the review to confirm |

Reference counts are unchanged.

## 2. The acceptance

| Check | Measured | Bound | Result |
|---|---|---|---|
| 1. transfer run, no vertical shear, against the registered product | largest dT 2.3e-6 K (level 29), altitude 1.1e-3 m, `r0` 0, relative N 1.4e-10; identity 5.790e-7 (5.777e-7) | 0.1 K, 10 m, 1 m | pass |
| 2. closure product | 9 of 9 array-equal | exact | pass |
| 3. run 7f at 5e4, 998.7 mbar isobar, largest off-node step of `S/g` | 0.0276 per rad against SPEC_05's 3.611, 131 times smaller | ten times | pass |
| 4. run 8 at 5e4, dT at 998.7 mbar from the no-shear run | -0.0000 K (SPEC_05 -5.872) | 0.1 K | pass |
| 5. run 5, T at 998.7 mbar, 5e4 against 2.5e4 | 116.398 and 116.412 K, 0.013 K apart; dT -16.77 and -16.75 K, the zone's value (SPEC_05 -8.75 and -9.06) | 0.5 K | pass |
| 6. largest identity, 5e4 against 2.5e4, every run of the SPEC_05 set | within 2.0 percent in run 7, 0.0 to 0.4 percent in the others | 10 percent | pass |
| 7. the identity at the stop pressure explained | section 4 | as stated there | pass |

Figures: `reports/step06_1/F8_sawtooth.png` (old and new `S/g` along the isobar); run 7f's new F8
in `forward/shear_r7f_decay_linp_fine/output/figures/`, its SPEC_05 F8 in
`reports/step06_1/spec05/r7f_figures/`, for the author's viewing.

## 3. The identity: no mesh, and what it is at `p_s`

The largest `|p/p_label - 1|`, SPEC_05 and now, with the 1 bar temperature change from run 2 (K):

| Run | identity 5e4 / 2.5e4, SPEC_05 | identity 5e4 / 2.5e4, now (level) | dT at 1 bar, SPEC_05 | dT at 1 bar, now |
|---|---|---|---|---|
| 2, 3a, 3b, 3c | 5.4e-7 to 5.8e-7 | 5.4e-7 to 5.8e-7 | | unchanged |
| 4 decay 12 | 2.07e-3 / 2.82e-3 | 1.64e-3 / 1.64e-3 (60) | -5.02 / -5.28 | -9.90 / -9.89 |
| 5 decay 20 | 2.60e-3 / 3.79e-3 | 1.51e-2 / 1.51e-2 (16) | -8.75 / -9.06 | -16.77 / -16.75 |
| 6 decay 40 | 2.23e-3 / 3.14e-3 | 6.41e-3 / 6.40e-3 (60) | -17.91 / -18.27 | -32.08 / -32.02 |
| 7 decay linear in p | 8.67e-3 / 2.82e-3 | 1.10e-2 / 1.08e-2 (60) | -40.56 / -42.64 | -68.65 / -68.61 |
| 7f, the finer wind grid | 9.10e-3 | 1.64e-2 (60) | -42.81 | -71.25 |
| 8 increase 25 | 5.46e-4 / 4.34e-5 | 1.45e-3 / 1.45e-3 (60) | -5.87 / -5.38 | -0.000 / -0.000 |
| 9 increase 50 | 1.12e-3 / 1.04e-4 | 2.89e-3 / 2.89e-3 (60) | -11.66 / -10.72 | -0.000 / -0.000 |

The residual is zero above the break and steps at it; below it the delivered pressure carries a
constant absolute offset, because the production integrates absolute layer masses from the top
(`lib.hydrostatic.pressure_from_top`), each layer's density taken log-linear across it. At `p_s`
the wind file's node at 1000 hPa sits just below level 59 (998.7 mbar), near the top of the layer
from 59 to 60, and the residual's step there has the size and sign check 7's law gives for a break
at `s` near zero (the second interim report's trapezoid estimate matched it within 15 percent in
all seven sheared runs). Runs 4 to 9 now read the 1 bar level at the zone's value, the same at
both spacings; SPEC_05's smaller identities came from the mesh blend spreading each break over
several levels, which also moved them with the spacing.

## 4. Check 7: the identity at the stop pressure

`tests/step06_1/stop_position.py` places the stop pressure of run 5's case (layer 15 to 16, 632.3
to 795.4 Pa) and run 6's case (layer 27 to 28, 7243 to 8728 Pa) at five positions `s` across the
layer, `p_stop = p_k (p_{k+1} / p_k)^s`, each added to the wind file's grid as a node so that the
file's break sits exactly there. The offset `p - p_label` gained across the layer (Pa):

| `s` | 0.02 | 0.25 | 0.50 | 0.75 | 0.98 |
|---|---|---|---|---|---|
| run 5's case | -10.53 | -5.10 | +1.15 | +7.75 | +14.22 |
| run 6's case | -265.79 | -132.39 | +31.68 | +217.43 | +410.28 |

- It arises in that layer alone: the neighbouring layers gain at most 0.09 Pa (run 5) and 0.56 Pa
  (run 6).
- It is linear in `s`, within 1.4 and 3.0 percent of its range, zero at `s` = 0.442 and 0.426.
- Its slope, 25.8 and 703.3 Pa, is the layer's mass (163.1 and 1485.0 Pa) times the jump in `ln T`
  across it relative to run 2 (0.1619 and 0.4716): ratios 0.976 and 1.004.

**The mechanism.** Between its nodes the wind file is linear in `ln p`, so its shear is a constant
in each interval and changes at a node (decision L). By the thermal wind, the delivered temperature
at 10 N then breaks at that node: the change of shear, integrated in latitude from the anchor,
arrives as a jump in temperature between the two levels on either side of it. The production takes
the density log-linear across each anchor layer, which spreads the jump over the whole layer. If
the break lies near the bottom of the layer (`s` near 1), the true layer is mostly the side above
the break, and the production gives the side below too much of it; where the side below is denser
(colder), the pressure beneath comes out high. Near the top (`s` near 0), the opposite. The offset
is the layer's mass times the jump in `ln T` times `s - s0`, `s0` a little under one half; its sign
is set by which half of the layer the break falls in and by the sign of the jump.

**The original runs.** Run 5's stop, 700 Pa, is not a node: the file breaks at its nodes 630.957 Pa
(layer 14, `s` = 0.991) and 794.328 Pa (layer 15, `s` = 0.994), both near the bottom of their
layers, so both offsets are positive. The law gives +5.06 and +6.66 Pa, +11.72 Pa in all, against
+12.02 Pa measured at level 16 (0.30 Pa apart). Run 6's break is the node 7943.28 Pa (layer 27, `s`
= 0.495), just above the zero crossing: the straight line gives +47.45 Pa against +27.19 Pa
measured, 20.3 Pa apart, 3.0 percent of the measured range of 676.1 Pa; the measured points
interpolated at `s` = 0.495 give +28.62 Pa. At `p_s` the node sits at `s` near zero, which is why
the earlier estimate fitted there and had the opposite sign at the stop pressure. The Lindal levels
from 1.6 to 16 mbar (levels 9 to 19) lie within 0.3 percent of the wind grid's ten-per-decade
nodes, most of them just above, so a break at one of those nodes falls near one end of a layer and
its sign follows that end.

**A bound restated.** The script first required the law to reproduce the original runs' offsets
within 15 percent; run 6's break lies at the zero crossing, where a relative bound measures nothing
(a 3 percent departure of the line is 75 percent of a small number). The bound is now within 5
percent of the measured range, stated in the script. The specification sets no number for check 7.

## 5. The two runs that cycled

At 1e-8 with this step (`tests/step06_1/cycled_runs.py`): run 6 at 5e4 converges in 9 passes, the
residual 6.4e-1, 9.6e-2, 1.3e-2, 1.4e-3, 1.2e-4, 9.1e-6, 5.7e-7, 3.2e-8, 1.6e-9; run 7 at 2.5e4 in
11 passes, to 2.1e-9. Under SPEC_05, and under the first interim stage, which changed only the
transfer, both held a constant residual for forty passes. The cycle was the mesh blend in the
tracing.

## 6. The SPEC_05 experiments rerun

All 21 complete: `reports/step06_1/table_final.md` beside `reports/step06_1/spec05/table.md`; the
comparison figure `reports/figures/step05_4_temperature_minus_run2.png` redrawn. Runs 2 and 3a to
3c are unchanged in temperature and `r0` (run 2 moves by 2.3e-6 K). In every sheared run `r0` is
unchanged and the outer loop takes the same or fewer passes. The largest temperature change moves
by at most 4.4 K (run 7f, -70.0 to -74.4 K); the 1 bar level is section 3's table.

**F8 of run 7f.** The stepping along the isobar is gone. The delivered temperature still steps
through the lowest 50 km, alternating level by level from 525 to 950 mbar: §6 ruling 3 places it as
decision L's staircase in the vertical, for SPEC_07's interpolant. Check 7's law is that staircase
seen in the pressure.

## 7. The regression

The full set through the driver, 30 suites:

| Suite | Checks | Time (s) | | Suite | Checks | Time (s) |
|---|---|---|---|---|---|---|
| `step1/accept_step1` | 6 of 6 | 2 | | `step02_6/accept_step02_6` | 7 of 7 | 297 |
| `step1/verify_review_changes` | 14 of 14 | 2 | | `step03_1/accept_step03_1` | 9 of 9 | 9 |
| `step2/accept_step2` | 8 of 8 | 2 | | `step03_2/accept_step03_2` | 8 of 8 | 14 |
| `step3/accept_step3` | 5 of 5 | 3 | | `step03_3/accept_step03_3` | 16 of 16 | 154 |
| `step4/accept_step4` | 7 of 7 | 0 | | `step03_4/accept_step03_4` | 13 of 13 | 75 |
| `step5/accept_step5` | 7 of 7 | 40 | | `step04_0/accept_step04_0` | 15 of 15 | 306 |
| `step6/accept_step6` | 6 of 6 | 28 | | `step04_1/accept_step04_1` | 13 of 13 | 27 |
| `step7/accept_step7` | 8 of 8 | 90 | | `step04_2/accept_step04_2` | 13 of 13 | 525 |
| `step8/accept_step8` | 9 of 9 | 3 | | `step04_3/accept_step04_3` | 9 of 9 | 1484 |
| `step9/accept_step9` | 6 of 6 | 5 | | `step04_4/accept_step04_4` | 11 of 11 | 975 |
| `step02_1/accept_step02_1` | 6 of 6 | 11 | | `step04_5/accept_step04_5` | 11 of 11 | 957 |
| `step02_2/accept_step02_2` | 7 of 7 | 44 | | `step05_1/accept_step05_1` | 11 of 11 | 20 |
| `step02_3/accept_step02_3` | 9 of 9 | 85 | | `step05_2/accept_step05_2` | **2 of 3** | 196 |
| `step02_4/accept_step02_4` | 9 of 9 | 165 | | `step05_3/accept_step05_3` | **5 of 7** | 193 |
| `step02_5/accept_step02_5` | 7 of 7 | 181 | | `step05_4/accept_step05_4` | **3 of 4** | 1039 |

`step04_4` and `step05_2` are their reruns after the two fixes of section 1 (in the full run they
crashed). 27 of 30 at their reference counts; the registered kind N, the closure inputs and product,
and the transfer inputs and product restored exactly; afterwards `git diff --name-only -- '*.nc'` is
empty and no file carries a `skip-worktree` mark.

**The three that fall short compare against the registered transfer product,** which this step
moves (check 1: 2.3e-6 K at most):

- `step05_2` check 2: of the 18 registered files, 17 bit-identical (the closure product's 136
  variables among them) and the transfer product 162 of 177;
- `step05_3` checks 5 and 6: the migrated transfer run and the new-run copy against the registered
  product, array-equal in `r0` only;
- `step05_4` check 1: run 2 against the registered product, array-equal in `r0` only.

Each is restored by the sweep the specification orders after the acceptance commit, which rebuilds
the two runs' registered files on the clean tree; these three suites are rerun then, and the
comparison added here.

## 8. What is committed on acceptance

`src/casspian/lib/kernel.py`, `src/casspian/forward/transfer.py`, `tests/step06_1/`, the five
accepted suites of section 1, the redrawn comparison figure, this report and the review; then,
after `python tests/step04_0/sweep.py --runs`, the rebuilt registered run files.

## 9. History: the two interim stages

**v0.2, the transfer alone** (first interim report). The transfer read the shear on the isobar; the
tracing still integrated the node kernel. Checks 1, 2, 3 and 5 passed; check 4 missed at -0.28 K,
the arrival geopotential (+1119 m2/s2, unchanged) carrying the tracing's blend, while the carried N
no longer changed (7.8e-8 relative, from 4.6e-2). The identity grew up to 40 times at every break,
since the tracing and the transfer now disagreed there, and the two runs that cycled still cycled.
Ruled (§5): extend the reading to the tracing, as a principle: the mesh never resamples the wind.

**v0.3, the tracing too** (second interim report). Checks 1 to 5 passed, check 4 at 0.0000 K; the
identity became the same at both spacings; the runs that cycled converged. The v0.3 check 6, the
identity below 1e-2 and at or below SPEC_05's, failed: at a break the identity is the production's
integration over the anchor's levels, which the trapezoid estimate matched at `p_s` but not, in
sign, at the stop pressure. Ruled (§6): check 6 restated as mesh independence; the stop pressure
explained before acceptance (check 7); F8's vertical stepping to SPEC_07; the standing rule that
the code must work on any profile and an unexplained behavior in an extreme case is explained, not
deferred.
