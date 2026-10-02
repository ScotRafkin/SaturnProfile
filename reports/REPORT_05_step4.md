# REPORT 05, Step 4. The named experiments

Coding agent, 2 October 2026. Specification: `docs/specs/SPEC_05_Shear_Experiments.md` v0.7, §4.
Working tree on `main` at `17b1539`; nothing of this step is committed.

**Result.** `tests/step05_4/accept_step05_4.py` passes 4 of 4. Run 2 reproduces the registered
transfer product exactly; runs 3a and 3b, the exact null, agree in N, p and T with no difference at
all; run 3a's `r0(10 N)` is the reviewing agent's no-wind surface to 0.004 m. Of the 21 runs (ten at
two spacings, 7f at one), 19 completed and two stopped with the outer loop's named failure, run 6
at 5e4 and run 7 at 2.5e4. The outcomes read as v0.7's expected outcomes say they should (section
3). Two things need the review: the two failures (finding 2), and a pressure identity of 2e-3 to
9e-3 in the decay runs that does not fall when the spacing is halved (finding 3). One defect was
found and fixed: F9 could not draw a run targeted at the anchor's own latitude (finding 1).

## 1. What was made

- **Eleven run directories under `forward/`**, each made by `casspian-new-run shear_r<N>_<tag>
  --from forward/lindal_transfer` and then edited in its `[shear]` section only: `case` and the
  parameters of the §4 table, `p_s = 1e5` Pa. Run 3b's namelist also takes the target at the
  anchor's `phi_c`, `30.80556842739218` deg; run 7f's build file also takes `[wind]
  pressure_grid_Pa` at twenty levels per decade, 121 values from 1 Pa to 1 MPa, `10^(k/20)` written
  to six significant digits. `[run] description` is left as copied, the person's to edit (§3).
- **`tests/step05_4/run_experiments.py`**: builds every run's inputs and runs it at 5e4 in its own
  directory (the run of record, figures F5 to F9 on), then at 2.5e4 in a copy under
  `reports/step05_4/spacing_2p5e4/<run>/` (the namelist's spacing halved, the anchor path rewritten
  for the copy's depth, figures off), 7f at 5e4 only. A run that stops is recorded with its failure,
  not raised. The working tree was not clean, so the `-dirty` refusal is relaxed in that process
  only, named in its output; every product records its stamp in `reports/step05_4/results.json`.
  Running it after the acceptance commit gives clean products. Runs named on its command line
  replace their own entries in `results.json` and keep the others.
- **`tests/step05_4/accept_step05_4.py`**: runs the experiments (unless `--no-run`), the four
  checks, the table below (`reports/step05_4/table.md`) and the comparison figure
  `reports/figures/step05_4_temperature_minus_run2.png`.
- **`src/casspian/tools/plots/figures_profile.py`**: F9's right panel draws the wind-grid nodes that
  bracket the transfer's range (finding 1).
- `tests/run_regression.sh`: the row `step05_4/accept_step05_4 4`.

## 2. The acceptance

| Check | Measured | Bound | Result |
|---|---|---|---|
| 1. run 2 against the registered transfer product: N, p, T, altitude, r0 | all five array-equal | exact | pass |
| 2. runs 3a and 3b: N, p, T level by level | largest relative difference 0.000e+00 in each | 1e-12 | pass |
| 3. run 3a's `r0(10 N)` | 60,092,307.686 m, -0.004 m from 60,092,307.69 | 0.1 m | pass |
| 4. the pressure identity of every completed run at both spacings | up to 9.10e-3 (run 7f); 5.4e-7 to 5.8e-7 in the unsheared runs | 1e-2 (REVIEW_05_step4) | pass |

## 3. The results

Every run delivers the anchor's 66 levels; differences from run 2 are taken level by level, by the
anchor's level index, since a run's isobar labels are its own anchor production under its own wind.
"dT" is the delivered temperature minus run 2's; "dz" the altitude minus run 2's; "identity" the
largest `|p/p_label - 1|`; "gauge" the 100 mbar level, "top" level 0 (0.2 mbar).

| Run | Spacing | Largest dT from run 2 (K) | at p (Pa) | dT at 1 bar | dT at gauge | dT at top | dz at 1 bar (m) | dz at gauge (m) | dz at top (m) | r0 (m) | identity | passes | wall (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| shear_r2_uniform | 5e4 | +0.000 | 19.95 | +0.000 | +0.000 | +0.000 | +0.0 | +0.0 | +0.0 | 60,128,612.97 | 5.78e-07 | 2 | 35.2 |
| shear_r2_uniform | 2.5e4 | +0.000 | 2882 | +0.000 | +0.000 | +0.000 | -0.0 | -0.0 | +0.0 | 60,128,612.97 | 5.77e-07 | 2 | 20.2 |
| shear_r3a_nowind | 5e4 | +1.626 | 31.65 | +1.467 | +0.912 | +1.529 | -3.4 | -227.1 | -949.8 | 60,092,307.69 | 5.41e-07 | 1 | 32.6 |
| shear_r3a_nowind | 2.5e4 | +1.626 | 31.65 | +1.467 | +0.912 | +1.529 | -3.4 | -227.1 | -949.8 | 60,092,307.69 | 5.41e-07 | 1 | 19.0 |
| shear_r3b_nowind_anchor | 5e4 | +1.626 | 31.65 | +1.467 | +0.912 | +1.529 | -9.2 | -8122.3 | -34355.3 | 58,516,188.29 | 5.41e-07 | 1 | 32.9 |
| shear_r3b_nowind_anchor | 2.5e4 | +1.626 | 31.65 | +1.467 | +0.912 | +1.529 | -9.2 | -8122.3 | -34355.3 | 58,516,188.29 | 5.41e-07 | 1 | 19.8 |
| shear_r3c_half | 5e4 | +0.821 | 31.65 | +0.741 | +0.460 | +0.772 | -1.7 | -114.5 | -478.7 | 60,110,315.85 | 5.61e-07 | 2 | 36.7 |
| shear_r3c_half | 2.5e4 | +0.821 | 31.65 | +0.741 | +0.460 | +0.772 | -1.7 | -114.5 | -478.7 | 60,110,315.85 | 5.59e-07 | 2 | 19.9 |
| shear_r4_decay12 | 5e4 | -10.031 | 8.314e+04 | -5.021 | -9.766 | -3.436 | -29.8 | -9834.4 | -36815.1 | 60,118,683.65 | 2.07e-03 | 6 | 39.2 |
| shear_r4_decay12 | 2.5e4 | -10.043 | 8.708e+04 | -5.281 | -9.773 | -3.471 | -4.5 | -9866.5 | -36898.0 | 60,118,683.65 | 2.82e-03 | 6 | 25.4 |
| shear_r5_decay20 | 5e4 | -17.200 | 7.576e+04 | -8.746 | -16.698 | +1.529 | -49.8 | -16823.7 | -36817.3 | 60,111,621.40 | 2.60e-03 | 7 | 40.6 |
| shear_r5_decay20 | 2.5e4 | -17.234 | 8.708e+04 | -9.058 | -16.713 | +1.529 | -9.7 | -16882.5 | -36891.9 | 60,111,621.40 | 3.79e-03 | 7 | 27.0 |
| shear_r6_decay40 | 5e4 | stopped: ValueError | | | | | | | | | | | 76.4 |
| shear_r6_decay40 | 2.5e4 | -33.796 | 8.71e+04 | -18.271 | -31.377 | +1.529 | -39.9 | -32914.6 | -36778.5 | 60,095,468.90 | 3.14e-03 | 9 | 30.4 |
| shear_r7_decay_linp | 5e4 | -70.306 | 8.72e+04 | -40.568 | -6.732 | +1.512 | -328.7 | -31669.9 | -35947.2 | 60,095,886.48 | 8.67e-03 | 10 | 45.9 |
| shear_r7_decay_linp | 2.5e4 | stopped: ValueError | | | | | | | | | | | 77.7 |
| shear_r7f_decay_linp_fine | 5e4 | -70.015 | 8.715e+04 | -42.805 | -7.138 | +1.512 | -336.8 | -31676.4 | -35962.6 | 60,095,886.48 | 9.10e-03 | 10 | 44.2 |
| shear_r8_increase25 | 5e4 | -9.444 | 1.294e+05 | -5.872 | +0.000 | +0.000 | +3.6 | -119.2 | -119.2 | 60,128,612.97 | 5.46e-04 | 6 | 38.6 |
| shear_r8_increase25 | 2.5e4 | -9.408 | 1.294e+05 | -5.384 | +0.000 | +0.000 | -1.7 | -58.0 | -58.0 | 60,128,612.97 | 4.34e-05 | 5 | 23.9 |
| shear_r9_increase50 | 5e4 | -18.905 | 1.294e+05 | -11.667 | +0.000 | +0.000 | +8.9 | -240.7 | -240.7 | 60,128,612.97 | 1.13e-03 | 6 | 38.7 |
| shear_r9_increase50 | 2.5e4 | -18.830 | 1.294e+05 | -10.721 | +0.000 | +0.000 | -2.8 | -117.3 | -117.3 | 60,128,612.97 | 1.04e-04 | 6 | 25.3 |

Read against v0.7's expected outcomes:

- **Run 2** is the accepted transfer product (check 1).
- **Runs 3a and 3b** are the exact null (check 2): one outer-loop pass, and the delivered N, p and T
  of 3a and 3b identical. Relative to run 2 the null is warmer by 1.47 K at 1 bar and 1.63 K at its
  largest: run 2's own change, undone. Their altitudes differ (3b at the anchor's latitude, -34 km
  at the top), which is geometry, as §4 says.
- **Run 3c** changes the temperature by 49.5 percent of run 2's change (0.805 of 1.626 K at the
  largest level), a little under half, as §4 expects from `Omega_abs` containing `u`.
- **Runs 4 to 7** change the delivered temperature between `p_s` and `p_stop`, colder toward the
  equator, and the more so the faster the decay: about -10 K (run 4), -17 K (run 5), -34 K (run 6,
  at 2.5e4) and -70 K (run 7), against temperatures of about 130 to 140 K near 1 bar. Above `p_stop`
  the delivered temperature is the anchor's: runs 5 and 6 read +1.529 K from run 2 at the top, the
  null's value. Run 4's `p_stop` is the column's top (20 Pa), so its whole column is in the zone.
  The reference surface moves with the decayed wind at the gauge: run 5's `r0(10 N)` is 60,111,621
  m, about half of the 36 km between the closure and no-wind surfaces, as §4 estimated.
- **Run 7f**, the same case on the finer wind grid, delivers within 0.3 K of run 7 at its largest
  and 2.2 K at 1 bar: the ten per decade grid's representation of the curved shape moves the answer
  by a few percent of the change.
- **Runs 8 and 9** change the temperature at the bottom of the column: -9.4 K and -18.9 K at its
  bottom level (1.294 bar), -5.9 K and -11.7 K at the 1 bar level (998.7 mbar, just above `p_s`),
  nothing at the gauge or the top, and the reference surface unchanged. Their altitudes above the
  zone shift by a constant (-119 m and -241 m at the gauge and the top, at 5e4), halving at 2.5e4.
  That the level just above `p_s` changes is a measurement, not yet explained.

The figure shows runs 3a, 3c, 4, 5, 7, 7f, 8 and 9 at 5e4. Run 3b is at another latitude; run 6
did not complete at 5e4.

## 4. Findings

1. **F9 could not draw a run targeted at the anchor's own latitude** (run 3b at 5e4: `zero-size
   array to reduction operation fmin`). The transfer's latitude range was 30.76 to 30.86 deg, which
   holds no node of the 0.5 deg wind grid, so the right panel's field was empty. The product was
   written before the figures and was intact; only the run's figures failed, and the 2.5e4 copy,
   figures off, completed. Fixed at the author's direction within this step: the panel draws the
   wind-grid nodes that bracket the transfer's range in latitude and in pressure, so it is never
   empty; for every other run the panel gains at most one node at each edge. Run 3b was rerun and
   completed at both spacings. The figure checks were rerun (section 5).
2. **Two runs stop with the outer loop's named failure, at opposite spacings.** Run 6 (40 percent per
   scale height) at 5e4 and run 7 (linear in p) at 2.5e4: `the outer loop did not converge in 50
   passes: the largest |d ln p| on the mesh is 1.509e-04` (run 6) and `1.392e-04` (run 7), against a
   tolerance of 1e-8. Each converges at the other spacing (9 and 10 passes). Every run that converged
   did so in 1 to 10 passes, its residual falling by a factor of about 10 to 30 a pass (run 5: 0.39,
   3.4e-2, 2.2e-3, 1.1e-4, 4.8e-6, 1.7e-7, 5.1e-9). A residual still at 1.5e-4 after 50 passes
   suggests the loop is cycling rather than converging slowly, but the per-pass history is not
   recorded when a run fails, so the mechanism is not established. Reported as results; nothing was
   tuned.
3. **The pressure identity of the decay runs does not fall with the spacing.** Run 4: 2.07e-3 at 5e4
   and 2.82e-3 at 2.5e4; run 5: 2.60e-3 and 3.79e-3; runs 6, 7 and 7f: 3.1e-3, 8.7e-3 and 9.1e-3.
   The increase runs fall about tenfold (run 8: 5.5e-4 to 4.3e-5; run 9: 1.1e-3 to 1.0e-4), as the
   SPEC_04 Step 5 synthetic sheared run did (3.3e-4 to 1.5e-4), and the uniform runs stay at the
   closure's 5.4e-7 to 5.8e-7. The finer wind grid of run 7f does not reduce it either. So in the
   decay runs something other than the mesh or the wind grid sets the identity, at 0.2 to 0.9 percent
   of the pressure. The largest sits at level 1 (run 4), 15 (run 5), 28 (run 6) and 59 (runs 7, 7f),
   inside each case's decay zone. The altitudes show the same sensitivity: run 4's at 1 bar differs
   by 25 m between the spacings. Check 4 is to be bounded from these measurements; this one should be
   understood first.
4. **The 2.5e4 runs are faster than the 5e4 ones** (20 s against 35 s for run 2) because the copies
   run with figures off.

## 5. The regression

Step 4 adds run directories and scripts, and the F9 fix reaches every suite that draws F9: its own
suite and the figure checks, `step02_5`, `step03_4`, `step04_5`, `step05_3`. Through the driver:

| Suite | Checks | Time (s) |
|---|---|---|
| `step02_5/accept_step02_5` | 7 of 7 | 190 |
| `step03_4/accept_step03_4` | 13 of 13 | 72 |
| `step04_5/accept_step04_5` | 11 of 11 | 905 |
| `step05_3/accept_step05_3` | 7 of 7 | 187 |
| `step05_4/accept_step05_4` | 4 of 4 | 1060 |

5 of 5 at their reference counts; the registered kind N (`b79c60cb...`), the closure inputs and
product, and the transfer inputs and product restored exactly. Afterwards `git diff --name-only --
'*.nc'` is empty and no file carries a `skip-worktree` mark. The `step05_4` suite reran all 21 runs
from scratch; the table of section 3 is that rerun's.

## 6. What is committed on acceptance

The eleven run directories' control files, `tests/step05_4/run_experiments.py`,
`tests/step05_4/accept_step05_4.py`, `tests/step05_4/review_additions.py`,
`src/casspian/forward/transfer.py`, `src/casspian/tools/plots/figures_profile.py`,
`tests/run_regression.sh`, the comparison figure `reports/figures/step05_4_temperature_minus_run2.png`
and this report. The runs' inputs and products are ignored, like every run's.

## 7. The review's two additions (REVIEW_05_step4)

`tests/step05_4/review_additions.py` makes both, with the `-dirty` refusal relaxed in its process
as in `run_experiments.py`; its output is `reports/step05_4/review_additions.txt` and `.json`.
Check 4 is bounded at 1e-2 relative, as the review set it, and passes.

**1. The outer loop's residual, pass by pass.** `forward/transfer.py` prints, on every pass, the
largest `|d ln p|` and the mesh's size, and the named failure's message ends with the whole history.
Nothing is stored; no product gains a variable. The two runs that stopped were rerun and stop as
before:

| Pass | Run 6 at 5e4 | Run 7 at 2.5e4 |
|---|---|---|
| 1 | 6.542e-01 | 7.263e-01 |
| 2 | 1.000e-01 | 1.810e-01 |
| 3 | 1.337e-02 | 5.363e-02 |
| 4 | 1.439e-03 | 1.256e-02 |
| 5 | 1.347e-04 | 2.005e-03 |
| 6 | 1.520e-04 | 3.005e-04 |
| 7 | 1.508e-04 | 1.799e-04 |
| 8 | 1.509e-04 | 1.391e-04 |
| 10 to 50 | 1.508572e-04, every pass | 1.391742e-04, every pass |

Neither alternates in the printed value, and neither decays: each converges normally for five to
eight passes, then settles to a constant held to seven digits for forty passes. The mesh does not
change (421 x 82 and 421 x 160 throughout; no extension). The residual is the distance between
successive maps, so a constant non-zero value is not a fixed point: the maps cycle with a step of
fixed size, most simply between two states. The printed scalar cannot tell a two-state cycle from
a longer one; the sign of the difference, or the map itself, would. Nothing was tuned.

**2. The kink at `p_s`.** Run 8's case with `shear_reference_pressure_Pa` at 1.05e5 Pa (the
review's run) and at 1.3e5 Pa (at the author's request, below), each at 5e4 in a copy under
`reports/step05_4/` (`r8_ps105/`, `r8_ps130/`): 5 passes each, identity 4.37e-4 and 5.95e-4. The
delivered temperature from run 2 (K), at the 1 bar level and the six levels below it:

| `p_s` | 59 (998.7 mbar) | 60 (1.045 bar) | 61 (1.095) | 62 (1.147) | 63 (1.201) | 64 (1.258) | 65 (1.294) |
|---|---|---|---|---|---|---|---|
| 1e5 Pa (run 8) | -5.872 | -7.815 | -9.224 | -9.318 | -9.351 | -9.424 | -9.444 |
| 1.05e5 Pa | -4.737 | -6.299 | -7.439 | -7.513 | -7.543 | -8.357 | -8.988 |
| 1.3e5 Pa | 0.000 (5.7e-14) | 0.000 | 0.000 | 0.000 | -0.005 | -4.033 | -6.828 |

**Confirmed.** With `p_s` at 1.3e5 Pa the 998.7 mbar level's change is zero to round-off, and so
are the next three levels'. The change reaches upward from the kink by about one mesh cell: with
`R/M` about 3780 J/(kg K) and `T` about 134 K, one unit of `ln p` is about 5.1e5 m2/s2, so level 64
(1.258 bar, -4.0 K) lies about 1.7e4 m2/s2 above the kink, a third of a 5e4 cell, and level 63
(1.201 bar, -0.005 K) about 4.0e4 m2/s2 above it, most of a cell. The review's run at 1.05e5 Pa
moved the kink only about 2.5e4 m2/s2 below 1 bar, half a cell, so the 998.7 mbar level was still
within reach, which is why its change fell by a fifth and not to zero. The one-mesh-cell
explanation of the review's ruling 2 holds.

**Regression.** The change to `forward/transfer.py` only prints. `step04_5` was rerun through the
driver as the evidence:

| Suite | Checks | Time (s) |
|---|---|---|
| `step04_5/accept_step04_5` | 11 of 11 | 934 |

1 of 1 at its reference count; its check 11, the closure production against the registered
product, passes. The registered kind N, the closure inputs and product, and the transfer inputs and
product restored exactly; afterwards `git diff --name-only -- '*.nc'` is empty and no file carries a
`skip-worktree` mark.
