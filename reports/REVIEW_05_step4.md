# REVIEW 05, Step 4. The named experiments

Reviewing agent, 2 October 2026. Report: `reports/REPORT_05_step4.md`. Specification:
`docs/specs/SPEC_05_Shear_Experiments.md` v0.7, §4. Read with the comparison figure
`reports/figures/step05_4_temperature_minus_run2.png`.

## Verdict

**Accepted.** The four checks pass, and the physics reads the way v0.7's expected outcomes say it
should. The three open findings are results of the kind SPEC_05 was written to produce; they are
recorded below as open items for the next specification and do not hold this step.

## The physics, read against the expected outcomes

- **The null is exact** (3a and 3b identical in N, p, T; `r0(10 N)` 0.004 m from the independent
  march), and relative to run 2 it is warm by run 2's own closure-wind change, 1.5 to 1.6 K.
- **Run 3c** gives 49.5 percent of run 2's change, a little under half, as `Omega_abs` containing
  `u` predicts.
- **The decay runs behave as thermal wind should.** Under a decay linear in ln p the shear is
  constant through the zone, so the temperature difference between the anchor and 10 N is constant
  with height there: about -10 K for run 4, -17 K for run 5 and -34 K for run 6, in proportion to
  the rates (11.7, 20.2, 39.6 percent per scale height). Above `p_stop` the change is exactly the
  null's +1.529 K, that is, the delivered temperature is the anchor's. Run 7's shear is strongest
  at 1 bar and falls off upward, and so does its change. The equator side is colder in every case.
  The reference surface moves with the wind at the gauge, run 5 by about half of the 36 km.
- **Runs 8 and 9** change only the levels below 1 bar, colder toward the equator, at the size the
  v0.3 correction predicted (comparable to run 5 at its levels).

## Rulings and explanations on the findings

1. **F9 at the anchor's latitude (finding 1).** The fix is right and was made at the author's
   direction; the figure checks were rerun. Accepted.
2. **The change at the level just above `p_s` (runs 8 and 9, section 3).** This is most likely the
   kink at `p_s` smeared over one geopotential cell of the mesh. The anchor's level at 998.7 mbar is
   about 700 m2/s2 above 1 bar, far inside one mesh cell (5e4 or 2.5e4), so the kernel interpolated
   between the mesh nodes on either side of the kink carries part of the shear onto that level. The
   decay runs show the same thing from the other side: run 5's change at that level is -8.7 K,
   about half its -17 K in the zone. It does not fall with the mesh spacing because both spacings are
   much coarser than 700 m2/s2. One run would confirm it cheaply: run 8 with `p_s` at 1.05e5 Pa,
   which puts the kink more than a mesh cell below the 998.7 mbar level; that level's change should
   then be close to zero. Optional; recorded as an open item.
3. **The outer loop's two failures (finding 2).** Recorded as results, as §4 asked, and nothing
   tuned. A residual parked at 1.4e-4 to 1.5e-4 for 50 passes, when converging runs fall by 10 to
   30 per pass, points to a limit cycle, most likely a discrete switch on the mesh (a traced curve
   moving between cells, or an extension decision) flipping back and forth between passes. Open
   item: record the per-pass residual history in the transfer record on every run, failed or not,
   so the next look can see whether it alternates. Under-relaxation would likely cure it, but that
   is a numerical change for a later step, not this one.
4. **The pressure identity of the decay runs (finding 3).** 2e-3 to 9e-3, not falling with the mesh
   spacing or the wind grid, largest inside each decay zone. In temperature that is a few tenths
   of a kelvin, 2 to 3 percent of the signal. Since neither the mesh nor the wind grid sets it, the
   candidates are the anchor's own fixed level spacing (the production integrates over the anchor's
   66 levels, and a profile whose temperature changes by 15 to 50 percent in the vertical stresses
   that integration in a way the closure wind never did) or the composition part of the identity
   that SPEC_04 §18 separated from the tracing part. Open item: report the identity's two parts
   separately for runs 5 and 7, which says which of the two it is.
   **Check 4's bound** (author, 2 October: fractions of a kelvin are meaningless here): 1e-2
   relative for every run, about 1 K at Saturn's temperatures. Measured up to 9.1e-3 in the sheared
   runs and 5.4e-7 to 5.8e-7 in the others. No further diagnosis of the identity is asked for.

## Two additions to this step (author, 2 October)

1. **The outer loop's residual, pass by pass.** `forward/transfer.py` prints the largest
   `|d ln p|` of every pass to the run's output, and the named failure's message carries the whole
   history. Printed, not stored: no product gains a variable, so no registered file changes. Then
   rerun the two runs that failed (run 6 at 5e4, run 7 at 2.5e4) and report the history and whether
   the residual alternates or settles. Nothing is tuned.
2. **The kink at `p_s`, confirmed or not.** One run: run 8's case with `shear_reference_pressure_Pa`
   at 1.05e5 Pa, at 5e4, in a copy under `reports/step05_4/` (no new run directory). Report the
   delivered temperature change from run 2 at the 998.7 mbar level. Near zero confirms the one
   mesh cell explanation of finding 2 above; if it is not, say what it is.

Both go into this report as a section of its own. **Regression:** the change to
`forward/transfer.py` only prints, so it cannot move a number; `step04_5` is rerun as the evidence
(its check 11 compares the closure production with the registered product), and `step05_4` is not
rerun beyond the two failed runs and the one new run.

The shear cases of this specification flex the code; they are not representative of Saturn. The
next step is a case built to be realistic, before any comparison with observations (author,
2 October).

## Order of work

1. The two additions above, reported in this step's report.
2. The author's go; the acceptance commit as section 6 lists plus `forward/transfer.py`, with check 4
   bounded as above. No AI attribution.
3. `STATE.md`: Step 4 accepted; SPEC_05's shear part complete.
