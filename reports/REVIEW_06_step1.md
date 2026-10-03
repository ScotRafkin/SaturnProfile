# REVIEW 06, Step 1. The shear read at the wind file's resolution, in the transfer and the tracing

Reviewing agent, 3 October 2026. Report: `reports/REPORT_06_step1.md`. Specification:
`docs/specs/SPEC_06_Kernel_On_The_Isobar.md` v0.4.

## Verdict

**Accepted.** 7 of 7, and the model no longer lets its mesh resample the wind anywhere the shear
enters.

## What was checked

- **The principle is implemented in both places.** The transfer reads the shear on the isobar at the
  curve's label; the tracing integrates it column by column with the wind file's pressure nodes as
  breakpoints and a two-point Gauss rule on each piece. The node kernel remains only where nothing
  reads it at a break.
- **Mesh independence**, the property the step was for: every run's largest identity and its 1 bar
  temperature agree between 5e4 and 2.5e4 (identity within 2 percent in run 7, 0.4 percent or less
  in the others; run 5's 998.7 mbar temperature within 13 mK).
- **The leak across `p_s` is gone** (run 8, 0.0000 K from -5.87 K), and the decay runs now read their
  zone's value at 1 bar instead of about half of it.
- **The two runs that cycled converge at 1e-8** in 9 and 11 passes. The cycle was the mesh blend in
  the tracing: removed at its cause, without loosening the tolerance.
- **No vertical shear, no change:** 2.3e-6 K at most, the closure product array-equal.
- **Check 7 is explained, and convincingly.** The stop pressure placed at five positions across a
  layer gives an offset confined to that layer, linear in the position, crossing zero a little
  under half way, with a slope equal to the layer's mass times the jump in `ln T` (ratios 0.976 and
  1.004). The law then reproduces run 5's measured offset within 0.3 Pa and explains why the sign at
  the stop pressure is opposite to the one at `p_s`: the wind file's node sits near the bottom of
  its anchor layer at the stop pressure and near the top at `p_s`. Restating the script's bound
  from 15 percent of the value to 5 percent of the measured range is right, since run 6's break lies
  near the zero crossing.
- **The accepted suites edited** are the signature change, two crash fixes in report lines (check
  conditions unchanged), and `step05_4` check 4 restated as SPEC_06 v0.4 restates check 6. All
  confirmed.

## The mechanism, and what it means for SPEC_07

The remaining identity at a break is the production's log-linear density across an anchor layer
smoothing a temperature jump that the wind's shear jump puts there by the thermal wind. Its size is
set by the anchor's levels, not by any numerical setting. A wind interpolant with a continuous first
derivative (SPEC_07) has no shear jumps, so the delivered temperature has no jumps, and both this
identity and the stepping in run 7f's lowest 50 km should fall together. SPEC_07's acceptance will
measure both against this step's values.

## Order of work

1. The author's go; the acceptance commit as section 8 lists. No AI attribution.
2. `python tests/step04_0/sweep.py --runs` on the clean tree; the rebuilt registered run files
   compared by value and committed; `step05_2`, `step05_3` and `step05_4` rerun, all at their
   reference counts expected; the comparison added to the report.
3. `STATE.md`: SPEC_06 Step 1 accepted, SPEC_06 closed.
