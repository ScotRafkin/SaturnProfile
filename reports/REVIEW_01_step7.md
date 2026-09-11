# REVIEW 01, Step 7. The wind tool for Lindal, kind W

Review of `reports/REPORT_01_step7.md` against SPEC_01 v0.8 Step 7. Reviewer: S. Rafkin,
11 September 2026. **Disposition: returned for rework under SPEC_01 v0.10.** The code did
what the specification said; the specification was wrong in one place that matters (the pole)
and weaker than it should be in another (the wind as bin means). Both are fixed in v0.10, and
Step 7 is to be rebuilt to it. SPEC_00 goes to v0.8 for the polar rule. Nothing is committed
from the current tree except the `control.py` change (item 7).

## The two rulings that change the design

**1. The wind is zero at both poles. Not a rule of thumb, a requirement.** Finding 5 is
correct in every particular, and the Step 4 review said this would be where it was caught.
The v0.8 extension rule (hold the last bin mean to the pole) was the defect, not the tool.
SPEC_00 §6.6 now enforces at load time that both poles are latitude nodes and `u_total_ms` is
exactly zero there, for every kind W file, reduction or forward. The tool brings the source to
zero under a declared, flagged rule; the model never patches it.

**2. The wind is a fitted function of latitude; the bins give the uncertainty.** Bin means
as the product were the wrong choice: they are discontinuous, they cannot be sampled off the
bin centers without a second rule, and they have no natural way to reach zero at the pole.
v0.10 replaces them with a cubic penalized B-spline in planetographic latitude, fitted to all
323 points with `u(±90°) = 0` as equality constraints in the solve, knots at the declared
spacing, and the smoothing parameter set by a scatter-matching rule: the RMS residual of the
points about the fit equals the pooled within-bin scatter, so the fit follows the jets and does
not absorb digitization noise. The bins remain in the file (`bin_mean_ms`, `bin_count`) and
their job is the uncertainty: the RMS of point residuals about the fit per populated bin, NaN
where the count is below the minimum, NaN in the polar caps. The polar caps (poleward of 81 N
and 73 S) are the constrained spline's decay to zero and are flagged `extrapolated`; the
report is to quantify what that choice costs at the equator against a linear taper, since the
current report already shows it is worth about 3 km. The fit is evaluated on a declared
planetocentric grid by the direct conversion `planetographic_from_planetocentric`, so the tool
no longer needs the fixed point at all.

## Rulings on the remaining findings

3. **Equatorial bins (finding 1).** Agreed; the criterion applies to bins at or above the
   minimum count. The single-point bin at 3 N is exactly what the flag is for.

4. **38 N (finding 2).** Agreed; the count criterion is dropped. What matters at the anchor
   is the fitted wind at `φ_g` = 36.3°, which enters the frozen `φ_c` and `r0`; v0.10 asks for
   that value, expects it between −20 and +25 m/s, and wants it quoted in the report.

5. **36.3° conversion (finding 3).** Agreed; reworded to the conversion of 36.3° itself,
   which the tool reproduces to 0.0 and which is the only version of the check that means
   anything.

6. **Zero shear (finding 4).** The report is right that the shear cannot be nonzero from this
   file's inputs, and right about why: the decomposition needs `r(φ, p)`, the model owns that
   mapping, and the wind tool is not to grow a hydrostatic integration to compute it. The v0.5
   expectation is withdrawn. For an altitude-independent field on a pressure coordinate the
   decomposition is trivial by declaration, `u_cylindrical = u_total`, `u_shear = 0`, and the
   file says so in a `decomposition` attribute (SPEC_00 v0.8). The decomposition becomes real in
   the generic wind tool of the forward-inputs spec, where a sheared component is supplied.

7. **`control.py` path resolution by named keys.** Correct fix, and the right lesson. Commit
   it with Step 7.

8. **Decision 4, numbers in `[wind]`.** Bin width, minimum count, knot spacing, grids: choices,
   and a control section may carry them; the principle 2 refusal is for keys naming physical
   quantities, and v0.10 says so. The observation level and its justification are a property
   of the Smith data, not a choice of this build, so they move to
   `data_static/winds/smith1982_fig4.toml` beside the digitization and the control file points
   to it.

9. **Decision 3, polar radius from the raw bundle.** Acceptable until kind D exists at Step 9;
   the pointer is declared and hashed, which is what matters.

10. **Figures.** Agreed that a report's figure must outlive its step directory. §0 now sends
    them to `reports/figures/step<N>_<name>.png`, committed with the report. The v0.9
    two-panel dynamical height figure (fitted wind above, `h(φ)` from the march against
    Lindal's Eq. 18 below) was not produced because the step ran against v0.8; it is part of
    the v0.10 acceptance.

## On the 7.86 km departure

The report weighed three candidates and missed the fourth, which is the largest. The march is
anchored at Lindal's fitted polar radius, 54,438 ± 10 km, and Eq. B3 carries that anchor
straight to the equator, so the comparison band against 60,367 ± 4 km is about ±11 km, not
±4. A departure of 7.86 km is inside it. That does not make the other candidates wrong (the
polar cap alone is worth about 3 km, and the wind Lindal actually used is not stated), but it
does mean the number is not evidence of a defect anywhere. v0.10 states the band. With the
fitted wind the number will change; report it, do not tune it.

## Actions

- Rebuild Step 7 to SPEC_01 v0.10 and SPEC_00 v0.8: the fit, the constraints, the
  scatter-matching rule, the bins as uncertainty, the polar zero check in `read`, the static
  `smith1982_fig4.toml`, the `decomposition` attribute, both figures under `reports/figures/`.
- Report the fitted wind at 36.3° N, the smoothing parameter with its sensitivity, the wind
  geoid equatorial radius with grid convergence and the polar cap sensitivity, and the
  Eq. 18 comparison.
- Then a new REPORT_01_step7 (replace, do not append); Step 8 waits on it.

## Addendum, 11 September 2026: rulings since v0.10, and what to build

The review above was written against v0.10. The coding agent then raised two questions on the
v0.11 spline rules and stopped, correctly. Since then the design has moved twice and the
build is to follow **SPEC_01 v0.15** and **SPEC_00 v0.8**. In order:

1. **The two v0.11 questions.** Both findings were right: check (a)'s local clause is wrong
   where the points within the window lie on one side of the latitude, and the smoothing
   ladder ran the wrong way for the over-smoothing detectors. Both are corrected in the spec
   (v0.13) for the record. They no longer bear on the build, because the spline is not the
   method.

2. **The method is no longer a fit at all (v0.14, v0.15).** The wind is the solid curve of
   Ingersoll and Pollard (1982) Fig. 5, digitized and placed at
   `data_static/winds/ingersoll_pollard1982_fig5_curve.csv` with its `.note.md`, the
   digitization script, and the overlay. Lindal cites that paper alongside Smith as his wind
   source; the curve is the authors' own smoothing of the Smith data, and the Smith points
   scatter about it with RMS 21.3 m/s, equal to their pooled within-bin scatter. The Smith
   points keep exactly one job: the uncertainty, as the RMS of the points about the curve per
   2° bin. `fit.py` stays in the module as the non-default alternative for a data set with no
   published curve (`wind_source_kind = "points"`), with the v0.13 corrections; it is not
   exercised by this build and needs no further work now.

3. **Assembly rules, all declared in `[wind]` and all flagged in `value_provenance`.**
   PCHIP through the digitized samples per segment. Ring gap: the northern segment reflected
   about the equator (their dashed curve), the ±1.3° strip bridged by a cubic Hermite (which
   reduces to `a + b φ²` by symmetry), and the reflection blended linearly into the southern
   segment over `join_window_deg = [-10.9, -15.0]` because the two do not meet (377 against
   346 m/s at −10.9°). Polar caps: PCHIP from the last two samples to exactly 0.0 at ±90°,
   flagged `extrapolated`. The §6.6 polar zero check goes into `lib.schema` and `read`
   refuses a file that fails it. The `.note.md` records that the ten Smith points on the
   southern flank lie 80 to 100 m/s below the reflection; that is not to be "fixed", it is a
   declared alternative (`gap_rule = "reflect_north_then_bins"`) for the sensitivity run.

4. **What survives from the queued v0.10 rework unchanged:** the static
   `smith1982_fig4.toml` for the observation level and its justification; the `[wind]`
   section carrying choices only; grid evaluation on the declared 0.5° planetocentric grid
   through the direct `planetographic_from_planetocentric` conversion (no fixed point in the
   tool); the trivial decomposition with its `decomposition` attribute; provenance code 4 on
   the extended levels; the raw bundle as the declared pointer for the polar radius; the
   `control.py` named-key change; both figures under `reports/figures/`; the entry point
   renamed `casspian-wind-from-curve`.

5. **Acceptance and report.** As v0.15 Step 7 states: curve through every sample, zero at
   the poles and the refusal test; the values at 36.3° N (1.9 ± 3 m/s, to be quoted, it
   enters the frozen anchor) and 30.8° N; the bin table and pooled RMS; the wind geoid
   equatorial radius against 60,367 km with the ±11 km band, grid convergence, and the two
   sensitivities (polar rule, gap rule) each as a difference in km; the Eq. 18 comparison
   figure. A new REPORT_01_step7, replacing the current one.

Nothing from the current working tree is to be committed until that report is accepted,
except that the `control.py` change may go in with the step.

## Second review, 11 September 2026: the v0.15 report

Report of eight of eight against v0.15. **Disposition: accepted with one change to `lib.geoid`
and a rerun of the geoid checks under SPEC_01 v0.16.** The wind file itself is accepted as
built; nothing in `curve.py` or `build_wind.py` changes.

**The equatorial step in Figure 2 is a construction error, not a wind result.** `wind_geoid`
marches Eq. B3 from each pole to the equator with the same polar radius at both ends. Eq. B3
is a first order equation with one constant, so that is two surfaces, and with a wind that is
not symmetric about the equator they cannot meet: the northern march arrives at 146 km of
dynamical height, the southern at 108, and the figure shows the 38 km step. The Step 5 spec
said "from each pole inward" and the author's own Eq. 18 curve in the same figure, integrated
continuously, already shows the answer (27 km at the south pole). An independent march of the
same wind with Null's set confirms every number: from the north pole at 54,438 km the equator
is 60,390.29 km (report: 60,390.31) and the south pole arrives at 54,466.7 km; from the south
pole the equator is 60,351.9 km. The surface has one constant and it has to be fixed once.

**Ruling (v0.16).** `wind_geoid` becomes one continuous march from the north pole to the south
pole with a declared `anchor_rule`. Default `"mean_polar_radius"`: the north polar start is
found by a one-parameter root find so that the two polar radii average to `r_anchor`, which is
what Lindal's phrase "mean polar radius" says he did. Alternatives `"north_pole"`,
`"south_pole"`, and `"latitude"`. The function returns both polar radii and their difference,
the polar asymmetry the wind produces, which is itself a result. The Eq. 18 curve is anchored
the same way (offset so its polar values average to zero) so the two curves in Figure 2 are
comparable.

**What that does to finding 4.** With mean polar anchoring the independent march gives an
equatorial radius of 60,371.1 km, north and south polar radii 54,423.7 and 54,452.3 km, a
bulge of 127 km over the no-wind 60,244 (Lindal's own two surfaces imply 123), and a departure
from 60,367 of **+4.1 km**, inside the ±11 km band. The 23 km departure of the report was
mostly the anchoring, not the wind representation; the report's inference that the wind
representation dominates by an order of magnitude was drawn from two builds that differed in
anchoring as well as wind and does not survive. The polar asymmetry of 28.7 km is the number
to carry forward: it says the wind Lindal fit was either more symmetric than this curve or
that his mean polar radius absorbed the difference. These are expected figures for the rerun,
to be reproduced by the code, not prescribed.

**Other rulings.** Finding 1: agreed, the sample check excludes the join window; spec reworded.
Finding 2: agreed, the `.note.md` now says 95 to 160 m/s against the assembled curve; the gap
bins carrying that as their uncertainty is right. Finding 3: correct and useful; `u sin φ`
kills the gap at the equator, and the same weighting is why the polar caps and the mid
latitude jets matter more than they look. Finding 5: the sign and size of the Eq. 18
difference are as expected; the observation that Eq. 18 lands on Lindal's implied 123 km is
noted but not to be read as evidence of his method, since with the corrected anchoring the
full march lands within 4 km of it too. Decisions 1 to 3: fine; the NumPy 2 scalar note is
worth keeping in the docstring.

**Actions.** Change `wind_geoid` per v0.16; rerun the Step 5 acceptance (the grid convergence
and the 45° checks are unaffected, the closure diagnostic now against the anchored surface)
and the Step 7 geoid checks with the anchor rule sensitivity added; regenerate Figure 2; a
short addendum to REPORT_01_step7 with the new table (no new report). Then commit Step 7,
rebuild `lindal_wind.nc` clean, and proceed to Step 8.
