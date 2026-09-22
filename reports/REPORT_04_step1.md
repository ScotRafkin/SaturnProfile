# REPORT 04, Step 1. The reference surface through the anchors, and the wind on the model's geometry

Coding agent, 21 September 2026; refreshed 22 September under REVIEW_04_step1 and SPEC_04 v0.10.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.10, sections 3 and 12.
Acceptance script `reports/step04_1/accept_step04_1.py`, output `reports/step04_1/output.txt`.
Regression `reports/step04_1/run_regression.sh`, output `reports/step04_1/regression.txt`.

**Thirteen of thirteen checks pass.**

The first filing was 12 of 14, the two failures both on the cylinder-extended wind. The review
verified the finding and withdrew the specification's values: the reviewing agent could not
reproduce 9.371 and -0.133 either, under gauge pinning or reference-level pinning, with a linear or
a PCHIP wind. Decision P at v0.10 fixes the construction (pinned on the file's reference level, the
library's radial columns under the flat-isobar map, a fixed point with the reference surface) and
places the file at Step 2, where the columns it needs exist. **Check 9's cylinder half and the
first filing's check 14 moved to Step 2 with the construction they test** (SPEC_04 section 12
ruling 2). Nothing was loosened: the closure half of check 9 is unchanged and still bounded at
1e-3 m/s, and check 13 now exercises `wind_on_mesh` and the pressure derivative on a minimal
synthetic sheared field rather than on the withdrawn construction.

Ruling 4 was applied: `step04_0`'s check 0 now takes a `-dirty` copy it writes itself as its
subject. Both suites were rerun. This section, section 3, finding 2 and section 5a are the refresh
the review's order of work item 1 asked for.

The sweep of item 4 then turned up one more, finding 5 in section 5c: `step04_0`'s check 6 counted
this step's two new product variables as changed values, because its list of admitted additions
names only Step 0's. The list was extended to name them, with the citation, and the finding is
recorded for ratification because it edits a closed step's acceptance record.

## 0. Before this step

HEAD `bd98c60`, Step 0 accepted at `1c9b310` with the sweep recorded at `f38c58b` and the deferred
suites at `bd98c60`. Every product on disk is the swept one, clean at `1c9b310`. This step changes
no input file, so nothing is relaxed anywhere: every check and every regression suite reads the
registered products as they stand.

It does change the kind `profile` product, which gains `u_column_ms` and its NaN companion, so by
SPEC_03 section 0 the acceptance runs `produce` in memory and compares against the registered
closure product on disk rather than rewriting it; the registered product is rebuilt on the clean
tree at the sweep, after the review.

## 1. What was built

**Deliverable 1, `lib.geoid.through_anchor` and `AnchorSurface`.** The Eq. B3 march with the
constant fixed by a stated point `(phi_i, r0_i)`: the march starts there and runs outward to each
pole, on the same dense grid `wind_geoid` uses, unioned with the caller's nodes, the equator and
the anchor latitude, so every latitude the caller asked for and the anchor itself are marched and
never interpolated. No root is found and no tolerance is declared. `AnchorSurface` carries the
radius on the caller's nodes, the march's own nodes and radii, the two polar radii, the asymmetry
and the equatorial radius, and two readers: `radius_at_latitude`, by the module's `radius_at` rule,
and `residual(phi, r)`, which is `r_measured - r0(phi)`, what the specification records as
`reference_surface_residual_<slug>` when a second occultation anchor is listed.

**Deliverable 2, `lib.windfield.WindField`.** The one place the kind W grid is read off. Built once
from the dataset, it sorts both axes and holds `u_total` on them, and answers:

* `wind_at(phi, p)`, linear in planetocentric latitude and linear in `ln p` (decision L), the two
  arguments broadcasting against each other;
* `wind_derivatives(phi, p)`, `(du/dphi)_p` and `(du/dln p)_phi` of that same interpolant, each
  constant across the cell in the direction it differences, the cell at a node being the one to the
  north or to higher pressure;
* `wind_on_mesh(latitude, pressure)`, `u` at every node of a set of columns;
* `reference_wind(phi)`, `wind_at` at the file's reference level, which is what the reduction
  anchored on and is provided so that the two read one field.

No coverage attribute and no declaration is read: the refusal is on the grid's own extent, and a
point outside it raises `ValueError` as `lib.composition.column_at` does for the same reason.

**Deliverable 3, the anchors' geopotential under the run's wind.** `produce` gains `u_column`, the
zonal wind in m/s on the profile's levels; omitted, it is the reference-level wind at the column's
latitude at every level, which is what closure mode means and what SPEC_03 formed. `Production`
gains `u_column_ms`, and the kind `profile` product gains `u_column_ms(level)` with its NaN
companion `u_column_uncertainty_ms` beside `u_at_phi_c_ms`, which keeps its meaning: the
reference-level value at `phi_c`. Nothing below `produce` changed; `geopotential_along_profile`
already carried `u` through element-wise arithmetic.

## 2. Decisions

Every choice the specification did not make.

1. **`through_anchor` is a function beside `wind_geoid`, not a new `anchor_rule` string inside it.**
   The specification calls it an anchor rule, and it is one, but `wind_geoid`'s body is organized
   around a secant on the north polar start, which this rule does not need: the constant is known
   before anything is marched. The two are compared node for node in check 4, which is what ties
   them to one surface.
2. **`lib.windfield` is one class, `WindField`, not free functions on the dataset.** The outer loop
   of SPEC_00 section 7.2 asks for the whole mesh once per pass and the file's grid does not change
   between passes, so the sorting and the grid are prepared once. The three methods carry the names
   the deliverable gives them.
3. **`wind_on_mesh` takes the column latitudes and the pressure map, not a mesh object**, which
   arrives at Step 2. What the deliverable needs of a mesh is those two things, and taking them
   directly keeps this step from implementing ahead.
4. **Coverage is the grid's own extent.** The deliverable says no declaration is read, and the
   `coverage_latitude_planetocentric_deg` and `coverage_pressure_Pa` attributes are declarations.
   Measured on both wind files: the declared coverage is `[-90, 90]` deg and `[1, 1e6]` Pa and the
   grid reaches exactly those, so the two agree here and the choice is visible only on a file where
   they would not.
5. **The refusal is `ValueError`**, the class `lib.composition.column_at` raises for the same
   refusal, extrapolation under decision N. It is not a schema fault and not a namelist fault.
6. **`wind_derivatives` returns the interpolant's exact partials**, not a finite difference of
   `wind_at` and not a separate fit: `(du/dphi)_p` is the difference of the two `ln p` interpolated
   columns over the latitude interval and `(du/dln p)_phi` the difference of the two latitude
   interpolated rows over the `ln p` interval. Each is constant across the cell in the direction it
   differences, which is what "piecewise constant" fixes, and each is exactly the derivative of the
   field `wind_at` returns. Check 12 measures the first against a central difference of `wind_at`.
7. **The node rule is `searchsorted(side="right") - 1`, clamped at the last node.** A point exactly
   on a node takes the interval to its north in latitude and to higher pressure in pressure, which
   is the rule the deliverable fixes. At the last node of either axis there is no such interval and
   the previous one is used, which is the only one there; the value `wind_at` returns is the node's
   own either way, so only the derivative sees it.
8. **`produce`'s new argument is `u_column`, defaulting to `None`.** It is checked against the
   profile's level count and otherwise taken as given. The specification's wording, "the scalar it
   forms today being the default", is implemented as the default being the scalar broadcast to the
   levels rather than a separate code path, so closure and transfer go through one line.
9. **`u_column_ms` carries a NaN uncertainty companion.** Kind `profile` requires one for every
   variable whose provenance is `modeled` (SPEC_03 decision 5, enforced in `lib.schema`), so
   `u_column_uncertainty_ms` is written NaN like the rest. `schema.PROFILE_MODELED_VARIABLES`,
   which is the written record of that set, gains `u_column_ms` with it.
10. **`through_anchor` snaps a latitude node a rounding step beyond a pole to the pole, and refuses
    one genuinely off the planet** (more than 1e-9 rad, about 6e-8 degrees, outside). Forced by
    finding 1.
11. **Superseded by decision P.** As filed, the extension was pinned on the gauge isobar, where
    the reference surface sits, with `u` constant on cylinders and equal to `u_reference` there.
    Ruling 2 fixes the pin on the file's reference level instead, on the ground that a kind W file's
    three parts make `u_total(phi, p_ref) = u_reference(phi)` an identity, so the gauge pinning put
    the reference-level wind at 0.917 m/s at the anchor and was not "the same reference-level wind".
    The construction is decision P's and is built at Step 2.
12. **The march is supplied `u_total` at the gauge isobar's pressure**, as the deliverable says. On
    the closure wind, which has no shear, that value at the anchor is bit-identical to the
    reference-level wind the reduction used, `u_at_anchor_ms = 2.16708095051181` m/s, so the
    surface this step marches and the reduction's are the same surface (checks 4 and 9).

## 3. Acceptance results

| # | Check | Result |
|---|---|---|
| 1 | the equatorial radius to 0.01 m | PASS, difference -1.05e-4 m |
| 2 | the polar radii and the asymmetry (28,740.37 m) to 0.1 m | PASS, -0.036, +0.035, +0.0003 m |
| 3 | the values at 10, 20, 45 and 60 N to 1 m | PASS, largest 0.034 m |
| 4 | agrees with the reduction's march at every dense node to 0.01 m | PASS, largest 1.05e-4 m |
| 5 | a second anchor entry on the same file gives a residual of zero | PASS, exactly 0 |
| 6 | a synthetic anchor moved by 1 km gives a residual of 1 km | PASS, exactly 1000 m |
| 7 | `wind_at` reproduces the file's nodes exactly | PASS, bit for bit, 22,021 nodes |
| 8 | `wind_at` refuses a point outside the file's coverage | PASS, four cases refused |
| 9 | the wind along the anchor's column under the closure wind to 1e-3 m/s | PASS, largest 8.1e-5 m/s |
| 10 | the closure rerun through `produce` with the wind array is bit-identical | PASS, bit for bit |
| 11 | beyond the specification: the surface passes through the stated point exactly | PASS, exactly 0 |
| 12 | beyond the specification: `wind_derivatives` is the interpolant's, and the node rule | PASS |
| 13 | beyond the specification: `wind_on_mesh`, and the pressure derivative on a sheared field | PASS, bit for bit |

Check 9's cylinder half and the first filing's check 14 are Step 2's at v0.10 and are not in this
table. What they measured is in finding 2, kept as the record of why they moved.

**The surface (checks 1 to 6, 11).** Through `(30.805568 deg, 58,516,188.288885 m)`: equator
60,366,999.9999 m, north pole 54,420,643.4644, south pole 54,449,383.8347, asymmetry 28,740.3703 m,
and 60,128,612.9664, 59,492,075.6754, 57,053,249.5161 and 55,675,128.9158 m at 10, 20, 45 and
60 N. Every one reproduces the specification's value and the reduction's own record. Against
`wind_geoid` under the `equatorial_radius` rule on the same 3,607 dense nodes the difference is one
sign throughout, +7.84e-5 to +1.05e-4 m, and at the equator it is +1.051947e-4 m, which is the
anchoring residual the reduction recorded, +1.050457e-4 m: the anchor radius this surface is
marched through was read off that march and carries that residual, while a fresh `wind_geoid` run
converges to +1.4e-7 m. So the two marches differ in the constant and in nothing else, which is
what check 4 was for, and `through_anchor` passes through its stated point exactly (check 11) where
the rule with a root find reports a residual.

**The wind (checks 7, 8, 12, 13).** `wind_at` returns all 22,021 of the closure file's nodes bit for
bit, and at the anchor its reference-level value is bit-identical to
`refrac.anchor.wind_of_latitude`. `(du/dphi)_p` in the cell bracketing the anchor is
-691.638087 m/s per radian, -12.07136 m/s per degree, which is the 30.5 to 31.0 degree interval
where the anchor sits; at the 30.5 degree node itself the value is that interval's and just south
of it the value is the southern interval's, -889.062216 per radian, so the rule of deliverable 2 is
the one in force. `(du/dln p)_phi` is exactly zero on both wind files, which carry no shear.

**The geopotential (check 10).** `produce` with the wind array, the run's own wind along the column
at 2.16708095051181 m/s at every level, returns the seven level arrays and the geopotential, its
increments, `|g_eff|` and `psi` bit-identical to the scalar path, and `pressure_Pa`,
`temperature_K`, `number_density_m3`, `mean_refractivity_m3`, `mean_molar_mass_kg_mol` and
`geopotential_m2s2` bit-identical to the registered closure product on disk at `1c9b310`. Closure
mode is unchanged, as the Appendix requires.

## 4. Findings

Findings 1 to 3 come from the acceptance. Finding 4 came from the regression and finding 5 from the
post-sweep rerun; each is recorded with the run that found it, in sections 5a and 5c.

**1. The dense march grid overshoots the pole, and the wind is then asked for a latitude kind W
does not cover.** `wind_geoid` builds its march grid as
`np.radians(np.arange(-90.0, 90.0 + 0.5 * march_step_deg, march_step_deg))`, and `arange`
accumulates: at a 0.01 degree step the last node comes out at 90.00000000004604 degrees. RK4 then
evaluates its midpoints against that node and asks `u_of_phi` for a latitude past the pole. Under
`refrac.anchor.wind_of_latitude` this has never shown, because that callable is `np.interp`, which
clamps silently; under `lib.windfield`, which refuses to extrapolate, it raises at once, which is
how it was found. `through_anchor` now snaps such a node to the pole and refuses one genuinely off
the planet (decision 10). **`wind_geoid` is left as it stands**: it is accepted code, the defect
changes no value it has ever produced, and touching it is the author's call. The same applies to
any caller that builds its own latitude grid with `arange`.

**2. The cylinder-extended wind's stated column values could not be reproduced, and the
specification did not state the construction. Ruled: the values are withdrawn, the construction is
decision P's, and the two checks move to Step 2.** Kept here as the record of why.

SPEC_04 v0.9 section 1 said only "the same reference-level wind extended along cylinders by the
acceptance script with the library's geometry", and Step 1 deliverable 3 stated the result: "it runs
from 2.167 m/s at the reference level to 9.371 m/s at the top level and -0.133 m/s at the bottom,
moving `Phi` at the top by -271.6 m2/s2".

Two conditions fix the field. `S = 0`, which is what the state is for, requires `u` to be constant
on cylinders coaxial with the rotation axis, so `u = U(s)` with `s` the distance from the axis; and
"the same reference-level wind" requires `U(s_ref(phi)) = u_ref(phi)`, which determines `U` once the
surface `u_ref` sits on is named. The specification did not name it. As filed this step pinned it on
the gauge isobar, where the reference surface sits, and measured on the anchor's column 2.1671 m/s
at the gauge, **7.6364** at the top, **0.7404** at the bottom and a top geopotential shift of
**-206.65** m2/s2, against the stated 2.167, 9.371, -0.133 and -271.6. The largest wind departure
was 1.73 m/s against a bound of 1e-3.

The ratios were the useful part, and the review confirmed it: 1.304 at the top and 1.297 at the
bottom on the cylindrical displacement, 1.314 on the geopotential shift, one factor near 1.30 rather
than a difference in kind, and neither `1 / cos phi_c` (1.164) nor `1 / cos^2 phi_c` (1.356). The
review identifies it as, to a few percent, the height above 1 bar over the height above the gauge at
the top of the profile, which is what pinning on the 1 bar level does to the displacement.

**Ruled (section 12.2).** The reviewing agent rebuilt the construction from scratch and could not
reach 9.371 or -0.133 either, under gauge pinning with the field-line height (7.636, this step's
number), gauge pinning with the radial column (7.465), or 1 bar pinning with the radial column
(9.490), nor under any of the three with a PCHIP wind, and no longer holds the code that produced
them. The v0.9 values are withdrawn. Decision P pins the extension on the file's reference level,
because a kind W file's three parts make `u_total(phi, p_ref) = u_reference(phi)` an identity; the
geometry is the library's radial columns under the flat-isobar map; and the reference surface and
the file are a fixed point, converged in the script. Since the columns are Step 2's, the file is
built there and these two checks are exercised there, against 9.490 m/s at the top, 3.717 at the
gauge, 2.168 nearest 1 bar and 1.926 at the bottom, and `Phi` lower at the top by 322.3 m2/s2 and
higher at the bottom by 18.0.

**3. Not a defect, recorded because check 10's strength rests on it.** `wind_at` reaches
2.1670809505118100446 m/s at the anchor, bit-identical to `refrac.anchor.wind_of_latitude`, through
two blends that are not obliged to land there. The pressure blend is `a (1 - t) + a t` over two
nodes holding the same value, which returns `a` exactly here; and the latitude blend is
`south (1 - t) + north t` against `np.interp`'s `y0 + (y1 - y0) t`, a different arrangement of the
same line, which agrees to the last bit here. Both were measured, neither is guaranteed in general,
and a wind file with shear will not have the first. So check 10's "bit-identical" is a measurement
on this file and this anchor, which is exactly what the Appendix's "closure mode and its expected
values unchanged" asks for, and not a claim that any wind array reproduces the scalar path
bit for bit.

## 5. Regression

`reports/step04_1/run_regression.sh`, output `reports/step04_1/regression.txt`. Every accepted
suite runs directly on the registered swept products: this step changes no input file, so there is
no candidate, no relaxation and no group that has to wait for a sweep.

The working tree carries this step's four source paths throughout, so every product a suite
rebuilds is written `-dirty`, and `refrac` and `forward` refuse a `-dirty` input. Steps 02_4 to
02_6 rebuild the registered kind N and `step03_4` runs `forward` on it, so the registered product,
its figures, the closure run and the transfer run's inputs are restored from set-aside copies after
**every** suite and `reports/figures/` is restored from git the same way. That is the Step 0 lesson
of REPORT_04_step0 section 5 applied at the start rather than after two failed rounds.

| Suite | Result | Reference |
|---|---|---|
| `step1/accept_step1` | 6 of 6 | 6 |
| `step1/verify_review_changes` | 14 of 14 | 14 |
| `step2/accept_step2` | 8 of 8 | 8 |
| `step3/accept_step3` | 5 of 5 | 5 |
| `step4/accept_step4` | 7 of 7 | 7 |
| `step5/accept_step5` | 7 of 7 | 7 |
| `step6/accept_step6` | 6 of 6 | 6 |
| `step7/accept_step7` | 8 of 8 | 8 |
| `step8/accept_step8` | 9 of 9 | 9 |
| `step9/accept_step9` | 6 of 6 | 6 |
| `step02_1/accept_step02_1` | 6 of 6 | 6 |
| `step02_2/accept_step02_2` | 7 of 7 | 7 |
| `step02_3/accept_step02_3` | 9 of 9 | 9 |
| `step02_4/accept_step02_4` | 9 of 9 | 9 |
| `step02_5/accept_step02_5` | 7 of 7 | 7 |
| `step02_6/accept_step02_6` | 7 of 7 | 7 |
| `step03_1/accept_step03_1` | 9 of 9 | 9 |
| `step03_2/accept_step03_2` | 8 of 8 | 8 |
| `step03_3/accept_step03_3` | 16 of 16 | 16 |
| `step03_4/accept_step03_4` | 13 of 13 | 13 |
| `step04_0/accept_step04_0` | 15 of 15, after the change of ruling 4 | 15 |
| `step04_1/accept_step04_1` | 13 of 13, after the changes of rulings 2 and 4 | this step |

The first filing was twenty of twenty-two, `step04_0` at 14 of 15 on finding 4 and `step04_1` at 12
of 14 on finding 2. Both were reworked under the rulings and rerun on the same tree; the other
twenty suites were untouched by either change and are not rerun. The registered kind N restored to SHA-256
`920304440ee4554648c3aa6321b713741ee26eb0eb1ff223f5a81a8f9b20fd37`, the swept value, and the
closure run's inputs and product and the transfer run's inputs all restored equal. Porcelain at the
end holds this step's four source paths and this report, and nothing else; the "paths dirty" counter
in the output rises from 4 to 5 partway down because this report was written while the run was in
progress, not because a restore failed.

## 5a. Finding 4, from the regression, and ruling 4 applied

**The finding.** `step04_0`'s check 0 asserted that the `-dirty` refusal fires on the working tree.
At Step 0 that held, because the step had rebuilt every product on the tree and each carried
`-dirty`. The Step 0 sweep then rebuilt them all clean at `1c9b310`, so there was nothing on disk
for the refusal to fire on and the check reported `refused before the relaxation: None`. The
relaxation itself was still announced and still patched only `control._refuse_dirty_commit`, so
nothing was hidden; the check had simply lost its subject. `step04_0` did not meet this in Step 0's
own regression because it ran in group 1, before the sweep, and was not one of the seven suites
deferred to after it.

**Ruled (section 12.4): the proposal is adopted, applied at this step's acceptance commit.** Check 0
now writes its own subject: a copy of the registered transfer wind under
`reports/step04_0/cases/dirty/`, with `casspian_git_commit` rewritten from the swept value to that
value plus `-dirty`. It then asserts two things, that `control._refuse_dirty_commit` refuses that
copy and that it does not refuse the same file with the clean commit, so the check means "the
refusal fires on a `-dirty` input and only on one" whatever state the tree is in. What the working
tree gives is still printed, labeled for information and not asserted, so the tree's state stays
visible without the check depending on it. The relaxation that follows is unchanged and still
named in the output.

**Measured after the change.** `step04_0` is back to 15 of 15. Check 0 refuses the `-dirty` copy
with `ControlFileError: lindal_transfer_wind.nc carries casspian_git_commit =
'1c9b310018979bb5343501a5a1dc6357199481f1-dirty'`, does not refuse the same file with the clean
commit, and prints that loading the run from this working tree refuses nothing, which is the state
the sweep left and the reason the old form of the check could not pass. Output
`reports/step04_1/rerun_accept_step04_0.txt`.

## 5c. The suites that read the closure product, on the swept product

`reports/step04_1/post_sweep_suites.sh`, output `reports/step04_1/post_sweep.txt`. The four suites
that read the closure product, rerun on the swept one, `fcc2e2c4`. `step02_1` and `step02_6` are not
here: they look up recorded hashes, but only the reduction chain's, which this step does not touch.

| Suite | Result | Reference |
|---|---|---|
| `step03_3/accept_step03_3` | 16 of 16 | 16 |
| `step03_4/accept_step03_4` | 13 of 13 | 13 |
| `step04_0/accept_step04_0` | 15 of 15, after the change of finding 5 | 15 |
| `step04_1/accept_step04_1` | 13 of 13 | this step |

Each was run from a clean tree and the tree restored after it, and the porcelain line after every
suite is empty. The registered kind N and the swept closure product both restored equal, the product
still `fcc2e2c4ae71554537d8aa53f07b247724b718a371abab73de31939bae9ae6af`.

**Finding 5, from this run.** `step04_0`'s check 6, "the cascade rebuilt every product and changed
no value, group by group", failed on the closure product with
`unexpected new variables ['u_column_ms', 'u_column_uncertainty_ms']`. The check compares each
rebuilt product against a copy taken before Step 0 and admits a named list of the additions Step 0
intended, `u_reference_ms` and `latitude_planetocentric_deg`; anything else it reports as a changed
value. Step 1's two variables are exactly what deliverable 3 required a later step to add, so the
list was extended to name them, with the reason and the citation in the code. Nothing else is
admitted, so a value that does change still fails, and every other product in the check reported
`unexpected differences none`. Measured after the change, `step04_0` is 15 of 15 with no failing
check, and the closure product's line reads `variables added ['latitude_planetocentric_deg',
'u_column_ms', 'u_column_uncertainty_ms', 'u_reference_ms'], removed ['u_cylindrical_ms'],
attributes removed ['decomposition', 'decomposition_geometry',
'latitude_planetocentric_absent_meaning']; unexpected differences none`. Output
`reports/step04_1/post_sweep_accept_step04_0.txt`. This is the same operation Step 0 ruling 3 performed on `step03_3`'s
baseline and the same class as Step 0 findings 5 to 8: it edits a closed step's acceptance record,
so it is recorded here for ratification rather than passed over.

The check also prints, and does not fail on, `closure_dropped_wind` losing `decomposition_geometry`
from the list of attributes the closure comparison dropped. That is Step 0's retirement of the
attribute, already tolerated by `provenance_attr` and printed rather than waved past, and it is
unchanged by this step.

## 5b. Cost, measured

Recorded because the regression's wall time invites a wrong inference about model cost. The 39
minutes the regression took is test scaffolding: `step04_0` alone is 5 min 15 s because it spawns a
full nested run of `step03_3`, which is itself 2 min 35 s and spawns a subprocess per refusal case,
and every one of those pays a fresh interpreter start and a fresh input load. Model cost, measured
on this machine:

| Quantity | Cost |
|---|---|
| `import casspian` | 1.00 s, once per process |
| `ctl.load_run_inputs`, every file with its schema checks, hashes and closure comparison | 7.60 s, once per run |
| `WindField(wind)` | 0.1 ms, once per run |
| `through_anchor`, 0.05 deg march, 3,607 nodes | 2.13 s, once per run |
| `wind_on_mesh` over 417 columns by 780 levels, 325,260 nodes | 34 ms |
| `wind_derivatives` over the same mesh | 53 ms |
| `produce`, one column of 66 levels | 1.7 ms |

Two of these bear on Step 1 directly. **`through_anchor` is 6.4 times cheaper than the rule it
replaces**: 2.13 s against 13.60 s for `wind_geoid` under `equatorial_radius`, because the secant
on the polar start costs seven or eight marches and this needs one. And **`wind_at` costs 17.9 us
per scalar call against 1.8 us for the `np.interp` callable the reduction uses**, ten times more,
from the broadcast, the `searchsorted`, the clip and the coverage refusal; the RK4 loop makes 14,424
such calls, so the march is 2.13 s with `wind_at` against 1.79 s with the cheap callable. That is
the price of refusing to extrapolate rather than clamping silently, which is finding 1's whole
point, and it is 16 percent of one march. A scalar fast path would recover it and is not proposed
here.

The march is the per-draw cost that matters later: 1.79 s of it is the RK4 Python loop, 124 us per
slope evaluation, which is `lib.gravity` on scalars and not the wind. It cannot be vectorized over
nodes, since an initial value problem needs the previous result, and SPEC_00 section 3.1 already
allows the loop for that reason. The 7.60 s input load is per process and hoists out of any wrapper
that loads once and draws in memory. Nothing is proposed; the numbers are recorded so the Monte
Carlo wrapper is specified against measurements.

## 6. Hashes, the Step 1 sweep

Step 1 changed one product: the kind `profile` closure product, which gains `u_column_ms` and its
NaN companion. Nothing else it touched enters a file, so only that product was rebuilt, by
`reports/step04_1/sweep.py` on the clean tree at the Step 1 acceptance commit `780a2de`, output
`reports/step04_1/sweep.txt`. The script refuses to start on a tree whose porcelain output is not
empty and reads `casspian_git_commit` back from the file it wrote. Every other product of the Step 0
sweep is untouched by this step and stays as that sweep left it, clean at `1c9b310`, so
**REPORT_04_step0 section 6 remains the record for the other sixteen** and the row below is the one
that changes. The rows are named by file, the form the `step02_1` and `step02_6` suites read; those
two read the reduction chain's hashes, which this step does not touch, so neither needed rerunning.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_closure_profile.nc` | profile | `fcc2e2c4ae71554537d8aa53f07b247724b718a371abab73de31939bae9ae6af` |

It carries `casspian_git_commit = 780a2de0d91b90e06f6906b504a65bc2308261e0`, not `-dirty`. It holds
24 variables where the Step 0 sweep's held 22: `u_column_ms`, 66 levels all at
2.16708095051 m/s, which is the closure wind, and `u_column_uncertainty_ms`, all NaN as SPEC_03
decision 5 requires. The Step 0 sweep's row for this file,
`543129c4a17e4155858bee3b4f6d66eb6332c3290d220e08c0a7c1087c2f70f0`, is superseded, and so is
REPORT_03_step4's earlier `e6693173...`. The forward runs' inputs and products are not committed
(`.gitignore`); they are regenerable from their build control files.

## 7. Next step

Step 2: `lib.mesh`, the columns, `z_lv`, and by the review's order of work item 5 the
cylinder-extended wind built by decision P as the first deliverable beyond the mesh, with check 9's
cylinder half and the first filing's check 14 exercised there against the values of section 12
ruling 1: on the anchor's column 9.490 m/s at the top, 3.717 at the gauge, 2.168 nearest 1 bar and
1.926 at the bottom, `Phi` lower at the top by 322.3 m2/s2 and higher at the bottom by 18.0, and at
10 N 427.09, 345.98 and 326.40 m/s at the top, gauge and bottom against 329.46 on the reference
level. The extension is a fixed point with the reference surface, so `through_anchor` is called
inside that loop and deliverable 1 is what it converges on.

One thing else carries forward, for the Step 3 review, as the Step 0 review asked. Deliverable 2 settles which
one-sided slope a mesh node that coincides with a wind-file node takes: the interval to the north,
or to higher pressure. The rule is implemented and check 12 measures it in force. What it costs on
the one state in hand, measured at the two nodes that bracket the anchor: the anchor sits in the
30.5 to 31.0 degree interval, whose slope is -691.638087 m/s per radian, -12.071 m/s per degree;
the interval immediately south of the 30.5 node is steeper, -889.062216 per radian, -15.517 per
degree. A mesh node landing exactly on 30.5 takes the northern, gentler of those two by this rule.
The two differ by 29 percent, so on a mesh whose nodes can coincide with the file's the choice is
visible in the kernel, and Step 3 is where that is bounded.
