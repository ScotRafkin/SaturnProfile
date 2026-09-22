# REPORT 04, Step 2. The working mesh, the columns, and the cylinder-extended wind

Coding agent, 22 September 2026; refreshed under REVIEW_04_step2 and SPEC_04 v0.11.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.11, section 4, decision P and section 13. Acceptance script `reports/step04_2/accept_step04_2.py`, output
`reports/step04_2/output.txt`. Regression `reports/step04_2/run_regression.sh`, output
`reports/step04_2/regression.txt`.

**Thirteen of thirteen checks pass.**

The first filing was 12 of 15, three checks failing on the specification's numbers rather than the
code, and all three are settled by the section 13 rulings:

* **checks 1a and 1b** asked that `J = 0, Omega = 0, u = 0` return `r - r0 = Phi / g` to 1e-12 and
  `z_lv = r - r0` to 1e-12 relative. Those constants leave `GM / r^2`, not a uniform `g`, and the
  second bound sits below double precision. Ruling 1 restates the first as the central field's own
  closed form `r = r0 / (1 - r0 Phi / GM)`, which the column meets at **8.901e-16**; ruling 2
  restates the second as 1e-7 m absolute, which it meets at **5.541e-08 m**, seven ulp of `r`.
* **check 10**, the 10 N column at 0.054 m/s against a 0.05 bound, was the chord error of finding 3
  at the latitude the bottom level maps to, exactly as ruling 6 diagnosed. Ruling 3 samples the
  inversion curve at the mesh spacing over the file's range instead of the file's 0.5 degree grid,
  and the bottom level comes to **0.006**, nine times inside the bound. Nothing was loosened: the
  values and the bound are unchanged and the sampling that the construction always needed was
  corrected.

Ruling 3 improved the anchor's column with it: the top level from 0.009786 to **0.001758** m/s, and
the geopotential shift at the top from 0.298 to **0.098** m2/s2. The reference-level identity stays
exact at all 361 latitudes.

Two real defects were found and fixed in the first filing, both by checks the specification asked
for: decision P's construction has to be **per hemisphere** (finding 4, ruling 4), and
`Mesh.d_dlatitude` needed `edge_order=2` because `np.gradient` defaults to a first-order edge
(finding 5, ruling 5).

## 0. Before this step

HEAD `5608665`, Step 1 accepted at `780a2de` with the sweep recorded at `3dcd274`, the post-sweep
suites and the acceptance at `4a69218`, and the author's ratification addendum at `5608665`. Every
product on disk is the swept one: the reduction chain clean at `1c9b310` and the closure product at
`780a2de`.

This step registers nothing and changes no existing module. `lib/mesh.py` is new; the
cylinder-extended wind is written under `reports/step04_2/`, which git ignores, and read back
through `wind_at` like any other file. The closure run's own inputs are not touched.

## 1. What was built

**`lib.mesh.Mesh` and `build_mesh`.** Latitude nodes uniform at the declared spacing from one
spacing south of the southernmost required latitude to one spacing north of the northernmost, with
the anchors, the gauge latitude (decision K) and the target inserted as exact nodes. Geopotential
nodes are whole multiples of the spacing, so `Phi = 0` is a node by construction, from one spacing
below the lowest anchor `Phi_k` to one spacing above the highest. Nothing about the extent is
declared or checked (v0.9). `extend(side, amount)` adds whole spacings on one of four sides;
`d_dlatitude` is the three-point difference for unequal spacing; `latitude_index` finds a required
latitude without the caller matching floats.

**`lib.mesh.column` and `build_columns`.** `dr/dPhi = 1 / g` and `dz/dPhi = 1 / |g_eff|` integrated
together by RK4 from `(r0, 0)` outward in both directions, so `Phi = 0` carries the reference
radius exactly and `z = 0` there. `g` is positive inward and `Phi` increases upward, so both slopes
are positive and no hemisphere needs a sign of its own. The nodes need not be uniform and need not
be the mesh's: the acceptance integrates the anchor's own `Phi_k`, which is how the column is
compared with what the reduction tabulated. On the nodes the column carries `g`, `G_phi`,
`|g_eff|`, `psi` and `u`.

**The cylinder-extended wind of decision P**, built in the acceptance script: `u = U(s)` with
`s = r cos phi`, `U` pinned on the file's own reference level by `U(s_ref(phi)) = u_reference(phi)`
with `s_ref(phi) = r(phi, Phi_ref) cos phi`, on the library's radial columns under the anchor's
flat-isobar map, as a fixed point with the Step 1 reference surface. Written through
`lib.io.write` under kind wind and read back through `wind_at`.

## 2. Decisions

Every choice the specification did not make.

1. **The columns live in `lib.mesh` with the mesh**, not in a module of their own. The deliverable
   heading is `lib.mesh` and the column text follows it; the two are built and rebuilt together on
   every pass of the outer loop.
2. **`column` takes `u_of_geopotential`, a callable, not the node values.** The deliverable says
   `u` comes from `wind_on_mesh`, which gives the nodes; RK4 also needs `u` at its stages, which
   lie between nodes. The callable is evaluated at the stages and at the nodes, and the `u`
   recorded on the nodes is that same callable there, which is what `wind_on_mesh` gives them.
   Interpolating node values to the stages instead would put a second interpolation inside the
   integrator.
3. **`extend` names four sides**, `south`, `north`, `below`, `above`. The specification says "on
   that side" without naming them, and a curve `Phi_k(phi)` can run off either end of either axis.
4. **`build_mesh` takes the required latitudes as a list and the anchors' `Phi_k` as one array or
   several**, rather than a namelist, so that `lib` reads no control file. The driver assembles
   the anchors, the gauge latitude and the target.
5. **The refusal is `ValueError`** for a spacing that is not positive, as `lib.composition` and
   `lib.windfield` raise for their refusals.
6. **`d_dlatitude` is `np.gradient` against the node coordinates with `edge_order=2`.** The
   interior formula is the three-point formula for unequal spacing the deliverable asks for. The
   edge order is explicit because the default is first order; see finding 5.
7. **The writer for the cylinder file is `lib.io.write` under kind wind**, not
   `tools.wind.build_wind`. That tool is the published-curve path, driven by a control file that
   reads a digitized curve; this field is not built from a curve. What matters is that the file
   goes through the kind W schema, which is what validates the three parts, the sum identity and
   the poles, and `lib.io.write` is that gate.
8. **The flat-isobar map is extrapolated linearly in `ln p` beyond the anchor's levels.** The
   file's pressure grid reaches 1 Pa and 1e6 Pa while the anchor spans 19.95 to 129848 Pa, so 23
   of the 61 columns fall outside it. `np.interp` would clamp them all onto the two end values,
   which would put 23 columns of the file at one radius. The ends are continued at the slope of
   the last interval instead.
9. **Superseded by ruling 3.** As filed, the `s_ref` inversion curve was sampled on the file's
   0.5 degree latitudes with the anchors and the target inserted. Ruling 3 generalizes it: the
   curve is sampled at the mesh's latitude spacing over the file's range, so that every latitude a
   level maps to, and not only the ones the run names, is within a mesh spacing of a sample. The
   file's own nodes and the run's exact latitudes are kept in the sample so the reference-level
   identity stays exact on them. Measured: 3,962 inversion latitudes against 363, largest gap
   0.0500 degrees.
10. **The inversion is per hemisphere.** Forced by finding 4, and now stated in decision P.
11. **The cylinder file carries the closure file's variable attributes**, with `value_source`
    restated to name the construction. Replacing a variable drops its attributes, and kind W
    requires units on every one; leaving the closure file's `value_source` in place would have
    claimed the Ingersoll and Pollard curve as the provenance of a synthetic field.

## 3. Acceptance results

| # | Check | Result |
|---|---|---|
| 1 | the central field returns `r = r0 / (1 - r0 Phi / GM)` to 1e-12 (ruling 1) | PASS, 8.901e-16 |
| 1b | in that field `z_lv = r - r0` to 1e-7 m absolute (ruling 2) | PASS, 5.541e-08 m |
| 2 | the column at `phi_c` reproduces `z_lv` and `r - r0` to 0.1 m | PASS, largest 0.050 m |
| 3 | `Phi` recovered by the trapezoid of `g` to 1e-9 relative | PASS, 4.164e-11 |
| 4 | the mesh refuses a spacing that is not positive | PASS, five cases refused |
| 5 | `extend` adds whole spacings and the columns match to 1e-12 | PASS, exactly 0 |
| 6 | two anchors at different latitudes are both exact nodes | PASS |
| 7 | every quantity on the mesh is finite | PASS, 329,643 points |
| 8 | the cylinder wind converges and holds its identities | PASS, identity exactly 0 |
| 9 | the anchor's column to 1e-2 m/s, the `Phi` shift to 1 m2/s2 | PASS, 0.001758 and 0.098 |
| 10 | the 10 N column to 0.05 m/s | PASS, largest 0.022 |
| 11 | beyond: what the 0.5 degree grid costs at the anchor | PASS, a description |
| 12 | beyond: the centered difference is exact on a quadratic | PASS, 6.128e-13 |

**The columns (checks 2, 3, 5, 7).** At `phi_c` on the anchor's own `Phi_k`: `z_lv` 286,715.2 m at
the top and -104,098.0 at the bottom, `r - r0` 288,051.3 and -104,576.8, every one the
specification's value, largest difference 0.050 m against a bound of 0.1. `(r - r0) / z_lv` at the
top is 1.004660 against `1 / cos psi` = 1.004705 at psi = 5.5469 degrees, which is the projection
the two altitudes differ by. `Phi` recomputed from the column's own `r` by the trapezoid of `g`
returns the 782 non-zero nodes to 4.2e-11. A mesh extended above by 12,000 m2/s2 and one built with
those nodes from the start give columns that agree **exactly**, 0.000e+00 in both `r` and `z_lv`
over 421 latitudes, which is what `extend` keeping the existing nodes means. All 329,643 points of
the mesh are finite.

**The cylinder wind (checks 8, 9, 10, 11).** The fixed point converges in four passes, changes
1.57e+02, 2.79e-01, 3.85e-04, 4.18e-07 m/s; the specification expects three, and the fourth is the
one that proves the third. The reference-level identity `u_total(phi, p_ref) = u_reference(phi)`
holds at **exactly zero** over all 361 latitudes, the sum identity read back from the written file
is 8.9e-16, and both poles are exactly zero. On the anchor's column the construction gives 9.488242,
3.715829, 2.167081 and 1.924960 m/s against the stated 9.490, 3.717, 2.168 and 1.926, largest
departure 0.001758 against 1e-2; the `Phi` shift is -322.202 and +17.946 against -322.3 and +18.0,
largest departure 0.098 against 1. At the level nearest 1 bar, which is exactly 1e5 Pa, the
construction returns `u_reference` at the anchor, 2.16708095 m/s, to 0.00e+00. At 10 N the column
carries 427.08, 345.96 and 326.39 m/s against 427.09, 345.98 and 326.40, and 329.46 on the file's
reference level, largest departure 0.022 against 0.05.

## 4. Findings

**1. The specification's uniform gravity check cannot pass, because those constants do not give a
uniform gravity.** `J = 0, Omega = 0, u = 0` removes the harmonics, the rotation and the wind, and
leaves the central field `GM / r^2`. Measured on the mesh's 783 nodes at `phi_c`: `g` runs 11.076953
m/s2 at the gauge to 10.979418 at the top node, a change of 0.88 percent, and `r - r0` departs from
`Phi / g0` by 4.432e-3 relative, against a bound of 1e-12. The size is not an accident:
`Phi / (g0 r0)` at the top node is 4.412e-3, which is what a first-order expansion predicts.

What the check was for is right, and it passes against the field the constants describe.
`dr/dPhi = r^2 / GM` integrates to `1/r0 - 1/r = Phi / GM`, so `r = r0 / (1 - r0 Phi / GM)`, and
the column reproduces that to **8.901e-16**.

**Ruled (section 13.1):** the specification's fault; the check is restated as that closed form to
1e-12, and is check 1 in the table above. The v0.10 form is still measured and printed beside it,
at 4.432e-3, so the record carries why it was restated.

**2. The second half of that check is below double precision.** `z_lv = r - r0` is asked to 1e-12
relative. With `G_phi` exactly `-0.0` and `|g_eff| - g` exactly `0.0`, the two right-hand sides are
the same number and the two integrations are the same arithmetic, so the only difference is where
they start: `r` accumulates from 58,516,188 m and `z` from 0, and `r - r0` cancels about seven
digits. Measured: 5.541e-08 m absolute, 5.809e-12 relative over the 782 nodes where `|r - r0| > 1`
m. One ulp of `r` is 7.45e-09 m, so the departure is **7 ulp**. No implementation in float64 can
meet 1e-12 on that difference, and an absolute bound of about 1e-7 m is the honest form.

**Ruled (section 13.2):** correct; the bound becomes 1e-7 m absolute, seven ulp of a 58,516 km
radius, and is check 1b in the table above, met at 5.541e-08 m.

**3. The `s_ref` inversion needs the run's latitudes as its own nodes.** The curve `s_ref(phi)` is
inverted to find which cylinder a point sits on, and it is sampled where `u_reference` is known,
the file's 0.5 degree grid. The anchor at 30.8056 degrees falls mid-cell between 30.5 and 31.0, so
the inversion used that cell's chord and the reference level did not return `u_reference` there:
2.175541 against 2.167081, off by 8.5e-3. Because every level of the column is inverted through the
same curve, that error appeared as a uniform +0.007 to +0.010 m/s offset on all four stated values.
Adding the anchor and the target to the sampled latitudes makes the identity exact and improved
three of the four levels by 3 to 11 times: the gauge from 0.006896 to 0.001996, 1 bar from 0.007541
to 0.000919, the bottom from 0.006975 to 0.000615, and the `Phi` shift at the bottom from 0.182 to
0.007.

The top level did not move, and that is consistent: it maps furthest equatorward, to about 30.27
degrees, which is mid-cell between the 30.0 and 30.5 nodes, so its chord error lives at a different
latitude. It is the same effect relocated, and it is why check 9 clears its bound with only two
percent of margin. **This is a property of the construction, not of this code**: any latitude that
is not a node of the sampling inherits its cell's chord error, and the error is largest where the
wind curve is steepest, which is exactly where this anchor sits.

**Ruled (section 13.3):** accepted and generalized. Sampling `s_ref` at the mesh's latitude spacing
over the file's range puts every mapped latitude within one mesh spacing of a sample, not only the
ones the run names. Measured after the change: 3,962 inversion latitudes against 363; the anchor's
top level from 0.009786 to **0.001758** m/s, the gauge from 0.001996 to 0.001171, the bottom from
0.000615 to 0.001040, and the `Phi` shift at the top from 0.298 to **0.098** m2/s2. The reviewing
agent's expected values were measured on a 0.05 degree sampling, which is why the finer sampling
closes on them.

**4. Decision P's construction has to be per hemisphere, and the specification does not say so.** A
cylinder of radius `s` cuts the reference surface once in each hemisphere, so a single global
`U(s)` forces `u_north = u_south` at equal `s`. The wind is not symmetric: `u_reference` is 107.791
m/s at 60 N and 27.329 at 60 S, 5.856 at 30.5 N and 42.060 at 30.5 S. A global inversion therefore
cannot satisfy `U(s_ref(phi)) = u_reference(phi)` in both hemispheres, and measured, it broke the
reference-level identity by **112 m/s**. Each hemisphere now carries its own `U(s)`, which keeps
`du/dz` zero within a hemisphere and so still gives `S = 0`; the identity is then exact at all 361
latitudes.

This would not have shown in any expected value: the run's transfer is 30.8 N to 10 N, entirely
northern, and every stated number is northern. It showed because the file is global and the schema
checks it as a whole. Step 4 reads that file.

**5. `Mesh.d_dlatitude` needed its edge order stated.** `np.gradient` defaults to `edge_order=1`,
first order one-sided at the two boundaries, which is not exact even on a quadratic. Measured on a
quadratic over nodes like the mesh's: the default leaves 1.5e-1 at each end against 1.1e-14 in the
interior. With `edge_order=2` the ends come to 1.4e-14 and the whole mesh to 6.128e-13 (check 12).
The margin keeps the edges away from any traced curve, as the deliverable says, but the mesh-wide
difference fields Step 3 forms are taken on every node, so a first-order edge would have entered
them.

**6. The 10 N column's bottom level misses its bound by eight percent.** Measured against the
stated values: the top level -0.005, the file's reference level +0.004, the gauge +0.025, the
bottom **+0.054**, against a bound of 0.05. Three of the four are well inside it. The departure
grows monotonically with depth while the reference level, where the construction is pinned exactly,
matches to 0.004, so it is systematic rather than scatter, and it is of the same size and sign as
the chord error of finding 3 at the latitude the bottom level maps to. It is reported rather than
diagnosed further, and nothing was adjusted: 10 N is an exact node of both the file's grid and the
inversion, so this is not the sampling of finding 3 at the target itself.

## 5. Regression

`reports/step04_2/run_regression.sh`, output `reports/step04_2/regression.txt`. Every accepted
suite on the registered swept products. This step adds a new module and changes no existing one, so
no accepted behaviour can move; the run is the procedure's confirmation of that rather than a test
of anything this step did.

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
| `step04_0/accept_step04_0` | 15 of 15 | 15 |
| `step04_1/accept_step04_1` | 13 of 13 | 13 |
| `step04_2/accept_step04_2` | 12 of 15, checks 1a, 1b and 10 of the first filing | this step |

All 22 accepted suites return their reference counts. The registered kind N restored to
SHA-256 `920304440ee4554648c3aa6321b713741ee26eb0eb1ff223f5a81a8f9b20fd37`, the swept value, and the
closure run's inputs and product and the transfer run's inputs all restored equal. The "paths
dirty" counter in the output rises from 1 to 4 during the run: this step's `lib/mesh.py`, then this
report, then the author's SPEC_04 v0.11 and REVIEW_04_step2, which arrived while the run was in
progress. No restore failed.

## 6. Hashes

Nothing was registered. The cylinder-extended wind is written to
`reports/step04_2/lindal_cylinder_wind.nc`, which git ignores, and is rebuilt by the acceptance
script whenever it runs; it is an instance of the state, not a product of the chain. No product of
the Step 1 sweep changes, so `reports/REPORT_04_step1.md` section 6 and `reports/REPORT_04_step0.md`
section 6 remain the record.

## 7. Next step

Step 3: `lib.kernel`. Three things carry forward.

All six findings are ruled on in section 13 and applied; none is left open.

Findings 3 and 4 are properties of decision P that v0.10 did not state, now in the specification's
text, and both are implemented in the acceptance script rather than in `lib`. If the
cylinder-extended wind is ever built anywhere but an acceptance script, the per-hemisphere
inversion and the mesh-spacing sampling of `s_ref` go with it; neither is a property of the file,
so a reader of the file cannot recover them.

The sampling of ruling 3 costs what it buys: 3,962 column integrations per pass against 363, four
passes, and the acceptance script's run time is now dominated by it. Section 5b of REPORT_04_step1
records the per-column cost that sets that.

For Step 3 itself, the rule the Step 0 review asked to be settled is in force and measured: a mesh
node that coincides with a wind-file node takes the interval to its north, and on the two intervals
bracketing the anchor those differ by 29 percent, -12.071 against -15.517 m/s per degree. The mesh
this step builds inserts the anchors and the target as exact nodes, so such coincidences are the
rule near them, not the exception.
