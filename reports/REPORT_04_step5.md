# REPORT 04, Step 5. The production at the target, altitude and datum, the product, transfer mode

Coding agent, 28 September 2026, refreshed under the review and again under the author's viewing.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.22, section 7, decisions B, C, F, G, H, P, Q and
R, and the section 15, 16, 17 and 18 rulings.

The interim filing `reports/REPORT_04_step5_interim.md` asked for a ruling on five decisions before
the acceptance script was written, and section 17 accepted all five. The first filing of this report
carried eleven checks with one failing; `reports/REVIEW_04_step5.md` accepted it with changes and
section 18 ruled. **This refresh is what the review's order of work asks for**: check 4 restated
against the v0.18 values and bound, check 2's convergence condition moved to the tracing's part, F8
and F9 rendered and covered by check 10, the regression rerun because the figure module changed, and
the head, the acceptance table and findings 3 and 4 rewritten.

**Section 18 ruling 5 then asked for one more change on the author's viewing**: F8's altitude panel
restated on a common datum, each profile above its own 1 bar level, which v0.21 carries. That is a
figure change, and under the section 0 rule at v0.22 it reruns the figure checks and no suite that
computes. The acceptance script gained `--carry` for it, so `output.txt` is one file at the final
code state with the six computing checks carried from the run of the whole suite, each row marked as
carried and naming where it was measured.

Acceptance script `reports/step04_5/accept_step04_5.py`, output `reports/step04_5/output.txt`,
regression `reports/step04_5/regression_figures.txt` with the full run of the first filing at
`reports/step04_5/regression_v017.txt`. The suite is pinned to 0.05 degrees by 50,000 m2/s2 and
asserts that the production namelist still carries them, because run 1 is the namelist's own run
through the driver (section 0 at v0.17). Four runs are kept: `run_v017.txt`, the first filing;
`output_full_v020.txt` with `run.txt`, the whole suite rerun under the review, eleven of eleven; and
`output.txt` with `run_carry.txt`, the carry pass that reran check 10 against the restated F8.

## 1. What is built

The four deliverables are as the interim filing describes them, and section 1 of that report stands.
Two things changed while the acceptance script was written, both in `forward/transfer.py`:

`run` is now three parts. `chain` is the outer loop, every anchor carried to the gauge, the estimate
there, `C` carried to the target and the production at it; `write_product` assembles and writes kind
`profile`; `run` reads the namelist, calls the two and renders the figures. Nothing moved between
them, and the reason is measurement: four of the five runs cannot go through a namelist, and without
the split the suite would measure a hand copy of the driver's chain rather than the driver. Check 2
reports the chain at the pinned spacings giving a production equal to the driver's by `array_equal`.

A copied subtree's own vertical dimension is no longer called `level`. The delivered profile's
`level` is the union of the anchors' arrival levels, 131 at M = 2, and an `xarray` DataTree child
cannot carry a dimension of its parent's name with another size. So the kind N copy under
`anchors/<slug>` carries `anchor_level`, a copied input under `inputs/<kind>` carries `input_level`,
and `isobars` carries each anchor's curves on `level_<slug>`, which also allows two anchors with
different level counts. The variables, their values and their attributes are the copied file's
unchanged. This is finding 2 below.

**Under the review, two figures and one product variable.** F8 draws the transfer along the
isobars: `S/g` against latitude on five isobars between the gauge and the target, the temperature
change accumulated along the same five, and the delivered `T` against altitude above the datum and
against radius, each beside the anchor's own on its own heights and radii. F9 draws the wind the run
assumed from the product's `inputs/wind` group, `u_reference(phi)` over the file's whole range with
the anchors, the gauge and the target marked and the mesh's span shaded, and `u_total(phi, p)` as
filled contours over the mesh's latitudes and the profile's pressures, the file's `source` and
`vertical_structure` in the titles. Both are in `tools/plots/figures_profile.py` and in
`FIGURES_TRANSFER`, so the driver and `casspian-plots` render them alike.

F8 required one thing the product did not carry. Deliverable 4 says the kernel is the `isobars`
group's and that nothing is recomputed for the figure; deliverable 3's list of that group named only
`Phi_k(phi_i)` and `ln N_k(phi_i)`. So `isobars` now carries `shear_kernel_<slug>_per_rad`, `S/g` at
the nodes of each anchor's own curves, sampled exactly as `transfer` samples it there (Eq. A27,
bilinearly on the mesh). The composition term is not added, because F8 draws `S/g`. This is
decision 5 below.

## 2. Decisions taken since the interim filing

Five decisions are ruled in section 17 and the acceptance script is written against them. Four more
were taken at the first filing and section 18 accepts all four; two more were taken under the
review, and section 18 ruling 5 accepts both, with the specification carrying them at v0.21.

1. **`chain` and `write_product` are extracted from `run`** so that runs 2 to 5 and run 5's product
   go through the driver's own code and not a copy of it. The alternative was 25 lines of the
   driver's chain restated in the acceptance script, where a divergence between the two would be
   invisible. No behavior moved: check 2 measures the chain against the driver bit for bit.

2. **A copied subtree's vertical dimension is renamed**, to `anchor_level`, `input_level` or
   `level_<slug>`, in transfer mode and at M = 1 as well, so that the product's structure does not
   change with the number of anchors. Closure mode is untouched, because there the anchor's levels
   are the product's. The alternative, renaming only when the sizes differ, would make the shape of
   a product depend on M, which a reader cannot write code against.

3. **The suite asserts the namelist's spacings rather than overriding them.** Section 0's rule at
   v0.17 has an acceptance script carry the spacings its bounds were set for. Step 5's bounds are
   stated at decision R's 0.05 degrees by 50,000 m2/s2, which is what the production namelist now
   carries, and run 1 is the namelist's own run through `casspian-forward`. So the pin is an
   assertion that the two still agree, naming both, and runs 2 to 5 take the pinned values. A
   namelist change would stop this suite rather than quietly move it.

4. **The sheared run's pair of spacings is the pinned one and half of it**, where the Step 4 check
   of the same file used the namelist's and twice it. Twice decision R's 50,000 m2/s2 puts the
   mesh's lower edge deep enough that the isobar map, continued at its last slope, asks the wind
   above 1 MPa, and `lib.windfield` refuses the extrapolation. That refusal is correct and has
   nothing to do with the shear, so the pair moved to the finer side. It is recorded here because it
   is a consequence of decision R that the specification did not state: at 50,000 m2/s2 the mesh's
   lower edge is close to the wind file's deepest pressure, and a run at twice that spacing is
   outside it.

5. **The `isobars` group carries `shear_kernel_<slug>_per_rad`.** F8 draws `S/g` along the isobars
   and deliverable 4 says the kernel is the `isobars` group's with nothing recomputed for the
   figure, but deliverable 3's list of that group named only `Phi_k(phi_i)` and `ln N_k(phi_i)`. The
   figure could have recomputed the kernel from the wind and the mesh, which deliverable 4 forbids,
   or the product could carry what the transfer already formed, which is this. It is `S/g` at the
   nodes of each anchor's curves, by the same call `transfer` makes, on the anchor's own vertical
   dimension. No schema entry was needed: the `isobars` group's variables are named for their
   anchor and the companion rule is the root's. The rest of F8 is read off the group's `ln N_k`,
   which is the integral of that kernel, and the run's own composition, so the second panel is the
   first panel integrated and nothing else. **Accepted at section 18 ruling 5 and now deliverable
   3's text (v0.21).**

6. **F8's altitude panel puts both profiles above their own datum level**, the anchor's located on
   its kind N heights by the log-linear rule the datum itself is placed with. The first drawing
   followed the v0.18 wording, the delivered profile above the run's datum and the anchor above its
   own anchor isobar, and the author read the 90 km between those two reference levels as a shift in
   the delivered profile. The panel is a comparison of two columns' thickness, so it now measures
   both from the same isobar; the radius panel keeps the offset, because there it is the reference
   surface and is the point. Ruled at section 18 ruling 5 and carried by v0.21.

   The reviewing agent's nit on that panel is fixed in the same pass. The legend named 1,602 km,
   which is the gap at the bottom level, while the words called it the reference surface, whose gap
   is 1,612 km at the gauge isobar: there the delivered radius is `r0` at the target and the
   anchor's is its own anchor isobar radius, so the difference is the surface and nothing else,
   while at the bottom it is the surface plus the two columns' thicknesses. The legend now names the
   gauge value and says where it is measured.

**Eleven checks, eleven passing.** Every run is at the pinned 0.05 degrees by 50,000 m2/s2, and
run 1 is the namelist's own run through `casspian-forward`. The whole suite costs 3,802 s, of which
check 4's cylinder construction is 2,023 s; the carry pass that reran check 10 against the restated
F8 costs 417 s. Checks 1, 3, 5, 6, 8, 9 and 11 exercise nothing that changed under the review and
reproduce the first filing to the last digit; checks 2, 4 and 10 are the ones the review asked to be
rerun. Of those, 2 and 4 are carried into `output.txt` from the run of the whole suite, because the
change that followed them reaches only the figures.

| Check | What is measured | Measured | Stated, bound | |
|---|---|---|---|---|
| 1 | run 1, T below the closure product's at the top, gauge and bottom | 1.1022e-02, 1.0884e-02, 1.0833e-02 | 1.103e-2, 1.089e-2, 1.084e-2; largest departure 7.624e-06 of 1e-4 | pass |
| 1 | run 1, `r0` at the target | 60,128,613.0 m | 60,128,613.0, 0.034 m from it of 1 m | pass |
| 1 | run 1, altitude above the 1 bar datum | 411,135, 98,187, -15,290 m | 411,132, 98,186, -15,290; largest 2.88 m of 5 m | pass |
| 2 | the delivered pressure identity | 5.777e-07 | 1e-6 (v0.18 ruling 1) | pass |
| 2 | its tracing part, the composition held on the tabulated axis, at the pinned spacings and halved | 1.293e-07 and 1.200e-07 | must fall; ratio 1.08 | pass |
| 2 | the delivered identity halved, for the record | 5.770e-07, ratio 1.00 | the floor, not a convergence | reported |
| 2 | the chain against the driver's production | equal by `array_equal` | beyond the list | pass |
| 3 | run 2 at 60 N, T below the closure product's | 2.5089e-03, 2.4818e-03, 2.4713e-03 | 2.51e-3, 2.48e-3, 2.47e-3; largest 1.844e-06 of 1e-4 | pass |
| 3 | run 2, altitudes | 328,127, 78,533, -12,239 m | 328,127, 78,533, -12,239; largest 0.34 m of 5 m | pass |
| 4 | run 3, the cylinder wind, T, p and N against the closure product | 9.684e-04, 1.901e-04, 9.525e-04 | 2e-3, decision Q's floor | pass |
| 4 | run 3, altitudes | 416,354, 99,309, -15,450 m | 416,320, 99,295, -15,457 (v0.18); largest 33.79 m of 60 m | pass |
| 4 | run 3, `r0` at the target under the cylinder wind | 60,129,732.3 m | the reviewing agent's fixed point 60,129,726 m, 6.3 m from it | reported |
| 4 | run 3, the isobar shift that sets the bound | +197.9, +0.0, -67.9 m2/s2 at the top, gauge and bottom; largest 235.5 | decision Q's 300 m2/s2 | reported |
| 5 | run 4 at the anchor's latitude, N and Phi against the closure product | 1.776e-15 in `ln N`, 0.000e+00 m2/s2 in Phi | 1e-12 | pass |
| 5 | run 4, p and T against the closure product | 5.686e-07, 7.441e-07 relative | 1e-5 | pass |
| 5 | run 4, altitudes | 376,783, 90,067, -14,031 m | 376,780, 90,067, -14,031; largest 2.56 m of 5 m | pass |
| 6 | run 5, M = 2, against run 1's product on the 129 union levels inside its range | 3.242e-08 in `ln N`, 7.583e-08 in p, 4.376e-08 in T, 0.00 m in altitude | 1e-5 | pass |
| 6 | run 5, the anchors groups and `D_ij` | 2 groups, `D_lindal_synthetic60` largest 6.484e-08 | one group per anchor; the floor is 6.5e-08 | pass |
| 7 | run 1's product validates and reads back with deliverable 3's groups | 24 groups, all present, the tabulated columns absent, one record group, the residual `derived` | deliverable 3 | pass |
| 8 | the datum refused above and below the produced range, both records, neither | four refusals, each naming the fault | deliverable 2, v0.17 ruling 3 | pass |
| 8 | a datum that is exactly a produced pressure | lands on that level exactly | beyond the list | pass |
| 9 | the sheared run's pressure identity at the pinned spacings and halved | 3.342e-04 and 1.454e-04, ratio 2.30, 6 passes each | must fall and converge | pass |
| 10 | F5, F6, F7, F8 and F9 by the driver and by `casspian-plots` | identical above the 19 px footer band, all five present, none skipped | identical apart from the footer (v0.19) | pass |
| 11 | the closure production against the registered product | equal by `array_equal` on all nine columns | beyond the list; the step's safety property | pass |

The pressure identity is reported at every level of every run, as decision H requires: check 1
prints all 66 levels, and every other check prints the statistics and where the largest sits. The
largest values are 5.777e-07 at 10 N, 5.713e-07 at 60 N, 8.834e-05 under the cylinder wind,
5.686e-07 at the anchor's own latitude, 5.763e-07 at M = 2 and 3.342e-04 under the shear.

**The author views F5, F7, F8 and F9 by eye and accepts them**, which section 7's acceptance asks
for at v0.19. They are in `forward/lindal_transfer/output/figures/`, with F6 beside them and
`lindal_transfer_diag.pdf` carrying all five. F7 at M = 1 has three panels; the M = 2 product at
`reports/step04_5/m2/output/` is the case with the fourth.

Two things the author asked about F8 on viewing it, with the answers section 18 ruling 5 gives, are
now in the figure's own docstring so that the next reader has them. The five kernel curves coincide
because this wind does not vary along a column, which is the file's own statement, its
`vertical_structure` being `altitude independent`; the line widths step down so that every member of
a coincident set still shows. And the staircase in `S/g` is the wind file, not the mesh: under
decision L the wind is linear between the file's half-degree nodes, so `S = 2 Omega cos(phi)
du/dphi` is constant on each interval and steps at every node, a mesh at 0.05 degrees resolves those
steps exactly, and a finer one would draw the same staircase. The alternation on the flank is the
source curve's own node-to-node roughness and the kinks in the second panel are its integral, which
cancels, which is why the accumulated offset is smooth. Smoothing belongs to the wind tool, not to
the model, which reads what it is given.

The altitude panel now measures both profiles from their own 1 bar level (decision 6), so it
compares the thickness of two columns; the radius panel keeps the gap between them and names it at
the gauge isobar, 1,612 km, where it is the reference surface between the two latitudes.

## 4. Findings

Four. The first two are the coding agent's own, found while the acceptance script was written, and
section 18 accepts them as reported. The third and fourth were measurements of the specification's
own construction that the first filing put up for a ruling; section 18 ruled on both, and they are
rewritten here as what was measured and what was ruled. Under those rulings check 4 passes, so this
refresh carries no failing check.

1. **The acceptance script compared the closure product's columns on the wrong pressure axis.** A
   transfer product's `pressure_label_Pa` is not the anchor's tabulated pressure: it is the closure
   production's own pressure at that level, and the two differ by the closure residual, 3.3e-3 at
   the bottom row. Interpolating the closure product's temperature onto the tabulated axis put
   6.078e-04 into the comparison against a 1e-4 bound. The delivered levels are the anchor's levels
   carried to the target, so the comparison is level for level; the labels agree with the closure
   production's pressures to 2.220e-15, which check 1 prints, and that is the measurement that says
   the pairing is right. With it, the departures are 7.6e-06. The defect was in the script, not in
   the product, and the interim report's table carried a smaller version of it.

2. **An M = 2 product could be written and not read back.** Four places carried a `level` dimension
   of the anchors' level count while the root carried the union's: the kind N copy under
   `anchors/<slug>`, the copied inputs under `inputs/<kind>`, the anchor's own columns in
   `anchors/<slug>/transfer`, and each anchor's curves in `isobars`. `xarray` refuses to open such a
   tree, naming the group and the two sizes. This is a defect of the product and it could only
   appear at M >= 2, which the interim filing listed as the path with no test at all. The fix is
   decision 2 above. Nothing in the M = 1 product's values changed; its dimension names did.

3. **The pressure identity does not fall under halving, and what falls is a tenth of it. Ruled at
   section 18 ruling 1.** Measured 5.777e-07 at the pinned spacings and 5.770e-07 with both halved,
   a ratio of 1.00: the delivered identity has a floor near 5e-7 for this state and does not
   converge. Holding the composition read on the anchor's tabulated pressures, at the same
   latitude and with everything else the same, it is 1.293e-07 and 1.200e-07: that part is the
   tracing's own discretization and it falls.

   The mechanism, which the ruling accepts: a level's label is the closure production's own
   pressure, which that production formed with the composition read at the anchor's tabulated
   pressures, while the transfer production re-forms the pressure with the composition read at the
   labels, as section 15 ruling 4 requires at a target where the levels are the union's. The two
   reads differ by the closure residual in their argument, 3.3e-3 at the bottom row, and by
   1.313e-06 in `ln R_bar`. Excluding the bottom row does not remove it: 4.019e-07 and 4.113e-07.

   **The ruling takes the first of the two alternatives this report offered.** The convergence
   condition of section 7 now applies to the tracing's part and the delivered identity is reported
   at every level and bounded at 1e-6; the closure path is not changed, because reading the
   composition at the produced pressure would mean iterating the closure and would move a
   registered product by a part in a million for nothing. Check 2 is written that way and passes on
   both counts. Section 7 carried this mechanism for the fourth run and not for the first, which
   the ruling records as the specification's omission, not the code's.

4. **The cylinder-extended run's altitudes: the check failed on the specification's stated values
   and bound, both now restated, and the delivered altitudes were right. Ruled at section 18
   ruling 2.** At the first filing the measured 416,354, 99,309 and -15,450 m sat 53.79 m from the
   v0.10 values against a 10 m bound. Two things were wrong with the comparison, neither in the
   delivered numbers.

   **The stated values were measured with the wrong reference surface.** Decision P makes the
   surface and the cylinder file a fixed point, the surface marched under the file's own gauge
   wind; the v0.10 values were measured with it marched under the closure wind. At 10 N the fixed
   point surface is 1,119 m higher, and this run's `r0` there is 60,129,732.3 m against the
   reviewing agent's 60,129,725.8, 6.3 m apart. Under the exact hypothesis the reviewing agent's
   altitudes are 416,320, 99,295 and -15,457 m, and those are the stated values from v0.18.

   **The bound was decision Q's floor at the wind but not in altitude.** The gridded cylinder
   file's arrival isobars are shifted from the anchor's own by +197.9 m2/s2 at the top level, 0.0
   at the gauge, where the gauge isobar cannot move, and -67.9 at the bottom, largest 235.5 against
   decision Q's 300; the review measures -127 at the datum, which is the level the altitude is
   referred to. The altitude of a level above the datum moves by the difference of the two shifts
   over `g`, and the bound from v0.18 is 60 m. **Measured 33.79 m**, which is the 34 m the ruling
   states. Check 4 passes.

   **Two claims from the first filing are withdrawn.** The 3.494 m/s I inferred in `u` at the
   target column is not a wind error: the review's own column, given this run's arrival
   geopotential and datum, returns 416,353.9, 99,308.7 and -15,450.3 m against this run's
   416,353.8, 99,308.7 and -15,450.3 at the production spacing, and the same at 25,000 and 5,000,
   so the altitude path is verified to 0.1 m and what I read as a wind error was the shifted
   isobars read as one. Decision Q's file-grid truncation, which I named as the candidate, is the
   cause, but it acts through the isobars and not through `u`. What stands from that diagnosis is
   the 12 m that moves with the geopotential spacing when the wind varies along a column, which the
   ruling records for the performance step, and the two candidates it killed.

## 5. The regression

**The first filing ran the full set: twenty five suites, twenty five passing, every one at its
reference count**, in `reports/step04_5/regression_v017.txt`, with each suite's own output beside
it. That run is the one that matters for this step, because it covers every module Step 5 changes:
`forward/production.py`, `lib/schema.py`, `lib/io.py`, `forward/transfer.py` and
`tools/plots/figures_profile.py`. The rows that carry the most of it are `step03_4` 13 of 13, the
closure run and F5 and F6 through the refactored production and the changed figure module;
`step03_3` 16 of 16, the namelist and schema refusals; `step04_0` 15 of 15, which reads the transfer
namelist whose geopotential spacing this step moved; and `step04_4` 11 of 11, at the 5,000 m2/s2 it
is pinned to, which is the section 0 pin doing its work rather than being asserted. It ran from
00:22 to 07:03, and `step04_4` alone is four and a half hours of that.

**Since that run the changes are F8, F9 and F8's restated panel in `tools/plots/figures_profile.py`,
and one variable added to a transfer product's `isobars` group, so under the section 0 rule at v0.22
the regression is `step02_5` and `step03_4`**, which are the suites whose figures that module draws,
together with this step's own check 10; they give 7 of 7 and 13 of 13 in
`reports/step04_5/regression_figures.txt`, the registered kind N, the closure inputs and product and
the transfer inputs all restored, and the porcelain carrying only this step's own paths. That set
ran twice, once for the restated altitude panel and once for the legend correction on the radius
panel, and gave the same rows both times.

The rule asks for the reason the others were not rerun, in one sentence: none of `step1` to
`step9`, `step02_1` to `step02_4`, `step02_6`, `step03_1` to `step03_3`, `step04_0` to `step04_3` or
`step04_4` imports the changed figure module or reads a product made with it, and the one change
that is not a figure, `shear_kernel_<slug>_per_rad`, adds a key to the dictionary
`_transfer_groups` builds and touches no function that computes a delivered value, the only
accepted suite importing `forward.transfer` being `step04_4`, which writes no product and never
calls it. The evidence that nothing moved is inside this step's own acceptance rather than asserted:
after the change, checks 1, 3, 5, 6, 8, 9 and 11 reproduce the first filing digit for digit, and
check 11 compares the closure production with the registered product by `array_equal` on all nine
columns. If the reviewing agent reads that claim differently, `step04_4` is four and a half hours
and settles it.

**One thing about the first filing's regression is the coding agent's fault and is recorded rather
than tidied.** Its restore verification was made by a separate script,
`reports/step04_5/verify_restored.sh`, because `run_regression.sh` did not reach its own tail. While
bash was executing it the coding agent edited the file, replacing one word in a comment to correct a
spelling. Bash reads a script by byte offset, so removing a byte shifted every later offset and the
file ended, as bash saw it, inside a quoted string: the console carries `line 83: unexpected EOF
while looking for matching '"'`. The loop had already been read into memory and ran in full, which
is why all twenty five rows are present and why the restore after each suite happened, including
after the last; what was lost is the tail that compares the digests and prints the porcelain. That
comparison is between the copies the loop set aside before it started and the files in the tree
afterwards, so it is the same comparison made later rather than in the same process, and the script
says so. No suite's result depends on it. The lesson is narrow and worth writing down: do not edit a
script while a shell is running it.

## 6. Hashes

Nothing is registered and no registered product changed, so there is no sweep. The registered kind N
is `920304440ee4...` and the registered closure product `fcc2e2c4ae71...`, both as they were at the
Step 4 acceptance commit `4802983`; check 11 rebuilds the closure production in a run directory of
its own and finds it bit-identical, so the registered file was read and never rewritten.

The transfer product this step delivers is
`forward/lindal_transfer/output/lindal_transfer_profile.nc`, SHA-256
`3d73e13883be...`, written at the pinned spacings with the working tree carrying this step's
changes, so its `casspian_git_commit` ends `-dirty`. It differs from the first filing's
`620f2deef79c...` by one thing, the `shear_kernel_<slug>_per_rad` the `isobars` group now carries
for F8 (decision 5); no delivered value moves, and checks 1, 3, 5, 6, 8, 9 and 11 reproduce the
first filing to the last digit. A
forward run product is transient and regenerable from its run directory, which is why
`forward/*/output/` is ignored by git; it is recorded here so that the figures the author views can
be tied to a file.

The M = 2 product is `reports/step04_5/m2/output/lindal_transfer_profile.nc`, inside the gitignored
step directory, written from a synthetic anchor whose commit is `1c9b310018979bb...-dirty` with the
refusal relaxed in the acceptance process only.

## 7. The next step

Step 5 is the last step of SPEC_04. What follows it, in the order section 0 and decision R set out:

1. **This filing's dispositions.** Finding 3, the pressure identity's floor, asks whether the bound
   should carry the tracing's part with the composition floor reported beside it. Check 4's failure
   asks what to do about the cylinder-wind altitude, where the coding agent has narrowed the
   departure to 3.494 m/s in `u` at the target column and cannot settle the rest without an
   independent measurement.
2. **The performance step** of decision R, before the Monte Carlo wrapper. This step's costs are on
   the record for it: a transfer run at the pinned spacings is 70 to 100 s, a run with both spacings
   halved is about 230 s, the cylinder construction alone is 1,881 s over 3,962 latitudes, and the M
   = 2 run on a mesh reaching 60 N is 248 s. The outer loop rebuilds every column on every pass,
   which is where the time goes.
3. **The combination**, A35, the posterior wind and the reference-surface constant, then SPEC_05,
   the end-to-end tests.

`docs/specs/STATE.md` is set to `reported` for Step 5 and nothing is committed beyond the branch
`step04_5-wip`, which carries the work in progress and no acceptance meaning (section 17).
