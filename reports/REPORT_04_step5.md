# REPORT 04, Step 5. The production at the target, altitude and datum, the product, transfer mode

Coding agent, 27 September 2026. Specification: `docs/specs/SPEC_04_Transfer.md` v0.17, section 7,
decisions B, C, F, G, H and R, and the section 15, 16 and 17 rulings.

The interim filing `reports/REPORT_04_step5_interim.md` asked for a ruling on five decisions before
the acceptance script was written. Section 17 accepted all five. This is the step's filing: it
restates nothing the interim report already carries, and the interim report's section 3 table is
superseded by section 3 here, which is measured by the acceptance script.

Acceptance script `reports/step04_5/accept_step04_5.py`, output `reports/step04_5/output.txt`,
regression `reports/step04_5/regression.txt`. The suite is pinned to 0.05 degrees by 50,000 m2/s2
and asserts that the production namelist still carries them, because run 1 is the namelist's own run
through the driver (section 0 at v0.17).

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

## 2. Decisions taken since the interim filing

Five decisions are ruled in section 17 and the acceptance script is written against them. Four more
were taken here.

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

## 3. The acceptance table

Eleven checks, ten passing. Every run is at the pinned 0.05 degrees by 50,000 m2/s2, and run 1 is
the namelist's own run through `casspian-forward`. The whole suite costs 3,582 s, of which check 4's
cylinder construction is 1,881 s.

| Check | What is measured | Measured | Stated, bound | |
|---|---|---|---|---|
| 1 | run 1, T below the closure product's at the top, gauge and bottom | 1.1022e-02, 1.0884e-02, 1.0833e-02 | 1.103e-2, 1.089e-2, 1.084e-2; largest departure 7.624e-06 of 1e-4 | pass |
| 1 | run 1, `r0` at the target | 60,128,613.0 m | 60,128,613.0, 0.034 m from it of 1 m | pass |
| 1 | run 1, altitude above the 1 bar datum | 411,135, 98,187, -15,290 m | 411,132, 98,186, -15,290; largest 2.88 m of 5 m | pass |
| 2 | the pressure identity at the pinned spacings and halved | 5.777e-07 and 5.770e-07 | must fall; ratio 1.00 (finding 3) | pass |
| 2 | the chain against the driver's production | equal by `array_equal` | beyond the list | pass |
| 3 | run 2 at 60 N, T below the closure product's | 2.5089e-03, 2.4818e-03, 2.4713e-03 | 2.51e-3, 2.48e-3, 2.47e-3; largest 1.844e-06 of 1e-4 | pass |
| 3 | run 2, altitudes | 328,127, 78,533, -12,239 m | 328,127, 78,533, -12,239; largest 0.34 m of 5 m | pass |
| 4 | run 3, the cylinder wind, T, p and N against the closure product | 9.684e-04, 1.901e-04, 9.525e-04 | 2e-3, decision Q's floor | pass |
| 4 | run 3, altitudes | 416,354, 99,309, -15,450 m | 416,300, 99,290, -15,456; largest 53.79 m of 10 m | **fail** |
| 5 | run 4 at the anchor's latitude, N and Phi against the closure product | 1.776e-15 in `ln N`, 0.000e+00 m2/s2 in Phi | 1e-12 | pass |
| 5 | run 4, p and T against the closure product | 5.686e-07, 7.441e-07 relative | 1e-5 | pass |
| 5 | run 4, altitudes | 376,783, 90,067, -14,031 m | 376,780, 90,067, -14,031; largest 2.56 m of 5 m | pass |
| 6 | run 5, M = 2, against run 1's product on the 129 union levels inside its range | 3.242e-08 in `ln N`, 7.583e-08 in p, 4.376e-08 in T, 0.00 m in altitude | 1e-5 | pass |
| 6 | run 5, the anchors groups and `D_ij` | 2 groups, `D_lindal_synthetic60` largest 6.484e-08 | one group per anchor; the floor is 6.5e-08 | pass |
| 7 | run 1's product validates and reads back with deliverable 3's groups | 24 groups, all present, the tabulated columns absent, one record group, the residual `derived` | deliverable 3 | pass |
| 8 | the datum refused above and below the produced range, both records, neither | four refusals, each naming the fault | deliverable 2, v0.17 ruling 3 | pass |
| 8 | a datum that is exactly a produced pressure | lands on that level exactly | beyond the list | pass |
| 9 | the sheared run's pressure identity at the pinned spacings and halved | 3.342e-04 and 1.454e-04, ratio 2.30, 6 passes each | must fall and converge | pass |
| 10 | F5, F6 and F7 by the driver and by `casspian-plots` | identical above the 19 px footer band, all three present, none skipped | identical apart from the footer | pass |
| 11 | the closure production against the registered product | equal by `array_equal` on all nine columns | beyond the list; the step's safety property | pass |

The pressure identity is reported at every level of every run, as decision H requires: check 1
prints all 66 levels, and every other check prints the statistics and where the largest sits. The
largest values are 5.777e-07 at 10 N, 5.713e-07 at 60 N, 8.834e-05 under the cylinder wind,
5.686e-07 at the anchor's own latitude, 5.763e-07 at M = 2 and 3.342e-04 under the shear.

**The author views F5 and F7 by eye and accepts them**, which section 7's acceptance asks for. They
are `forward/lindal_transfer/output/figures/lindal_transfer_diag_F5_geopotential.png` and
`_F7_delivered.png`, with F6 beside them and `lindal_transfer_diag.pdf` carrying all three. F7 at M
= 1 has three panels; the M = 2 product at `reports/step04_5/m2/output/` is the case with the
fourth.

## 4. Findings

Four. The first two are the coding agent's own and were found in this stretch. The third and fourth
are measurements of the specification's own construction and each asks for a ruling; the fourth is
check 4's failure.

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

3. **The pressure identity does not fall under halving, and what falls is a tenth of it.** Section 7
   has the identity reported at every level and the report's value falling under halving. It is
   5.777e-07 at the pinned spacings and 5.770e-07 with both halved: a ratio of 1.00. The check's
   condition is met, by a tenth of a percent, and that is not convergence and should not be read as
   it. Holding the composition read on the anchor's tabulated pressures instead of on the produced
   labels, at the same latitude and with everything else the same, the identity is 1.293e-07 at the
   pinned spacings and 1.200e-07 halved: that part is the tracing's own discretization and it falls.
   The rest is the read onto the labels, which moves `ln R_bar` by up to 1.313e-06 and does not
   depend on the mesh at all.

   The mechanism is a difference of axes between the two productions. A level's label is the closure
   production's own pressure there, which that production formed with the composition read at the
   anchor's tabulated pressures; the transfer production re-forms the pressure with the composition
   read at the labels, which section 15 ruling 4 requires because at a target the levels are the
   union's. The two reads differ by the closure residual in their argument, 3.3e-3, and by 1.3e-6 in
   `ln R_bar`, so the identity cannot close below about 5e-7 however fine the mesh is. Excluding the
   bottom row does not remove it: 4.019e-07 at the pinned spacings and 4.113e-07 halved.

   Nothing is loosened: the check keeps the specification's condition and reports both numbers. What
   a ruling would settle is whether the quantity that carries the bound should be the tracing's
   part, with the composition floor reported beside it as `identity_floor` already reports decision
   G's placement floor, or whether the composition should be read at one axis by both productions,
   which is a change to the closure path and therefore to an accepted step.

4. **Check 4 fails: the cylinder-extended run's delivered altitudes are 53.79 m from their stated
   values, against a 10 m bound.** Measured 416,354, 99,309 and -15,450 m against the stated
   416,300, 99,290 and -15,456 (v0.10). The same run's T, p and N pass at 9.684e-04, 1.901e-04 and
   9.525e-04 against decision Q's 2e-3, and the same altitude path reproduces its stated values
   under the closure wind to 2.88 m at 10 N, 0.34 m at 60 N and 2.56 m at the anchor's own latitude.
   Nothing is loosened. Three things were measured before anything was proposed, and two candidates
   are dead.

   **Not the node set.** The target column is integrated on the 67 arrival levels (section 16 ruling
   1). Subdividing that set by 2, 4, 8 and 16 moves the closure-wind altitudes by less than 0.1 m:
   the RK4 is converged on the levels themselves.

   **Not the field.** Rebuilt here by the Step 2 construction, the cylinder wind gives u = 9.488,
   3.716, 2.167 and 1.925 m/s at the anchor's four named levels, which is what the Step 4 acceptance
   measured and REPORT_04_step4 records for the same construction.

   **Only a twelfth of it is the mesh.** `reports/step04_5/diagnose_cylinder_altitude.py` repeats
   the run on one cylinder field at 50,000, 25,000 and 5,000 m2/s2. The top altitude goes 416,353.8,
   416,344.0, 416,341.7: the departure falls from 53.79 to 41.68 m and stops. `r0` at the target is
   60,129,732.338 m at all three, unmoved to the millimetre. Under the closure wind the same sweep
   moves the altitude by 0.0 m, which is itself the mechanism of the 12 m: at 10 N the closure wind
   is 329.464 m/s at every level, so the isobar map never enters the column, while the cylinder wind
   runs from 326.4 to 427.1 m/s there and the column reads it at `exp(ln_p_at(Phi))`, whose
   resolution is the geopotential spacing. So decision R's coarsening costs 12 m of delivered
   altitude when the wind varies along a column, and nothing when it does not. That is worth
   recording for the performance step whatever the rest turns out to be.

   **What is left is u at the target column, and its size is measured.**
   `reports/step04_5/diagnose_cylinder_wind_sensitivity.py` integrates the same column with `u`
   offset by a constant: the top altitude is linear at 15.4 m per m/s over plus or minus 1 m/s. The
   41.7 m that does not converge is therefore 2.7 m/s in `u` at the target column, and the whole
   53.79 m is 3.494 m/s, against the 427.101 m/s the cylinder wind has there: eight parts in a
   thousand. For scale, the cylinder and closure winds differ by up to 97.637 m/s at that column and
   the delivered top altitudes differ by 5,219 m.

   The candidate the coding agent cannot settle from here is decision Q's own truncation. Step 4's
   check 4 measured and reported it at the wind: a field constant on cylinders, stored on the wind
   file's 361 latitudes by 61 pressures and read back by decision L's interpolant, is not constant
   on cylinders any more, and that suite reported the departure without bounding it. This run reads
   the cylinder field the same way. Eight parts in a thousand of `u` near a 427 m/s jet is the size
   that truncation plausibly has, but this filing does not assert it: the measurement that would
   settle it is the reviewing agent's own delivered altitude under the accepted per-hemisphere
   construction, or one run with the cylinder hypothesis carried on a grid fine enough that the
   round trip is not the limit. The stated altitudes are v0.10's, and decision P's construction was
   settled at Step 2, v0.11, where the review withdrew the v0.9 cylinder wind values for the same
   reason; whether the stated altitudes were restated after that is not something the coding agent
   can tell from the specification.

## 5. The regression

**Twenty five suites, twenty five passing, every one at its reference count.**
`reports/step04_5/regression.txt` carries the rows and each suite's own output is beside it. This
step changes `forward/production.py`, `lib/schema.py`, `lib/io.py`, `tools/plots/figures_profile.py`
and `forward/transfer.py`, all of which accepted suites exercise, so the run is evidence and not a
formality. The rows that carry the most of it: `step03_4` 13 of 13, the closure run and F5 and F6
through the refactored production and the changed figure module; `step03_3` 16 of 16, the namelist
and schema refusals; `step04_0` 15 of 15, which reads the transfer namelist whose geopotential
spacing this step moved; and `step04_4` 11 of 11, at the 5,000 m2/s2 it is pinned to, which is the
pin of section 0 doing its work rather than being asserted. It ran from 00:22 to 07:03, and
`step04_4` alone is four and a half hours of that.

The tree was put back as it was found: the registered kind N is `920304440ee4...` before and after,
the closure inputs and product agree as a digest of digests, the transfer inputs agree, the
registered figures agree, and `git status --porcelain` is empty.

**One thing about the run is the coding agent's fault and is recorded rather than tidied.** The
verification in the paragraph above was made by a separate script,
`reports/step04_5/verify_restored.sh`, because `run_regression.sh` did not reach its own tail. While
bash was executing it the coding agent edited the file, replacing one word in a comment to correct a
spelling. Bash reads a script by byte offset, so removing a byte shifted every later offset and the
file ended, as bash saw it, inside a quoted string: the console carries `line 83: unexpected EOF
while looking for matching '"'`. The loop had already been read into memory and ran in full, which
is why all twenty five rows are present and why the restore after each suite happened, including
after the last; what was lost is the tail that compares the digests and prints the porcelain. That
comparison is between the copies the loop set aside before it started and the files in the tree now,
so it is the same comparison made afterwards rather than in the same process, and the script says
so. No suite's result depends on it. The lesson is narrow and worth writing down: do not edit a
script while a shell is running it.

## 6. Hashes

Nothing is registered and no registered product changed, so there is no sweep. The registered kind N
is `920304440ee4...` and the registered closure product `fcc2e2c4ae71...`, both as they were at the
Step 4 acceptance commit `4802983`; check 11 rebuilds the closure production in a run directory of
its own and finds it bit-identical, so the registered file was read and never rewritten.

The transfer product this step delivers is
`forward/lindal_transfer/output/lindal_transfer_profile.nc`, SHA-256
`620f2deef79cfe4507b39628c12fb141a88f053c9d7e2d7a0d40b8a48604532d`, written at the pinned spacings
with the working tree carrying this step's changes, so its `casspian_git_commit` ends `-dirty`. A
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
