# SPEC 06. The transfer kernel read on the isobar

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.4, 2 October 2026. Author of record: S. Rafkin. **Status: accepted by the author on
2 October 2026 at v0.2, with the rulings of §4 on the coding agent's reading of v0.1; v0.3 extends
the step to the tracing after the interim report (§5); v0.4 rules on the second interim report
(§6).** Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.22, the closed `SPEC_04_Transfer.md`
v0.23 and `SPEC_05_Shear_Experiments.md` v0.7, which it does not repeat. It starts after SPEC_05
Step 4 is accepted and committed. (The number SPEC_06 was earlier pencilled in for end-to-end tests;
those move to a later number.)

---

## 0. Why, and how this specification is to be worked

**Why (author, 2 October 2026).** In SPEC_05 run 7f, F8 shows the shear term along the 999 mbar
isobar as a sawtooth across the whole path from 31 N to 10 N. The cause is in how the transfer reads
its kernel. `lib.kernel.shear_kernel` forms `S` (Eq. A15) once at each mesh node, using the wind's
derivatives at that node's own pressure, and `forward.transfer.transfer` then reads `S/g` along each
curve by bilinear interpolation between nodes (`lib.kernel.transfer_kernel`). Where the wind's shear
changes abruptly in pressure (a kink, here at `p_s` = 1 bar at every latitude), an isobar lying
within one geopotential cell of it gets a blend of the nodes on either side. As the isobar tilts
across the mesh, the blend slides and resets at every node crossing: the sawtooth. It is a mesh
artifact, not the hypothesis. The Monte Carlo will draw hypotheses with kinks at random pressures,
and nobody will inspect every draw, so the transfer must not depend on where a kink falls relative
to the mesh.

**The remedy.** On an isobar the pressure is fixed: it is the curve's own label. So the part of the
kernel that comes from the wind is read exactly at the curve's point and label, from the wind file's
interpolant, and only the smooth geometry is read from the mesh. A kink in the wind is then resolved
at the wind file's own resolution, which is the hypothesis as written, and a kink placed on a wind
file node is represented exactly. Nothing is smoothed and no input changes.

**Isolated on purpose.** This is its own specification so that it can be reverted, or tried further,
without touching SPEC_05. One step.

**Protocol.** SPEC_04 §0 and SPEC_05 §0 apply: report and review, nothing committed before the
author's go, no AI attribution, no dashes, the regression rule of SPEC_04 v0.22 with the author's
direction of 1 October (the long suites rerun when `lib`, `refrac`, `forward` or the registered
products change; this step changes `lib` and `forward`, so the full set). Tolerances are physical,
not numerical (author, 2 October): a fraction of a kelvin is noise.

---

## 1. Step 1: the shear term evaluated on the curve

**Deliverables.**

1. **`lib.kernel`: the shear term at a curve point.** A function that, given the curve points
   `(phi, Phi)` and each curve's label pressure `p_label`, returns `S/g` there by Eq. A15:
   - from the wind file, exactly at `(phi, p_label)`: `u`, `(du/dphi)_p` and `(du/dln p)_phi`,
     through `WindField.wind_at` and `WindField.wind_derivatives` (the same interpolant, decision L);
   - from the mesh, bilinear at `(phi, Phi)` as now: the radius `r`, the radial gravity `g`, and the
     two isobar slopes `(dln p/dphi)_r` and `(dln p/dr)_phi` that `isobar_slopes` already forms on
     the nodes (these are smooth fields of the pressure map, and reading them bilinearly is sound);
   - then `(du/dphi)_r`, `(du/dr)_phi`, `(du/dZ)_R` and `S = 2 Omega_abs r (du/dZ)_R` exactly as
     `shear_kernel` forms them now, `Omega_abs` from `u` at the point, and `S/g`.
2. **`forward.transfer.transfer` uses it** for the shear part of `K` (Eq. A27) at every node of
   every curve, with the curve's label. The composition term is unchanged.
3. **The tracing reads the same shear (v0.3).** The tracing, `dPhi/dphi = -I` (Eq. A24), needs
   `I`, the integral of `S/g` in geopotential from the reference surface up to the curve (Eq. A16).
   It is formed from the same exact shear as the transfer, not from the node kernel interpolated:
   - **The principle, general and not about any one kink.** The wind enters the model only through
     the wind file's interpolant. The mesh never resamples it. Wherever the wind is needed it is
     read from the interpolant at the point itself, and wherever it is integrated the wind file's
     own pressure nodes are breakpoints of the quadrature. Every change of slope, curvature or
     higher derivative the wind file can hold sits at those nodes (between them the interpolant is
     a single smooth piece), so any number of them, anywhere, is integrated at the wind file's own
     resolution, whatever the mesh spacing.
   - **How.** In each column, `I` at the node below the curve is the sum over the cells beneath it,
     each cell integrated piece by piece between the breakpoints it contains (the pressures where
     the column crosses a wind file pressure node, found from the pressure map); then the partial
     cell from that node up to the curve's own geopotential, the same way. On each piece, between
     breakpoints, `S/g` is smooth, and a two-point Gauss rule (exact for cubics) is used, so the
     rule stays adequate if the interpolant later becomes cubic (the next specification). `S/g` at
     each quadrature point is Eq. A15 with the wind read from the interpolant at that point's
     pressure, as in deliverable 1, and the geometry bilinear on the mesh. At the RK4 stages between
     mesh latitudes, `I` is interpolated linearly in latitude between the two neighboring columns,
     whose latitudes are wind file latitude nodes or lie inside one wind cell.
   - Structure finer than the wind file's own grid is not in the hypothesis, and nothing here
     recovers it; the wind file's resolution is set by the person making the hypothesis.
   F8 draws the shear term the transfer actually used.
4. **No new key, no new product variable.** Reverting is a git revert of this step's commit.
5. **The label, not the curve's pressure.** The wind is read at the curve's label. A curve meets its
   label only to within the pressure identity (up to 9e-3 in SPEC_05's decay runs, 5.8e-7 in the
   transfer run); the report states the difference. That is the intended choice: the label is the
   isobar.
6. **A label exactly on a wind file node** takes the existing convention of
   `WindField.wind_derivatives`: at a node, the cell toward higher pressure (and toward the north in
   latitude). It is stated, not changed. A Monte Carlo draw puts a kink on a node with probability
   zero, and a kink on a node is still represented exactly in the wind's values.
7. **`transfer()` takes the labels and what the new kernel needs** (the wind field, the rotation
   rate and the node geometry formed once per state). Every caller passes them, including accepted
   suites and diagnose scripts that call `forward.transfer.transfer` directly (`step04_4`,
   `step04_5`, two diagnose scripts). That is a mechanical signature change, named in the report as
   such and separate from any moved expectation. `shear_kernel`, `transfer_kernel` and `step04_3`
   are unchanged. F8's variable description is reworded to say where the shear term is read.

**Acceptance (`tests/step06_1/accept_step06_1.py`).** Runs are made in copies under
`reports/step06_1/`, from the working tree with the `-dirty` refusal relaxed for the script only, as
in SPEC_05.

1. **A run with no vertical shear barely moves.** The transfer run (`lindal_transfer`) rebuilt with
   the new kernel against the registered transfer product: delivered temperature within 0.1 K at
   every level, altitude within 10 m, `r0(10 N)` within 1 m. The measured differences are reported;
   they are expected to be far smaller, since without vertical shear the wind term of `S` does not
   depend on pressure.
2. **The closure product does not move** (it has no transfer): `step04_5` check 11 passes.
3. **The sawtooth is gone.** Run 7f's case at 5e4, along the 999 mbar isobar: no pattern tied to
   the mesh's geopotential nodes. Between consecutive curve nodes that are not wind file latitude
   nodes, `S/g` may still ramp, since the wind's `(du/dln p)` varies linearly in latitude across each
   0.5 deg wind cell (that is the hypothesis as written); the largest such step must be at least ten
   times smaller than the largest such step in SPEC_05's run 7f along the same isobar. The old and
   new `S/g` are drawn together, and the author views the new F8.
4. **The kink no longer leaks across `p_s`.** Run 8's case (`increase_below`, `p_s` = 1e5 Pa) at
   5e4: the delivered temperature change from the no-shear run at the 998.7 mbar level, which lies
   above `p_s` where this case has no shear, within 0.1 K of zero (SPEC_05 measured -5.87 K). The
   tracing still reads `I` from the mesh (deliverable 3), and within the cell straddling `p_s` it
   carries a small share of the shear below; if the check misses, the measurement and the tracing as
   its cause are reported for a ruling, and nothing is tuned.
5. **The kink level does not depend on the mesh.** Run 5's case at 5e4 and 2.5e4: the delivered
   temperature at the 998.7 mbar level agrees between the two spacings within 0.5 K; it is expected
   to sit at the decay zone's value (SPEC_05 measured about half of it, -8.7 K against -17.2 K).

**Also delivered, reported not checked.**

- **The SPEC_05 Step 4 experiments rerun** with the new kernel, through
  `tests/step05_4/run_experiments.py` at both spacings: a new results table beside the old, the
  comparison figure redrawn, and run 7f's F8 viewed by the author. Expected: runs 2 and 3a to 3c
  essentially unchanged, the decay and increase runs changed only near `p_s` and the stop pressure,
  and the outer loop at least as well behaved as before.
- **The two runs that cycled** (run 6 at 5e4, run 7 at 2.5e4), with the per-pass residual printed:
  whether the mesh blend at the kink was what set them cycling. Reported; nothing tuned.

**The registered files.** The transfer product's values move (check 1), so after the acceptance
commit `python tests/step04_0/sweep.py --runs` rebuilds the two runs' registered files on the clean
tree; the closure files come out with unchanged values (check 2); both are recommitted and the
comparison is added to the report.

6. **The pressure identity does not depend on the mesh (v0.4).** For every run of the SPEC_05
   Step 4 set, the largest identity at 5e4 and at 2.5e4 agree within 10 percent; its size and level
   are reported, not bounded (§6 ruling 1).
7. **The identity at the stop pressure is explained (v0.4).** Runs 5 and 6 show at their stop
   pressure an identity whose sign is opposite to the trapezoid estimate that explains the one at
   `p_s`. The report gives the mechanism, with a measurement that shows it (§6 ruling 2).

Check 4 is expected to pass under deliverable 3, since the leak it measured was the tracing's.
The two runs that cycled are rerun under the new tracing, at the experiments' tolerance and at 1e-8,
and reported.

**Regression.** `lib` and `forward` change: the full set, about 100 minutes. An accepted suite whose
expectation moves only because the kernel is now read on the isobar (the SPEC_04 transfer and
cylinder-wind values, bounded at decision Q's floor; the SPEC_05 Step 4 suite) is updated in this
step, each named in the report with its old and new value and the reason. Reference counts stay
unless a check is added or removed. No accepted suite is edited silently or left failing.

---

## 2. Decisions of v0.1 (confirmed by the author, §4)

1. The tracing keeps reading `I` from the mesh (deliverable 3). If the experiments still show a
   mesh-dependent signal at a kink after this step, the tracing is the next place to look.
2. Reverting by git, with no switch between the old and new kernel in the code (deliverable 4). A
   switch would let one tree compare both, at the cost of a key and a second code path.
3. The bounds of checks 1, 4 and 5 (0.1 K, 0.1 K, 0.5 K); check 3 restated at v0.2 (§4 ruling 3).

---

## 3. Revision history

- v0.1, 2 October 2026: first draft, from the run 7f F8 sawtooth and the author's direction to fix it
  in its own specification.
- v0.2, 2 October 2026: accepted by the author; the coding agent's reading of v0.1 ruled on (§4):
  check 3 restated as no pattern tied to mesh nodes; deliverables 5 to 7 added (the label, the node
  convention, the signature change); check 4's dependence on the tracing stated.
- v0.4, 2 October 2026: the second interim report's rulings (§6): check 6 restated as
  mesh-independence; check 7 added, the stop pressure explained before acceptance.
- v0.3, 2 October 2026: the interim report's ruling (§5): the tracing reads the same exact shear,
  integrated with the wind file's pressure nodes as breakpoints (deliverable 3 rewritten as the
  general principle); check 6 added.

---

## 4. Rulings on the coding agent's reading of v0.1 (2 October 2026)

1. **The signature change.** Accepted as deliverable 7: the accepted suites that call `transfer`
   directly are edited mechanically and named in the report as a signature change, apart from any
   moved value.
2. **Label against the curve's pressure.** Accepted as deliverable 5.
3. **Check 3's bound.** The coding agent is right that `S/g` ramps inside a wind cell as the
   hypothesis is written; check 3 now tests what the remedy promises, no pattern tied to mesh nodes,
   at least ten times below SPEC_05's largest off-node step.
4. **Check 4 and the tracing.** Accepted: measured, and a miss reported with the tracing as its cause.
5. **A label on a wind node.** The existing convention, stated (deliverable 6).
6. **The rest.** F8's description reworded; decisions 1 and 2 of §2 confirmed; the bounds of checks
   1, 4 and 5 kept.

---

## 5. Rulings on the interim report (2 October 2026)

`reports/REPORT_06_step1.md` (interim): the sawtooth along the isobar is gone (check 3, 135 times
smaller), the run with no vertical shear moves by 1.4e-8 K, the closure does not move, and run 5's
998.7 mbar level now sits at the zone's value at both spacings. But with the transfer exact and the
tracing still blended at a kink, the two disagree there: the pressure identity grows up to 40 times,
above 1e-2 in four runs, the delivered temperature steps in the lowest 50 km of run 7f, and run 8
leaks -0.28 K at 998.7 mbar through the tracing.

1. **Extend the reading to the tracing** (the coding agent's option 1), not revert. Deliverable 3
   rewritten.
2. **Stated as a principle, not a patch for one kink** (author): a real or invented profile can
   change slope, curvature or higher derivatives anywhere. The model therefore never lets its mesh
   resample the wind; it reads the wind file's interpolant where it needs it and integrates with the
   wind file's nodes as breakpoints, so whatever the file holds is used at the file's resolution.
3. **The steps that remain along an isobar at every 0.5 deg** (the author's eye on the new F8) are
   the derivative of a wind read linearly between nodes (decision L), integrated exactly; they are
   not numerical instability. Removing them is a change of interpolant, the next specification
   (SPEC_07, decision L revisited), not this one.

---

## 6. Rulings on the second interim report (2 October 2026)

Under v0.3 (REPORT_06_step1 §7): checks 1 to 5 pass; run 8's leak is 0.0000 K; the identity is the
same at both mesh spacings to three digits; and the two runs that cycled under SPEC_05 now converge
at 1e-8 in 9 and 11 passes, so the cycle was the mesh blend in the tracing.

**The author's standing rule (2 October):** the code must work on any profile, not only realistic
ones; there are no strong constraints on what realistic means. An unexplained behavior in an extreme
test case is explained, not deferred.

1. **The identity at `p_s`.** It is the production's integration over the anchor's own 66 levels
   across a temperature break that the hypothesis puts between two adjacent levels (the report's
   trapezoid estimate matches it within 15 percent in all seven sheared runs). No mesh setting reaches
   it, so check 6 is restated as mesh-independence, with the size reported.
2. **The identity at the stop pressure** (runs 5 and 6) has the opposite sign to the same estimate
   and is not explained. It is explained in this step, before acceptance (check 7): what sets its
   sign and size, shown by a measurement (for example, the same estimate taken with the stop pressure
   placed on an anchor level, between two levels, and on a wind file node).
3. **The stepping in run 7f's lowest 50 km** is most likely decision L's staircase in the vertical
   (the wind's nodes every 12 percent in pressure, the anchor's levels every 4.5 percent). It goes to
   SPEC_07, whose continuous interpolant should remove it; SPEC_07's acceptance includes that F8.
4. **Order.** Check 7 first; then the full regression, the registered transfer product rebuilt and
   recommitted, and the report filed for acceptance.
5. **Noted for later, not in this step:** several extreme cases deliver superadiabatic profiles (run
   7f falls about 6 K/km below 1 bar at 10 N, several times the dry adiabat). A stability diagnostic,
   and what to do with an unstable result (reject the draw, or adjust and re-derive the wind), is a
   later specification.

