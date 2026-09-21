# REPORT 04, pre-execution. The coding agent's review of SPEC_04 v0.5

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 21 September 2026. Specification reviewed: `SPEC_04_Transfer.md` v0.5, 17 September 2026,
all of it, as the handoff of 17 September section 4 directs. Rulings: SPEC_04 v0.6 section 10,
with the changes applied in the body; the specification is accepted at v0.6 and carried at v0.7.
Status: **filed**. No code was written for this review.

**Seven findings and seven shorter items were sent. Every one was verified by the reviewing
agent and ruled correct; every one is accepted here.** Two entries of the record are corrected
in section 4, neither of which changes a ruling: the interpolant that produced the value 1.93991
m/s is PCHIP and not a cubic (a cubic spline on the same nodes gives 1.99409 m/s, measured), and
the bound on the composition regrid quoted in section 10 is the profile-wide bound, larger than
and correctly superseding the gauge-level number sent. One point is carried forward to Step 3
(section 5): with a piecewise-constant derivative and a mesh spacing that divides the file's, the
mesh nodes coincide with the file's nodes, where the derivative jumps, and which one-sided slope
such a node takes is not yet stated.

## 0. Before this review

| Commit | What |
|---|---|
| `1504056` | SPEC_03 Step 4, `forward/production.py`, `casspian-forward`, the closure run. |
| `545b385` | The SPEC_03 Step 4 sweep record. |
| `e2e7eb4`, `2cd41bd`, `5f23a8b` | `STATE.md`, SPEC_03 v0.14 closed, the REVIEW_03_step4 addendum. |

HEAD at the review was `5f23a8b`, not `545b385` as the handoff of 17 September section 2 states.
The three commits after `545b385` are the SPEC_03 closure, the `STATE.md` update and the review
addendum; none is a step record and none changed a product, so nothing of substance differed.
Reported as the handoff asks and corrected by the author. The untracked and modified author
documents were as section 2 describes.

## 1. Method

Read: the handoff of 17 September; SPEC_04 v0.5 in full; `STATE.md`; SPEC_00 sections 6.2, 6.6,
7.2 and 7.3; SPEC_03 section 0 and the Step 3 and Step 4 deliverables; REPORT_03_step3 and
REPORT_03_step4; `docs/CASSPIAN_Seasonal_Design_Note.md` sections 5.4 and 10; and in the
repository `lib/schema.py`, `lib/io.py`, `lib/control.py`, `lib/geopotential.py`,
`forward/production.py`, `refrac/anchor.py`, `tools/wind/curve.py`, `tools/wind/fit.py` and
`tools/composition/build_composition.py`.

The equation labels were checked against `CASSPIAN_AtmosphericModel_Draft9_2.docx`, the draft of
record, read outside the repository. Labels A1 to A40 and B1 to B8 are present; A10, A15, A16,
A24, A25, A27 to A35, A39 and A40 and B3, B5, B6, B7.1 to B7.5 and B8 were read directly.

Measured with the repository's own products and with numpy and scipy: the closure wind file
(`forward/lindal_closure/inputs/lindal_closure_wind.nc`, 361 latitude nodes at 0.5 degrees, 61
pressure nodes, `reference_level_pressure_Pa` 1.0e5, no variation along the column); the closure
product (`forward/lindal_closure/output/lindal_closure_profile.nc`, 66 levels); the closure
composition; and `occul_data/lindal/lindal_refractivity.nc`. Every number below was measured, not
taken from a document.

## 2. The findings as sent, the rulings, and the disposition

### Finding 1. The PCHIP rule of Step 1 deliverable 2 contradicts four of the draft's own checks

**Sent.** Step 1 deliverable 2 read `u_total` in latitude by PCHIP, "the rule the wind tool
used". Measured: the wind file's nodes bracketing the anchor are 5.8557 m/s at 30.50 degrees and
−0.1800 m/s at 31.00 degrees, so the anchor at 30.805568 degrees sits on a zero crossing whose
slope is about 12 m/s per degree. Linear interpolation gives 2.167086 m/s, which is what the
reduction and the closure used (`u_at_phi_c_ms` = 2.16708095 in the closure product's
`production_record`); PCHIP gives 1.939910 m/s, lower by 0.2272 m/s, 10.5 percent. The rule
therefore broke Step 1 deliverable 3's expected 2.167 m/s, Step 1's "bit-identical to the closure
product", Step 1's agreement with the reduction's own march to 0.01 m (the reduction marches with
`refrac.anchor.wind_of_latitude`, which is linear), and Step 5's fourth expected value at `phi_c`
to 1e-12. The parenthetical was also wrong about the tool: `tools/wind/fit.py` fits a constrained
penalized B-spline and PCHIP appears only in `tools/wind/curve.py` as the pole extension. The
model cannot carry two interpolants of one file, and which one it carries is a physics choice
here, because the anchor's own wind is set by it at the 10 percent level.

**Ruling (section 10 finding 1, decision L).** Correct in every part, the fault the
specification's. One rule for data on a latitude grid: linear between the nodes, the reduction's
rule, with the derivatives those of the interpolant. Every expected value remeasured under it in
section 1 and Steps 1, 3, 4 and 5. The reference surface under the linear rule returns the
reduction's equatorial radius to 0.03 m against 0.55 m under the cubic, which settles which rule
registered the anchor.

**Disposition: accepted.** The ruling resolves all four checks by making the model's rule the
rule that produced the anchor's registration, and the equatorial-radius residual is an
independent confirmation the review did not have. See section 4 item 1 on the name of the
interpolant that gave 1.93991.

### Finding 2. Rebuilding kind W cascades through every registered product

**Sent.** Step 0 deliverable 1 said only that three wind files are rebuilt with no value of
`u_total` changing. But `occul_data/lindal/lindal_refractivity.nc` records
`lindal_wind.nc sha256:fd1d82af...` in both the global `input_hashes` and `reduction_record`, and
`lib.io._check_refractivity_hashes` refuses a kind N file whose recorded hash is not listed
(`_REFRACTIVITY_HASHED` is the six inputs and the manifest, and the reader also requires exactly
seven entries). Rebuilding the reduction's wind file therefore forces the whole reduction chain,
the closure run's inputs and the closure product to be rebuilt and rehashed, so REPORT_03_step4's
product hash `e6693173...` and the `step02_1` suite's hash table both go stale.

**Ruling (section 10 finding 2).** Correct. Step 0 deliverable 1 now states the cascade: the
chain, the closure inputs and the closure product rebuilt as SPEC_03 Step 3 did, by the in-memory
candidate rule and then the sweep on the clean tree, no value changing in any product and only
variables added, verified group by group; the new hashes recorded in the report in the `step02_1`
row format; REPORT_03_step4's product hash superseded there.

**Disposition: accepted.** This is the work of Step 0 as it now stands.

### Finding 3. `produce` does not take the wind as an array

**Sent.** Step 1 deliverable 3 asserted that SPEC_03's `produce` takes the wind on the anchor's
levels as an array. It does not: `forward/production.py` line 114 is
`produce(profile, inputs, gauge, p_b)`, and line 132 forms a scalar internally through
`wind_of_latitude`; `Production.u_ms` is a scalar and is recorded as `u_at_phi_c_ms`.
`lib.geopotential.effective_gravity_magnitude` already broadcasts an array to the levels, so the
change is small, but it is a change to closed SPEC_03 code and to a product record, and the
specification should say it and say what the record becomes.

**Ruling (section 10 finding 3).** Correct. Step 1 deliverable 3 states the amendment: an
optional wind array with the scalar it forms today as the default, and `u_column_ms` recorded.
The Appendix lists the amendment against SPEC_03 Step 4, with closure mode and its expected
values unchanged.

**Disposition: accepted.** The optional argument with the present scalar as the default keeps the
closure bit-identical, which is what the Step 1 acceptance needs.

### Finding 4. The M = 2 identity test could not load

**Sent.** Step 0 deliverable 4 requires `gauge_isobar_Pa` to be a tabulated level of every
occultation anchor, matched exactly, and Step 4 wrote the synthetic anchor's thermo pressures as
its labels. Measured in the closure product: the label at the gauge level, which is the anchor's
produced pressure there, is 9998.4655 Pa, not 1.0e4 (the residual `p / p_tab − 1` is
−1.5345e-4). The M = 2 namelist would have been refused.

**Ruling (section 10 finding 4, decision M).** Correct, and the fix belongs in the synthetic
anchor, not in the rule. The rule stands and an anchor without the gauge among its tabulated
levels is refused as not implemented. The synthetic anchor carries the Lindal anchor's tabulated
pressures in its thermo group, the same isobars by construction, with temperatures
`p_tab R_bar / (k_B N)` so its reduction identity returns the transferred `N` exactly; its gauge
level is then a tabulated level matching `gauge_isobar_Pa` and its anchor isobar is the gauge
level.

**Disposition: accepted.** Keeping the rule and fixing the test article is the stronger reading:
the synthetic anchor is meant to be what an occultation at 60 degrees N reduced by Lindal's
method would tabulate, and it now is. The temperature form is B6 (`T = p R_bar / (k_B N)`), so
the reduction identity closes by construction.

### Finding 5. A composition on a latitude grid had no route into `produce`

**Sent.** `mean_properties` returns `R_bar` and `m_bar` on the composition's own levels and
`produce` uses them positionally. With `latitude_grid_deg = [-90, 90, 1.0]` the file becomes
(level, latitude) and nothing said who interpolates in latitude, nor whether the target's `R_bar`
and `m_bar` are taken by level index or interpolated in pressure onto the isobar labels. Measured
on the closure composition: the largest `|d ln R_bar / d ln p|` between levels is 4.495e-4, so
interpolating onto the labels rather than taking them positionally moves `ln R_bar` by up to
6.9e-8 at the gauge level, five orders above the 1e-12 of Step 5's fourth expected value.

**Ruling (section 10 finding 5, decision L).** Correct. Decision L covers kind C as well: a field
on `(phi, ln p)`, linear in latitude and log-linear in pressure onto the labels in transfer mode,
closure mode unchanged. Step 5's fourth expected value is restated as `N` and `Phi` to 1e-12 and
`p` and `T` to 1e-5, with the bound 1.5e-6 from the profile-wide residual.

**Disposition: accepted.** The split tolerance is right: `N` and `Phi` do not pass through the
composition and stay at 1e-12, while `p` and `T` do and cannot. See section 4 item 2 on the
bound.

### Finding 6. The "3.5 to 4.5 under halving" convergence claims do not hold for a C1 interpolant

**Sent.** Step 3 claimed the centered difference across columns reproduces `u'(phi)` with the
error falling by 3.5 to 4.5 under halving of the latitude spacing, and Step 4 repeated the ratio
for the out-and-back identity. Measured with PCHIP on the file's 0.5 degree grid, at the file
node 30.50 degrees, the centered difference is −790.35, −791.58, −786.65 and −782.80 m/s per
radian at spacings 0.5, 0.25, 0.125 and 0.0625 degrees against the interpolant's own derivative
−778.02, so the errors are 12.3, 13.6, 8.6 and 4.8: non-monotone and first order at best,
because the second derivative of a C1 interpolant jumps at every node. Refining the mesh below
the file's grid can converge only on the interpolant's derivative, and the file's values are
themselves a penalized spline fit sampled at 0.5 degrees, so `u'(phi)` has no truth behind it.
The reviewing agent's expected values were measured at 0.05 degrees with a cubic, a third rule.
One mitigating fact was measured and sent with the finding: the transfer's `Delta ln N`
telescopes to the endpoint wind difference, about 398 m/s from `phi_c` to 10 degrees N, so a
0.227 m/s endpoint error is about 6e-6 in `Delta ln N`, inside the stated 1e-4.

**Ruling (section 10 finding 6).** Correct; the claims are withdrawn. Under decision L nothing is
differenced in latitude on the mesh, `S/g` returns the closed form to round-off at the nodes, and
the tracing and the transfer are first order across the file's nodes; the halving ratios are
reported and not bounded, and the mesh spacing divides the file's. The telescoping is recorded in
decision L as the reason the transfer values are insensitive to the rule.

**Disposition: accepted.** Removing the difference from the mesh removes the question rather than
answering it, which is the better outcome: Step 3 now forms `S` from the interpolant's own
derivatives converted to fixed `r` through the isobar map, so there is nothing left to converge.
Reporting the ratios without a bound is right, since no order can be claimed. One consequence is
carried to Step 3 in section 5.

### Finding 7. `phi_r` is needed before the mesh exists

**Sent.** Step 2 required the gauge latitude as an exact mesh node, Step 4 computed `phi_r`
inside the estimate, and Step 5's driver order was load, surface and wind, mesh, kernels, loop,
trace, estimate. The three cannot all hold.

**Ruling (section 10 finding 7, decision K).** Correct. `phi_r` is fixed by the driver before the
mesh, from the anchors as they arrive, and inserted as a node;
`forward.estimate.gauge_latitude(anchors)` is the one place it is computed and the estimate takes
it as an input. The driver order is restated.

**Disposition: accepted.**

## 3. The shorter items

| Item as sent | Ruling | Disposition |
|---|---|---|
| `weight` is required by `lib.control` (`_ANCHOR_KEYS`) and defined in SPEC_00 section 7.2 as the anchor's role, construction 1 and validation 0, and SPEC_04's estimate read only `1 / (sigma^2 + P)`, so a validation anchor would be silently treated as construction. | Honored: Step 0 deliverable 4 defines it, a validation anchor is propagated, placed and traced and its `C_i` and every `D_ij` are reported but it enters neither `phi_r` nor `C`, and any other value is refused. Step 4's acceptance exercises it with the synthetic anchor at `weight = 0`. | Accepted. The reading proposed as an alternative (refuse anything but 1) is the weaker one; the ruling keeps the diagnostic value of a validation anchor, which is what SPEC_00 section 7.2 intends. |
| The centroid weight `1 / mean(sigma_i^2)` drops `P_i` and collapses `Phi`, a departure from A34 whose weights are A30's; the departure is what keeps `phi_r` out of a fixed point, but it should be recorded as a departure. | Decision K records it, with a note for the manuscript that A34's weights are circular through `P_i` and the text should say which weights place the gauge. | Accepted. |
| Section 1 attributed the variance inflation rule to "A10's rule"; A10 is the meridional momentum balance `(1 / rho) grad p = g_eff`, and the inflation statement is in the discussion of A33. | Corrected to the rule stated with A33. | Accepted. |
| Step 0 deliverable 4's text listed `[numerics.geopotential]` and `[numerics.hydrostatic]` among the restricted schemes, its own namelist omitted them, and its acceptance made an unknown `[numerics]` table an error. | The list is corrected: the geopotential and hydrostatic rules are SPEC_03's closed ones and are not namelist keys, as in the closure namelist. | Accepted. The namelist of deliverable 3 and the text now agree, and an unknown `[numerics]` table stays a refusal. |
| The synthetic anchor's kind N requirements were unstated: `lib.io.read(path, "refractivity")` requires the eight groups of the kind, exactly seven `input_hashes` entries, an `input_sha256_<key>` in `reduction_record` for each of the six inputs and the manifest with every recorded hash listed, and the root scalars with their uncertainty companions. | They are the schema's; the acceptance script writes the anchor through `lib.io.write` with the groups, the seven self-consistent entries and the root scalars the writer validates, and the report lists them. | Accepted. |
| Step 2 inserts up to three latitudes into a uniform grid, so short intervals exist (`phi_c` = 30.805568 sits between 30.5 and 31.0), and "centered with the neighbors" is first order on them. | The inserted nodes take the three-point formula for unequal spacing wherever a difference is taken on the mesh in latitude (Step 2). | Accepted. Under decision L no difference of `u` is taken in latitude at all, so what remains on the mesh are the isobar map's slopes, and the unequal-interval formula covers them. |
| The Appendix cited `CASSPIAN_Manuscript_Retrievals_and_Wind_Note.md` v0.4, which is not in `docs/`; and the draft of record for the equation labels was not named. | That note and the manuscript stay outside the repository while the manuscript is in flux, and neither is a dependency of any step; the Appendix now says so, and `CASSPIAN_AtmosphericModel_Draft9_2.docx` is named as the draft of record in the status line and in `STATE.md`. | Accepted. |

## 4. Corrections to the record

1. **The interpolant that gave 1.93991 m/s is PCHIP, not a cubic.** Section 10's preamble and
   decision L both attribute that value to a cubic. It is the value sent with finding 1 and it
   was measured with `scipy.interpolate.PchipInterpolator`. A cubic spline on the same nodes at
   the same latitude gives a different number: measured now,
   `scipy.interpolate.CubicSpline` 1.994094 m/s and `make_interp_spline(k=3)` 1.994094 m/s,
   against PCHIP 1.939910 and linear 2.167086. So the departure from linear is 10.5 percent
   under PCHIP and 8.0 percent under a cubic. Nothing in the ruling turns on which: both are far
   outside any tolerance in the specification, and decision L fixes the linear rule. Recorded so
   that the number and the interpolant are not paired wrongly if either is quoted again, and so
   that the reviewing agent's own expected values, which section 10 says were measured with a
   cubic, are understood to have come from the 1.994 branch and not the 1.940 one.
2. **The composition bound quoted in section 10 is the profile-wide one and supersedes the
   number sent.** Finding 5 was sent with 6.9e-8, which is the largest
   `|d ln R_bar / d ln p|` of 4.495e-4 times the pressure identity residual at the gauge level,
   1.5345e-4. Section 10 uses the largest residual anywhere on the profile, 3.29e-3 at the bottom
   row, and gets 1.5e-6. That is the bound the tolerance must respect, and it is the right one;
   the number sent was scoped to the gauge level only and was too small by the ratio of the two
   residuals. Recorded because Step 5's restated tolerance rests on 1.5e-6 and not on 6.9e-8.

## 5. Carried forward

**The value of a derivative at a node the mesh and the file share.** Decision L makes
`(du/dphi)_p` the interpolant's, piecewise constant, and Step 4 notes that the mesh spacing
divides the file's so that the nodes coincide. At a coincident node the piecewise-constant
derivative has two values and the specification does not say which the node takes. It is not
academic: measured on the closure wind at the file node 31.0 degrees, the slope southward is
−12.07 m/s per degree and northward −8.39, a difference of 30 percent in `S/g` at that node,
and Step 3's expected value for the largest `|S/g|` sits on the interval holding the anchor,
immediately south of that node. The transfer integral telescopes, so `Delta ln N` is unaffected,
but the largest `|S/g|` is an acceptance value at Step 3 and `S/g` is drawn in F5 and F7. This
bears on Step 3 and not on Step 0, so it is recorded here and will be raised as a finding at
Step 3 unless the specification settles it first. No rule is assumed in the meantime.

## 6. Next step

Step 0, by the procedure of the handoff of 17 September section 6: kind W in three parts, kind C
on a latitude grid, `forward/lindal_transfer/`, the transfer namelist for M anchors, and the
loader with the anchor object and the propagation hook, with the rebuild cascade of deliverable 1
on an in-memory candidate and the sweep on the clean tree. `STATE.md` carries the SPEC_04 section
and the status line says Step 0 proceeds once this report is filed.

One matter for the author before the sweep: `docs/CASSPIAN_Seasonal_Design_Note.md` is modified in
the working tree and was not named among the documents to commit. The sweep needs a completely
clean tree, so it is committed, stashed or reverted before the Step 0 sweep, at the author's word.
