# REVIEW 02, Step 4. The kind N product and `casspian-refrac`

Review of `reports/REPORT_02_step4.md` against SPEC_02 v0.6 Step 4, SPEC_00 v0.13 §6.7 and
SPEC_01 v0.19. Reviewer: S. Rafkin, 12 September 2026. **Disposition: accepted with one
change**, conditional on the pending regression passing. SPEC_02 goes to v0.7 and SPEC_00 to
v0.14; both are small.

## Verified independently

The product was opened outside the acceptance script and the following were checked from the
file alone: `n = p / (k_B T)` reproduces `number_density_m3` to the last bit at all 66 levels;
`N` at 1 bar is 2.5985605e-4 with the tabulated 10.9 ppm ammonia; `refractivity_uncertainty /
refractivity` is 0.0233185 at every level, and the same fraction holds for the mean
refractivity companion; `radius_m` equals `r0` exactly at the anchor level and `r0` minus
90,000 m at 1 bar; the mean molar mass is 2.135083e-3 kg/mol dry and 2.136264e-3 at the bottom
level; the manifest text embedded in the `manifest` group is byte identical to
`lindal_reduction.toml` and its `sha256` attribute is the hash of that text; the three scalar
companions close in quadrature from the recorded terms (`r0` 16,376.8 m from 12,373.3 and
10,728.6; `phi_c` 0.1092° from 1.9057e-3 and 3.619e-5 rad; `psi` 0.00661° from 1.0963e-4 and
3.619e-5 rad); the anchor-rule spread of ±17.5 km equals the latitude-held partial 1.217 times
the 14.37 km by which either polar rule moves the effective anchor radius; and no em dash or
en dash appears in any attribute of any group. The report's numbers stand.

## Rulings on the findings

1. **The anchor-rule spread exceeds the propagated uncertainty.** Agreed on the substance and
   on keeping them apart. One change: the spread is to be measured by rerunning the whole
   Step 2 fixed point under each rule, not by two latitude-held marches. The v0.6 wording
   ("two extra marches") was the reviewer's shorthand and it leaves out the coupling through
   `dphi_c/dr_anchor`, which is exactly the term the Step 3 review ruled into the anchor
   term. A rule that moves the surface also moves the frozen latitude, and the recorded
   spread should be the whole effect of the choice. The cost is two more fixed points of five
   marches each. Record `r0` and `phi_c` under each rule. Expected: about ±17.8 km
   (58,537.7 and 58,502.1 km, each ±0.1 km) and `phi_c` shifted by ∓0.0030°, to about
   30.8020° under `north_pole` and 30.8079° under `south_pole`. Update the
   `anchor_rule_spread_note` to say what was done. SPEC_02 v0.7 Step 4 carries the expected
   values and decision 5 is restated.

2. **`_kg_mol` added to the unit suffixes.** Agreed; no specification text changes. The side
   effect on `lib.control`'s physical-value detection is in the right direction.

3. **The `input_hashes` warning for the other derived kinds.** Agreed that it is not this
   step's work. Deferred to the reader work SPEC_03 will need, recorded as SPEC_02 decision 8
   so it is not lost.

## On the decisions

All eight stand.

- Decision 1 (the hash recorded twice, refusal on a mismatch between the two copies, warning
  against the disk) is a clean reading of SPEC_00 §8 and the acceptance exercises both paths.
- Decision 2 (copied companions keep their declared kind; only combined companions are
  converted) is correct and is now SPEC_02 decision 7. The label's `0.2, range` in kind N is
  what kind T said, and the converted 0.115° appears where it was used, in the derived
  companions and in `reduction_record`.
- Decision 4: `index` for `anchor_isobar_pressure_Pa` is the right reading. SPEC_00 §5 said
  "coordinate" and now says "coordinate or scalar" (v0.14) so nobody has to reread the
  intent.
- Decision 5 is superseded by ruling 1.
- Decisions 3, 6, 7 and 8 are as specified.

## Order of work

1. Change the anchor-rule spread to the full fixed-point rerun under each rule; record `r0`
   and `phi_c` under each; update the note; rerun `accept_step02_4.py` with the v0.7 expected
   values added as a check; replace `REGRESSION_PENDING` in the report with the regression
   result and add a short addendum for the change.
2. Commit. Sweep: `lindal_refractivity.nc` is rebuilt clean; record its SHA-256 in the report
   addendum, superseding nothing (it is the first hash of that file).
3. Step 5, the standard diagnostics, against SPEC_02 v0.7 §5 and SPEC_00 v0.14 §3.5 and §7.1.
   Read the status lines and STATE.md first.
