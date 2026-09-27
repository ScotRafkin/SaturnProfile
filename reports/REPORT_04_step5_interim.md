# REPORT 04, Step 5, interim. The production at the target, altitude and datum, the product

Coding agent, 27 September 2026.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.16, section 7, decisions B, C, F, G, H and R,
and the section 15 and 16 rulings.

**This is not a step filing and asks for no acceptance disposition.** Step 5 is part built: the
four deliverables exist and run end to end, and the acceptance script does not exist yet. It is
filed for review in the manner of `reports/REPORT_04_preexecution.md`, because five decisions were
taken that the specification did not take, the acceptance script will harden around them, and it is
cheaper to rule on them now than after it is written. The coding agent's context filled at this
point; `docs/CASSPIAN_CodingAgent_Handoff_2026-09-27.md` hands the work on.

## 1. What is built

**Deliverable 1, the production at the target.** `forward.production.produce_on_geopotential(N,
Phi, R_bar, m_bar, p_b)` is the production both paths enter (B4, B5, B6); `produce` forms `Phi` and
calls it, so closure and transfer are one production and not two.
`mean_properties_on_labels` reads the composition log-linearly onto the isobar labels (section 15
ruling 4), and `produce` takes `label_pressure_Pa` to use it. On the file's own levels that read
returns bit-identical values to the level-by-level one, which is why closure mode is untouched.

**Deliverable 2, altitude and datum.** `datum_geopotential` places the datum log-linearly in
pressure, exactly on a level when the datum is one, and refuses one outside the produced range.
`produce_at_target` returns the production, the radius, the height along the local vertical and the
altitude above the datum, integrating the target column on the arrival levels themselves rather
than interpolating between the mesh's nodes.

**Deliverable 3, kind `profile` in transfer mode.** `_transfer_product` and `_transfer_groups` in
`forward/transfer.py`; ten transfer-mode variables added to the schema as optional. The product
carries `anchors/<slug>`, `reference_surface`, `isobars`, `estimate` and `transfer_record` beside
the inputs and the namelist.

**Deliverable 4, the driver and the figures.** `forward.transfer.run(namelist)` does the whole
chain: load, the hook, the gauge latitude, the outer loop, every anchor carried to the gauge, the
estimate, `C` carried to the target, the production, the datum, the product, the figures.
`casspian-forward` dispatches on the mode. F5 gains its across-latitude panel, F6 draws the
pressure identity in transfer mode, F7 is new.

## 2. Decisions taken, for a ruling

1. **`anchors/<slug>` stays the kind N file verbatim, and this run's per-anchor arrays go in an
   `anchors/<slug>/transfer` subgroup.** Section 7 lists both under the one group. Splitting them
   keeps "the kind N file verbatim" literally true: what the anchor carries is the anchor's, what
   the transfer made of it is the transfer's. The subgroup holds `C_i`, the two uncertainty columns
   with `season_term`, the uncertainty and weight on the union levels, the presence mask, the
   arrival geopotential, and as attributes the role and weight, `P_i`, the reference-surface
   residual, the seasons and the hook's record.

2. **`pressure_identity_residual` is `derived`, not `modeled`, and so carries no uncertainty
   companion.** It is arithmetic on two variables already in the file. Every `modeled` variable
   does carry a NaN companion, scalars included, which is what the schema requires.

3. **A profile carries exactly one record group**, `production_record` in closure mode or
   `transfer_record` in transfer mode. `groups_required` no longer names `production_record` and
   the profile check enforces the choice, naming the fault when neither or both are present.

4. **The production namelist now carries `geopotential_spacing_m2s2 = 5.0e4`** (decision R), and
   `reports/step04_4/accept_step04_4.py` is **pinned** to the 0.05 degrees by 5,000 it was accepted
   at, with an assertion on the latitude spacing so the pin cannot drift silently. Only that suite
   read the key. Without the pin the namelist change would have re-run an accepted suite at ten
   times the spacing its bounds were set for, and broken the regression for every later step.

5. **F5's panel draws the reference surface as the line `Phi = 0`**, which is what it is in the
   `(phi, Phi)` plane, with the radius range stated on the axis; the varying radius is in the
   `reference_surface` group. **F7's fourth panel, `D_ij`, appears only at M >= 2** and is left out
   at M = 1 rather than drawn empty.

## 3. What is measured

At the namelist's 0.05 degrees by decision R's 50,000 m2/s2. A run costs about 60 s.

| | measured | stated | bound |
|---|---|---|---|
| pressure identity, 10 N | 5.777e-07 | 4.1e-7 | must fall under halving |
| altitudes, 10 N | 411,135 / 98,187 / -15,290 m | 411,132 / 98,186 / -15,290 | 5 m |
| T below the anchor's, 10 N | 1.1084e-02 / 1.0943e-02 / 1.0892e-02 | 1.103e-2 / 1.089e-2 / 1.084e-2 | 1e-4 |
| r0(10 deg) | 60,128,613.0 m | 60,128,613.0 | 1 m |
| at phi_c: N, Phi | 1.776e-15, 0.000e+00 | the closure product | 1e-12 |
| at phi_c: p, T | 5.686e-07, 7.441e-07 relative | the closure product | 1e-5 |
| at phi_c: altitudes | 376,783 / 90,067 / -14,031 m | 376,780 / 90,067 / -14,031 | 5 m |

**The closure production is bit-identical to the registered product** after the refactor:
pressure, temperature, number density, mean refractivity, mean molar mass and geopotential all
compare equal by `array_equal`. That is the safety property of the step and it should be checked
again after any change to `forward/production.py`.

Two readings worth recording. The temperature comparison is against **the closure product's**
temperature, not the anchor's tabulated column: the closure residual is 3.3e-3 in both `p` and `T`,
far outside the 1e-4 bound, and A40 is about the model's own gradients. And the `p` and `T`
agreement at `phi_c`, 6e-7 rather than round-off, is the composition regridded onto the produced
labels, which moves `ln R_bar` by 1.313e-06 against the 1.5e-6 section 10 finding 5 states.

## 4. Defects found and fixed, all the coding agent's

1. **The gauge geopotential was passed to the production where the arrival at the target
   belonged.** The pressure identity read 1.100e-02 against a stated 4.1e-7, which is exactly the
   transfer's own `d ln N`: the isobar shift is what keeps the hydrostatic integral returning each
   label, and dropping it misses by the whole of the refractivity change. The parameter is now
   `arrival_geopotential_m2s2` with the mechanism in the docstring. Decision H's internal check
   found this within a minute of first running.
2. **The temperature was compared against the anchor's tabulated column** rather than the closure
   product's. A misreading of section 1, not a code fault.
3. **Three schema omissions**, each caught by the schema: the record group requirement, the
   companions for `radius_m` and the height, which are modeled at the target and copied in closure
   mode, and the companions for `refractivity_gauge`, `geopotential_gauge_m2s2` and the five
   modeled scalars.
4. **F5's panel drew the reference surface by multiplying a radius by zero.**

## 5. What is not done, and what has had no test at all

The acceptance script `reports/step04_5/accept_step04_5.py` does not exist. Runs 2 (60 N), 3 (the
cylinder wind) and 5 (Step 4's M = 2 run carried through production) are not measured. The
regression has not been run since the step's changes.

**Unexercised code**, which matters for how much weight section 3's table can carry: the product
with two anchors, which is the only path that writes more than one `anchors/<slug>` group and the
only one where `D_ij`, the reduced chi-square and F7's fourth panel are not empty; the datum
refusal; the sheared run; and `casspian-plots` rendering the figures by hand rather than through
the driver. The measured values above come from two runs of five.

The values in section 3 are **targets for the acceptance script to reproduce, not results to be
carried into a filing**. If the script's own run of 10 N does not give 5.777e-07 and 411,135 m,
something changed between now and then and that is worth knowing.

## 6. What a ruling would most help with

Whether the five decisions of section 2 stand, before the acceptance script is written against
them. In particular decision 1, the per-anchor split, which shapes every `anchors/<slug>` group the
product will ever carry, and decision 4, the pin, which is the coding agent's own guard on an
accepted step rather than anything the specification asked for.
