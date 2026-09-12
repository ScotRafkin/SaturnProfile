# REPORT 02, Step 3. Number density, absolute radius and refractivity (B3.1, B1, B3.3)

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.5, Step 3,
against `SPEC_00_Architecture_and_Data_Files.md` v0.12 and the closed
`SPEC_01_Lindal_Tool_Chain.md` v0.18 (commit `ece58d2`). Status: reported, awaiting review.
**Nine of nine acceptance checks pass.**

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product, so the directory sweep after the acceptance commit will find nothing to rebuild.

---

## 1. What was built

**`src/casspian/lib/reduction.py`**, pure functions, NumPy in and out:

- `number_density(p, T)`, B3.1, with `k_B` from `lib.constants` (CODATA 2018).
- `absolute_radius(h, h_ref, r0)`, B1. The difference is taken first, so the anchor level
  returns `r0` exactly.
- `mean_over_species(x, q)`, and `refractivity(n, R_bar)`, B3.3.
- `temperature_from_refractivity(p, R_bar, N)`, the inverse, for the closure check.
- `share_of_remainder_partials(x, q, share, partner)`: `d q_bar/d s` and `d q_bar/d x_k` for the
  species assigned first, with the closure applied.
- `quadrature(terms)`: the stated terms in quadrature, NaN terms left out, and the labels
  included and unstated returned with the total.

**`src/casspian/refrac/reduce.py`**, the orchestration:

- `reduce_profile(inputs, manifest, anchor)` returns a frozen `Reduction`. It holds every per
  level quantity, eight companions (each value with `included` and `unstated` label tuples),
  the partials, the terms, the closure it found, and the CODATA release.
- `anchor_partials(inputs, manifest, anchor)`: `dr0/dr_anchor` by central difference on the
  anchored march (plus and minus 1 km, latitude held). `dphi_c/dphi_g = 1/(1 + dpsi/dphi_c)`
  with `dpsi/dphi_c` by central difference along one march with `phi_c +- 1e-4 rad` as nodes.
  `dr0/dphi_c = r0 G_phi / g` from Eq. B3, with the march difference kept beside it.
  `dphi_c/dr_anchor` and the `psi` partials follow from `phi_c + psi = phi_g`. Three marches.
- `composition_closure(composition)`: identifies the share, partner and species assigned first
  from the file's own values (decision 1).

**`src/casspian/refrac/anchor.py`**, two changes to the accepted Step 2 code:

- `fixed_point_closure_rad` is added to the record, per REVIEW_02_step2 ruling 4. It measures
  -5.50e-10 rad (-3.15e-8 deg).
- The shared setup (constants, anchor radius, seed flattening, wind, march) moves into
  `geoid_setup`, so the partials march exactly the surface `freeze_anchor` marched.

The Step 2 suite passes 7 of 7 after the change.

## 2. Decisions

1. **The share of remainder structure is established from the composition file's values.**
   Kind C does not declare its closure in a machine readable form (`remainder_split` is prose).
   The share species is taken as the one with a stated `x_<species>_uncertainty` whose ratio to
   exactly one other species is constant at every level (to 1e-9). That species is the partner,
   and every other species is assigned first. For Lindal this gives H2, He and NH3. Any other
   shape is refused, because Step 3 specifies the constrained derivative for this declaration
   only. See finding 4.

2. **Term labels.**
   - `refractivity`: `composition`, `pressure`, `temperature`. The composition label covers the
     share term, any stated `x_k` of a species assigned first, and the per molecule terms
     `x_i dR_i`, all from the one file.
   - `radius_uncertainty_m`: `anchor_radius`, `label_latitude`, `height`, `anchor_height`.
   - The three scalars: `anchor_radius`, `label_latitude`.

   Inside the composition term the finer labels are kept on `mean_refractivity_uncertainty_m3`
   and `mean_molar_mass_uncertainty_kg_mol`. Those two companions are computed here because the
   spec asks for `m_bar`'s propagation too. Whether kind N writes them is a Step 4 question:
   SPEC_00 section 6.7 lists only the three.

3. **A label partly stated appears in both lists.** When a term is finite at some levels and
   NaN at others, `quadrature` lists it as included and as unstated. That is what the file can
   truthfully say with two attributes per variable. Lindal has no such term.

4. **At the anchor level, `h - h_ref` carries no height uncertainty.** `h` and `h_ref` are the
   same table entry there, so their difference is exactly zero. Adding the same `dh` twice in
   quadrature would propagate an uncertainty onto a quantity that is identically zero. Elsewhere
   the two are different entries and are added in quadrature as specified. Lindal's heights
   state no uncertainty, so no number here depends on it.

5. **`dr0/dr_anchor` is the latitude held partial, as specified.** Finding 2 gives the size of
   the coupling that leaves out.

6. **The label uncertainty is converted from degrees in `refrac.reduce`**, at the boundary, and
   the declared uncertainty kinds (`1sigma` for the anchor radius, `range` for the label) are
   recorded beside the partials.

## 3. Acceptance results

Run by `reports/step02_3/accept_step02_3.py`; full output in `reports/step02_3/output.txt`.

| Check | Measured |
|---|---|
| 1 bar | **Pass.** T 134.8 K; `n` = 5.3731235e25 m-3. x_NH3 = 1.09e-5, `R_bar` = 4.8362195e-30 m3. **`N` = 2.5985605e-4** (2.59856e-4). |
| 794.33 mbar | **Pass.** x_NH3 exactly 0. `N / n` = 4.836272215504724e-30 m3, identical to 0.94 R_H2 + 0.06 R_He from the species group (difference 0). |
| Top and bottom | **Pass.** Top: 20 Pa, 138.7 K, `n` = 1.044408e22. Bottom: 129,848 Pa, 146.2 K, `n` = 6.432868e25, `R_bar` = dry x (1 - 7.9253e-5), `N` = 3.1108633e-4. |
| Radius at anchor and 1 bar | **Pass**, both exact. `radius_m[29]` = 58,519,883.28568849 m = `r0`. At 1 bar `h` = 0 and `radius_m` = `r0` - 90,000 m. |
| Fractional uncertainty of `N` | **Pass.** 0.0233184547 at all 66 levels, spread 1.0e-17. It matches the closed form `0.03 (1 - x_NH3)(R_H2 - R_He)/R_bar` to 6.9e-18. Included `('composition',)`, unstated `('pressure', 'temperature')`. |
| Radius uncertainty | **Pass.** **22,212.692 m** at every level. Included `('anchor_radius', 'label_latitude')`, unstated `('height', 'anchor_height')`. Anchor term: `dr0/dr_anchor` = **1.216954** x 10 km = **12,169.5 m**. Label term: `dr0/dphi_c` = -5,629.722 km/rad (march difference -5,629.722) x `dphi_c/dphi_g` = **0.945601** x 0.2 deg = **18,582.4 m**. |
| Inverse closes | **Pass.** max relative departure 4.0e-16, at 100 Pa. |
| Beyond the specification: term lists and the NaN rule | **Pass.** Every companion listed with its terms. `number_density_uncertainty_m3` has no stated term and is NaN. A unit case confirms partly stated terms appear in both lists. |
| Beyond the specification: closure refusal | **Pass.** A composition with `x_H2_uncertainty` dropped is refused, naming zero share pairs. |

Scalar companions: `r0` 22,212.7 m, `phi_c` 0.189132 deg, `psi` 0.011076 deg.

**Both anchor partials were cross checked by rerunning the whole fixed point** on in memory
copies of the inputs.
- Label moved to 36.3 +- 0.05 deg: `dphi_c/dphi_g` = 0.945601, identical to the propagation
  to six places. `dr0/dphi_g` = -5,323.467 km/rad, against -5,323.469 by the chain.
- Anchor radius moved +- 1 km: the total `dr0/dr_anchor` is 1.237328 (finding 2).

## 4. Findings

1. **`dr0/dr_anchor` is 1.217, not "near unity".** SPEC_02 Step 3 expects it near unity, and
   the v0.3 draft said the Step 7 march showed unit sensitivity. The march says otherwise, and
   the value is physical. Raising the polar radius raises `r0` in proportion
   (`r0 / r_anchor` = 1.075) and also makes the planet more oblate, because the centrifugal to
   gravity ratio grows as `r^3`. A first order ellipse estimate gives 1.075 + 0.18 cos^2(phi_c),
   about 1.21. The anchor term is therefore 12.2 km, not 10.

2. **The latitude held partial omits a coupling of 1.7 percent.** Moving the anchor radius also
   moves `phi_c`, because `psi` depends on radius: `dphi_c/dr_anchor` = -3.62e-9 rad/m. Through
   the Eq. B3 slope that adds 0.0204 to `dr0/dr_anchor`. The full rerun gives the total 1.237328
   and the chain gives 1.216954 + 0.020374 = 1.237328. So the first order anchor term is
   12,373 m, not 12,170, and the radius uncertainty 22,322 m, not 22,213. Built as specified
   (latitude held). **Proposing** the total derivative, which is the first order propagation the
   spec intends; it is one term in `anchor_partials`, and both values are in the output.

3. **`dphi_c/dphi_g` is 0.9456, not 0.923.** The v0.5 figure 0.923 was my ellipse estimate in
   REPORT_02_step2. Three surfaces give three values:

   | surface | `dpsi/dphi_c` | `dphi_c/dphi_g` |
   |---|---|---|
   | ellipse through the fitted 100 mbar radii | | 0.923 |
   | no wind reference geoid, at its own `phi_c` | 0.07075 | **0.9339**, v0.4's "about 0.93" |
   | wind geoid, at the frozen `phi_c` | 0.05753 | **0.9456**, confirmed by the full fixed point rerun to six places |

   So v0.4's 0.93 was right for the no wind geoid, and v0.5's 0.923 replaced it with a cruder
   number. The wind moves the value from 0.934 to 0.946. The kind W column falls from 76 m/s at
   26 deg planetocentric to 2.2 m/s at 30.8 deg, and that gradient changes how the tilt varies
   with latitude near the anchor. The label term is 18.58 km, which is the spec's "about 18 km".

4. **For SPEC_00 section 6.2: kind C should declare its closure.** Decision 1 infers the share
   structure from values, which works for Lindal but is an inference. **Proposing** two
   attributes, for example `closure_rule = "share_of_remainder"` and
   `closure_species = "H2 He"`, written by the composition tool and read by `refrac`. The check
   on values would then confirm the declaration instead of standing in for it.

5. **The per molecule refractivity uncertainties are 0.0, not NaN.** The `lindal1985` species
   group writes `refractivity_uncertainty_m3 = [0, 0, 0]` with
   `uncertainty_method = "zero where the source states none"`. SPEC_00 v0.7 makes an unstated
   uncertainty NaN, never 0.0. As written they enter the composition term and contribute nothing,
   so no number changes, but the file cannot distinguish "stated zero" from "not stated". This
   belongs to the SPEC_01 Step 8 tool; a rebuild would change `lindal_composition.nc` and its
   hash. Not changed here.

6. **The two anchor terms combine a `1sigma` with a `range`.** Kind D states the anchor radius
   uncertainty as `1sigma`. Kind T states the label uncertainty (0.2 deg) as `range`. First order
   propagation adds them in quadrature as if both were one standard deviation. If the range is a
   half width of a uniform distribution the label term would be 18.58 / sqrt(3) = 10.7 km. Both
   kinds are recorded beside the partials; the rule for combining them is the author's.

## 5. Regression

The ten SPEC_01 suites and the SPEC_02 Step 1 and Step 2 suites were rerun after the
`anchor.py` change (`reports/step02_3/regression.txt`): all pass, 6, 14, 8, 5, 7, 7, 6, 8, 7 and
6 for SPEC_01, 6 of 6 for Step 1 and 7 of 7 for Step 2.

No product under `occul_data/` was rebuilt, so all still carry `ece58d2`. No em dash or en dash
appears in any file written in this step.

## 6. Next step

Step 4, the kind N product and `casspian-refrac`, is not started. It waits for the review of this
report.

---

## 7. Addendum after review: the SPEC_02 v0.6 changes

Added 12 September 2026, per `REVIEW_02_step3.md` (accepted with changes). The order of work
there was followed. First came the SPEC_01 v0.19 Step 8 amendment: commits `6ee8690` and
`31c30ef`, recorded in REPORT_01_step8 section 7. `lindal_composition.nc` was rebuilt clean at
`6ee8690` with SHA-256 `a1250b60...`. Then came these changes.

**What changed in the code.**

1. **Total derivative** (finding 2). `anchor_partials` returns `dr0_dr_anchor` as the latitude
   held central difference plus the coupling `dr0/dphi_c x dphi_c/dr_anchor`, and keeps both
   parts (`dr0_dr_anchor_latitude_held`, `dr0_dr_anchor_latitude_coupling`).
2. **Kinds converted before quadrature** (finding 6, SPEC_00 section 5 v0.13).
   - `lib.reduction.to_one_sigma` applies the declared rule: `1sigma` unchanged, `stated` taken
     as `1sigma`, `range` divided by sqrt(3).
   - An unknown kind is refused. So is a finite uncertainty with no declared kind.
   - Every `Companion` carries `conversions`, the texts of the conversions applied to the terms
     that entered. Step 4 writes it as `uncertainty_kind_conversions`.
3. **Closure read from the declaration** (finding 4, SPEC_00 section 6.2 v0.13).
   - `composition_closure` reads `closure_rule` and `closure_species`.
   - It refuses a missing or unknown rule, a rule other than `share_of_remainder`, and a species
     pair not among the file's species.
   - It then runs the v0.5 inference as the check, refusing if the values do not show exactly
     the declared pair.
4. **The mean refractivity and molar mass companions** are unchanged in form. Kind N now
   carries them (SPEC_00 section 6.7 v0.13).

**Acceptance, rerun.** `reports/step02_3/accept_step02_3.py`, updated to v0.6, nine of nine.
Checks 1 to 4 and 7 are unchanged to the last digit. The changed checks:

| Check | Measured |
|---|---|
| 5. Fractional uncertainty of `N` | **Pass.** 0.0233184547 at all 66 levels, as before. Included `('composition',)`, unstated `('pressure', 'temperature')`. `uncertainty_kind_conversions = ('share: stated taken as 1sigma',)`. Inside the composition term only `share` entered: the per molecule uncertainties are NaN since the rebuild, and `x_NH3` states none. |
| 6. Radius uncertainty | **Pass.** **16,376.812 m** at every level (16.38 km). Anchor term: total `dr0/dr_anchor` 1.237328 (held 1.216954 plus coupling 0.020374) x 10 km `1sigma` = **12,373.3 m** (12.37). Label term: -5,629.722 km/rad x 0.945601 x 0.115470 deg (0.2 deg `range` / sqrt(3)) = **10,728.6 m** (10.73). Conversions: the label range only. The full fixed point rerun gives total `dr0/dr_anchor` 1.237328, differing from the propagation by 1.3e-8. `phi_c` carries **0.109208 deg** (0.109), `psi` 0.006615 deg. |
| 8. Term lists | **Pass.** As in section 3, with the new values above. |
| 9. Closure refusals | **Pass**, four cases. `x_H2_uncertainty` dropped. `closure_rule` removed. `closure_species = "He H2"`: the values show `('H2', 'He')`. `closure_species = "H2 NH3"`: likewise. |

Findings 1 to 6 of section 4 are all resolved by the rulings; none is open. Regression after
these changes: see `reports/step02_3/regression.txt`, summarized in the commit.
