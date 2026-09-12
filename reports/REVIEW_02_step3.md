# REVIEW 02, Step 3. Number density, absolute radius and refractivity

Review of `reports/REPORT_02_step3.md` against SPEC_02 v0.5 Step 3. Reviewer: S. Rafkin,
12 September 2026. **Disposition: accepted with changes.** All six findings are agreed. They
touch three specifications: SPEC_00 goes to v0.13, SPEC_01 (closed) takes a Step 8 amendment
at v0.19, SPEC_02 goes to v0.6. Order of work is at the end.

## Rulings on the findings

1. **`dr0/dr_anchor` = 1.217, not near unity.** Agreed, and the explanation is right: raising
   the anchor radius raises the oblateness as well as the scale. The spec no longer says "near
   unity".

2. **Total derivative, 1.237.** Agreed. `phi_c` is not an independent input; it is a function
   of the anchor radius through `psi`, so the first-order propagation the spec intends is the
   total derivative. Implement it as the chain the report already computed (latitude-held
   partial plus the coupling through `dphi_c/dr_anchor` and the B3 slope), confirmed by the
   full fixed-point rerun. The anchor term becomes 12.37 km.

3. **`dphi_c/dphi_g` = 0.9456 on the wind geoid.** Agreed; the spec's 0.923 was an ellipse
   estimate and 0.934 was the no-wind value. The measured wind-geoid value is the one to use
   and the spec now records all three so nobody chases the difference again.

4. **Kind C declares its closure.** Agreed. SPEC_00 §6.2 v0.13 adds `closure_rule` and
   `closure_species`, written by the composition tool; `refrac` reads them and checks the
   values against the declaration instead of inferring the structure. The inference code stays
   as the check.

5. **Per-molecule uncertainties NaN, not 0.0.** Agreed; that is the v0.7 rule and the Step 8
   tool missed it. Findings 4 and 5 together are a small amendment to the closed SPEC_01
   (v0.19) with a rebuild of `lindal_composition.nc` only; nothing else under
   `occul_data/lindal/` depends on that file until kind N exists, which is why now is the
   cheap moment.

6. **Combining a `range` with a `1sigma`.** The author's rule, now SPEC_00 §5 v0.13: each
   declared kind is converted to one standard deviation before quadrature under a declared rule
   (`range` is the half-width of a uniform distribution, divided by √3; `stated` is taken as
   `1sigma`), the result is `1sigma`, and the conversions applied are listed in an attribute on
   the derived companion. The label term becomes 0.115° × 0.9456 × 5,630 km/rad = 10.73 km,
   the anchor term 12.37 km, and `radius_uncertainty_m` 16.38 km, with `phi_c` carrying 0.109°.
   Both readings of the swath are defensible; the uniform one is the honest default, and the
   transcription can declare otherwise if the source's words ever justify it.

## On the decisions

All six stand. Decision 2's `mean_refractivity` and `mean_molar_mass` companions are now in
the kind N schema (SPEC_00 §6.7 v0.13), so Step 4 writes them. Decision 4 (no height term at
the anchor level, where `h − h_ref` is identically zero) is correct and is now in the spec.

## Order of work

1. SPEC_01 Step 8 amendment (v0.19): the composition tool writes the closure attributes and NaN
   per-molecule uncertainties; rebuild `lindal_composition.nc`; rerun the Step 8 acceptance plus
   the two new checks; a short addendum to REPORT_01_step8; commit; record the new hash.
2. Step 3 changes (SPEC_02 v0.6): total derivative; kind conversion before quadrature with the
   conversions attribute; read the closure declaration; rerun; commit.
3. Step 4.
