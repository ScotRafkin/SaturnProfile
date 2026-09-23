# SPEC 04. The transfer of refractivity in latitude

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.12, 23 September 2026. Author of record: S. Rafkin. **Status: accepted by the author
at v0.6 on 21 September 2026; v0.12 states the truncation floor of a gridded wind hypothesis
(decision Q) and restates the cylinder-wind checks of Steps 3, 4 and 5 as the line integral the
transfer applies, at that floor, §14; v0.11 corrected the Step 2 column checks (the central-field closed
form in place of a uniform gravity that the constants cannot produce, an absolute bound on the
cancellation) and completes decision P (per hemisphere, the inversion curve sampled at the mesh
spacing), §13; v0.10 fixed the construction of the cylinder-extended wind
(decision P) and moves it to Step 2, with the v0.9 values for it withdrawn and remeasured
(§12); v0.9 made kind C one structure, a field on (level, latitude),
for every use (decision O), and the mesh self-extending with no declared extent or margin;
v0.8 applied the author's direction of 21 September on
consistency checks (decision N) and the rulings of REVIEW_04_step0 (§11), Step 0 to be
completed under it (decision L confirmed; the manuscript draft of record for the
equation labels is `CASSPIAN_AtmosphericModel_Draft9_2.docx`, outside the repository). Step 0
proceeds once `reports/REPORT_04_preexecution.md` is filed; Steps 1 to 5 in order, each after
the review of the one before.** v0.6 ruled on the coding agent's pre-execution
review of v0.5 (§10): one interpolation rule for data on a latitude grid, linear between nodes,
the rule the reduction used, with the derivatives those of the interpolant (decision L), and
every expected value remeasured under it; the rebuild of kind W cascading through the chain at
Step 0; `produce` taking the wind along the column; the gauge latitude fixed before the mesh
(decision K); the synthetic anchor's thermo group (decision M); the composition as a field on
latitude and pressure; the validation role of `weight` honored. v0.5 added, at the author's direction of 17 September,
the propagation hook and the seasonal term of every anchor's uncertainty as a slot the code
carries from the start, at every M including M = 1, for the propagator specification to fill
(decision J); the estimate takes each anchor's uncertainty as it arrives rather than reading it
from the file. v0.4 applied the author's markup of 16 September on v0.3: the code is structured for M anchors from the first step, with the estimate of A29 to A34
implemented now and the kernel correction of A35 left to the combination specification; kind W
is data along the local vertical in three human-readable parts and the model reads the total
only, with no declarations or parameterizations in the model; kind C is data on a latitude grid;
the tags of v0.3 are withdrawn; the seasonal order (propagate, then transfer) is stated; every
expected value is labeled as an instance for the one anchor in hand. Depends on
`SPEC_00_Architecture_and_Data_Files.md` v0.20 and the closed `SPEC_03_Forward_Production.md`
v0.14, which it does not repeat, on `docs/CASSPIAN_Seasonal_Design_Note.md` v0.5, and on the
manuscript's Appendix A (A13 to A34, A39) and B (B3 to B7). See §9 for the revision history.

---

## 0. How this specification is to be worked

The protocol of SPEC_03 §0 applies unchanged: one step at a time; the status line and
`STATE.md` read before any step; acceptance scripts under `reports/step04_<N>/`; figures a
report cites under `reports/figures/`; `REPORT_04_step<N>.md` and `REVIEW_04_step<N>.md`;
commit on acceptance, then the sweep over `occul_data/` and `forward/`; the full regression
after any change to `lib`, `refrac`, `forward` or `tools`; radians in `lib`; no dashes; no
silent choices. The in-memory candidate rule for a step that changes an input file (SPEC_03 §0)
applies to Step 0, which rebuilds input files.

**What is checked (decision N, author, 21 September 2026).** The pipeline assumes its inputs
came from the prior steps. The code refuses only what would silently produce a wrong number or
what the schema catches for free: a file of the wrong kind or missing a required variable
(the schema); an unknown, missing or malformed namelist key, or a scheme not implemented (a
typo guard); a point used outside the coverage of the data it needs (an interpolant refuses to
extrapolate); a numerical failure the step names (a curve leaving the mesh, curves crossing,
the outer loop not converging). Everything else is recorded in the product, with a warning
where the record alone might be missed: a season that differs from the run's, a lineage hash
that differs from the one recorded, a role or scale combination. No refusal exists to police
what a user might do, and no acceptance test exercises a refusal that this paragraph does not
name. This rule overrides "else refused" wherever earlier text says it.

**What this specification covers and what it does not.** The transfer half of the forward
model, Appendix B4 to B6 and the parts of B7 the closure did not need: the reference surface at
every latitude, the wind on the model's own geometry, the working mesh in latitude and
geopotential, the shear kernel and the transfer kernel, the isobar tracing, the transfer of each
anchor's refractivity along the isobars to the gauge latitude, the estimate of the anchor
constant from M anchors (A29, A30, A33, A34; an identity at M = 1, exercised at M = 2 in the
acceptance), the transfer from the gauge to the target, the production there, the delivered
altitude and datum, and the `transfer` mode of `casspian-forward`. **Nothing in the code assumes
one anchor or names a profile; the anchor list is the primary object from Step 0.** What is left
to later specifications: the kernel correction η of A35 and the posterior wind (the combination
specification, next after this one, which slots into Step 4's estimate without restructuring);
the seasonal propagator (the specification that fills the hook and the seasonal uncertainty
term this specification carries empty, decision J); the retrieval instance of kind N with
geopotential as its coordinate
(the retrieval leg); the forward-inputs tools (the composition generator, the wind tool for
parameterized shear, the diagnostic tool that decomposes a kind W file on a given geometry into
its barotropic and baroclinic parts), which follow SPEC_05 so that the transfer can then be
exercised under other assumed states; the end-to-end tests (SPEC_05); the Monte Carlo wrapper.

**The state the steps are coded against.** Every step codes against the closure inputs, the
one state in hand: the swept anchor `occul_data/lindal/lindal_refractivity.nc` and the inputs of
`forward/lindal_closure/`, rebuilt at Step 0 in the forms this specification fixes. Two wind
files are used throughout: the closure wind, which is the cloud-level wind at every level (what
Lindal assumed), under which the transfer is a small baroclinic one measured independently; and
the same reference-level wind extended along cylinders by the acceptance script with the
library's geometry (decision P: pinned on the file's own reference level, built at Step 2),
under which the kernel vanishes and the transfer is an identity. The general
paths, a wind with axial shear, a composition varying with latitude, and M = 2, are exercised
by synthetic files built inside the acceptance scripts, as the linear-T column was for SPEC_03.
Every expected value below is an instance for the Lindal anchor under these inputs and is
labeled as such; the code paths are generic.

**Expected values.** Those for the closure inputs were measured by the reviewing agent with
code independent of the repository (§1 and the steps) and are to be reproduced. Those for the
synthetic sheared and latitude-varying states will be measured before the step that needs them
proceeds; until then those checks rest on the identities and convergence stated in each step.

---

## 1. The formulation as it is implemented, and what the closure inputs give

**The field and its coordinates.** The working coordinates are planetocentric latitude `φ` and
the effective geopotential `Φ`, with `Φ = 0` on the gauge isobar (100 mbar for the Lindal
runs). At each latitude the column is the radial line through the reference surface `r0(φ)`,
and `Φ` on it is Eq. A13, `Φ(r, φ) = ∫_{r0(φ)}^{r} g(r', φ) dr'` with `g` the radial component
of the effective gravity built with the run's wind (A4). The reference surface is the gauge
isobar, whose slope is Eq. B3, `g dr0/dφ = r0 G_φ`, marched with the run's gravity, rotation
and wind from a measured radius (§3). An occultation anchor's levels sit at `(φ_i, Φ_k)`, with
`Φ_k` the line integral of `|g_eff|` along the profile's own vertical (SPEC_03 Step 1) under the
run's fields, the latitude drift along that vertical neglected in `φ` and recorded (0.027° at
the top of the Lindal profile). The two constructions of `Φ`, radial on the column and along
the field line for an anchor, are values of one scalar field wherever the effective gravity is
conservative and differ by the shear's circulation around the drift loop otherwise, second
order in the drift and not carried. One consequence is stated once because it will be noticed:
the radius of the column at `(φ_i, Φ_k)` is not the anchor's `radius_m` at level `k`; the
anchor's level lies on the tilted field line, so the same equipotential meets the column at a
larger radius, by `tan²ψ (h_k − h_ref)` (2,668.5 m at the top of the Lindal profile, −955.1 m
at its bottom, measured). That is geometry, and it is why the delivered altitude is defined
along the local vertical (decision B).

**The wind (decision A).** Kind W is data: the zonal wind against latitude and pressure along
the local vertical, `u_total(φ, p)`, relative to a declared rotation system, in three parts a
person can read: the reference-level wind `u_reference(φ)` at `p_ref`, the total, and the shear
along the local vertical `u_shear = u_total − u_reference` (SPEC_00 §6.6 amended, Appendix).
How the file was built, from cloud tracking, a published table, a thermal-wind analysis, a
parameterized decay, a cylinder extension, or an invented field, is the wind tool's business,
recorded in the file's provenance as prose the model does not read. The model reads
`u_total`, checks that the three parts sum, and does the one thing it must: it places the
field on its own geometry. At a mesh point `(φ, Φ)` it needs `u` at the radius of that point,
and the file gives `u` at a pressure, so the model evaluates `u_total(φ, p(φ, Φ))` with
`p(φ, Φ)` from the isobar map. The map depends on the wind through the kernel, so the mapping
is a fixed point (the outer loop of SPEC_00 §7.2): start from flat isobars (`p(Φ)` from the
anchors' productions, the barotropic guess), build the columns and kernels, trace, form the new
map, repeat to a declared tolerance in `ln p`. A file whose columns do not vary with pressure
converges on the first pass; the code does not special-case it. From the field on the mesh the
axial derivative of A15 is formed by centered differences. Two facts about the one-anchor
state are worth having in mind: the closure wind, the cloud-level wind at every level, is not
the barotropic state (a wind that varies with latitude and not with height has an axial
derivative through the meridional term of A15, and by A39 a temperature gradient on isobars),
and the same reference-level wind extended along cylinders is, with the kernel vanishing. Both
are hypotheses; the tools make them; the Monte Carlo draws them.

**Seasons (decision I).** The run declares one season, and its W and C are the hypothesis for
that season, checked to carry the run's season or to be declared uniform in season. An anchor
at another season is moved to the run's season by the propagator before it is transferred, in
`N` at fixed pressure at its own latitude; only then is it transferred under the run's W and C.
The wind of the anchor's own season never enters the forward run. Until the propagator exists
(a later specification), a season mismatch between an anchor and the run is not refused
(design note decision 2); the product records the run's season beside each anchor's, and the
missing propagation is a recorded, unmodeled term. The one anchor in hand and its transfer run
are at the same season.

**The uncertainty of an anchor as it arrives (decision J).** Three terms carry uncertainty
into the estimate, and the code carries all three from the start. The measurement term
`σ_i^meas(Φ)` is the anchor's own refractivity uncertainty as reduced, which already holds
whatever composition and wind its source assumed, scaled by its `measurement_uncertainty_scale`.
The season term `σ_i^season(Φ)` is the propagator's, for the seasonal distance the anchor was
moved (design note §5.4: the fit's parameter covariance, the interannual floor, the residual of
the rule against the record); it is not a declared number and it is not zero at the anchor's
own season, and the propagator specification computes it from the record. The transfer term
is `P_i`, from the run's W and C through the kernel over the latitude distance. The anchor
arrives at the estimate with `σ_i²(Φ) = σ_i^meas² + σ_i^season²` on its levels and `P_i` is
added there. This specification implements the first and the placeholder of the third, and
carries the second as a column that is zero and recorded absent until the propagator exists;
the propagation acts at every M, so the slot exists at M = 1 and the M = 1 product carries
it.

**The kernels.** The shear kernel is Eq. A15, `S = 2 Ω_abs r (∂u/∂Z)_R` with
`(∂u/∂Z)_R = sin φ (∂u/∂r)_φ + (cos φ / r)(∂u/∂φ)_r`, with the derivatives of `u` those of
the wind file's interpolant in its own coordinates, `(φ, ln p)`, linear between nodes
(decision L), converted to fixed `r` with the isobar map's slopes on the mesh. `I(φ, Φ) = ∫_0^Φ (S / g) dΦ'` (A16) is the vertical integral from the reference
surface at each latitude. The transfer kernel is Eq. A27, `K = S / g + (∂ ln(ℛ̄ / m̄) / ∂φ)_p`,
the second term from the run's composition on the isobar, by centered differences in latitude
on the composition file's own grid, zero when its columns are identical.

**The tracing, the transfer, the estimate.** Each level `k` of anchor `i` is an isobar of
label pressure `p_k^(i)`, the anchor's own produced pressure at that level (SPEC_03 Step 4 at
the anchor's latitude). Its path in the `(φ, Φ)` plane is the characteristic
`dΦ/dφ = −I(φ, Φ)` (A24, B6.2); along it Eq. A28, `ln N = ln N_k^(i) + ∫ K dφ'`. The curves
cannot cross; a crossing is refused. Every anchor is traced to the gauge latitude `φ_r`, the
weighted centroid of the anchors' latitudes (A34; `φ_r = φ_1` at M = 1), where each supplies
its estimate `C_i(Φ)` of the anchor constant (A29). The estimate `C(Φ)` is the inverse-variance
mean of A30 on the union of the anchors' arrival levels at `φ_r`, each anchor's `ln N`
interpolated log-linearly in `Φ` onto the union within its own span, with weights
`w_i = 1 / (σ_i² + P_i)`: `σ_i` the anchor's own refractivity uncertainty scaled by its
`measurement_uncertainty_scale`, and `P_i = (σ_K (φ_i − φ_r))²` from a declared per-run kernel
uncertainty `σ_K`, the placeholder for A31 that the combination specification replaces with
the covariance of A35. `D_ij = C_i − C_j` (A33) on the levels both cover, and the reduced
chi-square of the `C_i` about `C` per level, which inflates the variance of `C` where it
exceeds one (the rule stated with A33). At M = 1 every one of these is the identity: `φ_r = φ_1`, the union
is the anchor's own levels, `C = ln N_1`, no interpolation, no `D`. From `C(Φ)` at `φ_r` the
isobars are traced to the target and `C` transferred along them (A28); the production of SPEC_03
on the arrival levels, with `p_b` the label of the topmost isobar of the highest anchor
(decision F), gives `p` and `T`. Because pressure is what the characteristics carry, the
produced `p` on each level must return its label: exactly for a single anchor in the
barotropic state, to the discretization error of the tracing and the kernel integral
otherwise, and at M ≥ 2 the labels of different anchors on the same level must agree in the
same sense. This **pressure identity** is the internal check of the whole chain and is
reported at every level of every run.

**What the closure inputs give (measured independently; the Lindal anchor; the reviewing
agent's mesh at 0.05°, the wind linear between the file's nodes with its derivative piecewise
constant as decision L rules, `I` bilinear, columns and curves integrated to 1e-12; v0.6, the
v0.5 values were measured with PCHIP and are withdrawn).**
Under the closure wind `S / g = 2 Ω_abs cos φ u'(φ) / g`, nearly uniform in the vertical, so
the transfer is close to a uniform scaling of `N` with a linearly growing isobar shift. From
`φ_c` = 30.805568° to 10° N: largest `|S/g|` 0.0977 per radian; isobar shift −31,258.1 m²/s²
at the top level, zero at the gauge, +11,339.7 at the bottom; `ln N` up by 1.1090e-2 (top),
1.0988e-2 (10 mbar), 1.0950e-2 (gauge), 1.0899e-2 (bottom); the produced temperature at the
target lower than the anchor's on every isobar by 1.103e-2 (1.53 K) at the top, 1.089e-2
(0.91 K) at the gauge, 1.084e-2 (1.58 K) at the bottom, A40's equal and opposite gradients;
the pressure identity to 4.1e-7. To 60° N: largest `|S/g|` 0.0540; shift −7,117.0 and
+2,585.1 m²/s²; `ln N` up by 2.5119e-3, 2.4922e-3, 2.4840e-3, 2.4749e-3; temperature down by
0.35, 0.21 and 0.36 K at the top, gauge and bottom. Under the cylinder-extended wind: `S ≡ 0`
analytically, isobars flat, `N`, `p`, `T` unchanged at every latitude (1e-15 in the reviewing
agent's computation, where the extension was analytic; the model's, from a file on a grid,
is bounded in Step 4). The geometry, common to both: the reference surface through
`(φ_c, 58,516,188.29 m)` has radius 60,367,000.0 m at the equator, 60,128,613.0 m at 10° N,
59,492,075.7 m at 20° N, 57,053,249.5 m at 45° N, 55,675,128.9 m at 60° N, 54,420,643.5 and
54,449,383.8 m at the poles (an independent march that reproduces the reduction's record).
Delivered altitudes above the 1 bar datum along the local vertical (decision B): at `φ_c`
376,780 m at the top level, 90,067 m at the gauge level, −14,031 m at the bottom, against
Lindal's 376,700, 90,000 and −14,100 m (the closure residual at 1 bar moves the datum by 68 m;
the drift geometry gives 15 m at the top); at 10° N under the closure wind 411,132, 98,186 and
−15,290 m and under the cylinder wind 416,300, 99,290 and −15,456 m (v0.10, decision P: the
cylinders through the 10° N column reach the 1 bar surface nearer the equator, where the wind
is stronger, 427 m/s at the top level against 329 m/s on the file's reference level at 10° N,
and that wind is in the column's gravity); at 60° N under the closure wind 328,127, 78,533 and
−12,239 m.

---

## 2. Step 0: kind W in three parts, kind C on a latitude grid, the transfer run directory, the namelist for M anchors, the loader

**Deliverable 1: kind W in three parts (SPEC_00 §6.6 amended).** `casspian-wind-from-curve`
writes `u_reference_ms(latitude)` at `reference_level_pressure_Pa`, `u_total_ms(latitude,
pressure)` and `u_shear_ms(latitude, pressure) = u_total − u_reference`, with
`u_total_uncertainty_ms`; the `decomposition`, `u_cylindrical_ms` and `decomposition_geometry`
of SPEC_00 v0.8 are retired from the file (the baroclinic decomposition on a given geometry is
a later diagnostic tool). The reader checks the sum identity and the poles (zero at ±90°) and
otherwise validates the schema; it reads no declaration of how the field was built. The
reduction's and the closure's wind files are rebuilt by the §0 procedure with no value of
`u_total` changing. **The rebuild cascades (§10, finding 2):** kind N records the wind file's
hash in `input_hashes` and `reduction_record` and the reader refuses a mismatch, so the whole
reduction chain, the closure run's inputs and the closure product are rebuilt as SPEC_03 Step 3
did (the in-memory candidate rule, then the sweep on the clean tree), with no value changing in
any product, only variables added, verified group by group; the new hashes are recorded in the
report in the row format the `step02_1` suite reads, and REPORT_03_step4's product hash is
superseded there.

**Deliverable 2: kind C is one structure, a field on (level, latitude) (decision O, author,
21 September 2026).** `casspian-composition-lindal` takes a required control-file key
`latitude_grid_deg = [south, north, step]` and writes every level variable on
`(level, latitude_planetocentric)`, the source's column at every node (the Lindal hypothesis,
uniform in latitude), `pressure_Pa` one dimensional; the one-latitude form and its
`latitude_planetocentric_absent_meaning` marker are retired from the kind (SPEC_00 §6.2 and
SPEC_01 Step 8 amended, Appendix). Every reader takes the composition as a field on
`(φ, ln p)`: the reduction and closure mode read the column at the anchor's latitude by the
interpolant of decision L (for a uniform field the identical column, so no value changes) and
use it on the anchor's levels by level as SPEC_03 does; transfer mode interpolates log-linearly
in pressure onto each isobar's label (§10, finding 5). The interpolant refuses a latitude
outside the file's grid. Whether a field is uniform is a property of its values, not a tag.
The reduction's composition control file and the closure's build file declare the grid
`[-90, 90, 1.0]`; the rebuild of deliverable 1 carries this change through the same cascade.

**Deliverable 3: `forward/lindal_transfer/`.** The namelist `lindal_transfer.toml` (below),
`lindal_transfer_build.toml` with the four sections under the prefix `lindal_transfer_` and
the composition on the grid `[-90, 90, 1.0]`, `inputs/` by `casspian-run-inputs`, `output/`.

**Deliverable 4: the namelist in transfer mode.** SPEC_00 §7.2 with these rules: `mode =
"transfer"`; `[[anchors]]` one or more entries, each a kind N file with its slug, `weight` (its role,
SPEC_00 §7.2: 1 construction, 0 validation; a validation anchor is propagated, placed and
traced to the gauge like any other and its `C_i` and every `D_ij` are reported, but it enters
neither `φ_r` nor `C`; any other value refused) and `measurement_uncertainty_scale`; `[target] latitude_planetocentric_deg`, one value; `[inputs]`
as SPEC_03; `[hydrostatic_boundary]` exactly one of `p_b_rule = "anchor_profile_top"` (the
label of the topmost isobar of the highest anchor) and `p_b_Pa`; `[isobars]` `gauge_isobar_Pa`
(a tabulated level of every occultation anchor, matched exactly) and `datum_isobar_Pa`;
`[grid]` `geopotential_spacing_m2s2`, `latitude_spacing_deg` (the mesh's extent is not
declared: Step 2 builds it from the anchors' levels and extends it when a traced curve needs
more, v0.9);
`[numerics]` with `scheme` restricted to what is implemented (`reference_surface`,
`isobar_tracing`: `"rk4"`; `shear_integral`, `transfer`, `altitude`: `"trapezoid"`;
`outer_loop`: `relative_tolerance_ln_p`, `max_iterations`; the geopotential and hydrostatic
rules are SPEC_03's closed ones and are not namelist keys, as in the closure namelist); `[estimation]` `gauge_latitude_rule =
"weighted_centroid"`, `kernel_uncertainty_per_rad` (`σ_K`, required, no default in code; 0.02
confirmed by the author, §8), `model_error_correlation_length_deg = 0.0` (accepted and recorded; A35 later);
`[output]`, `[diagnostics]`. An unknown key is an error; a §7.2 key not listed is refused as
not implemented in SPEC_03's message form. The closure-mode refusals are unchanged.

```toml
[run]
name        = "lindal_transfer"
description = "Lindal refractivity transferred to 10 N under the closure inputs"
mode        = "transfer"
solar_longitude_deg = 18.2
date        = "1981-08-26"

[[anchors]]
slug   = "lindal"
path   = "../../occul_data/lindal/lindal_refractivity.nc"
weight = 1.0
measurement_uncertainty_scale = 1.0

[inputs]
composition = "inputs/lindal_transfer_composition.nc"
gravity     = "inputs/lindal_transfer_gravity.nc"
rotation    = "inputs/lindal_transfer_rotation.nc"
wind        = "inputs/lindal_transfer_wind.nc"

[hydrostatic_boundary]
p_b_rule     = "anchor_profile_top"
p_b_location = "top_of_anchor_profile"

[isobars]
gauge_isobar_Pa = 1.0e4
datum_isobar_Pa = 1.0e5

[target]
latitude_planetocentric_deg = 10.0

[grid]
geopotential_spacing_m2s2 = 5.0e3
latitude_spacing_deg      = 0.05

[numerics.reference_surface]
scheme = "rk4"
[numerics.shear_integral]
scheme = "trapezoid"
[numerics.isobar_tracing]
scheme = "rk4"
[numerics.transfer]
scheme = "trapezoid"
[numerics.altitude]
scheme = "trapezoid"
[numerics.outer_loop]
relative_tolerance_ln_p = 1.0e-8
max_iterations          = 50

[estimation]
gauge_latitude_rule                = "weighted_centroid"
kernel_uncertainty_per_rad         = 0.02
model_error_correlation_length_deg = 0.0

[output]
directory = "output"
product   = "lindal_transfer_profile.nc"

[diagnostics]
figures = true
format  = "png"
dpi     = 150
```

**Deliverable 5: `load_run_inputs` in transfer mode.** As SPEC_03 (existence, `-dirty`, prefix,
role, the wind against the rotation, the sum identity and poles of W) for every anchor in the
list and the four inputs, with no closure comparison. Coverage is the interpolants' business
(Step 1 deliverable 2 and the composition field refuse a point outside their data), not the
loader's. The seasons of W, C and every anchor are recorded beside the run's, a difference
warned and recorded, never refused (decision N; this amends decision I's "else refused"). The
retrieval instance of kind N does not pass the schema today, so the loader carries no check of
its own for it; the retrieval leg adds what it needs.

**Deliverable 6: the anchor as it arrives.** The loader returns every anchor as one object
carrying, on the anchor's own levels, its `ln N`, its coordinate, its label pressures, and its
uncertainty in two named columns, `sigma_ln_N_measurement` (the file's `refractivity_uncertainty`
divided by its `refractivity` and multiplied by the scale: an uncertainty of `ln N`, which is
what A30's weights act on; REPORT_04_step0 finding 1) and `sigma_ln_N_season` (zero, with the attribute `season_term = "absent"`),
beside the anchor's season and the run's. Every later step reads the anchor's uncertainty from
this object and never from the file.

**Deliverable 7: the propagation hook.** `forward.propagate(anchor, season)` is the one place
an anchor is moved in season, called for every anchor between loading and Step 1's placement.
In this specification it is the identity: it returns the anchor with its `N` unchanged, its
season term zero, and a record of the anchor's season, the run's, and `propagation = "none:
propagator not implemented"`. The propagator specification replaces its body ((S8), (S9),
(S9b) and §5.4 of the design note) and nothing that calls it. A run whose anchors' seasons all
equal the run's passes through the hook the same way, so the code path is the same at every M
and every season.

**Acceptance.** The hook called for every anchor the loader returns, its record on the
object, the object's two uncertainty columns present with the season column zero and
`season_term = "absent"`; the three wind files (reduction, closure, transfer) rebuilt with the
three parts, `u_total` unchanged to the bit, the sum identity holding, a copy with the identity
broken refused; the transfer composition on 181 latitude nodes with every column equal to the
point file's, a copy with one node edited refused nowhere (it is data) but showing a nonzero
composition term at Step 3; `casspian-run-inputs` writes the four inputs; the transfer
namelist loads and resolves; refusal cases, each message quoted (decision N: these and no
others): a missing `[target]`, an unknown `[numerics]` table, a `scheme` not implemented, a
`weight` other than 0 or 1; recorded and warned, not refused: a wind or composition at another
season (a copy with the attribute edited); a namelist with two `[[anchors]]` entries pointing at the
same file loads (the loader does not refuse M > 1); the closure namelist still loads and its
thirteen refusals hold; the rebuilt files differ from their swept copies only in the new
variables and the writer globals.

---

## 3. Step 1: the reference surface through the anchors, and the wind on the model's geometry

**Deliverable 1: `lib.geoid` anchor rule `through_anchor`.** The march of Eq. B3 (SPEC_01
Step 5, unchanged in method) with the constant fixed by a stated point `(φ_i, r0_i)`, starting
there and running to both poles, on the caller's nodes; the wind supplied to the march is
`u_total` at the gauge isobar's pressure. With several occultation anchors the surface is
marched from the first listed, every other anchor's measured radius is compared with the
surface at its latitude, and the differences are recorded as `reference_surface_residual_<slug>`
(the weighted estimate of the constant from several measured radii belongs to the combination
specification; retrieval anchors, when they exist, contribute no radius).

**Deliverable 2: `lib.windfield`.** `wind_at(phi, p)`: `u_total` interpolated linearly in
latitude between the file's nodes and linearly in `ln p` (decision L: the rule the reduction
used in `refrac.anchor.wind_of_latitude`, which is what registered the anchor; the wind tool's
fit is its own business and is sampled on the file's grid), refusing a point outside the file's
coverage; `wind_derivatives(phi, p)`: `(∂u/∂φ)_p` and `(∂u/∂ln p)_φ` of the same interpolant,
piecewise constant, the value at a node being that of the interval to its north or to higher
pressure; `wind_on_mesh(mesh, pressure_map)`: `u` at every mesh node from
its latitude and the pressure the map gives it. No declaration is read and no field is
extended: what the file does not cover, the run refuses.

**Deliverable 3: the anchors' geopotential under the run's wind.** SPEC_03's `produce`
(`produce(profile, inputs, gauge, p_b)`, which today forms one scalar wind at the anchor's
latitude) gains an optional wind on the anchor's levels as an array, the scalar it forms today
being the default so that closure mode is unchanged; the product records the column as
`u_column_ms(level)` beside `u_at_phi_c_ms`, which keeps its meaning (the reference-level
value); the amendment to SPEC_03 is listed in the Appendix; the transfer supplies `wind_at(φ_i, p_k^(i))`
with `p_k^(i)` the anchor's tabulated pressures on the first pass and its produced pressures
after. For the closure wind this is 2.167 m/s at every level of the Lindal anchor, which is what
SPEC_03 used, so the closure product is unchanged. For the cylinder-extended wind, built at
Step 2 by decision P and exercised there, it runs from 2.168 m/s at the 1 bar level (the file's
reference level, where the extension is pinned) to 3.717 m/s at the gauge level, 9.490 m/s at
the top level and 1.926 m/s at the bottom, moving `Φ` at the top by −322.3 m²/s² and at the
bottom by +18.0 (an instance, measured by the reviewing agent at v0.10; the v0.9 values 9.371,
−0.133 and −271.6 were measured under a construction the reviewing agent could not identify
afterward and are withdrawn, §12).

**Expected values (instance: the Lindal anchor, the closure inputs).** The surface through
`(30.805568°, 58,516,188.29 m)`: 60,367,000.0 m at the equator (the reduction anchored there;
residual 1e-4 m), 54,420,643.5 and 54,449,383.8 m at the poles (the reduction's record),
60,128,613.0 m at 10° N, 59,492,075.7 m at 20° N, 57,053,249.5 m at 45° N, 55,675,128.9 m at
60° N. The wind along the anchor's column under each wind file as above.

**Acceptance.** The surface reproduces the equatorial radius to 0.01 m, the polar radii and
the asymmetry (28,740.37 m) to 0.1 m, the values at 10, 20, 45 and 60° N to 1 m, and agrees
with the reduction's own march (`equatorial_radius` rule) at every dense node to 0.01 m; a
second anchor entry pointing at the same file gives a residual of zero at its latitude; a
synthetic anchor entry with its radius moved by 1 km gives a residual of 1 km; `wind_at`
reproduces the file's nodes exactly and refuses a point outside coverage; the wind along the
anchor's column under both wind files matches the values above to 1e-3 m/s; the SPEC_03
closure rerun through `produce` with the wind array is bit-identical to the closure product.

---

## 4. Step 2: the working mesh and the columns

**Deliverable: `lib.mesh`.** Latitude nodes uniform at the declared spacing over the interval
from the southernmost to the northernmost of the anchors, the gauge latitude and the target,
with one spacing of margin each side and every anchor latitude, the gauge latitude (known
before the mesh, decision K) and the target inserted as exact nodes; any centered difference
taken on the mesh in latitude uses the three-point formula for unequal spacing at the inserted
nodes; geopotential nodes uniform at the declared spacing from the lowest to the highest
`Φ_k` of any anchor (Step 1 deliverable 3), one spacing beyond each, with `Φ = 0` a node.
Nothing about the extent is declared or checked (v0.9): the mesh is built from the data, and
when Step 4's tracing finds a curve reaching the mesh's edge the mesh is extended on that side
by what the curve needs, rounded up to whole spacings, the columns and kernels tabulated on the
new nodes, and the pass repeated; the extensions are recorded in `transfer_record`. The mesh
`extend(side, amount)` is this deliverable's. At each latitude node the column: `r(φ, Φ_j)` by
`dr/dΦ = 1 / g(r, φ, u)` from `(r0(φ), 0)` in both directions (RK4 on the `Φ` nodes, `u` from
`wind_on_mesh`), `z_lv(φ, Φ_j)` by `dz/dΦ = 1 / |g_eff|` on the same integration (decision B),
and on the nodes `g`, `G_φ`, `|g_eff|`, `ψ`, `u`. The columns are rebuilt on every pass of the
outer loop.

**Deliverable: the cylinder-extended wind file (decision P).** The acceptance script builds
the second wind file of §1 here, where the columns exist: `u_total(φ, p) = U(s)` with
`s = r cos φ` the distance from the rotation axis and `U` fixed on the file's own reference
level, `U(s_ref(φ)) = u_reference(φ)` with `s_ref(φ) = r(φ, Φ_ref) cos φ` the 1 bar isobar's
distance from the axis at each latitude, so that the file's three parts are the closure
file's `u_reference`, this `u_total`, and their difference, and the sum identity holds. The
geometry is the library's: `r(φ, Φ)` the radial columns of this step from the reference
surface of Step 1, under the flat-isobar map `Φ(p)` of the anchor (interpolated in `ln p`),
which is exact under `S = 0`; the map is continued beyond the anchor's levels at the slope of
its last interval in `ln p`. `U` is one function per hemisphere (v0.11, §13): a cylinder cuts
the reference surface once in each hemisphere and the wind is not symmetric about the equator,
so a single `U(s)` cannot return `u_reference` in both; within a hemisphere `du/dz` on a
cylinder is still zero, so `S = 0` holds wherever a run's span lies in one hemisphere, and a
symmetric deep wind is a different hypothesis for a tool to make. The inversion curve
`s_ref(φ)` is sampled at the mesh's latitude spacing over the file's whole latitude range, not
on the file's 0.5° grid (v0.11: the chord of a 0.5° cell misplaces a mid-cell latitude by
several thousandths of a degree, which on the steep flank of a jet is several hundredths of a
metre per second, the size of the departures REPORT_04_step2 measured). Because the reference surface is marched with the wind at the
gauge isobar, and that wind is the file's own, the construction is a fixed point: the closure
file's surface and the closure anchor's `Φ_k` on the first pass, the file's own on the next,
repeated until the file changes by less than 1e-6 m/s (three passes on the state in hand). The
file is written on the closure wind's grid with the writer, carries the run's season, and is
read back through `wind_at` like any other. Expected values (instance, the reviewing agent's
radial columns, the reference surface marched under the closure wind, a difference below the
tolerances): on the Lindal anchor's column `u` is 9.490 m/s at the top level, 3.717 at the
gauge level, 2.168 at the level nearest 1 bar and 1.926 at the bottom; the anchor's `Φ`
under this wind is lower at the top by 322.3 m²/s² and higher at the bottom by 18.0 than
under the closure wind (to 1e-2 m/s and 1 m²/s²); at 10° N the column carries 427.09 m/s at
the top level, 345.98 at the gauge level and 326.40 at the bottom against 329.46 on the file's
reference level (to 0.05 m/s).

**The staggering (SPEC_00 §7.3 rule ii).** Every quantity lives on the `(φ_i, Φ_j)` nodes
except: `S` and `S/g`, on the nodes, from the interpolant's derivatives (Step 1 deliverable
2) and the isobar map's slopes by centered differences (one-sided at the mesh edges, which the
margin keeps away from any traced curve); `I`, on the
nodes, the trapezoid of `S/g` in `Φ` from the `Φ = 0` node outward; the composition term of
`K`, on the isobar labels and the composition file's latitude grid, interpolated to the curve;
the curves `Φ_k(φ)`, on the latitude nodes, with `I` bilinear between nodes.

**Expected values (instance).** At `φ_c`, on the Lindal anchor's `Φ_k`: `z_lv` 286,715.2 m at
the top level and −104,098.0 m at the bottom (against the tabulated `h − h_ref` of 286,700 and
−104,100: the 15 m and 2 m are the second-order drift terms); `r − r0` 288,051.3 m at the top
and −104,576.8 m at the bottom, `tan²ψ` more than the projected altitude, as §1 says.

**Acceptance.** The central field (`J = 0`, `Ω = 0`, `u = 0`, which leaves `GM / r²`, not a
uniform gravity; v0.11, §13) returns `r = r0 / (1 − r0 Φ / GM)` to 1e-12 relative and
`z_lv = r − r0` to 1e-7 m absolute (the difference cancels seven digits of `r`, so a relative
bound below 1e-10 is below double precision); the column at `φ_c` reproduces the values above to 0.1 m;
`Φ` recomputed from the column's own `r` by the trapezoid of `g` returns the nodes to 1e-9
relative; the mesh refuses a spacing that is not positive; `extend` adds whole spacings on the named side and the columns on the new nodes match a mesh built with them from the start to 1e-12; a namelist with two anchors at different latitudes gives a mesh with both as exact
nodes; every quantity on the mesh is finite.

---

## 5. Step 3: the kernels

**Deliverable: `lib.kernel`.** `shear_kernel(mesh)`: `S` on the nodes by A15 from the
interpolant's derivatives at the node's `(φ, p)` (Step 1 deliverable 2, decision L), converted
to fixed `r`: `(∂u/∂φ)_r = (∂u/∂φ)_p + (∂u/∂ln p)_φ (∂ln p/∂φ)_r` and `(∂u/∂r)_φ =
(∂u/∂ln p)_φ (∂ln p/∂r)_φ`, the slopes of the isobar map by centered differences on the mesh
(the pressure map is smooth); nothing differentiates `u` on the mesh; `Ω_abs` from `u` at the
node. `shear_integral(mesh, S)`: `I` on the nodes. `composition_term(...)`:
`(∂ ln(ℛ̄ / m̄) / ∂φ)_p` on the composition file's latitude grid at each isobar label, centered
differences, zero where the columns are identical. `transfer_kernel(...)`: `K` at a point of a
curve from `S/g` (bilinear on the mesh) and the composition term (linear in latitude).

**Expected values (instance).** Closure wind: `(∂u/∂Z)_R = (cos φ / r) u'(φ)` exactly, since
`u` does not vary along the column, so `S / g = 2 Ω_abs cos φ u'(φ) / g` at every node with
`u'` the derivative of the file's interpolant; the largest `|S/g|` on the mesh from `φ_c` to
10° N is 0.0977 per radian and to 60° N 0.0540 (the reviewing agent's values on a 0.05°
mesh, `u'` the file's piecewise-constant slopes: 12.07 m/s per degree southward of the node at
31.0°, so the largest `|S/g|` sits on the interval holding the anchor); the model returns
`S/g = 2 Ω_abs cos φ u'/g` at every node to round-off, since no difference is taken (the v0.5
convergence claim for a centered difference is withdrawn, §10 finding 6). Cylinder-extended wind (v0.12, decision Q): `S ≡ 0` analytically (`u = f(R) r cos φ`, for
which the two terms of `(∂u/∂Z)_R` cancel for any `f`), and on the mesh the residual is the
truncation of the file's linear interpolant, whose derivatives jump at the file's nodes, so
the cancellation fails inside the cells holding a jump and nowhere else: `S/g` is a spiky,
sign-alternating field there whose largest value does not fall under refinement while its
integral does. What is measured is therefore the line integral of `S/g` along each of the
anchor's isobars over the run's span (flat isobars, the first-order estimate of `Δ ln N`),
which must be below 2e-3 at every level, and the vertical integral `I` and the isobar shift
`−∫ I dφ` it implies, below 300 m²/s²; the largest `|S/g|` is reported, not bounded. Expected
(instance, the reviewing agent's kernel on the coding agent's Step 2 file, the map continued
per Step 2 decision 8, a 0.05° mesh): `Δ ln N` +1.2e-4 at the top level, −5.2e-4 at 10 mbar,
−1.7e-4 at the gauge, −1.1e-3 at the bottom, largest 1.1e-3; largest `|S/g|` 0.050 per
radian on the span; shift +107 at the top, −10 at the bottom, largest 224 m²/s². The line
integral is remeasured with the file's latitude sampling at 0.25° and 0.1° (the cylinder
construction resampled, the mesh unchanged) and the ratios reported, no order claimed. The
map is continued at the slope of its last interval beyond the anchor's levels, as Step 2
decision 8 states; a clamped map puts a slope discontinuity at the ends and a spurious
residual an order of magnitude larger there (§14). Synthetic: `u = β r sin φ` applied at
every node (β = 1e-6 per second) has `(∂u/∂Z)_R = β` exactly, so `S = 2 Ω_abs r β` to the
truncation error of the interpolant, reported; a solid-body offset `u = ΔΩ r cos φ` has
`S = 0` analytically and on the mesh `|S/g|` below 1e-4 per radian (v0.12; the truncation of
the same interpolant, measured 3.4e-5); a synthetic
composition file whose `x_He` varies as `0.06 + 0.02 sin φ` with H2 as the closure gives the
composition term equal to its closed form (evaluated analytically in the script from the
species table) to 1e-8, and the transfer composition (identical columns) gives zero.

**Acceptance.** The values above; `I = 0` at the `Φ = 0` node of every column exactly.

---

## 6. Step 4: the tracing, the transfer, the outer loop, the estimate

**Deliverable 1: `forward.transfer`.** `trace(mesh, I, Phi_levels, phi_from, phi_to)`: for
each level, RK4 of `dΦ/dφ = −I(φ, Φ)` on the latitude nodes with `I` bilinear, returning
`Φ_k(φ_i)` on every node between; a curve reaching the mesh's `Φ` edge stops the pass and
returns the side and the excess, the caller extends the mesh (Step 2) and repeats, recorded;
refuses a pair of curves that cross ("isobars k and k + 1 cross at φ; the wind field is not balanced there").
`transfer(mesh, K, curves, lnN)`: the trapezoid of `K` along each curve on the same nodes.
`outer_loop(...)`: from the barotropic guess (flat isobars, `p(φ, Φ)` from the anchors'
productions at their own latitudes, log-linear between labels and, between anchors, the map
of the nearest anchor in latitude), form `u` on the mesh, rebuild the columns and kernels,
trace every anchor's levels to every node, form the new `p(φ, Φ)` from the curves, and repeat
until the largest `|Δ ln p|` on the mesh is below `relative_tolerance_ln_p` or
`max_iterations` is reached (refuse, naming the residual). The loop runs for every wind file;
a file whose columns do not vary with pressure converges on the first pass, and the record
says how many passes were taken.

**Deliverable 2: `forward.estimate`.** The estimate takes each anchor as it arrives (Step 0
deliverable 6, after the hook of deliverable 7), with `σ_i²(Φ) = sigma_ln_N_measurement² +
sigma_ln_N_season²` on its levels; it reads no uncertainty from any file. The gauge latitude
`φ_r` is an input, fixed by the driver before the mesh exists (decision K: the weighted
centroid of the construction anchors' latitudes with `1 / mean(σ_i²)` over each anchor's
levels as the weight; `φ_r = φ_1` at M = 1; `forward.estimate.gauge_latitude(anchors)` is the
one place it is computed). Each anchor traced from its latitude to `φ_r` with `ln N`
transferred along its curves: `C_i(Φ)` on its arrival levels (A29). The union of arrival levels at `φ_r`; each
anchor's `C_i` interpolated log-linearly in `Φ` onto the union within its own span, absent
outside it. Weights `w_i(Φ) = 1 / (σ_i² + P_i)`, `σ_i` as it arrives, `P_i = (σ_K (φ_i −
φ_r))²`; `C(Φ)` by A30
over the anchors present at each level; `D_ij(Φ)` (A33) for every pair on their common levels;
the reduced chi-square per level and the inflated variance of `C` where it exceeds one; all of
it in the product's `estimate` group. At M = 1 the union is the anchor's levels, `C = ln N_1`
with no interpolation, and no `D` exists. The isobar labels at `φ_r`: each union level's
pressure from the anchors that reach it (their labels, which must agree to the pressure
identity's tolerance; the disagreement recorded per level), the label carried to the target.
From `C(Φ)` at `φ_r` the curves are traced to the target and `C` transferred along them.

**Expected values.** Instance, closure wind, M = 1 (measured, §1): to 10° N the isobar shift
−31,258.1 m²/s² at the top level, 0 at the gauge, +11,339.7 at the bottom, and `Δ ln N`
+1.1090e-2, +1.0988e-2, +1.0950e-2, +1.0899e-2 at the top, 10 mbar, gauge and bottom levels;
to 60° N the shift −7,117.0 and +2,585.1 m²/s² and `Δ ln N` +2.5119e-3, +2.4922e-3, +2.4840e-3,
+2.4749e-3. To be reproduced to 1e-4 in `Δ ln N` and 1 percent in the shift at the default
spacings, and to 2e-5 and 0.2 percent with both spacings halved, the ratio reported. Instance,
cylinder-extended wind: `Φ_k(φ) = Φ_k` to 300 m²/s² and `ln N_k(φ) = ln N_k` to 2e-3 at every
node (v0.12, decision Q: the truncation floor of the file's grid, measured at Step 3 as 224
m²/s² and 1.1e-3 by the first-order line integrals; the traced values reported beside them,
and their ratio to the Step 3 estimates). **M = 2, the identity test:**
the closure-wind run to 60° N is written by the acceptance script as a synthetic kind N anchor
at 60° N (`radius_m` from the column at its arrival levels, `height_above_anchor_isobar_m` the
`z_lv` there, `refractivity` the transferred `N`, the thermo group's pressures the Lindal anchor's
tabulated pressures (the same isobars, by construction; decision M, so the gauge level matches
`gauge_isobar_Pa` exactly and the anchor isobar is the gauge level) and its temperatures
`p_tab ℛ̄ / (k_B N)` so that the reduction identity returns the transferred `N` exactly, its
anchor radius `r0(60°)` at that level, its `ψ` at 60°, the run's four
inputs embedded, the Lindal companions copied, slug `synthetic60`); a namelist with both
anchors and the target at 10° N: `φ_r` at the weighted centroid (the two `σ_i` equal, so the
midpoint 45.40°), `D_12` below 1e-5 at every common level, the reduced chi-square below one,
the target profile equal to the M = 1 product's `ln N`, `p` and `T` to 1e-5, and the
reference-surface residual of the synthetic anchor below 1 m. Sheared synthetic state: a wind
file built in the script with `u_total(φ, p) = u_reference(φ) [1 + β ln(p_ref / p)]` for
`p < p_ref` and `u_reference` below, `β = 0.1`, `p_ref` 1 bar, the three parts written:
expected values to be measured before this step proceeds; until then the identities: out and
back (`φ_c` to 10° N and back) returns `Φ_k` and `ln N_k` to a tolerance that falls under halving of
both spacings, the ratio reported (no order is claimed: `u'` is piecewise constant, so the
trapezoid and RK4 are first order across the file's nodes and the mesh spacing divides the
file's so that the nodes coincide); the curves stay ordered; the outer loop converges in under
ten passes, its residuals reported.

**Acceptance.** The values and identities above; the crossing and range refusals on
constructed inputs; the loop's pass count reported for every file (one for the closure wind);
a validation anchor (the synthetic anchor with `weight = 0`) leaving `φ_r` and `C` at the M = 1
values and its `C_i` and `D_12` reported;
the M = 2 run repeated with a nonzero `sigma_ln_N_season` column set on the synthetic anchor by
the acceptance script after the hook (not by the hook), which moves `φ_r` and the weights and
appears in the product's anchor group, and changes nothing else in the identity test.

---

## 7. Step 5: the production at the target, altitude and datum, the product, `casspian-forward` in transfer mode

**Deliverable 1: production at the target.** `produce_on_geopotential(N, Phi, R_bar, m_bar,
p_b)`, the SPEC_03 production entered with `Φ` given (the closure path computes `Φ` first and
calls the same function). At the target the arrival levels and `exp C` there, `ℛ̄` and `m̄`
from the run's composition at the target latitude on the isobar labels, `p_b` by decision F.
The pressure identity `p_produced / p_label − 1` at every level, its largest magnitude,
recorded and drawn.

**Deliverable 2: altitude and datum (B7.3, B7.4 under decision B).** `z_lv(φ_t, Φ)` from the
column interpolated to the arrival levels; the datum `Φ_d` where the produced `p` equals
`datum_isobar_Pa`, by log-linear interpolation between the bracketing levels, refused if
outside the produced range; `altitude_m = z_lv(Φ) − z_lv(Φ_d)`; `radius_m` on the levels from
the column. B8 is evaluated with `|g_eff|`, not the radial component; the radial distance is
carried as `radius_m`.

**Deliverable 3: kind `profile` in transfer mode.** As SPEC_03 deliverable 3 with:
`pressure_tabulated_Pa` and `temperature_tabulated_K` absent; `pressure_label_Pa(level)`
(`index`) and `pressure_identity_residual(level)` (`modeled`); `refractivity` `modeled`
(transferred) beside `refractivity_gauge(level)` (`exp C` at `φ_r`, `modeled`);
`geopotential_gauge_m2s2(level)` beside the coordinate (the arrival); `altitude_m(level)`
above the datum along the local vertical, `radius_m(level)` at the target;
`latitude_planetocentric_deg` the target, `gauge_latitude_planetocentric_deg`; scalars
`datum_isobar_Pa`, `datum_geopotential_m2s2`, `datum_radius_m`, `reference_surface_radius_m`
(`r0(φ_t)`); groups `anchors/<slug>` for every anchor (the kind N file verbatim, with its
season and the run's, the hook's record, its reference-surface residual, its `C_i`, its
`sigma_ln_N_measurement`, `sigma_ln_N_season` (with `season_term`) and `P_i` on the union
levels, and its weights); `reference_surface` (`r0(φ)` on the latitude nodes); `isobars` (every anchor's
`Φ_k(φ_i)` and `ln N_k(φ_i)` on the nodes, and the gauge-to-target curves); `estimate`
(`φ_r`, the union levels, `C`, its variance, `D_ij`, the chi-square, the label agreement);
`transfer_record` (the mesh spacings, range and node counts, the outer loop's passes and
residuals, the largest `|S/g|`, the largest isobar shift and its level, the pressure identity
statistics, the estimation keys, the wind and composition seasons against the run's, the
neglected drift). Uncertainty companions NaN with the terms unstated, as SPEC_03, `kernel_error` and
`propagation` added to the list; the estimate's variance is a diagnostic in its group, not the product's
uncertainty.

**Deliverable 4: the driver and the figures.** `casspian-forward` dispatches `mode =
"transfer"` to `forward.transfer.run(namelist)`: load, propagate (the hook), the gauge
latitude (decision K), surface and wind, mesh, kernels, loop, trace every anchor to the gauge,
estimate, trace to the target, produce, altitude, datum, write, figures. F5 gains its across-latitude panel: the reference surface and the traced
isobars in the `(φ, Φ)` plane with every anchor, the gauge and the target marked; F6 in
transfer mode draws the pressure identity residual against the label pressure; a new F7 draws
the delivered `T(p)` and `N(Φ)` beside each anchor's, the isobar shift against `p`, and at
M ≥ 2 `D_ij` against `p`.

**Expected values (instance, measured, §1).** `casspian-forward
forward/lindal_transfer/lindal_transfer.toml` (10° N, closure wind, M = 1): `p` equal to the
labels to the tracing's discretization error, reported at every level (4.1e-7 on the reviewing
agent's mesh; the report's value must fall under halving); `T` on every isobar below the
anchor's by 1.103e-2 at the top level, 1.089e-2 at the gauge, 1.084e-2 at the bottom (to
1e-4); `r0(10°)` 60,128,613.0 m (to 1 m); `altitude_m` 411,132 m at the top level, 98,186 m at
the gauge level, −15,290 m at the bottom (to 5 m). A second run at 60° N: `T` lower by
2.51e-3, 2.48e-3 and 2.47e-3; altitudes 328,127, 78,533 and −12,239 m. A third at 10° N with
the cylinder-extended wind: `T`, `p`, `N` equal to the anchor's to 2e-3 (v0.12, decision Q;
the Step 4 values), altitudes 416,300,
99,290 and −15,456 m to 10 m (v0.10). A fourth at `φ_c` itself: `N` and `Φ` equal to the closure
product's to 1e-12, and `p` and `T` to 1e-5 (the composition regridded from the file's
tabulated levels onto the produced labels moves `ln ℛ̄` by at most 1.5e-6, at the bottom row;
§10 finding 5), and altitudes 376,780, 90,067 and −14,031 m. A fifth, the M = 2 run of
Step 4, carried through production: its product equal to the first run's to 1e-5.

**Acceptance.** The five runs; the product validates as kind `profile` and reads back as a
`DataTree` with the groups listed, one `anchors/<slug>` group per anchor; the datum refused
when placed outside the produced range; the sheared synthetic run carried through production
with its pressure identity reported and converging; F5 with its panel, F6 and F7 rendered by
the driver and by `casspian-plots` by hand, identical apart from the footer; the author views
F5 and F7 and accepts them by eye; the full regression passes.

---

## 8. Decisions this draft proposes (the author's markup)

A. **(Author, 16 September 2026) Kind W is data along the local vertical, in three parts a
person can read, and the model reads the total only.** The reference-level wind, the total, and
their difference along the local vertical; how the field was built is the wind tool's business,
recorded as prose the model does not read; no declaration and no parameterization exists in
the model. The model places the field on its geometry through the isobar map and forms the
axial shear by finite differences; the mapping is a fixed point run for every file. The closure
wind (the cloud-level wind at every level, what Lindal assumed) and the same wind extended
along cylinders are two hypotheses; both are made by tools, the second by the acceptance script
until the wind tool has the option. A diagnostic tool that decomposes a kind W file on a given
geometry into barotropic and baroclinic parts is a later tool. SPEC_00 §6.6 amended (Appendix).

B. **The delivered altitude is measured along the local vertical**, `z = ∫ dΦ / |g_eff|` on
the column, the source's convention (Lindal's altitude column is recovered to 15 m at the top
level) and what a descending probe traverses; B8's radial form differs by `tan²ψ`, 2,668.5 m
at the top of the Lindal profile, and is carried as `radius_m`.

C. **No regridding of the refractivity at M = 1.** Each anchor's levels are the isobars
transferred; the uniform `Φ` grid is the working mesh for the kernels and the columns only. At
M ≥ 2 the anchors' arrival levels at the gauge are combined on their union, each anchor's
`ln N` interpolated log-linearly onto it within its own span, the one interpolation of `N` in
the chain, recorded.

D. **(Author, 16 September 2026) The code is structured for M anchors from Step 0, and the
estimate of A29, A30, A33 and A34 is implemented now.** At M = 1 it is the identity; the
acceptance exercises M = 2 with a synthetic second anchor the code itself produced. The kernel
correction η of A35 and the posterior wind, the reconciliation of the anchors with the wind
hypothesis, are the combination specification, next after this one, which replaces the
placeholder `P_i` with the covariance of A35 and adds the correction inside Step 4's estimate;
the weighted estimate of the reference surface's constant from several measured radii is
there too. Kind N's retrieval instance (geopotential as coordinate, no radius) is the
retrieval leg's; this specification refuses it as not implemented.

E. **The development state is the closure inputs**, rebuilt in the forms of Step 0, with the
composition carried to every latitude; target 10° N; the general paths on synthetic files. The
forward-inputs tools follow SPEC_05, after which the transfer is exercised under other states.

F. **The pressure boundary is planet-wide.** The pressure differential is exact (B7.2), so the
top isobar is one surface and `p_b` one scalar: the label of the topmost isobar of the highest
anchor, or a declared `p_b_Pa`, exactly one. The author's reservation is recorded: the behavior
of the delivered top of the column is to be watched in SPEC_05.

G. **An occultation anchor's levels are placed at its latitude** with `Φ_k` from the
field-line integral under the run's fields; the drift is neglected in latitude and recorded;
the column radius at `(φ_i, Φ_k)` differs from `radius_m` by the slope geometry.

H. **The pressure identity is reported at every level of every run**, and at M ≥ 2 the
agreement of the anchors' labels on shared levels beside it.

I. **(Author, 16 September 2026) One season per run; propagate first, then transfer.** W and C
are the hypothesis for the run's season and must carry it or be uniform in season; an anchor
at another season is moved to the run's season by the propagator in `N` before transfer, and
the wind of its own season never enters; until the propagator exists the mismatch is recorded,
not refused, and (v0.8, decision N) a W or C at another season is likewise recorded and warned,
not refused. The wind's seasonal cycle is not propagated by any machinery of its own: it is in
the delivered wind where the propagated anchors put it through the estimate, and in the
hypothesis elsewhere.

J. **(Author, 17 September 2026) The seasonal term of every anchor's uncertainty, and the
propagation itself, are slots the code carries from Step 0.** The anchor arrives at the
estimate as one object with its uncertainty in named columns, measurement and season, and
every later step reads it there; `forward.propagate` is called for every anchor at every M and
is the identity until the propagator specification fills it, with the absence recorded in the
product. No number is declared for the seasonal term in the namelist: it is computed from the
record by the propagator specification (design note §5.4) and is not zero even at an observed
season, so a declared per-degree rate would have the wrong shape and would be trusted once
forgotten. The estimate's interface does not change when the propagator exists.

K. **The gauge latitude is fixed before the mesh.** `φ_r` is the centroid of the construction
anchors' latitudes weighted by `1 / mean(σ_i²)` over each anchor's levels, computed once from
the anchors as they arrive, before the mesh, and inserted as a node. This is a departure from
A34 as written, whose weights carry `P_i` and depend on `Φ` and on `φ_r` itself; the
departure keeps `φ_r` out of a fixed point and out of the vertical, and is recorded here and
in the product. A note for the manuscript: A34's weights are circular through `P_i`, and the
text should say which weights place the gauge.

L. **One interpolation rule for data on a latitude grid: linear between the nodes, with the
derivatives those of the interpolant.** Kind W is read linearly in latitude and in `ln p`, the
rule the reduction used to register the anchor (`refrac.anchor.wind_of_latitude`), so the
anchor's own wind is the closure's (2.16709 m/s at `φ_c`; PCHIP gives 1.940, 10.5 percent
lower, and a cubic spline 1.994, because the anchor sits on a zero crossing of the curve between nodes
0.5° apart). The model never differentiates `u` on its mesh; `(∂u/∂φ)_p` and `(∂u/∂ln p)_φ`
are the interpolant's, piecewise constant, and the conversion to fixed `r` uses the isobar
map. Kind C is read the same way in latitude and log-linearly in pressure onto the isobar
labels. What a hypothesis looks like between the nodes is the tool's business: a tool that
wants a smoother field writes a finer grid, and the model's rule does not change. The
transfer's `Δ ln N` telescopes to the endpoint wind difference under a wind uniform along the
column, so the interpolant moves the transfer by a few 1e-6 in `ln N`; it moves the anchor's
registration, the reference surface (6.5 m at 10° N between the two interpolants) and the
kernel's local values, which is why one rule is fixed here.

N. **(Author, 21 September 2026) The code is not idiot-proof; the pipeline assumes its inputs
came from the prior steps.** What is refused is stated once in §0 and nowhere else grows:
schema, namelist typo guards, extrapolation, named numerical failures. Everything else is
recorded, warned where useful, never refused. Applied at v0.8 to Step 0 (the refusal list cut
from fifteen to four, the season mismatch recorded and warned, the `Φ` extent no longer
declared) and to every later step as written. A number the code can compute from the data is
never asked of the user: the mesh extent is built from the anchors and grown when a curve needs
it (v0.9). The kind N reader's refusal (SPEC_02) is of a file whose two hash records disagree,
a file edited after it was written; a changed input beside it only warns, as for every derived
kind. It stays; the Step 0 cascade was the closure's own definition (its inputs must be the
reduction's), not the reader's.

Q. **A wind hypothesis on a grid carries a truncation floor, and the model reports it rather
than hiding it.** Under decision L the kernel is formed from the derivatives of the file's
linear interpolant, which jump at the file's nodes. For a wind with vertical shear the two
terms of A15 cancel only to that truncation, so a field that is barotropic analytically (the
cylinder-extended wind) is transferred as if it carried a small shear: on the Lindal state
with the 0.5° grid, `Δ ln N` up to 1.1e-3 over 21° of latitude (0.1 K) and isobar shifts up to
224 m²/s², against 1.1e-2 and 31,000 m²/s² under the closure wind. This is a property of the
data's resolution, not of the code, and it is the floor a wind hypothesis on that grid can
claim; a tool wanting a smaller floor writes a finer grid (the integral falls with the
sampling; the largest local value does not). The checks that call the cylinder state an
"identity" are stated at that floor, the local maximum is reported and not bounded, and the
manuscript carries the floor once where the wind hypothesis is described. A wind with no
vertical shear is exempt: the vertical derivative is exactly zero and the kernel is the closed
form to round-off, which is what the closure wind shows.

P. **The cylinder-extended wind is pinned on the file's reference level and built with the
library's columns at Step 2.** A kind W file's three parts make `u_total(φ, p_ref) =
u_reference(φ)` an identity, so "the same reference-level wind extended along cylinders" can
only mean `U(s)` fixed on the 1 bar isobar; pinning it on the gauge isobar, which is what the
words left open, gives a file whose reference-level wind is not the closure's. The geometry is
the model's own radial columns under the flat-isobar map, exact under `S = 0`, so the file is
only as approximate as the wind grid's interpolation, which is what Steps 3 and 4 bound; the
reference surface and the file are a fixed point, converged in the script; `U` is one
function per hemisphere and the inversion curve is sampled at the mesh spacing (v0.11). The v0.9 values
for this file (9.371, −0.133 m/s, −271.6 m²/s²) could not be reproduced by the coding agent
under any construction the words determine, nor by the reviewing agent afterward; they are
withdrawn and the v0.10 values stand.

O. **(Author, 21 September 2026) Kind C is one structure for every use:** a field on
(level, latitude), pressure the vertical coordinate, uniform in latitude when that is the
hypothesis, the grid declared to the tool that writes it. Every input kind is then a global
field, read by one interpolation rule (decision L), and the one-latitude form of kind C is
retired. Applied at Step 0 inside the same rebuild.

M. **The synthetic anchor of the M = 2 test carries the Lindal anchor's tabulated pressures
in its thermo group** and temperatures `p_tab ℛ̄ / (k_B N)`, so that it is what an
occultation at 60° N reduced by Lindal's method would tabulate on the same isobars: its
reduction identity returns the transferred `N` exactly, its gauge level is a tabulated level
matching `gauge_isobar_Pa`, and its labels, produced at 60° N, are the Lindal anchor's labels
to the pressure identity. An occultation anchor whose tabulated levels do not include the
run's gauge isobar is refused as not implemented in this specification (the regauge by
log-linear interpolation is a later rule).

---

## 9. Revision history

| Version | Date | Change | Cause |
|---|---|---|---|
| 0.12 | 2026-09-23 | Decision Q (the truncation floor of a gridded wind); the cylinder-wind checks of Steps 3, 4 and 5 restated as line integrals at that floor with the maximum reported; the solid-body bound; §14 rulings on REPORT_04_step3 | REVIEW_04_step3 |
| 0.11 | 2026-09-22 | Step 2's column checks restated (central-field closed form; absolute cancellation bound); decision P completed (per hemisphere; `s_ref` sampled at the mesh spacing; the map continued in `ln p`); §13 rulings on REPORT_04_step2 | REVIEW_04_step2 |
| 0.10 | 2026-09-22 | Decision P (the cylinder-extended wind pinned on the file's reference level, built with the columns at Step 2 as a fixed point); its v0.9 values withdrawn and remeasured, the Step 5 altitudes under it restated; §12 rulings on REPORT_04_step1 | REVIEW_04_step1 |
| 0.9 | 2026-09-21 | Decision O (kind C one structure); the mesh self-extending, `geopotential_margin_m2s2` withdrawn; decision N's last sentence corrected (the kind N reader warns on a changed input; the cascade was the closure's definition) | author's direction of 21 September |
| 0.8 | 2026-09-21 | Decision N (what is checked), applied to Step 0 and the mesh; the `Φ` range derived from the anchors with a declared margin; deliverable 6's uncertainty as `σ_N / N`; §10 corrected (PCHIP); §11 rulings on REPORT_04_step0 | author's direction of 21 September; REVIEW_04_step0 |
| 0.7 | 2026-09-21 | Accepted by the author at v0.6; Step 0 proceeds; the manuscript draft of record named | author's acceptance of 21 September 2026 |
| 0.6 | 2026-09-17 | Rulings on the coding agent's pre-execution review (§10): decision L (one linear interpolation rule, derivatives of the interpolant, every expected value remeasured), the Step 0 cascade, `produce` with the wind along the column, decision K (gauge latitude before the mesh), decision M (the synthetic anchor's thermo group), the composition as a field on latitude and pressure, `weight` honored, the `[numerics]` list and the A33 citation corrected, convergence claims withdrawn | coding agent's review of v0.5, 17 September 2026 |
| 0.5 | 2026-09-17 | The anchor as it arrives with named uncertainty columns; the propagation hook as the identity; decision J; the estimate reading no uncertainty from files; `kernel_uncertainty_per_rad` 0.02 confirmed, no open question | author's direction of 17 September 2026 |
| 0.4 | 2026-09-17 | M anchors from Step 0 with the estimate (A29, A30, A33, A34) implemented and an M = 2 identity test; kind W as data in three human-readable parts, the model reading the total and running the isobar map as a fixed point for every file, no declarations or parameterizations in the model; kind C on a latitude grid; the v0.3 tags withdrawn; the seasonal order stated (decision I); every expected value labeled an instance; the forward-inputs tools placed after SPEC_05; decisions A, C, D, E, F rewritten | author's markup of 16 September 2026 |
| 0.3 | 2026-09-16 | The wind hypothesis in the user's coordinates with declarations; write-back; decision D naming the combination specification | author's decisions of 16 September |
| 0.2 | 2026-09-16 | The meaning of wind shear; closure-state expected values measured under the human reading | author's decision of 16 September |
| 0.1 | 2026-09-16 | First draft | SPEC_03 closed |

---

## 10. Rulings on the coding agent's pre-execution review of v0.5 (17 September 2026)

The findings were verified by the reviewing agent from the files on disk before ruling: the
wind file's nodes bracketing the anchor (5.8557 m/s at 30.5°, −0.1800 at 31.0°), linear
2.16709 against PCHIP 1.93991 m/s (a piecewise cubic Hermite; a cubic spline gives 1.99409, as
the coding agent's record corrects); the closure product's produced pressure at the gauge level,
9998.4655 Pa; the largest `|d ln ℛ̄ / d ln p|` of the Lindal composition, 4.5e-4, and the
largest `|ln(p_produced / p_tabulated)|`, 3.29e-3 at the bottom row, so the regrid of finding
5 moves `ln ℛ̄` by at most 1.5e-6.

1. **The interpolant.** Correct in every part, and the fault was the specification's: the
   parenthetical about the wind tool was wrong (the tool fits a constrained penalized B-spline;
   PCHIP is its pole extension), and the reviewing agent's expected values were measured with
   PCHIP, a third rule (the reviewing agent's code used `PchipInterpolator`, so its v0.5 values
came from the 1.940 branch, not a cubic spline's; corrected at v0.8). Decision L: linear between nodes, the reduction's rule,
   the derivatives those of the interpolant. Every expected value is remeasured under it (§1,
   Steps 1, 3, 4, 5); the reference surface under the linear rule returns the reduction's
   equatorial radius to 0.03 m, where the cubic returned it to 0.55 m, which confirms which
   rule registered the anchor.
2. **The cascade.** Correct. Step 0 deliverable 1 now says it: the chain, the closure inputs
   and the closure product are rebuilt as SPEC_03 Step 3 did, no value changing, hashes
   re-recorded in the `step02_1` row format, REPORT_03_step4's product hash superseded.
3. **`produce`.** Correct. Step 1 deliverable 3 states the amendment (optional wind array,
   scalar default, `u_column_ms` recorded) and the Appendix lists it against SPEC_03.
4. **The gauge at M = 2.** Correct, and the fix is in the synthetic anchor, not the rule:
   decision M. The rule "a tabulated level of every occultation anchor, matched exactly"
   stands; an anchor without the gauge among its tabulated levels is refused as not
   implemented.
5. **The composition.** Correct. Decision L covers kind C: a field on `(φ, ln p)`,
   log-linear in pressure onto the labels in transfer mode, closure mode unchanged; the fourth
   expected value of Step 5 is restated (`N`, `Φ` to 1e-12; `p`, `T` to 1e-5, the bound
   1.5e-6).
6. **The convergence claims.** Correct; withdrawn. Under decision L nothing is differenced
   in latitude on the mesh, `S/g` returns the closed form to round-off at the nodes, and the
   tracing and transfer are first order across the file's nodes; the halving ratios are
   reported, not bounded, and the mesh spacing divides the file's. The telescoping the coding
   agent measured is recorded in decision L as the reason the transfer values are insensitive
   to the rule.
7. **`φ_r` before the mesh.** Correct: decision K, the driver order restated, the estimate
   takes `φ_r` as an input.

Shorter items: `weight` is honored as the anchor's role (Step 0 deliverable 4; the acceptance
of Step 4 exercises it); the centroid weights are decision K, recorded as a departure from A34
with a note for the manuscript; "A10's rule" corrected to the rule stated with A33; the
`[numerics]` list corrected (the geopotential and hydrostatic rules are SPEC_03's closed ones
and not namelist keys); the synthetic anchor's kind N requirements are the schema's, and the
script writes it through `lib.io.write` with the groups, the seven self-consistent
`input_hashes` entries and the root scalars the writer validates, listed in the report; the
inserted nodes take the three-point formula for unequal spacing wherever a difference is taken
on the mesh in latitude (Step 2). The measured instance values the coding agent checked for
mutual consistency (the altitudes, the datum shift, A40 to half a percent) are the v0.5
values; the v0.6 values differ from them by the interpolant alone (a few 1e-5 in `Δ ln N`,
6.5 m in `r0(10°)`, 13 m in the altitude at the top level) and are consistent in the same way.

---

## 11. Rulings on REPORT_04_step0 (21 September 2026)

1. **`sigma_ln_N_measurement`.** Correct; the specification's parenthetical was dimensionally
   wrong and deliverable 6 now reads `σ_N / N` times the scale. Verified on the rebuilt kind N:
   `σ_N / N` is 2.331845e-2 at all 66 levels while `σ_N` spans 1.18e-9 to 7.25e-6.
2. **The `[grid]` refusal.** The check is removed rather than moved. The author's direction
   (decision N) and the finding point the same way: the range is not something a user should
   declare and the code police; Step 2 builds it from the anchors' levels and extends it when
   a traced curve needs more (v0.9), and no refusal is left.
3. **Check 15 of `accept_step03_3`.** The `before/` baseline of `step03_3` is refreshed at the
   Step 0 acceptance commit, so the check means "no value changed since SPEC_04 Step 0"; check
   13 of Step 0 then passes on the numbers already measured. Nothing is loosened: this step's
   own check 6 is the same comparison against the swept copies and passes.
4. **The unreachable retrieval refusal.** Removed from the loader (decision N); the retrieval
   leg adds its own check when the instance exists.
5. **`accept_step7`.** The change is the Appendix's amendment applied to a closed step's
   acceptance, recorded in the report; accepted.

Decisions 1, 3, 4 and 6 to 13 are accepted as reported. Decision 2 is superseded by ruling 2.
Decision 5 is superseded by decision O: kind C is a field on (level, latitude) for every use,
and the one-latitude form is retired. The fifteen refusal cases of check 11 reduce to the four the Step 0 acceptance now
names, with the recorded-and-warned season case beside them; the closure namelist's thirteen
refusals are SPEC_03's and stand.

---

## 12. Rulings on REPORT_04_step1 (22 September 2026)

1. **The march grid overshooting the pole.** `through_anchor` snapping a node within 1e-9 rad
   of the pole and refusing one beyond is right; `wind_geoid` is left as it stands, its defect
   changing no value it produced. Accepted.
2. **The cylinder-extended wind.** Correct that the specification did not state the
   construction, and correct that the stated values could not be reproduced: the reviewing
   agent could not reproduce them either, under gauge pinning or reference-level pinning,
   with a linear or a PCHIP wind, and withdraws them. Decision P fixes the construction:
   pinned on the file's reference level (the file's own definition decides between the two
   readings the report names), the library's radial columns under the flat-isobar map, a fixed
   point with the reference surface. Since the columns are Step 2's, the file is built and
   its values exercised there; Step 1 deliverable 3 keeps the sentence for the record and its
   checks 9 (the cylinder half) and 14 are moved to Step 2, not loosened. Decision 11 of the
   report is superseded. The report's measurement that the ratio of the stated to the built
   displacement is one factor of 1.30 was the right observation: it is, to a few percent, the
   height above 1 bar over the height above the gauge at the top of the profile.
3. **`wind_at` bit-identical to the reduction's callable.** Recorded as measured, not
   guaranteed; accepted.
4. **Check 0 of `accept_step04_0` on a clean tree.** The proposal is adopted: the check takes a
   copy of one registered input with a `-dirty` commit written under `reports/step04_0/` as
   its subject, so it means "the refusal fires on a `-dirty` input" whatever the tree's state.
   Applied at the Step 1 acceptance commit, recorded in the report.

Decisions 1 to 10 and 12 are accepted as reported; `WindField` as one class, the coverage as
the grid's own extent, `ValueError` for extrapolation, the exact partials with the north and
higher-pressure rule at a node, `u_column` defaulting to the broadcast scalar, the NaN
companion. Section 5b's costs are noted for the Monte Carlo wrapper; nothing is asked.

---

## 13. Rulings on REPORT_04_step2 (22 September 2026)

1. **The uniform gravity check.** Correct: with `GM / r²` in `lib.gravity` those constants give
   a central field, and the specification's closed form was wrong. The check is restated as the
   central field's closed form, `r = r0 / (1 − r0 Φ / GM)`, to 1e-12, which the column meets at
   8.9e-16.
2. **The cancellation bound.** Correct; `z_lv = r − r0` to 1e-7 m absolute, seven ulp of `r`.
3. **The inversion curve's sampling.** Accepted and generalized: `s_ref` is sampled at the mesh's
   latitude spacing over the file's range, so every latitude a level maps to, not only the
   anchors and the target, is within a mesh spacing of a sample. The reviewing agent's expected
   values were measured with a 0.05° sampling, which is why the finer sampling should close on
   them.
4. **Per hemisphere.** Correct, and a property of the hypothesis the specification should have
   stated; decision P now does. Nothing in the run's span crosses the equator.
5. **`edge_order = 2`.** Correct; recorded in Step 2's staggering as the rule.
6. **The 10° N bottom level at 0.054 m/s.** Not loosened. The departure is the size and sign of
   the chord error of finding 3 at the latitude the bottom level maps to (between file nodes at
   10.0° and 10.5°), and ruling 3 removes it; the check is rerun under the finer sampling
   against the same values and bound. If it still misses, the report says so with the measured
   values and the reviewing agent remeasures under the same sampling before any bound moves.

Decisions 1 to 11 are accepted as reported; decision 8 (the flat-isobar map continued at the
slope of its last interval) is now stated in Step 2, and decisions 9 and 10 are superseded by
rulings 3 and 4 as generalized.

---

## 14. Rulings on REPORT_04_step3 (23 September 2026)

The reviewing agent rebuilt the kernel independently on the coding agent's own Step 2 cylinder
file (the bilinear interpolant, its piecewise-constant partials, the conversion to fixed `r`
with the column map, the columns marched under the file's wind) before ruling.

1. **Finding 1, the solid-body check.** Correct: "to round-off" was written for analytic
   derivatives; under decision L the residual is the interpolant's truncation. Restated as
   `|S/g|` below 1e-4 per radian, measured 3.4e-5.
2. **Finding 2, the map's continuation.** The Step 3 script clamps the flat-isobar map at the
   anchor's ends; Step 2 decision 8 continues it at the slope of its last interval, which is
   continuous in slope and leaves no kink. The eightfold degradation at the ends, the NaN of
   finding 4 (two pressure nodes coinciding) and most of check 9's maximum are that clamp:
   the reviewing agent measures the largest `|S/g|` on the span at 0.050 per radian with the
   map continued and 0.136 with it clamped, the excess at the clamped rows and their
   neighbors through the two-sided difference. The script follows decision 8; no change to
   the specification's text.
3. **Finding 3, the exact zero.** The difference form is right and is the rule for the
   composition term. Accepted.
4. **Finding 4.** The synthetic states posed on the mesh's own geometry is right. The NaN is
   the clamp's (ruling 2), and a valid kind W file has strictly increasing nodes, which is the
   schema's business; no guard in `lib.windfield` (decision N). Accepted.
5. **Finding 5, the cylinder-extended wind.** Right about the mechanism, wrong about the
   measure and the size. The kernel of a piecewise-linear wind does have a spiky,
   sign-alternating residual in the cells holding a slope jump, and its maximum does not fall
   under refinement. But the transfer applies the integral of the kernel along an isobar, not
   its maximum, and `max |S/g|` times the span is not `Δ ln N`: the reviewing agent's line
   integrals along the anchor's isobars give `Δ ln N` up to 1.1e-3 over the span (1.2e-4 at
   the top level, 5.2e-4 at 10 mbar, 1.7e-4 at the gauge, 1.1e-3 at the bottom) and isobar
   shifts up to 224 m²/s², where the report's bound implied 7.6e-2. Decision L and A15 are
   not in tension; a gridded wind has a truncation floor, decision Q states it, and the
   cylinder checks of Steps 3, 4 and 5 are restated at that floor as line integrals, with the
   maximum reported. Neither decision L nor decision P changes. The report's correction of its
   own `d ln N` column is noted; the replacement column was still the wrong quantity.

Decisions 1 to 6 are accepted as reported.

---

## Appendix. Amendments to SPEC_00 and SPEC_03, to be applied at acceptance

- §3.4: "The model does not decompose the wind" stands, restated: the model reads `u_total`,
  places it on its own geometry through the isobar map, and forms the shear kernel from that
  field; it reads no declaration of how the field was built.
- §6.6: kind W carries `u_reference_ms(latitude)` at `reference_level_pressure_Pa`,
  `u_total_ms(latitude, pressure)`, `u_shear_ms = u_total − u_reference` (along the local
  vertical), `u_total_uncertainty_ms`, `value_provenance`, the rotation system, the season, the
  coverage attributes, and prose provenance (`method`, `vertical_structure`,
  `parameterization`) the model does not read; `decomposition`, `u_cylindrical_ms`,
  `decomposition_geometry` and the `Ω_abs` cylinder check retired (a diagnostic tool on a given
  geometry later). The reader checks the sum identity and the poles.
- §6.2: kind C is a field on (level, latitude) for every use, `latitude_grid_deg` required by
  the tool, the one-latitude form and `latitude_planetocentric_absent_meaning` retired
  (decision O). SPEC_01 Step 8 amended the same way; the reduction reads the column at the
  anchor's latitude by decision L.
- §6.8 (kind `profile`): the transfer-mode variables, scalars and groups of Step 5.
- §7.2: transfer mode as Step 0 states; `[[anchors]]` one or more, `weight` as the anchor's
  role honored by the estimate; `[estimation]` gains `kernel_uncertainty_per_rad`; the
  `[numerics]` schemes as implemented.
- §6.2 and §6.6: the interpolation rule of decision L stated once, for kinds W and C.
- SPEC_03 Step 4 (closed; amended at acceptance): `produce` gains an optional wind on the
  anchor's levels, the scalar it forms today the default; the product gains
  `u_column_ms(level)`; closure mode and its expected values unchanged.
- §7.3: closed. The vertical grid is the working mesh only; the latitude grid uniform with the
  anchors, the gauge latitude and the target as nodes; the integrations and the outer loop as
  the steps state; the estimate at M anchors as Step 4 states.
- Manuscript notes: B8 restated with `|g_eff|` (decision B); B4, B6.4 and B8 per
  `CASSPIAN_Manuscript_Retrievals_and_Wind_Note.md` v0.4. That note and the manuscript are the
  author's, kept outside the repository while the manuscript is in flux; neither is a
  dependency of any step, and the equation labels this specification cites are those of the
  manuscript draft the author supplies to the coding agent.

**Ruled (author, 17 September 2026):** `kernel_uncertainty_per_rad = 0.02` per radian, the
declared uncertainty of the kernel that sets `P_i` until A35 replaces it (a fifth of the
closure wind's largest kernel). It acts on nothing at M = 1 and on the M = 2 identity test only
through the gauge latitude. It is a namelist value, required with no default in code, recorded
in the product, and changed by editing the namelist and rerunning; the combination
specification replaces it with the covariance formed from the declared W and C uncertainties.
No question remains open in this draft.
