# SPEC 02. `refrac`: the reduction to refractivity, the kind N product, and the standard diagnostics

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.7, 12 September 2026. Author of record: S. Rafkin. Status: accepted by the author
at v0.4 (12 September 2026); Steps 1 to 4 accepted (Step 3 with the v0.6 changes, Step 4 with
the v0.7 change); Step 5 proceeds. See §8 for the revision history.
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.14 and on the closed
`SPEC_01_Lindal_Tool_Chain.md` v0.19 (closed at commit `ece58d2`, Step 8 amended at
`31c30ef`), which it does not repeat.

---

## 0. How this specification is to be worked

The protocol of SPEC_01 §0 applies unchanged: one step at a time; acceptance script under
`reports/step02_<N>/` (ignored; SPEC_01's `reports/step<N>/` directories remain and its suites
still run, v0.3), figures a report cites under `reports/figures/` (committed);
`reports/REPORT_02_step<N>.md` and `reports/REVIEW_02_step<N>.md`; `STATE.md` updated; commit
on acceptance, then the directory sweep of SPEC_00 §8 (every `-dirty` product under
`occul_data/` rebuilt, whichever step wrote it); do not implement ahead. Every angle in `lib` is
in radians; the manifest's `fixed_point_tolerance_deg` is converted once, at the `refrac`
boundary, and nowhere else. No em dash or en dash anywhere.

**What this specification covers and what it does not.** `refrac` is Appendix B3 plus the two
constants frozen per profile (SPEC_00 §3.3): number density from the tabulated pressure and
temperature, the anchor latitude and anchor radius on the wind-included geoid, absolute radius
per level, refractivity from number density with the reduction composition, and the kind N
product with everything it consumed embedded. Nothing in it uses a chosen input; every input is
the source's own. The geopotential coordinate (B5), the hydrostatic production (B7) and the
transfer (B6) are the forward model and belong to SPEC_03, whose first step is the round trip:
the forward model run at the source latitude with the source's own inputs, which must return
the tabulated pressure and temperature. The standard diagnostics module (Step 5 here) is built
against kind N; its geopotential and hydrostatic closure figures are specified here so the
module has its full shape from the start, and they render when the file offered carries those
fields, which is from SPEC_03 on.

---

## 1. Step 1: `lib.control` reduction manifest and input loading

**Purpose.** Read `<profile>_reduction.toml` (SPEC_00 §7.1) and load the six inputs under
their kinds, refusing everything SPEC_00 §8 says to refuse, before any arithmetic.

**Deliverable: `src/casspian/lib/control.py`** (the general control-file reader SPEC_00 §3.1
names; `tools/gravity/control.py` is absorbed into it, its named-key path resolution kept):

- `read_reduction_manifest(path)`: parses the file, checks every section and key against the
  §7.1 vocabulary, refuses an unknown key, refuses a bare physical value, resolves paths
  relative to the manifest's directory, and returns a frozen record.
- `load_reduction_inputs(manifest)`: reads the six files with `lib.io.read` under their kinds
  and returns them with their SHA-256 hashes. Refuses when: any input carries `-dirty`; the
  wind file's rotation system or rate differs from the rotation file; the wind file's
  `reference_level_pressure_Pa` is not on its pressure grid; the anchor isobar of the manifest
  is not a surface of the geodesy file; the anchor isobar is not a tabulated level of the
  thermo file (the reduction needs `h_ref` at that exact level, not an interpolation); the
  composition file's levels are not the thermo file's levels.

**Acceptance.** `lindal_reduction.toml` loads with all six inputs and the hashes of
REPORT_01_step9 §6. A copy of the manifest naming a `-dirty` product is refused with the
file named. A copy of the manifest with `pressure_Pa = 9000` under `[anchor_isobar]` is
refused because 9000 Pa is neither a geodesy surface nor a tabulated level. A copy with an
extra key `polar_radius_m = 54438000` under `[geoid]` is refused as a physical value in a
control file.

---

## 2. Step 2: the frozen constants, `phi_c` and `r0`

**Purpose.** The two numbers that fix the profile in space: the planetocentric latitude of the
anchor and the absolute radius of the anchor isobar there, on the wind-included geoid, from
the source's own G, R, W and D. Handoff §9A; SPEC_01 Steps 5 to 7.

**Deliverable: `src/casspian/refrac/anchor.py`**, one function `freeze_anchor(inputs,
manifest)` returning a record, built from `lib` only:

1. `u_of_phi`: the kind W column at `reference_level_pressure_Pa`, interpolated linearly in
   planetocentric latitude on the file's grid (0.5° for Lindal), zero at the poles by the
   file's own rule, as a callable in radians returning arrays shaped like its argument.
2. The wind geoid: `lib.geoid.wind_geoid` with `anchor_rule` and `r_anchor` from the manifest
   and the geodesy file (`radius_polar_m` on the anchor surface), the source's G and R, and
   `u_of_phi`. The north polar start, both polar radii, the asymmetry and the anchoring
   residual are kept for the record.
3. `phi_c`: `lib.latitude.planetocentric_fixed_point` for the profile's planetographic
   latitude (from kind T), with `surface` a callable that marches the wind geoid with the
   iterate's latitude inserted as a node (SPEC_01 v0.7: the anchor radius is never
   interpolated), `u_of_phi` as above, the ellipsoid seed from the anchor surface's oblateness
   (kind D), tolerance from the manifest converted to radians, and the iterates kept.
4. `r0 = r_wind(phi_c)`, read from the final march at the node `phi_c`; `psi` at `phi_c`.
5. Two reference values computed alongside and recorded, not used: `phi_c` and `r0` on the
   **no-wind** reference geoid (`lib.geoid.reference_geoid` solved at the iterate), which are
   the handoff §9A.8 numbers, so the effect of the wind on the frozen pair is visible in the
   file.

**Acceptance.** No-wind reference values: `phi_c` = 30.8182 ± 0.001° and `r0` = 58,453.1 ±
0.5 km with Null's GM (SPEC_01 Steps 5 and 6). Wind-included values reported: an independent
march with the SPEC_01 v0.16 wind and anchoring gives the wind surface near 58,518.5 km at the
no-wind latitude (about 65 km of dynamical height); the frozen `phi_c` itself moves about
0.013° equatorward under the wind and the surface rises about 1.3 km over that shift, so the
frozen `r0` lands near 58,519.9 km (REPORT_02_step2: 58,519.883 km, `phi_c` = 30.804949°);
these are expected figures, to be reproduced, not prescribed. Convergence of the fixed point in at most eight iterations to the
manifest tolerance; the iterates listed. The anchoring residual below 1e-3 m. Polar radii and
asymmetry as REPORT_01_step7 §6 (54,423.6 and 54,452.4 km; 28.7 km) to 0.1 km, since the
inputs are the same files.

---

## 3. Step 3: number density, absolute radius, refractivity (B3.1, B1, B3.3)

**Purpose.** The reduction proper. Pure `lib` functions plus their orchestration.

**Deliverable: `src/casspian/lib/reduction.py`** (pure) and `src/casspian/refrac/reduce.py`
(orchestration):

1. **B3.1.** `number_density(p_Pa, T_K)` = `p / (k_B T)` with `k_B` from `lib.constants`
   (CODATA release recorded).
2. **B1.** `absolute_radius(h_m, h_ref_m, r0_m)` = `r0 + (h − h_ref)`, with `h_ref` the
   tabulated height at the anchor isobar (Step 1 guarantees the level exists).
3. **B3.3.** From the composition file and its species group: `mean_refractivity_m3(level)`
   = Σ x_i ℛ_i and `mean_molar_mass_kg_mol(level)` = Σ x_i M_i, with the per-molecule ℛ_i in
   m³ as the species group carries them; `refractivity(level)` = `n ℛ̄`.

**Uncertainty companions (v0.4): first-order propagation of the declared input uncertainties,
and nothing else.** The companions of kind N are the linearized (first-order Taylor)
propagation of the uncertainties the input files declare, through the algebra above, added in
quadrature, with no covariance assumed between different inputs. This is bookkeeping: the
partial derivatives are properties of the algebra at the solution and involve no choice.
Everything beyond it (covariance between inputs, asymmetry, nonlinearity, alternative rules,
Monte Carlo) is the job of the generic Monte Carlo wrapper of §7, which perturbs input files
and reruns the pipeline; nothing of that kind lives in `refrac`.

- `n`: `(δn/n)² = (δp/p)² + (δT/T)²` from the thermo companions.
- `ℛ̄` and `m̄`: the derivative is taken **with the closure Σ x_i = 1 applied**, never by
  treating the mole fractions as independent. For a composition declared as a share `s` of a
  remainder (Lindal: `x_H2 = s (1 − x_NH3)`, `x_He = (1 − s)(1 − x_NH3)`, `x_H2_uncertainty`
  = δs = 0.03), `δℛ̄ = (1 − x_NH3)(ℛ_H2 − ℛ_He) δs` and likewise for `m̄` with the molar
  masses. A stated per-molecule uncertainty `δℛ_i` from the species group enters as `x_i δℛ_i`
  in quadrature (the `lindal1985` set states none). A stated uncertainty on a species mole
  fraction that is not the closure species enters through the same constrained derivative.
- `N`: `(δN/N)² = (δn/n)² + (δℛ̄/ℛ̄)²`.
- `r`: `δr² = δr0² + δh² + δh_ref²` (at the anchor level `h − h_ref` is identically zero and
  carries no height term), with `δr0` from two independent declared sources added in
  quadrature after conversion to one standard deviation (SPEC_00 §5 v0.13): the anchor radius
  uncertainty of kind D (`radius_polar_uncertainty_m`, 10 km `1sigma` for Lindal) through the
  **total** derivative `dr0/dr_anchor`, which includes the coupling through `phi_c` since
  `psi` depends on radius (v0.6; the latitude-held partial 1.217 plus the coupling 0.020 gives
  1.237 for Lindal, REPORT_02_step3 findings 1 and 2; it is not near unity because raising the
  anchor also raises the oblateness), evaluated by a central difference on the anchored march
  with the fixed point rerun; and the label latitude uncertainty of kind T (0.2° declared
  `range`, so 0.115° as one standard deviation) through `dφ_c/dφ_g` (0.9456 on the wind geoid
  at the frozen pair; 0.934 on the no-wind geoid) and `∂r0/∂φ_c = r0 G_φ / g` from Eq. B3 at
  the solution (−5,630 km per radian). The scalars `phi_c`, `psi` and `r0` carry the
  corresponding companions, all `1sigma`, with `uncertainty_kind_conversions` listing the
  range conversion.
- **Unstated terms are left out, not set to zero, and the file says so.** Each companion is
  the quadrature of the stated terms only; two attributes on each, `uncertainty_terms_included`
  and `uncertainty_terms_unstated`, list which inputs entered and which were NaN, so a reader
  knows that Lindal's `N` uncertainty is composition and anchor and nothing else, not that the
  temperature contributes zero. A companion with no stated term at all is NaN.
- The anchor-rule spread (about ±19 km) is a choice, not a declared uncertainty, and is not
  in the companions; it is recorded as a number in `reduction_record`.

**Acceptance.** At the 1 bar level (p = 100,000.0 Pa, T = 134.8 K): `n` = 5.3731e25 m⁻³
(Phase 1 check; 5.373124e25 with CODATA 2018 k_B), `N` = 2.59856e-4 with the tabulated 10.9
ppm of NH3 at that level (ℛ̄ = 4.836272e-30 × (1 − 10.9e-6) m³; v0.2 quoted the dry value
2.59859e-4 by mistake). At the 794.33 mbar level, where NH3 is exactly zero by the SPEC_01
v0.17 rule, `N / n` equals the dry ℛ̄ = 4.836272e-30 m³ to round-off. Top level (20 Pa, 138.7 K): `n` = 1.04441e22 m⁻³. Bottom level (129,848 Pa, 146.2
K): `n` = 6.43287e25 m⁻³, `N` reported with the NH3-reduced ℛ̄ (a factor 1 − 79.3e-6 below the
dry value). `radius_m` at the anchor level equals `r0` exactly and at the 1 bar level equals
`r0 − h_ref` (the 1 bar height is zero by the table's datum). Fractional uncertainty of `N` at
every level equals the composition term, 0.0233 ± 0.0002 (0.03 × (136 − 35) / 129.94; v0.2
said 0.0227, an arithmetic slip), and is constant across levels because ammonia dilutes H2 and
He by the same factor; `uncertainty_terms_included` on `refractivity` reads composition only
and `uncertainty_terms_unstated` reads pressure and temperature. `radius_uncertainty_m` is the
same at every level (the height companions are NaN) and equals the quadrature of the anchor
term (10 km × 1.237 = 12.37 km) and the label term (0.115° × 0.9456 × 5,630 km/rad = 10.73
km), 16.38 km, with both partials, both terms and the range conversion reported; `phi_c`
carries 0.109°. Step 3 also computes `mean_refractivity_uncertainty_m3` and
`mean_molar_mass_uncertainty_kg_mol` for kind N to carry (SPEC_00 §6.7 v0.13), and reads the
closure declaration of kind C (SPEC_00 §6.2 v0.13) rather than inferring it, checking the
declaration against the values. Recovering
`T = p ℛ̄ / (k_B N)` from the products returns the tabulated temperature to 1e-12 relative at
every level (the inverse closes).

---

## 4. Step 4: the kind N product

**Purpose.** Write `<profile>_refractivity.nc` per SPEC_00 §6.7 (v0.10), the only file the
forward model reads from the reduction.

**Deliverable: `src/casspian/refrac/product.py`** and the entry point `casspian-refrac`
(thin: takes the manifest path, calls one function, returns), writing:

- Variables on `level`: `radius_m` (coordinate), `height_above_anchor_isobar_m`,
  `number_density_m3`, `refractivity`, `mean_refractivity_m3`, `mean_molar_mass_kg_mol`,
  and the uncertainty companions of the first four by the SPEC_00 naming rule.
- Scalars with uncertainties: `latitude_planetocentric_deg`, `psi_deg`,
  `latitude_planetographic_deg` (copied from kind T with its swath and uncertainty),
  `anchor_isobar_pressure_Pa`, `anchor_isobar_radius_m`, `anchor_isobar_height_m`.
- Groups: `inputs/thermo`, `inputs/composition`, `inputs/geodesy`, `inputs/gravity`,
  `inputs/rotation`, `inputs/wind`, verbatim copies; `manifest` with `text` and `sha256`;
  `reduction_record` with attributes: the fixed-point iterates, `anchor_rule`, north polar
  start, both polar radii, polar asymmetry, anchoring residual, no-wind `phi_c` and `r0`,
  the anchor-rule spread (v0.7: the whole Step 2 fixed point rerun under `north_pole` and
  under `south_pole` anchoring, recording `r0` and `phi_c` under each rule; not the
  latitude-held march, for the same reason the Step 3 anchor term is the total derivative;
  recorded as a choice), the partial derivatives of the Step 3 propagation, the CODATA
  release, `casspian_version`, `casspian_git_commit`.

**Expected values (v0.7).** With the anchor radius effectively moved by half the polar
asymmetry, 14.37 km, and the total derivative 1.237, the spread is about ±17.8 km at the
anchor: `r0` near 58,537.7 km under `north_pole` and 58,502.1 km under `south_pole` (each
±0.1 km), against 58,519.883 km under `mean_polar_radius`; `phi_c` moves by
∓0.0030° (the coupling `dphi_c/dr_anchor` = −3.62e-9 rad/m), to about 30.8020° and 30.8079°.
The latitude-held values the v0.6 build recorded, ±17.5 km, differ by the coupling term.
- Global `input_hashes` over the six inputs and the manifest.

**Acceptance.** The file reads back as kind N (a DataTree); `read` refuses a copy missing the
`inputs/wind` group and a copy whose `input_hashes` entry for the thermo file has been edited;
every embedded input group is byte-identical in variables and attributes to the file it came
from (compare through `xarray.Dataset.identical`); the six scalars match Step 2 and Step 3;
`lindal_refractivity.nc` appears in `occul_data/lindal/` and the SPEC_00 §2.2 listing is now
complete for the profile; the anchor-rule spread in `reduction_record` matches the expected
values above (v0.7).

---

## 5. Step 5: the standard diagnostics, `casspian.tools.plots`

**Purpose.** A pipeline feature, not a report artifact: a generic tool under `tools/`, tied to
no profile, that renders a fixed set of figures from a file alone, dispatching on the file's
`casspian_kind`, called by `refrac` after it writes kind N when the manifest asks for it, and
callable by hand on any CASSPIAN file (v0.2: placed under `tools/plots/`, not a top-level
subpackage, at the author's direction). Figures are generated from
the file, never from the run's memory, so a figure always shows what the file contains and
can be regenerated for any product ever written. The same module serves the forward model's
product later; it renders the figures whose fields the offered file carries and says which it
skipped.

**Manifest.** SPEC_00 §7.1 (v0.10) gains

```toml
[diagnostics]
figures = true
format  = "png"        # "png" or "pdf"; a combined multi-page PDF is written in either case
dpi     = 150
```

**Output.** `occul_data/<profile>/figures/<prefix>_diag_F<k>_<name>.<format>` and
`<prefix>_diag.pdf` (all figures, one per page). Regenerable, therefore ignored by git (a
`.gitignore` rule beside `reports/step*/`). Every figure carries a footer with the product
file name, its `casspian_git_commit`, the first twelve characters of its SHA-256, and the
generation time, so a figure is traceable to a file.

**Deliverable: `src/casspian/tools/plots/`** with `render(path, out_dir, format, dpi)` and
the entry point `casspian-plots <file.nc> [--out DIR]`. `render` reads `casspian_kind` and
dispatches: kind N renders the figure set below; a kind W, T or C file offered on its own
renders a single figure of that input (the F1 panel that concerns it, plus for W the `u(φ)`
curve with its uncertainty band and provenance shading), which is the quick look a new input
file gets before it enters a manifest; the forward product renders its own set when SPEC_03
defines it; any other kind is refused by name. matplotlib only, one style module
(`tools/plots/style.py`: fonts, line widths, the color per quantity, no seaborn, no external
stylesheet). Pressure is the vertical axis wherever a profile is shown, log scale, decreasing
upward, in mbar on the axis label with Pa in the file. The anchor isobar is drawn as a
horizontal line on every profile panel and the profile latitude as a vertical line on every
latitude panel. Nothing in `tools/plots` computes a number that reaches a product; where a
panel needs a derived quantity (gravity along the profile), it calls `lib` on the embedded
inputs and says so in the panel title.

**The figures.**

- **F1, inputs.** Four panels: `T(p)` with its uncertainty band where finite; `u(p)` at
  `phi_c` from the embedded kind W (a vertical line for an altitude-independent wind; the
  panel exists because forward winds carry shear) with an inset of `u(φ)` on the reference
  level and `phi_c` marked; `x_i(p)` on a log axis for every species in the set, with the
  per-level provenance of NH3 shown by marker; `m̄(p)` and `ℛ̄(p)` on twin axes.
- **F2, gravity along the profile.** From `lib.gravity` on the embedded G, R, W at `phi_c`
  and `radius_m(level)`: `g_N(p)` and `g_eff(p)`; their difference split into the rigid
  centrifugal part and the wind part; `G_φ(p)`; `ψ(p)`.
- **F3, gravity and shape across latitude.** On the anchor isobar: `g_N(φ)` and `g_eff(φ)`;
  the no-wind reference geoid and the wind geoid `r₀(φ)` with the dynamical height between
  them (the SPEC_01 Step 7 figure promoted to a standard output, with Lindal's Eq. 18 curve
  omitted here); `ψ(φ)`; `phi_c` and `r0` marked, both polar radii annotated.
- **F4, the product.** `n(p)` and `N(p)` with their uncertainty bands; `N` against
  `radius_m`; the recovered `T = p ℛ̄ / (k_B N)` against the tabulated `T` as a fractional
  difference, which must be at round-off and is drawn so that it would be visible if it were
  not.
- **F5, geopotential** (renders when the file carries `geopotential_m2s2(level)`): `Φ(p)` with
  `h(p)` on a twin axis; `Φ` on the anchor isobar against latitude.
- **F6, hydrostatic closure** (renders when the file carries `pressure_hydrostatic_Pa(level)`):
  `p_hydro / p_tab − 1` against `p`, with the declared boundary pressure and its level marked.
  This is the panel that says whether the source profile, the gravity and the composition are
  mutually consistent, and it is the first honest estimate of the profile's internal error.
**Acceptance.** `casspian-refrac` with `figures = true` writes F1 to F4 and the combined PDF
for `lindal_refractivity.nc` and reports F5 and F6 skipped with the field each needs;
`casspian-plots` on the same file by hand writes identical figures (byte-identical PNGs apart
from the generation time in the footer, which the acceptance strips before comparing); every
figure carries the footer; `casspian-plots lindal_wind.nc` writes the single-figure wind
view and `casspian-plots lindal_gravity.nc` is refused by kind; F4's recovered-temperature panel shows a fractional difference
below 1e-12 at every level; F2's `g_eff` at the 1 bar level agrees with `lib.gravity` called
directly to round-off; the author views the figures (attached to the report) and accepts
them by eye, which is the only check a figure can finally have.

---

## 6. What comes after

SPEC_03 opens the forward model with the round trip: B5 (geopotential from the effective
gravity along the profile, the numerics and the level-versus-layer accounting stated before
coding), B7.1 and B7.2 at the source latitude with the source's own inputs and `p_b` at the
top tabulated level, recovering the tabulated pressure and temperature; F5 and F6 render for
the first time on that product. Then the forward-inputs specification (generic wind tool,
composition scenario generator) and the transfer. **The Monte Carlo wrapper (v0.4)** is a
generic tool, specified after the forward model exists so that it covers the whole chain from
input files to delivered temperature and pressure: it draws perturbed input files from their
declared uncertainties and from declared alternative rules into a scratch profile or run
directory, reruns the pipeline, and collects the products. Nothing in `refrac` or `forward`
carries a sensitivity study; the pipeline's job is to be deterministic, file-driven and fast
enough to be wrapped.

---

## 7. Decisions made in this document

1. `refrac` ends at kind N; geopotential and hydrostatic pressure are forward-model quantities
   and arrive with the round trip in SPEC_03. The diagnostics module is specified in full here
   so that it does not grow piecemeal, and renders what the offered file carries.
2. Figures are generated from product files only, never from run state, and are regenerable
   and uncommitted.
3. Kind N gains `mean_refractivity_m3` and `mean_molar_mass_kg_mol` (SPEC_00 v0.10) because
   every consumer of the file needs them and they are cheap to carry.
4. Kind N's uncertainty companions are first-order propagation of declared input
   uncertainties only, uncorrelated between inputs, with the closure constraint applied inside
   the composition; unstated terms are left out, not zeroed, and the file lists which entered.
5. The anchor-rule spread is recorded as a choice, not folded into the radius uncertainty.
   It is measured by rerunning the fixed point under each rule (v0.7), so that the recorded
   spread is the whole effect of the choice, latitude shift included. At the anchor it is
   about ±17.8 km, larger than the 16.4 km the declared uncertainties propagate to; the
   anchor rule is the largest single term in where the profile sits in radius, and the Monte
   Carlo wrapper is where that choice is exercised (REPORT_02_step4 finding 1).
6. No sensitivity study or Monte Carlo lives in the pipeline; a generic wrapper over input
   files does that, specified after the forward model (v0.4; the v0.1 Step 6 is withdrawn).
7. A companion copied from an input (the planetographic label and its `range`) keeps the kind
   its source declared; the SPEC_00 §5 conversion applies to companions that combine terms
   (REPORT_02_step4 decision 2).
8. The `input_hashes` consistency warning of SPEC_00 §8 is implemented for kind N, whose
   inputs are embedded. Extending it to the other derived kinds is deferred to the reader work
   SPEC_03 will need (REPORT_02_step4 finding 3).

---

## 8. Revision history

| Version | Date | Change | Cause |
|---|---|---|---|
| 0.1 | 2026-09-12 | First draft: Steps 1 to 6 with acceptance numbers; diagnostics as a pipeline feature after the kind N product | SPEC_01 closed; author direction on standard output figures |
| 0.2 | 2026-09-12 | Diagnostics placed under `tools/plots/` as a generic, kind-dispatching tool (entry point `casspian-plots`), with single-figure views of input kinds | author direction |
| 0.3 | 2026-09-12 | Status remains draft (Step 1 built ahead and accepted on its own; no further step before acceptance); step directories `reports/step02_<N>/`; Step 3 acceptance corrected (composition term 0.0233; N at 1 bar with the tabulated NH3, dry check moved to 794.33 mbar) | REPORT_02_step1 §3 |
| 0.4 | 2026-09-12 | Accepted by the author. Step 6 (sensitivity and wind Monte Carlo) withdrawn in favor of a generic Monte Carlo wrapper after the forward model; Step 3 uncertainty companions specified as first-order propagation with the closure constraint and the label and anchor terms on `r0`; F7 dropped; `[sensitivity]` removed from the manifest | author discussion of sensitivity philosophy |
| 0.5 | 2026-09-12 | Step 2 accepted; Step 2 expected values split into the fixed-latitude dynamical height and the latitude shift; Step 3 slope and label term corrected to −5,630 km/rad and about 18 km | REPORT_02_step2 findings 1 and 2 |
| 0.6 | 2026-09-12 | Step 3 accepted with changes: total derivative `dr0/dr_anchor`; uncertainty kinds converted before quadrature (label `range` divided by √3); measured partials recorded; closure read from kind C's declaration; mean refractivity and molar mass companions carried into kind N | REPORT_02_step3 findings 1 to 6 |
| 0.7 | 2026-09-12 | Step 4 accepted with one change: the anchor-rule spread by a full fixed-point rerun under each rule, `r0` and `phi_c` recorded; decisions 7 and 8; dependency lines updated to SPEC_00 v0.14 and SPEC_01 v0.19 | REPORT_02_step4 findings 1 and 3, decisions 2 and 5 |
