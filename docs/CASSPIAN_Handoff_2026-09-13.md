# CASSPIAN handoff, 13 September 2026

Saturn atmosphere reference model. Handoff at the close of SPEC_02, for a fresh reviewing
agent (the role that writes specifications and reviews reports) and a fresh coding agent (the
role that builds to the specifications and writes reports). This document is meant to stand on
its own without the chat that produced it. Everything it summarizes is also in the repository;
where the two disagree, the repository is right and this document is stale.

Author of record: S. Rafkin, Southwest Research Institute.

---

## 1. What CASSPIAN is, and where it stands

CASSPIAN is a one-dimensional Saturn atmospheric structure model for probe mission design,
with quantified uncertainty. The physics is in the manuscript
`C:\Users\srafkin\OneDrive - SWRI\Documents\Proposals\SaturnProbe\Environment\CASSPIAN_AtmosphericModel.docx`
and its earlier handoff `CASSPIAN_Handoff_Coding_Phase.md` (same folder). The manuscript is the
authority on the physics; the specifications in the repository are the authority on the code.

The pipeline has two halves. The **reduction** (`refrac`) takes a published radio occultation
profile (so far: Lindal, Sweetnam and Eshleman 1985, Voyager 2 ingress at 36.3 N planetographic)
and reduces it to refractivity N on a registered absolute radius grid, with the inputs that went
into it (gravity, rotation, wind, composition, thermo profile, geodesy) embedded in the product.
The **forward model** (`forward`, not yet built) takes N plus a chosen composition and gravity
and produces temperature and pressure by hydrostatic integration, at the source latitude or
transferred along isobars to another latitude.

As of commit `2cf4457` the reduction is complete and closed for the Lindal profile.
`occul_data/lindal/` holds the raw bundle, the six inputs, the manifest, the kind N product
`lindal_refractivity.nc`, and the standard diagnostic figures. Every product carries a clean
commit and the SHA-256 of everything it was made from. The next specification to write is
SPEC_03, the forward model round trip.

---

## 2. Roles, and how the work is done

There are three parties. The **author** (Scot) makes every decision of substance, accepts or
rejects each step, and is the only one who commits by his own hand outside the agents. The
**reviewing agent** writes specifications into `docs/specs/`, reviews each coding report
against the specification, and writes the review into `reports/`. The **coding agent** builds
exactly what the current specification step says, runs the acceptance, writes the report, and
does not commit until the step is accepted.

The loop for every step: the specification states the step, its deliverables, its expected
values and its acceptance; the coding agent builds it and writes
`reports/REPORT_0<spec>_step<N>.md` with the acceptance results, its decisions and its findings;
the reviewing agent verifies the numbers independently where it can, writes
`reports/REVIEW_0<spec>_step<N>.md` with a disposition and rulings, revises the specifications
if needed, and updates `docs/specs/STATE.md`; the author reads both and says go; the coding
agent commits (the acceptance commit), then runs the **sweep**: every product under
`occul_data/` carrying a `-dirty` commit is rebuilt on the clean tree and its new hash is
recorded in the report. Only then does the next step start.

The coding agent **must read the status line at the top of each specification and
`docs/specs/STATE.md` before starting any step**. It once read "ready" as acceptance and
started early; the status lines now say exactly what is accepted and what proceeds.

Steps are small and modular, one at a time, each verified before the next. Architecture before
code. Tools before models. All data in netCDF with units, provenance, and citations. Ad hoc
testing (acceptance scripts under `reports/step<dir>/`), not a test framework. Every
specification revision is recorded in its revision history so that a clean v1.0 "as built" can
be produced at the end. Simplest method that works; nothing fancy. Accuracy, precision and
exactness are the priorities, and a fix is made regardless of cost when something is found to be
wrong.

**Standing rules for every document, comment, commit message and chat reply:** no em dashes or
en dashes used as parentheticals anywhere (a minus sign in a number is fine); American
spellings; no first person in the manuscript; no AI attribution in any document, code or
commit; do not explain elementary physics to the author; verify, do not assert; call balls and
strikes; do not refer to text by paragraph or line number. The work is SwRI funded and
proprietary: implementation details do not flow to external collaborators or into
team-facing material, and questions to outside experts (for example Julie Moses) are phrased
so as not to disclose the method.

The author works in Word and Windows. The repository is `C:\Users\srafkin\SaturnProfile`,
mirrored at GitHub `ScotRafkin/SaturnProfile`. The reviewing agent reaches it through the
Claude desktop app's device bridge (staging files to read, committing files to write).

---

## 3. Repository layout and the specifications

```
docs/specs/
  SPEC_00_Architecture_and_Data_Files.md   v0.15  architecture, layout, file kinds, control files, rules
  SPEC_01_Lindal_Tool_Chain.md             v0.20  closed at ece58d2; Steps 0 to 9; amendments v0.19, v0.20
  SPEC_02_Refrac_and_Diagnostics.md        v0.9   closed at 2cf4457; Steps 1 to 6
  STATE.md                                        one row per step, status, report, date
  CASSPIAN_Fig1_RepositoryLayout.png/.svg, fig1_repository_layout.py
src/casspian/
  lib/        constants, schema, io, control, gravity, geoid, latitude, reduction   (pure functions, radians)
  tools/      lindal/ (build_raw, build_inputs), gravity/, wind/, composition/, plots/
  refrac/     anchor.py, reduce.py, product.py                                     (entry point casspian-refrac)
  forward/    not yet built (SPEC_03)
data_static/  species_master.toml; winds/ (Ingersoll and Pollard digitization, Sanchez-Lavega table, notes)
occul_data/lindal/
  raw/        lindal_table1.csv, lindal_scalars.toml, lindal_raw.nc, notes.md
  lindal_build.toml            control file for the Lindal tools (pointers and declared choices only)
  lindal_gravity.nc (G)  lindal_rotation.nc (R)  lindal_wind.nc (W)  lindal_composition.nc (C)
  lindal_thermo.nc (T)   lindal_geodesy.nc (D)
  lindal_reduction.toml        the manifest refrac reads
  lindal_refractivity.nc (N)   the product
  figures/                     standard diagnostics, regenerable, not committed
reports/
  REPORT_01_step0..9.md, REVIEW_01_step0..9.md, REPORT_02_step1..6.md, REVIEW_02_step1..6.md
  step1..9/ and step02_1..6/   acceptance scripts, output.txt, regression.txt
  figures/                     committed copies of figures cited by reports
```

**SPEC_00** defines the principles (every number the model consumes comes from a
self-describing netCDF file; control files hold pointers and declared choices only and the
parser refuses physical values and unknown keys; NaN means no information and is never filled
silently; unstated uncertainties are NaN, never 0.0), the file kinds and their required
variables and attributes (§6), the provenance vocabulary (`measured`, `derived`, `assumed`,
`inferred`, `interpolated`, `parameterized`, `extrapolated`, `modeled`, `index`), the
uncertainty naming rule (`_uncertainty` inserted before the unit suffix), the rule for combining
uncertainty kinds (`range` is a uniform half-width divided by √3, `stated` is taken as 1σ,
conversions recorded in an attribute), the manifest vocabulary (§7.1), the `-dirty` rule and
the sweep (§8), the `input_hashes` warning, and the standard diagnostics in `tools/plots`
(§3.5). §7.3 (grid and numerics for the forward model) is still open and is SPEC_03's job.

**SPEC_01** built the library and the Lindal tool chain. **SPEC_02** built `refrac` and the
diagnostics. Both carry a decisions section and a revision history; read those before anything
else in them.

---

## 4. The physics as implemented, and the conventions that were verified

Effective gravity (manuscript Eq. A1 to A5): g_eff = −g r̂ + G_φ φ̂ with g positive inward,
Ω_abs = Ω + u/(r cos φ), g = g_N − Ω_abs² r cos²φ, G_φ = −(1/r)∂V/∂φ − Ω_abs² r cos φ sin φ,
taken along increasing planetocentric latitude and therefore negative in the north. The tilt of
the local vertical is ψ = arctan2(−G_φ, g) and the planetographic latitude is φ_g = φ_c + ψ. The
manuscript, Lindal's Eqs. 2, 3, 5 and 8, and the code all agree on this convention; it was
checked three ways. (The old handoff §9A.2 shorthand line drops the sign of g_r; it is a defect
in that document, not the code.) Lindal's Eq. 11 U = V − ½ω²r²cos²φ with g = −grad U is identical
to the code's U_rigid.

Gravity potential V = −(GM/r)[1 − Σ J_l (R/r)^l P_l(sin φ_c)], g_N = ∂V/∂r positive inward. The
trap to remember: the radial derivative gives a factor (l + 1) on J_l, not (2l + 1); the wrong
factor gives 9.29 m/s² at the Lindal latitude instead of 9.10. Lindal's Table II values for
checking: 8.951 m/s² at the equator with u = 450, 12.137 at the pole; Jupiter 23.116 and 27.015.

The geoid (Eq. B3): g dr₀/dφ = r₀ G_φ. `lib.geoid.wind_geoid` marches it once, continuously,
from the north pole to the south pole (RK4, a dense 0.05° grid unioned with the caller's nodes,
both poles and the equator always exact nodes), with a latitude-dependent wind supplied as a
callable. The slope equation has one free constant, fixed by a declared `anchor_rule`:
`mean_polar_radius`, `north_pole`, `south_pole`, or `equatorial_radius` (added at SPEC_02 Step 6).
The no-wind reference geoid is a Newton solve on the rigid-rotation potential, anchored where
the rule anchors. `radius_at` interpolates 1/r² linearly in sin²φ and refuses intervals that
straddle the equator. Two things the march taught us: a surface in balance with an asymmetric
wind has different polar radii (28.7 km, south higher, for the Ingersoll and Pollard wind), and
the asymmetry scales as r³ with the size of the surface because the slope carries accelerations
proportional to r over gravity proportional to 1/r².

The frozen constants of the reduction: the fixed point φ_c + ψ(φ_c) = φ_g with φ_g = 36.3° the
Fig. 4 label; ψ evaluated at φ_c on the final march; r₀ the radius of the anchor isobar (100 mbar,
the surface Lindal fitted his geoid on) at φ_c. Absolute radius r = r₀ + (h − h_ref) (Eq. B1).
Number density n = p/(k_B T) (B3.1). Refractivity N = n ℛ̄ with ℛ̄ = Σ x_i ℛ_i (B3.3), stored
unscaled (not times 10⁶). Temperature recovers from N to machine precision as the closure check.

Uncertainty companions in kind N are first-order propagation of declared input uncertainties
only, uncorrelated between inputs, with the composition closure Σx = 1 applied inside the
derivative, each input kind converted to 1σ before quadrature, unstated terms left out and
listed in `uncertainty_terms_unstated`. Discrete choices (the anchor rule) are recorded in
`reduction_record`, never folded into a companion. No sensitivity study lives in the pipeline;
a generic Monte Carlo wrapper over perturbed input files is specified after the forward model.

---

## 5. The Lindal reduction as built: inputs and decisions

**Thermo (kind T).** Lindal Table I, 66 levels from 0.2 mbar to 1298 mbar, T and altitude
relative to 1 bar; no uncertainties stated (NaN). Label 36.3 N planetographic from the Fig. 4
label, swath 36.3 to 36.7, uncertainty 0.2° declared as `range`.

**Geodesy (kind D).** Lindal's fitted 100 mbar surface: equatorial radius 60,367 ± 4 km, mean
polar radius 54,438 ± 10 km, and the 1 bar pair 60,268 / 54,364. Fit residual stated as a range
3 to 4 km (per-surface NaN, `fit_residual_range_m`). His Fig. 9 caption says the south polar
radius may be of the order of 10 km greater than the north, possibly seasonal.

**Gravity and rotation (G, R).** Null et al. 1981 GM and J2, J4, J6 with Campbell 1984 personal
communication as Lindal used them (`[gravity_used_by_source]` in `lindal_scalars.toml`);
System III rotation.

**Wind (kind W).** The digitized solid curve of Ingersoll and Pollard 1982 Fig. 5 (Lindal's cited
wind source), 708 samples at 0.2°, PCHIP between samples, peak 490.5 m/s at 7.4°, 75 m/s near
30.8°, 2 m/s near 36°. Uncertainty from the Smith et al. 1982 Fig. 4 points as RMS about the
curve per 2° bin (21 m/s). The ring-obscured southern equatorial flank is filled by reflecting
the north (their dashed curve), bridged across the equator with a cubic, blended into the
southern observed segment over −10.9 to −15°, and both polar caps go to exactly zero at ±90
(PCHIP to zero). Provenance codes per node: 0 observed, 1 interpolated, 2 parameterized,
3 extrapolated, 4 extended by source assumption. Altitude independent (trivial decomposition).
The digitization script, overlay and note are in `data_static/winds/`. A constrained spline
fitter (`fit.py`) is retained as a non-default alternative. The Sanchez-Lavega et al. 2000
tabulated Voyager profile (`vasavada_saturn_winds.txt`, via Vasavada via Moses; not Cassini;
the ±90 rows are added, not data) is registered with a note as a second declared wind source
and is **not** used yet; its epoch and the meaning of its u_rms column are to be confirmed from
the paper. It is for the wind uncertainty study later.

**Composition (kind C).** H2 and He split 0.94/0.06 of the remainder after NH3, with
`x_H2_uncertainty` 0.03 as Lindal states; NH3 from Table I at nine levels (831.76 to 1258.93
mbar), interior gap interpolated, one level extrapolated downward, and **zero above 831.76 mbar
by our assumption** (saturation reasoning recorded; Lindal does not say). Closure declared in
the file (`closure_rule = "share_of_remainder"`, `closure_species = "H2 He"`). Refractivities
from the `lindal1985` set in `species_master.toml`, in which NH3 is zero (his Fig. 3 caption
reading); per-molecule uncertainties NaN. Mean molar mass 2.135083 amu dry.

**The numbers** (equatorial anchoring, the current default): φ_c = 30.805568°, ψ = 5.494432°
(φ_c + ψ = 36.3 to 3e-8°), r₀ = 58,516.188 km, radius uncertainty 11.35 km (anchor term 3.70 km
from ±4 km times dr₀/dr_anchor 0.924; label term 10.73 km from 0.2°/√3 times dφ_c/dφ_g 0.9456
times 5,630 km/rad), φ_c uncertainty 0.109°, N fractional uncertainty 0.0233 at every level
(composition only). Predicted mean polar radius 54,435.01 km (Lindal 54,438 ± 10), north
54,420.6, south 54,449.4, asymmetry 28.740 km. Under the previous default `mean_polar_radius`:
φ_c = 30.804949°, r₀ = 58,519.883 km, uncertainty 16.38 km, equator 60,371.0 km (+4.0 over
Lindal). No-wind reference values (Null GM, polar anchoring): φ_c = 30.818189°, r₀ = 58,453.120 km.

**Why the equator (SPEC_02 decision 9).** Both of Lindal's radii come from one ellipse fit; the
equatorial radius is the better-constrained end (±4 against ±10) and it is a single point of
the surface, so anchoring there removes any reading of his polar number and makes the polar
radii predictions. The intuition that the pole is the safer anchor because the wind is zero
there did not survive measurement: r₀ at 30.8° is the integral of the slope from the anchor to
the profile under either rule, and the mean-polar rule needs the whole pole-to-pole march to
fix its mean, so polar anchoring is the more wind-sensitive (3.38 km against 2.60 km for a
±5 percent wind scaling). The polar rules remain in the code and the manifest vocabulary, the
spread under them is recorded in `reduction_record`, and the SPEC_02 Step 2 to 4 acceptance
values remain valid as values under `mean_polar_radius` (those suites pin that rule in memory).

**Open physics question carried forward:** the march's 28.7 km polar asymmetry against
Lindal's "order of 10 km". The asymmetry is a property of the wind field (weakest in the south:
reflected flank, observed only to −72.8°, extended caps), unchanged by the anchor. It is the
most wind-sensitive quantity produced and belongs to the wind uncertainty study, not to any
step before it.

---

## 6. Kind N, the product the forward model reads

`lindal_refractivity.nc`, SPEC_00 §6.7. On `level`: `radius_m` (coordinate, positive up),
`height_above_anchor_isobar_m`, `number_density_m3`, `refractivity`, `mean_refractivity_m3`,
`mean_molar_mass_kg_mol`, each with its uncertainty companion carrying `uncertainty_kind`,
`uncertainty_method`, `uncertainty_terms_included`, `uncertainty_terms_unstated`,
`uncertainty_kind_conversions`. Scalars with companions: `latitude_planetocentric_deg`,
`psi_deg`, `latitude_planetographic_deg` (copied from T with its `range` kind), 
`anchor_isobar_pressure_Pa` (provenance `index`), `anchor_isobar_radius_m`,
`anchor_isobar_height_m`. Groups: `inputs/thermo`, `inputs/composition` (with
`inputs/composition/species`), `inputs/geodesy`, `inputs/gravity`, `inputs/rotation`,
`inputs/wind`, verbatim copies; `manifest` (attributes `text`, `sha256`); `reduction_record`
(about 60 attributes: iterates, anchor rule and quantity, polar radii and mean, asymmetry,
equatorial radius marched, residual, no-wind pair, spread under the other rules with φ_c and r₀
under each, every partial derivative and term of the propagation, declared kinds, CODATA
release, commit, version, `input_sha256_<input>`). Global `input_hashes` lists the six inputs
and the manifest; `read` refuses a file whose recorded hash is missing from `input_hashes` and
warns when an input on disk no longer matches. Levels run from the top (20 Pa) to the bottom
(129,848 Pa), so radius decreases with index.

Reading it: `casspian.lib.io.read(path)` returns an `xarray.DataTree` for kinds N and C and a
`Dataset` for the others, validating against the schema. The author inspects files with
`ncvue` (installed on his machine at his request).

---

## 7. Standard diagnostics

`casspian.tools.plots`, entry point `casspian-plots <file.nc> [--out DIR] [--format png|pdf]
[--dpi N]`, generic and kind-dispatching, also called by `casspian-refrac` when the manifest's
`[diagnostics] figures = true` (the Lindal manifest now says so). Figures are generated from the
file only, never from run state, and carry a footer with the file name, commit, SHA-256 prefix
and generation time. For kind N: F1 inputs, F2 gravity along the profile, F3 gravity and shape
across latitude, F4 the product; F5 geopotential and F6 hydrostatic closure render only when the
file carries `geopotential_m2s2(level)` and `pressure_hydrostatic_Pa(level)` (an attribute
`boundary_pressure_Pa` marks the boundary), which SPEC_03 is to adopt or rename in both places.
Kind W, T or C offered alone gets a single-figure view. Output goes to
`occul_data/<profile>/figures/`, ignored by git; committed copies cited by reports live under
`reports/figures/`. The author has judged the figures adequate for now and will do an appearance
cleanup pass himself later, after F5 and F6 exist; the mechanism is what mattered.

---

## 8. What comes next, in order

**SPEC_03, the forward model round trip (to be written by the reviewing agent before any
code).** Purpose: at the source latitude, with the source's own inputs, recover Lindal's
tabulated T and p from N, which proves the forward machinery before any transfer or new
composition. Contents the specification must settle before coding: the geopotential (B5) from
the effective gravity along the profile, with the numerics and a level-versus-layer staggering
table stated explicitly (which quantities live on levels, which on layers, how the integral is
discretized, what order, and the expected discretization error); the gauge isobar where Φ = 0
(normally the anchor isobar, but a separate declared choice); the hydrostatic production (B7.1,
B7.2) with m̄/ℛ̄ inside the integral and the boundary pressure p_b declared at the top tabulated
level (about 0.2 mbar); the acceptance as the recovered T and p against Table I with a stated
tolerance and a stated reason for the tolerance; F5 and F6 rendering for the first time; the
forward product kind (SPEC_00 §6 gains it, with `modeled` provenance) and its manifest; and the
`input_hashes` warning extended to derived kinds generally, since the forward reader will need
it. SPEC_00 §7.3 (grid and numerics) is closed by this specification.

**The forward-inputs specification.** A generic wind tool (any declared wind source into kind W
with the same rules), and a composition scenario generator from the Moses photochemical files:
helium as a knob with H2 as the closure; NH4SH equilibrium in the Lewis form (K in atm², pressure
squared; the Fortran fragment from J. Moses uses pressure to the first power and her file was
built with it, which reproduces her NH4SH base near 2 bar where the squared form gives about
5.5 bar; to be raised with her, wording was drafted); NH3 and H2S saturation at a declared
relative humidity; min(saturation, photochemical) aloft; PH3 shape movable. A modern NH3
refractivity entry (Debye static value with a temperature law) is still to be entered in
`species_master.toml`. Composition for forward runs is a choice, not a retrieval; many
scenarios are expected.

**The Monte Carlo wrapper**, generic, after the forward model: draws perturbed input files from
their declared uncertainties and declared alternative rules (wind draws from the bin scatter,
anchor rule, gap rule, label latitude, composition scenario) into scratch directories, reruns
the pipeline, and collects the products. This is where the anchor rule and the wind are
exercised under both anchoring rules, and where the Sanchez-Lavega profile is used.

**Then the transfer** of N along isobars to other latitudes (the manuscript's central idea) and
the delivered profiles.

---

## 9. Open items and notes, none blocking

Figure appearance cleanup (author, later). Kinds T and D carry a `control_file` hash from the
build they were written under, not the current `lindal_build.toml` (kept by content at Step 6;
recorded and accepted). The `input_hashes` warning is implemented for kind N only. Confirm the
Sanchez-Lavega 2000 epoch and u_rms meaning from the paper before it is used. Raise the NH4SH
pressure power with Moses. Manuscript notes to carry: the geoid section should state the
equatorial anchoring and its reasons (decision 9); §B2 Campbell and Anderson J4 discrepancy;
the old handoff §9A.2 sign shorthand; "Lindal asserted a temperature at the top" is an
inference, not his statement. External context: Julie Moses supplied the photochemical
composition files (`kinsatfull_4casspian.pun`, 27 S, 198 levels, 127 species, He 0.100 above
versus 0.110 deep, the decline being real diffusion per her) and answered six questions;
attribution "personal communication" is fine with her; the assessment and email drafts are in
the author's notes, not the repository.

---

## 10. Instructions to give a fresh coding agent

Read, in this order: this document; `docs/specs/SPEC_00_Architecture_and_Data_Files.md` in
full; the status line, decisions and revision history of SPEC_01 and SPEC_02; `docs/specs/STATE.md`;
the most recent report and review (`reports/REPORT_02_step6.md`, `reports/REVIEW_02_step6.md`)
to see the expected shape and standard of a report. Then wait. Nothing is to be built until
SPEC_03 exists and its status line says a step proceeds. When it does: build only that step,
run its acceptance, write the report with the acceptance table, decisions and findings, state
plainly what fails and why, do not commit, and do not start the next step until the review is
in and the author says go. Run the full regression (all earlier suites) after any library
change. After the acceptance commit, run the sweep and record the hashes. Keep every standing
rule in §2, especially: no dashes as parentheticals in anything you write, no attribution, and
no silent choices (a choice the specification did not make is a decision in the report, or a
question to the author, never a default in the code).

## 11. Instructions to give a fresh reviewing agent

Read the same set, plus the manuscript's Appendix A and B (Eqs. A1 to A5, B1, B3, B5, B7) and
Lindal 1985 (Eqs. 2, 3, 5, 8, 11, 18, Table I, Table II, Fig. 4, Fig. 5, Fig. 9 caption), which
the author will supply from the paths in §1. Write SPEC_03 as a draft for the author's markup
first; do not let coding begin on a draft. When reviewing a report, open the product files and
check the numbers independently (xarray and numpy are enough) rather than trusting the
acceptance script; state what you verified and how. When the coding agent is right and the
specification was wrong, say so in the review and fix the specification, with the revision
recorded. Numbers in this document worth keeping at hand for sanity checks are in §4 and §5.
