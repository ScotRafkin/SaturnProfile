# SPEC 01. The Lindal tool chain, the library it needs, and the repository reset

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.8, 11 September 2026. Author of record: S. Rafkin. Status: in work; Steps 0 to 6 accepted. See §13 for the revision history.
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.6, which it does not repeat.

---

## 0. How this specification is to be worked

**One step at a time.** This document is a sequence of steps, each producing one deliverable
with its own acceptance checks. The coding agent implements one step, runs its acceptance
checks from a script kept under `reports/step<N>/` (ignored; SPEC_00 §2), writes `reports/REPORT_01_step<N>.md` stating what was built, what deviates from the
spec and why, the acceptance results with the actual numbers, and any question that blocked it,
updates `docs/specs/STATE.md`, and **stops**. The next step begins only after the report has
been reviewed. Do not implement ahead. Do not build two steps in one commit.

**Housekeeping at the start of Step 2 (one commit, before the step's code).** Move the
existing `REPORT_01_step0.md`, `REPORT_01_step1.md`, `REVIEW_01_step0.md`, and
`REVIEW_01_step1.md` from `docs/specs/` to `reports/` with `git mv`; change the `.gitignore`
rule from `reports/` to `reports/step*/`; commit. From then on reports and reviews are written
to `reports/` directly.

**Commit on acceptance, then rebuild.** A step's code is committed when its review says
`accepted`. Any data product the step wrote before that commit carries `-dirty` in
`casspian_git_commit`; after the acceptance commit the step's tool is rerun so that every
product on disk carries the clean commit, and only those products are inputs to the next step
(SPEC_00 §8).

**The spec does not describe tests to be written.** The acceptance checks below are numbers to
be reproduced and reported; how they are checked (a script, a notebook, a printed value) is the
coding agent's choice. No test framework is required at this stage.

**Rules that apply to every step.** SPEC_00 §1 (principles), §3.1 (array rule, naming rule),
§5 (file conventions), §8 (provenance and refusal). No physical value in code except the
constants of nature in `lib/constants.py`. No em or en dashes in code comments, docstrings,
commit messages, or reports. No AI attribution anywhere. American spelling. Every module names
in a comment which manuscript equations it implements.

**Numerical rules already fixed** (SPEC_00 §7.3): near-exponential quantities are log-linear
across a layer; other quantities linear; every variable has an assigned location on the grid.
This specification reaches only as far as the reduction inputs, so the only integrals in it are
the geoid construction and the planetographic fixed point, and both are stated explicitly.

---

## 1. Step 0: repository reset and skeleton

**Purpose.** Bring the local repository and the GitHub remote to a single clean state that
matches SPEC_00 Figure 1, without losing the Phase 1 record.

**Procedure.**

1. Commit the current uncommitted working tree exactly as it is, with the message
   `Phase 1 prototype as of September 2026, superseded by SPEC_00`. Before committing, add
   any file that the current `.gitignore` excludes but that is part of the Phase 1 record, with
   `git add -f`; as of this writing that is `data/lindal1985_voyager2.nc`, ignored by the
   `*.nc` rule. Confirm with `git status --ignored` that nothing else under the tree is left
   behind (caches and `__pycache__` are not part of the record and stay ignored). Tag the
   commit `phase1-prototype`. Push the commit and the tag. This is the only time Phase 1 code
   is committed; it is retrievable from the tag forever and is never merged forward. Note for
   the record: that netCDF carries the superseded 36.5° latitude attribute (handoff §9A.8); it
   is preserved as history, not as data.
2. In one further commit, remove `occultation_recon/`, `tests/`, `docs/this_spec.md`, and
   `docs/CASSPIAN_OccultationReconstruction_BuildSpec.md`. Move
   `data/lindal1985_table1.csv` to `occul_data/lindal/raw/lindal_table1.csv` and delete
   `data/`. Leave `docs/Lindal_et_al_1985_AJ90_1136.pdf`, `docs/specs/`, `data_static/`, and
   `occul_data/lindal/raw/lindal_scalars.toml` where they are. Do not create an `attic`
   directory; the tag is the attic.
3. Create the skeleton: `src/casspian/{lib,tools,refrac,forward}/__init__.py`,
   `src/casspian/tools/{lindal,wind,composition,gravity}/__init__.py`, empty `forward/` at the
   repository root with a `README.md` saying what goes there, and `occul_data/lindal/` with
   `raw/` as above.
4. Replace `pyproject.toml`: package name `casspian`, `src` layout, `requires-python >= 3.11`,
   dependencies `numpy`, `scipy`, `xarray`, `netCDF4`; `tomllib` is standard library. No
   console entry points yet; they are added as tools appear. Version `0.1.0`.
5. Replace `.gitignore`: ignore `__pycache__/`, `*.egg-info/`, `build/`, `dist/`,
   `.pytest_cache/`, and `forward/*/output/`. Do **not** ignore `*.nc` globally. Ignore
   `occul_data/*/*.nc` and `occul_data/*/raw/*.nc` for now, with a comment that they are
   un-ignored once a reduction is accepted (SPEC_00 §2.4).
6. Replace `README.md` with ten lines: what the repository is, that `docs/specs/SPEC_00` is the
   architecture, that `docs/specs/STATE.md` is the current state, and that the package is
   proprietary to Southwest Research Institute.
7. Create `docs/specs/STATE.md`: a table with one row per step of this specification, columns
   step, deliverable, status (`not started`, `in progress`, `reported`, `accepted`), report
   file, date. Step 0 is its first `reported` row.
8. Commit, push. Confirm that `git status` is clean and that the GitHub remote shows the same
   HEAD.

**Acceptance.** `git log --oneline` shows the two commits above the Phase 1 history; the tag
resolves; `pip install -e .` succeeds and `python -c "import casspian"` runs; the tree matches
SPEC_00 Figure 1 for every directory that exists at this stage; `STATE.md` exists.

**Out of scope.** Nothing from `occultation_recon/` is copied into `casspian`. Two Phase 1
values survive as acceptance numbers in later steps (2.1351 amu; 5.3731e25 m⁻³), and the
Table I CSV survives as data. Everything else is rewritten to this specification.

---

## 2. Step 1: `lib.constants` and the netCDF conventions layer

**Purpose.** The two pieces every tool needs before it can write a file.

**Item 0, before any code: line endings.** Add `.gitattributes` at the repository root pinning
text files to LF (`*.csv`, `*.toml`, `*.md`, `*.py`, `*.txt`, `.gitignore`, `.gitattributes`:
`text eol=lf`) and marking `*.nc`, `*.png`, `*.svg`, `*.pdf` as `binary`; run
`git add --renormalize .`; commit and push. Provenance hashes are of the bytes on disk, and this
makes those bytes the same on every platform. One commit, before Deliverable 1a.

**Version rule.** `pyproject.toml` and `casspian.__version__` must agree. `write()` obtains the
version through `importlib.metadata` and refuses if the two differ.

**Deliverable 1a: `src/casspian/lib/constants.py`.** CODATA 2018 values, each with a comment
giving the CODATA name: Boltzmann constant, Avogadro constant, atomic mass constant, molar gas
constant, Newtonian constant of gravitation, Loschmidt constant at 273.15 K and 101325 Pa, and
the two defining values 273.15 K and 101325 Pa. A module constant `CODATA_RELEASE = "2018"`.
Nothing else. No Saturn quantity of any kind.

**Deliverable 1b: `src/casspian/lib/schema.py` and `src/casspian/lib/io.py`.** The conventions
of SPEC_00 §5, and the kinds of SPEC_00 §6, as code:

- A registry of kinds (`thermo`, `composition`, `geodesy`, `gravity`, `rotation`, `wind`,
  `refractivity`, `raw`), each with its required dimensions, variables, per-variable
  attributes, and global attributes, and a schema version integer. The `raw` kind requires only
  the §5 globals.
- `write(path, dataset, kind)`: takes an xarray Dataset, validates it against the kind, fills
  in the §5 globals that are the writer's job (`Conventions`, `casspian_kind`,
  `casspian_schema_version`, `created_by`, `created_at`, `casspian_git_commit`,
  `codata_release`), and writes netCDF-4. Refuses on any missing required item with a message
  naming it. Supports groups (needed for the `species` group of kind C, the `inputs/*` groups of
  kind N, and the raw bundle).
- `read(path, kind)`: opens, checks `casspian_kind` against the kind asked for, checks the
  schema version is not newer than the reader's, validates, returns an xarray Dataset (with
  groups accessible). Refuses per SPEC_00 §8.
- `sha256(path)` and a helper that formats `input_hashes` entries.
- `history_append(dataset, line)`.

Latitude, radius, and pressure names are enforced by the registry: a variable or attribute
named `latitude`, `r`, `z`, or `p` with no convention or unit in its name is refused.

**Acceptance.** Write a two-level `raw` file with one variable and read it back with every §5
global present and correctly typed; `read` refuses it when asked for kind `thermo`; `read`
refuses a file whose `casspian_kind` attribute has been edited to `wind` but which lacks the
wind variables; a `sha256` of a file matches `sha256sum` on the command line; a dataset with a
variable named `latitude` is refused by `write`.

---

## 3. Step 2: stage one of the Lindal tool, the raw bundle

**Purpose.** Turn the two ASCII transcriptions into `occul_data/lindal/raw/lindal_raw.nc`
(SPEC_00 §2.2.1).

**Inputs.** `occul_data/lindal/raw/lindal_table1.csv` (trusted for its values: 66 rows;
columns pressure in mbar, temperature in K, NH3 in ppm with blanks where not measured, height
in km above the source's 1 bar level) and `occul_data/lindal/raw/lindal_scalars.toml`
(transcription of every scalar, with `value_source` per entry; its keys are the vocabulary).

**Deliverable: `src/casspian/tools/lindal/build_raw.py`**, with a console entry point
`casspian-lindal-raw`. It reads the two files and writes the raw bundle with:

- Group `table1`: dimension `level` ordered by increasing pressure; variables
  `pressure_Pa`, `temperature_K`, `nh3_mole_fraction` (NaN where the CSV is blank),
  `height_m`; units converted (mbar to Pa, ppm to mole fraction, km to m) and the conversions
  named in `history`; per-variable `provenance = "measured"` for height and `"derived"` for
  the rest, `value_source = "Table I"`.
- Group `scalars`: every table of the TOML reproduced as attributes on a sub-group of the same
  name, values as numeric attributes where numeric and strings otherwise, with each
  `value_source` beside its value. Nothing is interpreted, converted, or dropped; the bundle
  is the transcription in netCDF form.
- Globals per §5, plus `raw_sources` listing both ASCII files and `raw/notes.md` with hashes
  (the notes are documentation, not read by the tool, hashed so the record is complete).

**Acceptance.** 66 levels; the 1 bar row reads exactly 100000.0 Pa, 134.8 K, 0.0 m; the top
row 20.0 Pa, 138.7 K, 376700.0 m; the bottom row 129848.0 Pa, 146.2 K, −14100.0 m; exactly
nine finite NH3 values, from 831.76 mbar (2.6 ppm) to 1258.93 mbar (66.9 ppm), with NaN at
1047.13 and 1298.48 mbar (those two are interior and end gaps that the composition tool fills
later; no filled value may appear in the raw bundle); the `scalars/latitude` group carries
`planetographic_deg = 36.3` and the swath; `read(path, "raw")` succeeds and
`read(path, "thermo")` refuses. The report lists every attribute in the `scalars` groups so the
transcription can be checked against the paper by a second reader.

---

## 4. Step 3: the gravity and rotation tools, kinds G and R

**Purpose.** Standard harmonic-set and rotation-system files from the `data_static`
transcriptions, so that `lib.gravity` in the next step reads real files.

**Inputs.** `data_static/harmonics/null1981.toml`, `data_static/harmonics/iess2019.toml`,
`data_static/rotation/system_iii.toml`. The TOML keys are the vocabulary; the files carry a
`status` per value, and the tool copies it into the netCDF as a per-variable attribute.

**Deliverable: `src/casspian/tools/gravity/`**, two entry points, `casspian-gravity-file` and
`casspian-rotation-file`, each taking a TOML control file (SPEC_00 §2.2.2) that names the
`data_static` source, the output path, and the prefix. For this step the control files are
`occul_data/lindal/lindal_build.toml` (sections `[gravity]` and `[rotation]` only; the other
sections are added by later steps) and the outputs are `occul_data/lindal/lindal_gravity.nc`
(from `null1981`) and `occul_data/lindal/lindal_rotation.nc` (from `system_iii`). Kind G per
SPEC_00 §6.4 with `harmonic_convention = "CASSPIAN-J1"`, the code defined in `lib.schema`; kind R
per §6.5. The planet GM goes into kind G; the system GM and the derivation note go into an
attribute.

**Acceptance.** `lindal_gravity.nc` reads back with `J` = [16479e-6, −937e-6, 84e-6] on
`degree` = [2, 4, 6], `J_status` = [fitted, fitted, assumed], `normalization_radius_m` =
6.0e7, `GM_m3s2` = 3.7929085e16; a second run pointed at `iess2019` gives six degrees to 12 and
`normalization_radius_m` = 6.033e7; `lindal_rotation.nc` gives `period_s` = 38362.4 and
`angular_rate_rad_s` = 2π/38362.4 to double precision; the reader refuses a G file whose
`harmonic_convention` attribute is a code it does not know.

---

## 5. Step 4: `lib.gravity`

**Purpose.** The Newtonian potential and the effective gravity components. Implements
manuscript Eqs. A2, A3, A4, A5 and handoff §9A.5.

**Deliverable: `src/casspian/lib/gravity.py`.** Pure functions, NumPy in and out, every input
an argument:

- `legendre_even(degrees, x)` and its derivative with respect to `x`, for even degrees to 12,
  vectorized over `x`.
- `potential_V(r, phi_c, GM, J, degrees, R_norm)`:
  `V = −(GM/r) [1 − Σ J_ℓ (R/r)^ℓ P_ℓ(sin φ_c)]`, ℓ over the even degrees supplied.
- `g_newton(r, phi_c, GM, J, degrees, R_norm)`:
  `g_N = (GM/r²) [1 − Σ (ℓ+1) J_ℓ (R/r)^ℓ P_ℓ(sin φ_c)]`, positive inward.
  **The coefficient is ℓ + 1 with ℓ the degree: 3 on J2, 5 on J4, 7 on J6.** This was misread
  once already (data_static/harmonics/null1981.toml, `check_values`).
- `G_phi_newton(r, phi_c, ...)`: `−(1/r) ∂V/∂φ = −(GM/r²) Σ J_ℓ (R/r)^ℓ cos φ_c P'_ℓ(sin φ_c)`.
- `omega_abs(u, r, phi_c, Omega)`: Eq. A2.
- `g_eff_radial(u, r, phi_c, Omega, GM, J, degrees, R_norm)`: Eq. A3 (equivalently A4),
  positive inward.
- `G_phi_eff(u, r, phi_c, Omega, GM, J, degrees, R_norm)`: Eq. A5.
- `g_eff_vector(...)` returning both, the magnitude, and the angle ψ = arctan2(−G_φ, g).

**Sign convention for `G_φ` (ruled at the Step 4 review, v0.6).** `G_φ` is the component of
the effective gravity along increasing planetocentric latitude, exactly as the formula
`−(1/r) ∂V/∂φ` defines it. It is negative in the northern hemisphere and positive in the
southern, the bulge and the centrifugal term both pulling toward the equator. It is never
redefined as "equatorward positive": a component along a coordinate direction keeps one
meaning in both hemispheres, and Eq. B3, `g dr₀/dφ = r₀ G_φ`, holds as written with this
convention (dr₀/dφ < 0 north of the equator) and would need a hemisphere dependent sign with
the other. The angle ψ = φ_g − φ_c is the tilt of the local vertical from the radial
direction toward the pole, positive in the north, so ψ = arctan2(−G_φ, g); the minus sign is
part of the definition of ψ, not a correction. Handoff §9A.5 and the manuscript's statement
of Lindal Eq. 5 are to be read with `G_φ` in this convention; the manuscript is to be checked
for consistency (manuscript notes list).

All functions accept scalars or broadcastable arrays. No file I/O in this module.

**Acceptance.** With Null's set at 60,000 km, the planet GM read from `iess2019.toml`
(3.7931206234e16), System III, and the 1 bar radii from `lindal_scalars.toml` (60,268 and
54,364 km): `g_eff_radial` at the equator with u = 450 m/s is 8.951 ± 0.005 m/s² and with
u = 0 is 9.102; at the pole 12.137 (Lindal Table II: 8.96, 12.14). With Null's GM
(3.7929085e16) the same to within 1e-3 m/s² (the two GM values differ by 5.6e-5 relative, so
the shift is 6e-4 at the equator and 6.8e-4 at the pole; v0.5 said 5e-4, which was set without
allowing for the polar value). `G_phi_eff` is zero at the equator and the pole to round-off and
negative (equatorward, by the convention above) at 45° north on a sphere of the 1 bar
equatorial radius with u = 0, value reported; on that sphere it is −1.056 m/s² and ψ is 6.35°,
to be read as a sign and an order of magnitude, since the geoid radius at 45° is Step 5's
product and this check is to be repeated on the constructed surface once it exists. Jupiter
check optional but recommended: with J2 = 14736e-6, J4 = −587e-6, J6 = 31e-6 at 71,398 km,
GM = 1.26686534e17, period 9h 55m 29.7s, radii 71,492 and 66,854 km, u = 100 m/s at the
equator: 23.116 and 27.015 (Lindal: 23.12, 27.01). A deliberately wrong coefficient (5, 9, 13)
must give 9.29 at the equator with u = 0 (against the correct 9.102; with u = 450 the wrong
value is 9.14 against 8.951) and 11.68 at the pole (against 12.137), demonstrating that the
check discriminates. Every Saturn number in the acceptance script is read from the static
transcriptions; the Jupiter values, which appear in no static file, are quoted from this
paragraph and labeled as such.

---

## 6. Step 5: `lib.geoid`

**Purpose.** The reference surface: Lindal Appendix Eqs. 11 to 14, and manuscript Eq. B3.

**Deliverable: `src/casspian/lib/geoid.py`.** Pure functions:

- `U_rigid(r, phi_c, Omega, GM, J, degrees, R_norm)`: Lindal Eq. 11, the potential of a
  uniformly rotating fluid, `U = V − ½ Ω² r² cos² φ_c`, with `g_eff = −grad U` and U → 0 on
  the axis at infinity. This is Lindal's Eq. 11 term for term and in the same sign (v0.6 said
  his form carried the opposite overall sign; it does not, and the docstring is not to claim
  it). With `g` positive inward, `dU/dr = g`.
- `reference_geoid(phi_c_grid, r_polar, Omega, GM, J, degrees, R_norm, tol_m, max_iter)`:
  Lindal Eqs. 12 to 14. Set `U_ref = U_rigid(r_polar, π/2)`; at each latitude iterate
  `r ← r − (U(r, φ) − U_ref) / g` from an ellipsoid seed until `|Δr| < tol_m` (Lindal Eq. 14
  is `Δr = ΔU / g_r` with his outward `g_r = −g`; same update). Returns `r_ref(φ_c)`, the
  iteration count per latitude, and the final `|ΔU|`. This is the **no-wind** reference geoid.
  It is solved independently at every latitude supplied, so any latitude at which the radius
  is wanted is to be put in the grid and solved, never interpolated.
- `wind_geoid(phi_c_grid, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm)`: the surface a
  latitude dependent, altitude independent wind produces (Lindal's assumption), which is the
  manuscript's Eq. B3, `g dr₀/dφ = r₀ G_φ`, integrated from each pole inward with a
  fourth-order scheme on the supplied latitude grid using the Step 4 components with `u(φ)`.
  Since with wind the field is not conservative when u varies with latitude (§A6), this surface
  is the level set of nothing, and no potential closure is defined for it. The function returns
  the radius and, as a diagnostic, `U(r₀(φ), φ) − U_ref` against the **no-wind** potential
  through the same polar radius, which is the dynamical height of the wind surface above the
  no-wind geoid in potential units (divide by `g` for meters). No pseudo-potential built from
  `Ω_abs` is to be evaluated or reported; its variation along the surface restates the wind
  kinetic term and measures nothing. The march accepts any grid, so a latitude at which the
  radius is wanted (the anchor latitude above all) is to be inserted as a node, never
  interpolated. No `tol_m` or `max_iter` arguments: an initial value problem does not iterate.
- `radius_at(phi_c, phi_c_grid, radii)`: interpolation of a constructed geoid to a latitude,
  for the vectorized uses of Step 7 and later, never for the frozen anchor radius. Scheme
  (ruled at the Step 5 review, v0.7): bracket the target between its two neighboring grid
  nodes in `φ_c`, then interpolate `1/r²` linearly in `sin² φ_c` between them. For an ellipse
  `1/r² = sin² φ / b² + cos² φ / a²` holds exactly, so the scheme interpolates only the
  harmonic and wind residual; measured error at the anchor latitude is eight times smaller
  than linear in `sin φ_c` (v0.6) at equal spacing. Because `sin² φ_c` is not monotonic across
  the equator, the bracketing is done in `φ_c`, and an interval that straddles the equator is
  refused; a grid is to contain the equator as a node. Refuses to extrapolate.

**Acceptance.** No-wind reference geoid with Null's set and Null's GM (the reduction uses
the GM of `lindal_gravity.nc`, which is Null's; the handoff numbers below were evidently
computed with the modern GM, and the two differ by 0.2 km, inside the tolerance), System III,
`r_polar` = 54,438 km: equatorial radius 60,244 ± 2 km and `r(31.0°)` = 58,435 ± 2 km
(handoff §9A.6); `r(30.8185°)` = 58,452.9 ± 2 km (handoff §9A.8); the surface is an
equipotential to a stated residual. The 45° check of Step 4 repeated on this surface, value
reported. Wind geoid: grid convergence of the equatorial radius under refinement is the
acceptance; a stand-in `u(φ)` does not determine the bulge (profiles equally consistent with
"450 m/s at the equator falling to zero by 35°" give 66 to 215 km), so no bulge figure is
prescribed here. The comparison of the wind geoid equatorial radius with Lindal's fitted
60,367 ± 4 km belongs to Step 7, with the real Smith wind, where the number is to be reported,
not prescribed.

---

## 7. Step 6: `lib.latitude`

**Purpose.** The planetographic to planetocentric conversion, handoff §9A.2, Lindal Eq. 5.

**Deliverable: `src/casspian/lib/latitude.py`.**

- `planetocentric_from_ellipsoid(phi_g, flattening)`: the seed, `tan φ_c = (1−f)² tan φ_g`.
- `planetocentric_fixed_point(phi_g, surface, u_of_phi, Omega, GM, J, degrees, R_norm,
  tol_rad, max_iter, flattening)`: iterate `φ_c ← φ_g − ψ(φ_c)` with ψ from Step 4 evaluated at
  `r = surface(φ_c)`, where `surface` is a Step 5 construction solved at the current
  iterate's latitude (for the no-wind geoid, a direct Newton solve at that latitude; for the
  wind geoid, a march with that latitude inserted as a node), not an interpolated radius
  (Step 5 review, v0.7). Seeded from the ellipsoid `tan φ_c = (1−f)² tan φ_g` with `f` the
  100 mbar oblateness of the geodesy scalars (0.09822), the surface the Fig. 4 label refers
  to; the converged value does not depend on the seed, only the count does. Every angle in
  `lib` is in radians, the tolerance included (`tol_rad`; v0.7 said `tol_deg`, the only
  degree valued quantity in `lib`, corrected at the Step 6 review). Returns `φ_c`, ψ, the iterates, and the count. Vectorized over `phi_g` for the
  wind tool.
- The inverse, `planetographic_from_planetocentric`, which is a direct evaluation.

**Acceptance.** For `phi_g` = 36.3 on the no-wind reference geoid of Step 5 with Null's set
and Null's GM: `φ_c` = 30.8185 ± 0.001°, ψ = 5.4815 ± 0.001° (Null's GM gives 30.81819 and
5.48181, the modern GM 30.81846 and 5.48154; the handoff numbers were computed with the
modern GM, as at Step 5, and the reduction uses Null's), converging from the ellipsoid seed
30.8524 (handoff: 30.8526) in three iterations to 1e-3°, four to 1e-4°, and no more than
eight to 1e-9° (handoff §9A.8 lists the iterates). Round trip: the inverse applied to the
result returns 36.3 to 1e-9°. For `phi_g` = 36.5 the result is 31.005 (handoff §9A.2); the
0.2° change in the label moves `φ_c` by 0.187°, so the label's uncertainty passes through the
conversion almost one for one and the frozen anchor latitude is known no better than the
figure label (the swath in `lindal_scalars.toml` is the record of that). With the Step 7 wind
included, the value **moves**, by an amount to be reported (the Step 5 stand-in moves it
0.016° equatorward); that value, not 30.8185, is the one the reduction will freeze.

---

## 8. Step 7: the wind tool for Lindal, kind W

**Purpose.** `occul_data/lindal/lindal_wind.nc`: the cloud-top wind Lindal carried in his
geoid, on the model's axes, with its uncertainty.

**What this wind is for, and is not for.** In the reduction phase the wind enters exactly one
computation: the effective gravity on the anchor isobar, for `phi_c` (Step 6) and `r0`
(Step 5). Lindal states that he assumed the wind independent of altitude (Appendix p. 1144).
The reduction wind file therefore carries zero shear, not as a simplification but as a
transcription of the source's assumption, which a different assumption would change through
`r0`. Kind W is two-dimensional in general and the forward-phase wind will carry real vertical
structure (cloud-top profile at its level, sheared component from thermal retrievals through
Eq. A39, parameterized decay below); that is the generic wind tool and a later spec. This
step builds only the Lindal instance, using the same file layout.

**Inputs.** `data_static/winds/smith1982_fig4_points.csv` (323 digitized points, planetographic
latitude, System III; see the `.note.md` beside it), `lindal_gravity.nc`, `lindal_rotation.nc`,
the `[wind]` section of `lindal_build.toml`, which declares: the source file, the bin width in
degrees (2.0), the minimum count per bin below which a bin is flagged (3), the observation level
in Pa assigned to the cloud-top wind and its justification string, the interpolation rule for
empty bins (linear in latitude between the nearest populated bins, flagged `interpolated`),
the extension rule poleward of the last populated bin (hold the last bin mean, flagged
`extrapolated`, to ±90°), `vertical_structure = "altitude independent"`, and the pressure grid
on which the file is written (a declared list, for example ten levels per decade from 1 Pa to
1.0e6 Pa), which must cover every pressure any reduction or run will touch.

**Deliverable: `src/casspian/tools/wind/`** with entry point `casspian-wind-from-points`,
which:

1. Bins the points in planetographic latitude at the declared width; per bin computes count,
   mean, and sample standard deviation (NaN for count 1); flags bins below the minimum count.
2. Converts bin-center planetographic latitudes to planetocentric by `lib.latitude` on the
   no-wind reference geoid of Step 5 built from the named G and R files. (The wind itself is
   what will later make the geoid depart from no-wind; using the no-wind surface for the
   conversion is a second-order choice and is recorded in `latitude_conversion_inputs`.)
3. Fills empty bins and extends to the poles per the declared rules, flagging each value in
   `value_provenance` (0 observed, 1 interpolated, 2 parameterized, 3 extrapolated).
4. Extends the one-level profile onto the declared pressure grid as a full two-dimensional
   array, every column identical, with `value_provenance` at every level other than the
   observation level set to a fifth code, 4, `extended_by_source_assumption`, so the
   assumption is visible in the data and the file has the same shape as every other kind W.
5. Writes kind W (SPEC_00 §6.6): `u_total_ms(latitude, pressure)` as the binned mean,
   `u_total_uncertainty_ms` as the bin standard deviation with `uncertainty_kind = "1sigma"` and
   `uncertainty_method = "sample standard deviation of the digitized points in the bin"`,
   `bin_count` as an extra variable, `reference_level_pressure_Pa` equal to the observation
   level, `u_cylindrical_ms` computed from the field at the reference level extended along
   cylinders on the no-wind geoid geometry, and `u_shear_ms = u_total_ms − u_cylindrical_ms`.
   For an altitude-independent field on a spherical-shell geometry the shear component is
   not identically zero (a column at fixed latitude tilts away from a cylinder), and the small
   nonzero values are the correct decomposition of Lindal's assumption; report their maximum.
   Provenance attributes per §6.6, including `source_latitude_convention = "planetographic"`
   and the citation and caption quoted in the `.note.md`.

**Acceptance.** Bin table reported in full (latitude, count, mean, std, flag). Sum of counts
is 323. The equatorial bins (centers near 3 to 9° N) have means between 440 and 495 m/s. The
bin containing 38° N planetographic has a mean between −20 and +10 m/s and a count of at
least 10. Planetocentric conversion of the 36.3° bin center matches Step 6 to 0.01°. Every
column of `u_total_ms` is identical. `read` of the file as kind W succeeds and the SPEC_00
§6.6 load-time checks (sum identity to round-off; constancy of `Omega_abs` from
`u_cylindrical_ms` on cylinders to the declared tolerance) pass; the maximum of
`|u_shear_ms|` is reported. The report includes a plot of the binned mean with ±std against
latitude over the digitized points. With the binned wind as `u(φ)`, the Step 5 wind geoid
through `r_polar` = 54,438 km is run and its equatorial radius reported and compared with
Lindal's fitted 60,367 ± 4 km (moved here from Step 5 at v0.7); the number is to be reported,
not prescribed, and a departure is a finding about the wind, the anchoring, or Lindal's fit,
to be discussed rather than tuned away.

---

## 9. Step 8: the composition tool for Lindal, kind C

**Purpose.** `occul_data/lindal/lindal_composition.nc`: the reduction composition with the
species group carrying the `lindal1985` refractivity set.

**Inputs.** The raw bundle (Step 2) for the NH3 column and pressure levels;
`data_static/species_master.toml`; the `[composition]` section of `lindal_build.toml`, which
declares: the species set name (`lindal1985`), the H2 and He fractions of the remainder
(0.94 and 0.06 as the ratio to hold), the NH3 rule (interpolate and extrapolate linearly in
pressure with `scipy.interpolate.interp1d(fill_value="extrapolate")`, clamp at zero, then
assign NH3 first and split the remainder 0.94 to 0.06), and the vertical coordinate
(`pressure`).

**Deliverable: `src/casspian/tools/composition/`** with entry point
`casspian-composition-lindal` (the generic generator for forward compositions comes later and
is out of scope here), writing kind C per SPEC_00 §6.2 with `composition_role = "reduction"`,
`source_statement` quoting the Table I footnote, per-species `provenance` (`assumed` for H2
and He, `measured` where NH3 is tabulated and `interpolated` or `extrapolated` elsewhere),
`x_H2_uncertainty` = 0.03, `latitude_absent_meaning = "point"` with the profile latitude in
the source convention beside it (SPEC_00 §5), and the `species` group with molar masses and the `lindal1985`
refractivities (H2 136, He 35, NH3 0 with the caption note) converted to m³ per molecule by
the Loschmidt constant of `lib.constants`, with each value's `status` copied from the master
table and the master table's hash recorded.

**Acceptance.** Mole fractions sum to one at every level to 1e-12. NH3 is zero at and above
794.33 mbar (the clamp), 15.9 ± 0.1 ppm at the 1047.13 mbar interior gap, and the
extrapolated value at 1298.48 mbar is reported (Phase 1 found 79.3 ppm). Mean molar mass with NH3 set to zero is 2.1351 amu
(Phase 1 check) and the value with NH3 at the deepest level is reported. Mean molecular
refractivity at a level with zero NH3 is (0.94 × 136 + 0.06 × 35) × 1e-6 / n_Loschmidt to
six figures. `read` as kind C succeeds and refuses when one mole fraction is perturbed by
1e-6.

---

## 10. Step 9: stage two of the Lindal tool, kinds T and D and the manifest

**Purpose.** Complete `occul_data/lindal/` so that `refrac` (SPEC_02) can run.

**Inputs.** The raw bundle; `lindal_build.toml`, now with all sections.

**Deliverable: `src/casspian/tools/lindal/build_inputs.py`**, entry point
`casspian-lindal-inputs`, which:

1. Writes `lindal_thermo.nc`, kind T (SPEC_00 §6.1): `pressure_Pa`, `temperature_K`,
   `height_m` from the raw bundle's `table1` group, the three uncertainty companions
   (`pressure_uncertainty_Pa`, `temperature_uncertainty_K`, `height_uncertainty_m`) present and
   NaN, and
   every global attribute of §6.1 filled from the `scalars` groups (latitude with convention
   in the name, its `value_source`, swath, uncertainty, longitude, date, bands, datum,
   source top boundary statement, source gravity, rotation, and wind citations), with
   `latitude_absent_meaning = "point"` and `thermo_instance = "source_profile"` (SPEC_00 §5,
   §6.1). **No composition and no geodesy in this file.**
2. Writes `lindal_geodesy.nc`, kind D (SPEC_00 §6.3): the two surfaces from the `scalars`
   geodesy groups, with the fit residual, `fit_inputs`, and `fit_latitude_convention`.
3. Calls the Step 3, 7, and 8 tools (or verifies their outputs exist and match the hashes in
   the control file) so that G, R, W, C are present.
4. Writes `lindal_reduction.toml` (SPEC_00 §7.1) pointing at the six files, with the anchor
   isobar at 1.0e4 Pa, the fixed-point and geoid tolerances from the control file, and the
   product name.
5. Records the raw bundle hash and the control file hash in every file it writes.

**Acceptance.** All six files validate under their kinds; `lindal_thermo.nc` has exactly the
three profile variables and their three uncertainty companions and no `x_*` or geodesy
content; `latitude_planetographic_deg = 36.3` with `value_source` present;
`lindal_geodesy.nc` carries 60,367 / 54,438 / 0.09822 on the 1e4 Pa surface and 60,268 /
54,364 / 0.09796 on the 1e5 Pa surface with their uncertainties; the manifest parses and
every path in it exists; the directory listing matches SPEC_00 §2.2 exactly, with every file
prefixed `lindal_`.

---

## 11. What comes after

SPEC_02 specifies `refrac` (SPEC_00 §3.3): B3.1, the frozen `phi_c` and `r0` using the
wind-included geoid, B1, B3.3, and the kind N product with embedded inputs. Its first
acceptance number will be the wind-included `phi_c` reported in Step 6 and the wind-included
`r0` reported in Step 5. SPEC_03 opens the numerical methods discussion for the forward model.

---

## 12. Decisions made in this document

1. The Phase 1 code is committed once, tagged, and deleted from the working tree; no attic
   directory (§1).
2. The planetographic conversion in the wind tool uses the no-wind reference geoid, with the
   choice recorded in the file (§8, item 2).
3. The Lindal wind is binned at 2° with a minimum count of 3, empty bins interpolated and the
   poles extrapolated by holding the last populated bin, all flagged (§8).
4. The reduction wind is written as a full two-dimensional kind W on a declared pressure grid
   with identical columns flagged as extended by the source's assumption, so that no reader
   special-cases a one-level file (§8, item 4). Its decomposition is computed like any other
   wind file's; the shear component is small but not identically zero (§8, item 5).
5. The NH3 refractivity in the reduction species set is zero, per the Fig. 3 caption (§9).
6. The wind geoid of Step 5 is solved as the Eq. B3 ODE with a potential-closure diagnostic,
   rather than as Lindal's iteration with wind, because with a latitude-dependent wind the
   effective field is not conservative and the two do not agree exactly (§6). The difference
   is itself a useful number and is to be reported.

---

## 13. Revision history

| Version | Date | Change | Cause |
|---|---|---|---|
| 0.1 | 2026-09-10 | First draft: Step 0 reset plus nine steps with acceptance numbers | SPEC_00 v0.2 accepted |
| 0.2 | 2026-09-10 | Step 0 force-adds ignored Phase 1 files before tagging; header dependency corrected; Step 7 rewritten as a two-dimensional reduction wind with flagged extension | coding agent pre-execution review of Step 0; discussion of Steps 7 to 9 |
| 0.3 | 2026-09-10 | Step 1 item 0 (`.gitattributes`); version rule; `raw/notes.md` hashed into the raw bundle | REVIEW_01_step0 |
| 0.4 | 2026-09-10 | `harmonic_convention` code; `uncertainty_kind = "1sigma"` with `uncertainty_method`; `reports/` for acceptance scripts; commit-on-acceptance-then-rebuild rule; `thermo_instance` on the Lindal kind T | REVIEW_01_step1 |
| 0.5 | 2026-09-10 | Uncertainty companion names spelled out in Step 9; Step 2 housekeeping commit moving reports and reviews to `reports/`; this section added | REPORT_01_step1 §6; author request |
| 0.6 | 2026-09-10 | Step 4: `G_φ` sign convention ruled (along increasing latitude; ψ = arctan2(−G_φ, g)); Null GM tolerance 1e-3; wrong coefficient values paired with their no-wind counterparts; GM read from `iess2019.toml`; the 45° check placed on the equatorial 1 bar sphere with a repeat on the Step 5 surface | REPORT_01_step4 §3 |
| 0.7 | 2026-09-10 | Step 5: `radius_at` scheme changed to `1/r²` linear in `sin² φ_c` with bracketing in `φ_c`; anchor radius solved or marched, never interpolated; closure diagnostic fixed as the no-wind potential departure, pseudo-potential excluded; stand-in bulge figure dropped, 60,367 comparison moved to Step 7; Null GM stated for the reduction geoid; `U_rigid` sign claim about Lindal Eq. 11 withdrawn; unused arguments removed from `wind_geoid` | REPORT_01_step5 §3; author direction on the interpolant |
| 0.8 | 2026-09-11 | Step 6: `tol_rad` replaces `tol_deg`; seed flattening stated (100 mbar oblateness); iteration counts tied to tolerances; Null GM values recorded beside the handoff values; label uncertainty pass-through noted | REPORT_01_step6 §3 |
