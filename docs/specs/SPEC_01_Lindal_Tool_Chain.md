# SPEC 01. The Lindal tool chain, the library it needs, and the repository reset

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.22, 14 September 2026. Author of record: S. Rafkin. Status: closed at commit ece58d2; Steps 2, 8 and 9 amended at v0.21, the snap rule restated at v0.22 before Step 0 began (the Table I pressure grid: `pressure_printed_Pa` beside `pressure_Pa` on the declared `10^(k/100)` mbar grid, applied at stage one, read by every downstream tool) for SPEC_03 Step 0, with a rebuild of the whole chain; Step 8 amended at v0.19 (closure declaration, NaN per-molecule uncertainties) with a rebuild of the composition product; Steps 5 and 9 amended at v0.20 (the `equatorial_radius` anchor rule in `lib.geoid`; `[stage_two]` carries the anchor quantity and rule and the `[diagnostics]` keys into the manifest) for SPEC_02 Step 6. See §13 for the revision history.
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.15, which it does not repeat.

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

**Figures that a report refers to survive the step (v0.10).** `reports/step<N>/` is ignored
and deletable, so any figure a report cites is written to `reports/figures/step<N>_<name>.png`
and committed with the report; the step directory keeps only the script and its output.

**Commit on acceptance, then rebuild.** A step's code is committed when its review says
`accepted`. Any data product the step wrote before that commit carries `-dirty` in
`casspian_git_commit`; after the acceptance commit **every** product under `occul_data/` that
carries `-dirty`, whichever step wrote it, is rebuilt so that every product on disk carries
the clean commit, and only those products are inputs to the next step (SPEC_00 §8, v0.9: the
rule is a property of the directory, since a later step can rerun an earlier tool in
passing).

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
- **v0.21, the pressure grid (SPEC_03 Step 0 deliverable 1).** `lindal_scalars.toml` carries a
  `[pressure_grid]` table declaring that Table I's pressures lie on `10^(k/100)` mbar
  (`rule`, `denominator = 100`, `unit = "mbar"`, `excluded_printed_values_mbar = [1298.48]`,
  `basis`, `value_source` stating that this is an inference, not a statement of the source).
  `table1` then carries `pressure_printed_Pa` (the CSV value in Pa, `value_source = "Table
  I"`) and `pressure_Pa` on the grid, by the rule of SPEC_03 v0.4 Step 0 (v0.22): for each row
  the set of integers `k` whose grid value, rounded to the CSV's printed decimals, equals the
  printed value; one match fixes the row; several matches take the one equally spaced in `k`
  with the fixed neighbors (the top row from the two rows below it); still ambiguous, or no
  match outside the exclusion list, is a refusal with the row named; an excluded row keeps its
  printed value; the per-level flag `pressure_grid_applied` says which. No spacing numbers in
  the code; the acceptance checks the resulting `k` sequence. For Table I, 62 rows have one
  match and three (0.20, 0.25, 0.32 mbar) have two. Both columns stay in the bundle, so nothing
  printed is dropped.

**Acceptance.** 66 levels; the 1 bar row reads exactly 100000.0 Pa, 134.8 K, 0.0 m; the top
row 19.9526 Pa on the grid with `pressure_printed_Pa` 20.0 Pa (v0.21; then 25.1189, 31.6228,
39.8107, 50.1187, 63.0957, 79.4328, 100.0 Pa), 138.7 K, 376700.0 m; 65 rows snapped and the
bottom row as printed; the `k` sequence −70, −60, ... to 130 by 10, then by 8 to 194, by 6 to
230, by 4 to 270, by 2 to 310; a CSV copy with one printed value edited to something no grid
value rounds to is refused with the row named; the bottom row 129848.0 Pa, 146.2 K, −14100.0 m; exactly
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
- `wind_geoid(phi_c_grid, r_anchor, anchor_rule, u_of_phi, Omega, GM, J, degrees, R_norm)`:
  the surface a latitude dependent, altitude independent wind produces (Lindal's assumption),
  which is the manuscript's Eq. B3, `g dr₀/dφ = r₀ G_φ`, integrated as **one** march from
  the north pole through the equator to the south pole with a fourth-order scheme on the
  supplied latitude grid using the Step 4 components with `u(φ)`. Eq. B3 is first order and
  has one constant; a surface marched from each pole separately with the same polar radius
  is two surfaces, and with a wind that is not symmetric about the equator they do not meet
  (REPORT_01_step7 v0.15, Figure 2: a 38 km step at the equator). The constant is fixed by a
  declared `anchor_rule`: `"mean_polar_radius"` (default; Lindal states a **mean** polar
  radius, and the two polar radii of the marched surface are made to average to `r_anchor`
  by a one-parameter root find on the north polar start), `"north_pole"` or `"south_pole"`
  (`r_anchor` applied at that pole), or `"latitude"` with a declared latitude (for a surface
  anchored at an observed radius). The function returns the radius, the two polar radii and
  their difference (the polar asymmetry the wind produces), and the anchoring residual.
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

**Inputs.** `data_static/winds/ingersoll_pollard1982_fig5_curve.csv` (v0.14: the solid zonal
wind curve of Ingersoll and Pollard 1982, Fig. 5, digitized at 0.2° from +81.1° to +1.3° and
−10.9° to −72.8° planetographic; see its `.note.md`), which is the wind;
`data_static/winds/smith1982_fig4_points.csv` (323 digitized cloud-tracking points, the data
that curve was drawn through; see its `.note.md`), which is the uncertainty;
`data_static/winds/smith1982_fig4.toml` (the observation level in Pa assigned to the cloud-top
wind with its justification string, and the latitude convention of the sources; properties of
the data, so they live beside it); `lindal_gravity.nc`, `lindal_rotation.nc`; the raw bundle
(for the polar radius of the no-wind geoid until kind D exists at Step 9; declared as a
pointer); and the `[wind]` section of `lindal_build.toml`, which declares choices only: the
source pointers; the bin width in degrees (2.0) and the minimum count per bin (3) for the
uncertainty; the gap rule and the polar rule (below); the planetocentric latitude grid of the
product (a declared spacing, 0.5°, from −90 to +90 with both poles and the equator as nodes);
`vertical_structure = "altitude independent"`; and the pressure grid on which the file is
written (a declared list, ten levels per decade from 1 Pa to 1.0e6 Pa), which must cover every
pressure any reduction or run will touch. A control section may carry such choices; the
principle 2 refusal is for keys that name a physical quantity.

**Why this wind (v0.14).** Lindal's appendix cites Smith et al. (1982) and Ingersoll and
Pollard (1982) as the wind sources. The Ingersoll and Pollard curve is a published smoothing of
the Smith data by its authors' hands, and the Smith points scatter about it with RMS 21.3 m/s
and mean −4.3 m/s, equal to their pooled within-bin scatter (21.4 m/s), so it is the mean curve
of the data to the precision of the digitizations. Using it as the wind removes every fitting
choice from this step (the v0.10 to v0.13 designs, spline and then PCHIP through bin means,
are superseded; PCHIP survives only as the interpolant between digitized samples and as the
polar rule). The Smith bins keep their role as the uncertainty. **A zonal wind that is nonzero
at a pole is a defect, in this file and in every kind W file**; SPEC_00 §6.6 enforces at load
time that `u_total_ms` is exactly zero at ±90°, which must be grid nodes.

**Deliverable: `src/casspian/tools/wind/`** with entry point `casspian-wind-from-curve`
(replacing `casspian-wind-from-points`; the point-based path is retained in the module for a
data set that has no published curve, selected by `wind_source_kind = "points"`), which:

1. Reads the digitized curve and builds `u(φ_g)` as a shape-preserving cubic interpolant
   (PCHIP) through its samples, per segment.
2. Fills the ring gap (−10.9° to +1.3°) by the declared rule `gap_rule = "reflect_north"`:
   the northern segment reflected about the equator, which is Ingersoll and Pollard's own
   dashed curve (their Fig. 3 caption defines it). The strip between the northern end (+1.3°)
   and its reflection (−1.3°) is bridged by a cubic Hermite matched to value and slope at
   both ends, which by the symmetry reduces to `a + b φ²` (v0.15). The reflection does not
   meet the southern segment where the latter begins (at −10.9° the reflection reads about
   377 m/s and the southern solid 346), so the two are blended linearly over a declared
   `join_window_deg = [-10.9, -15.0]`, reflection weight one at the first edge and zero at
   the second, so that the assembled wind is continuous and lies on the southern data from
   −15° southward (v0.15; the window is part of the sensitivity study). Values in the gap
   and the window are flagged `parameterized` (code 2) with the rule named in
   `parameterization`. The `.note.md` records that the ten Smith points on the
   southern flank (−7.5° to −10.9°) lie 80 to 100 m/s below the reflection; the alternative
   `gap_rule = "reflect_north_then_bins"` (reflection only where no Smith bin exists, the
   populated Smith bins otherwise) is a declared choice for the sensitivity study, not the
   default.
3. Brings the wind to zero at each pole by the declared rule `polar_rule = "pchip_to_zero"`:
   PCHIP through the last two digitized samples of the segment and a node of exactly 0.0 at
   ±90°, which is monotone by construction; flagged `extrapolated` (code 3). The alternative
   `polar_rule = "linear_to_zero"` is a declared choice for the sensitivity study.
4. Evaluates the curve on the declared planetocentric grid. Each grid `φ_c` is converted to
   `φ_g` by the **direct** relation `planetographic_from_planetocentric` on the no-wind
   reference geoid built from the named G and R files (no fixed point is needed in this
   direction). The choice of the no-wind surface for the conversion is second order and is
   recorded in `latitude_conversion_inputs`.
5. Provenance per grid latitude: 0 `observed` where the value comes from the digitized curve
   (it is Ingersoll and Pollard's reading of observations; the `provenance` attribute of the
   variable says `derived` and the source); 2 `parameterized` in the gap; 3 `extrapolated` in
   the polar caps and at the poles (value exactly 0.0).
6. Uncertainty: the Smith points are binned at the declared width in planetographic latitude
   and, per populated bin with at least the minimum count, the RMS of the points about the
   curve is the 1σ (`uncertainty_kind = "1sigma"`, `uncertainty_method = "RMS of the Smith
   et al. 1982 Fig. 4 points about the Ingersoll and Pollard 1982 Fig. 5 curve within a 2
   degree planetographic bin; NaN where the bin has fewer than the minimum count"`); assigned
   to grid latitudes by linear interpolation in `φ_g` between populated bin centers; in the gap
   the uncertainty is the mirrored northern value where the reflection supplies the wind and
   the bin RMS where Smith bins exist; NaN in the polar caps and at the poles. The bins are
   written as `bin_rms_ms` and `bin_count` on a `bin_latitude` coordinate (planetographic,
   stated).
7. Extends the one-level profile onto the declared pressure grid as a full two-dimensional
   array, every column identical, with `value_provenance` at every level other than the
   observation level set to code 4, `extended_by_source_assumption`.
8. Writes kind W (SPEC_00 §6.6): `u_total_ms(latitude, pressure)`, `u_total_uncertainty_ms`,
   `value_provenance`, the bin variables, `reference_level_pressure_Pa`. Decomposition: for a
   field declared altitude independent on a pressure coordinate, `u_cylindrical_ms =
   u_total_ms` and `u_shear_ms = 0` identically, with `decomposition =
   "trivial_altitude_independent"`. Provenance attributes per §6.6, including
   `source_latitude_convention = "planetographic"`, `wind_source = "ingersoll_pollard1982_fig5"`,
   and the citations and captions quoted in both `.note.md` files.

**Acceptance.** The wind curve passes through every digitized sample outside the join window
to round-off (inside it the blend replaces the samples by declaration) and is exactly 0.0 at
both poles; the §6.6 polar check passes; a copy of the file with one polar value
perturbed to 1e-6 m/s is refused by `read`. Peak 490.5 ± 2 m/s at +7.4 ± 0.3° planetographic;
the reflected peak at −7.4°; the value at `φ_g` = 36.3° is 1.9 ± 3 m/s (this value enters the
frozen anchor and is to be quoted in the report); the value at +30.8° is 75 ± 3 m/s. The tool's
conversion of `φ_g` = 36.3° matches Step 6 to 1e-6°. Bin table reported in full (bin latitude,
count, RMS about the curve, flag); sum of counts is 323; the RMS over all populated bins,
pooled, is 21.3 ± 0.5 m/s. Every column of `u_total_ms` is identical. `read` as kind W succeeds
and the §6.6 load-time checks pass. Figure 1 (to `reports/figures/step7_wind_curve.png`, committed): the Smith points, the
Ingersoll and Pollard curve with its ±1σ band, the reflected gap fill and the polar
extrapolations distinguished, and the provenance ranges shaded, against planetographic
latitude. With the wind curve as `u(φ)`, the Step 5 wind geoid is run as a single march
with `anchor_rule = "mean_polar_radius"` and `r_anchor` = 54,438 km (v0.16); its equatorial
radius is reported and compared with Lindal's fitted 60,367 ± 4 km, remembering that the
anchor itself carries ±10 km which passes straight into the equatorial value, so the
comparison band is about ±11 km (an independent march with these rules gives 60,371 km,
north and south polar radii 54,424 and 54,452 km, a 28.7 km asymmetry; these are the
expected figures, to be reproduced, not prescribed); grid convergence across an eightfold
refinement is reported; the sensitivity of the equatorial radius to the anchor rule
(`north_pole` and `south_pole` against the default), to the polar rule (`pchip_to_zero`
against `linear_to_zero`) and to the gap rule (`reflect_north` against
`reflect_north_then_bins`) is reported, each as a difference in km. The number is to be reported, not prescribed, and
a departure beyond the band is a finding about the wind, the anchoring, or Lindal's fit, to
be discussed rather than tuned away. Figure 2 (`reports/figures/step7_dynamical_height.png`,
committed): the wind curve above, and below it the dynamical height
`h(φ) = r_wind(φ) − r_nowind(φ)` from the Eq. B3 march over all latitudes, overlaid with
Lindal's small-wind approximation, his Eq. 18,
`h(φ) ≈ (2 ω r_p / ⟨g⟩) ∫_φ^{π/2} u sin φ dφ`, integrated continuously from the north pole
to the south pole and offset so that its two polar values average to zero, the same
anchoring as the march (v0.16), evaluated from the same wind curve with
`⟨g⟩` taken as the no-wind gravity on the reference geoid at that latitude; the anchor
latitude is marked and the equatorial value annotated against the 123 km implied by
60,367 − 60,244. The march and Eq. 18 are expected to agree to within the cyclostrophic term
and the oblateness corrections Eq. 18 drops; a larger departure is to be explained.

---

## 9. Step 8: the composition tool for Lindal, kind C

(v0.21: the tool reads `table1/pressure_Pa`, the grid values, so kind C's levels are kind T's.)

**Purpose.** `occul_data/lindal/lindal_composition.nc`: the reduction composition with the
species group carrying the `lindal1985` refractivity set.

**Inputs.** The raw bundle (Step 2) for the NH3 column and pressure levels;
`data_static/species_master.toml`; the `[composition]` section of `lindal_build.toml`, which
declares: the species set name (`lindal1985`), the H2 and He fractions of the remainder
(0.94 and 0.06 as the ratio to hold), the NH3 rule, and the vertical coordinate
(`pressure`).

**The NH3 rule (v0.17; the v0.5 wording, "interpolate and extrapolate linearly in pressure,
clamp at zero", produced a value of 0.2 ppm at 794.33 mbar and a puzzle about where the
clamp fell, REPORT_01_step8 finding 1).** Lindal's Table I lists NH3 at nine levels from
831.76 to 1258.93 mbar, inferred from the S and X band absorption where absorption was
measurable; above 831.76 mbar the table is blank, which is "not measured", not "zero". The
rule has three parts. Interior gap (1047.13 mbar): linear interpolation in pressure between
the neighboring tabulated values, `interpolated`. Below the deepest tabulated level (1298.48
mbar): linear extrapolation in pressure from the two deepest tabulated values,
`extrapolated`. Above the highest tabulated level (every level at or above 794.33 mbar):
**zero, by assumption**, `assumed`, with the reason recorded in an attribute
(`nh3_aloft_assumption`): Lindal states the troposphere is saturated with ammonia at these
levels, and the saturation mixing ratio over NH3 ice on his own temperatures falls from a
few ppm near 0.8 bar to below 1e-7 by 0.5 bar and below 1e-12 at the tropopause, so the
assumption costs less than 1e-6 of the mean molar mass and nothing in the refractivity,
where the `lindal1985` set carries NH3 at zero. No upward extrapolation, no clamp. NH3 is
assigned first and the remainder split 0.94 to 0.06, so the three fractions sum to one
identically.

**Deliverable: `src/casspian/tools/composition/`** with entry point
`casspian-composition-lindal` (the generic generator for forward compositions comes later and
is out of scope here), writing kind C per SPEC_00 §6.2 with `composition_role = "reduction"`,
`source_statement` quoting the Table I footnote, per-species `provenance` (`assumed` for H2
and He, `measured` where NH3 is tabulated and `interpolated` or `extrapolated` elsewhere),
`x_H2_uncertainty` = 0.03 read from the raw bundle (a stated value is data, not a control
choice), `latitude_planetocentric_absent_meaning = "point"` (the SPEC_00 §5 form
`<dimension>_absent_meaning`; v0.5 wrote the shorthand) with the profile latitude in the
source convention beside it, a per-level `nh3_provenance` flag beside the variable-level
`provenance` (same pattern as kind W's `value_provenance`), and the `species` group with molar masses and the `lindal1985`
refractivities (H2 136, He 35, NH3 0 with the caption note) converted to m³ per molecule by
the Loschmidt constant of `lib.constants`, with each value's `status` copied from the master
table and the master table's hash recorded. `is_polar` is a property of the molecule and is
read from the `[is_polar]` table of the master file (v0.17; it was per refractivity set and
absent from the `lindal1985` NH3 entry, which would have written ammonia as nonpolar);
`temperature_dependence` stays per set, since it describes the value.

**Amendment after closure (v0.19, from REPORT_02_step3 findings 4 and 5).** The tool writes
the closure declaration of SPEC_00 §6.2 v0.13, `closure_rule = "share_of_remainder"` and
`closure_species = "H2 He"`, and writes the species group's `refractivity_uncertainty_m3` as
NaN where the master table states no uncertainty (the `lindal1985` set states none), with
`uncertainty_method` saying so, instead of 0.0. `lindal_composition.nc` is rebuilt and its new
hash recorded; nothing else under `occul_data/lindal/` depends on it until kind N exists.
Acceptance of the amendment: the two attributes read back; `refractivity_uncertainty_m3` is
NaN for all three species; every other Step 8 acceptance number is unchanged.

**Acceptance.** Mole fractions sum to one at every level to 1e-12. NH3 is exactly zero at
and above 794.33 mbar, the first grid level above the highest tabulated value, flagged
`assumed`; 2.6 ppm at 831.76 mbar, `measured`; 15.9 ± 0.1 ppm at the 1047.13 mbar interior gap, and the
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
   (v0.20: also the geoid `anchor_quantity` and `anchor_rule` from `geoid_anchor_quantity`
   and `geoid_anchor_rule`, and a `[diagnostics]` section from `diagnostics_figures`,
   `diagnostics_format` and `diagnostics_dpi` when `diagnostics_figures` is present; the
   tool refuses a rule and quantity that do not belong together, as `lib.control` does.)
5. Records the raw bundle hash and the control file hash in every file it writes.
6. (v0.21) Kind T carries `pressure_printed_Pa` and the global `pressure_grid_rule` quoting the
   `[pressure_grid]` declaration (SPEC_00 §6.1 v0.16). The keep-by-content rule of SPEC_02
   Step 6 decision 4 compares `input_hashes` as content, so a changed raw bundle rebuilds T
   and D.

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

Items carried to SPEC_02 from Step 7 (v0.12): a Monte Carlo on the reduction wind, drawing
bin means from their 1σ (independent per bin as the default, with a declared alternative of
correlated draws), rebuilding the wind geoid and the anchor for each draw, and reporting the
spread of `r0`, `phi_c` and the equatorial radius; the treatment of the polar cap, where the
uncertainty is NaN, as a declared alternative (PCHIP decay against linear taper) inside the
same Monte Carlo; and the pass-through of the 0.2° label uncertainty into `phi_c`.

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
| 0.9 | 2026-09-11 | Step 7: two-panel diagnostic figure of the wind and the dynamical height, with Lindal Eq. 18 overlay | author request |
| 0.10 | 2026-09-11 | Step 7 redesigned: fitted wind (constrained penalized spline, zero at both poles, scatter-matching smoothing) as the product, bins retained for the uncertainty; §6.6 polar zero check; observation level moved to a static file beside the digitization; trivial decomposition declared for altitude-independent files; acceptance rewritten (bin criteria restricted to populated bins, anchor wind quoted, ±11 km comparison band, polar cap sensitivity); `reports/figures/` rule in §0 | REPORT_01_step7 §3; author direction |
| 0.11 | 2026-09-11 | Step 7: fit quality rules (hard checks a to d, soft checks e to g, smoothing escalation, PCHIP fallback), stated as general rules for any latitude fit of point data | author request |
| 0.12 | 2026-09-11 | Step 7: default method changed to a shape-preserving interpolant (PCHIP) through the populated bin means with polar zeros; the v0.11 spline and its rules retained as a non-default alternative; uncertainty is the bin sample standard deviation; acceptance simplified accordingly | author direction: simplest method that works |
| 0.13 | 2026-09-11 | Header dependency corrected to SPEC_00 v0.8; alternative spline rules corrected: check (a) local clause applies only where points bracket the latitude, smoothing search runs in the direction the failing check indicates | coding agent findings on v0.11 |
| 0.14 | 2026-09-11 | Step 7: the wind is the digitized Ingersoll and Pollard (1982) Fig. 5 curve, with the Smith points as its uncertainty; gap rule (reflection, their dashed curve) and polar rule declared with alternatives for sensitivity; entry point `casspian-wind-from-curve`; acceptance numbers from the digitization | author direction; new static data set |
| 0.15 | 2026-09-11 | Step 7: equatorial bridge specified (cubic Hermite, `a + b φ²` by symmetry); linear blend of the reflection into the southern segment over a declared join window | author review of the assembled curve |
| 0.16 | 2026-09-11 | Step 5 `wind_geoid`: one continuous march from pole to pole with a declared `anchor_rule` (default mean polar radius), returning the polar asymmetry; Step 7 acceptance: anchor rule sensitivity, Eq. 18 anchored the same way, expected figures from an independent march; sample check excludes the join window | REPORT_01_step7 (v0.15) Figure 2 discontinuity |
| 0.17 | 2026-09-11 | Step 8: NH3 rule restated (zero above the tabulated range by assumption with the saturation reason, no upward extrapolation or clamp; interior interpolation and downward extrapolation kept); `is_polar` moved to a molecular table in `species_master.toml`; `x_H2_uncertainty` from the raw bundle; `<dimension>_absent_meaning` spelled out; per-level `nh3_provenance` | REPORT_01_step8 §3 and author question |
| 0.18 | 2026-09-11 | Step 9 accepted; §0 rebuild rule restated as a directory sweep; header dependency to SPEC_00 v0.9; `gravity_used_by_source` added to `lindal_scalars.toml` by the author (raw bundle and all downstream products to be rebuilt) | REPORT_01_step9 §3 |
| 0.19 | 2026-09-12 | Step 8 amendment: kind C closure declaration attributes; per-molecule refractivity uncertainty NaN where unstated; composition product rebuilt | REPORT_02_step3 findings 4 and 5 |
| 0.22 | 2026-09-14 | Step 2: the grid snap rule restated as SPEC_03 v0.4 states it (equal spacing with fixed neighbors; three ambiguous rows named) | coding agent's pre-execution review of SPEC_03 v0.3 |
| 0.21 | 2026-09-14 | Step 2 amendment: `pressure_printed_Pa` beside `pressure_Pa` on the declared `10^(k/100)` mbar grid, the `[pressure_grid]` transcription and its refusal rule; Step 8 reads the grid values; Step 9: kind T carries the printed values and the rule, keep-by-content includes `input_hashes`; whole chain rebuilt at SPEC_03 Step 0 | SPEC_03 v0.3 decision 2 |
| 0.20 | 2026-09-12 | Step 5 amendment: `lib.geoid.wind_geoid` gains the `equatorial_radius` rule (the `latitude` rule at an exact equator node, under its own name); Step 9 amendment: `[stage_two]` carries the anchor quantity and rule and the `[diagnostics]` keys into the manifest; manifest and product rebuilt at SPEC_02 Step 6 | SPEC_02 v0.8 Step 6 |
