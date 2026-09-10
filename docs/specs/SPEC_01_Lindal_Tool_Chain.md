# SPEC 01. The Lindal tool chain, the library it needs, and the repository reset

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.4, 10 September 2026. Author of record: S. Rafkin. Status: in work; Steps 0 and 1 accepted; v0.4 aligns with SPEC_00 v0.4 (harmonic code, uncertainty vocabulary, reports directory, commit and rebuild rule).
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.4, which it does not repeat.

---

## 0. How this specification is to be worked

**One step at a time.** This document is a sequence of steps, each producing one deliverable
with its own acceptance checks. The coding agent implements one step, runs its acceptance
checks from a script kept under `reports/<step>/` (not committed; SPEC_00 §2), writes `docs/specs/REPORT_01_<step>.md` stating what was built, what deviates from the
spec and why, the acceptance results with the actual numbers, and any question that blocked it,
updates `docs/specs/STATE.md`, and **stops**. The next step begins only after the report has
been reviewed. Do not implement ahead. Do not build two steps in one commit.

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
- `g_eff_vector(...)` returning both, and the magnitude and the angle ψ = arctan(G_φ / g)
  (Lindal Eq. 5 with the sign convention of handoff §9A.5: Lindal's outward `g_r` is `−g`).

All functions accept scalars or broadcastable arrays. No file I/O in this module.

**Acceptance.** With Null's set at 60,000 km, GM = 3.7931206e16, System III, and the 1 bar
radii from `lindal_scalars.toml` (60,268 and 54,364 km): `g_eff_radial` at the equator with
u = 450 m/s is 8.951 ± 0.005 m/s² and with u = 0 is 9.102; at the pole 12.137 (Lindal Table II:
8.96, 12.14). With Null's GM (3.7929085e16) the same to within 5e-4 m/s². `G_phi_eff` is zero
at the equator and the pole to round-off and positive (toward the equator, by the sign of
Eq. A5 with the J2 term) at 45° at 1 bar radius with u = 0, magnitude reported. Jupiter check
optional but recommended: with J2 = 14736e-6, J4 = −587e-6, J6 = 31e-6 at 71,398 km,
GM = 1.26686534e17, period 9h 55m 29.7s, radii 71,492 and 66,854 km, u = 100 m/s at the
equator: 23.116 and 27.015 (Lindal: 23.12, 27.01). A deliberately wrong coefficient (5, 9, 13)
must give 9.29 and 11.68, demonstrating that the check discriminates.

---

## 6. Step 5: `lib.geoid`

**Purpose.** The reference surface: Lindal Appendix Eqs. 11 to 14, and manuscript Eq. B3.

**Deliverable: `src/casspian/lib/geoid.py`.** Pure functions:

- `U_rigid(r, phi_c, Omega, GM, J, degrees, R_norm)`: Lindal Eq. 11, the potential of a
  uniformly rotating fluid, `U = V − ½ Ω² r² cos² φ_c` in the sign convention where `g_eff =
  −grad U` and U → 0 on the axis at infinity (Lindal's form is stated with the opposite overall
  sign; use one convention and say which in the docstring).
- `reference_geoid(phi_c_grid, r_polar, Omega, GM, J, degrees, R_norm, tol_m, max_iter)`:
  Lindal Eqs. 12 to 14. Set `U_ref = U_rigid(r_polar, π/2)`; at each latitude iterate
  `r ← r + (U(r, φ) − U_ref) / g_r(r, φ)` from an ellipsoid seed until `|Δr| < tol_m`. Returns
  `r_ref(φ_c)`, the iteration count per latitude, and the final `|ΔU|`. This is the **no-wind**
  reference geoid.
- `wind_geoid(phi_c_grid, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm, tol_m, max_iter)`:
  the same construction with the full effective gravity of Step 4, `u` altitude independent
  (Lindal's assumption), which is the manuscript's Eq. B3 solved as a potential closure rather
  than as an ODE in latitude. Since with wind the field is not exactly conservative when u
  varies with latitude (that is the whole point of §A6), this function integrates the slope
  equation `g dr₀/dφ = r₀ G_φ` (Eq. B3) from the pole inward with a fourth-order scheme on the
  supplied latitude grid, using the Step 4 components with `u(φ)`, and reports the closure
  `U(r₀(φ), φ) − U_ref` as a diagnostic of how far the wind takes the surface from an
  equipotential. Both the ODE result and the closure are returned.
- `radius_at(phi_c, geoid)`: interpolation of a constructed geoid to a latitude (linear in
  sin φ_c, stated).

**Acceptance.** No-wind reference geoid with Null's set, System III, `r_polar` = 54,438 km:
equatorial radius 60,244 ± 2 km and `r(31.0°)` = 58,435 ± 2 km (handoff §9A.6 numbers, which
were computed the same way); `r(30.8185°)` = 58,452.9 ± 2 km (handoff §9A.8). Wind geoid with
the Step 6 wind (or, if Step 6 is not yet built, a stand-in `u(φ)` with 450 m/s at the equator
falling to zero by 35°): equatorial radius larger than the no-wind value by roughly 120 km
(handoff: Lindal's fitted 60,367 against 60,244 no-wind, the wind bulge). The exact number with
the real Smith wind is to be reported, not prescribed, and compared to 60,367 ± 4.

---

## 7. Step 6: `lib.latitude`

**Purpose.** The planetographic to planetocentric conversion, handoff §9A.2, Lindal Eq. 5.

**Deliverable: `src/casspian/lib/latitude.py`.**

- `planetocentric_from_ellipsoid(phi_g, flattening)`: the seed, `tan φ_c = (1−f)² tan φ_g`.
- `planetocentric_fixed_point(phi_g, surface, u_of_phi, Omega, GM, J, degrees, R_norm,
  tol_deg, max_iter)`: iterate `φ_c ← φ_g − ψ(φ_c)` with ψ from Step 4 evaluated on the
  supplied surface (a constructed geoid from Step 5) at `r = surface(φ_c)`, seeded from the
  ellipsoid. Returns `φ_c`, ψ, the iterates, and the count. Vectorized over `phi_g` for the
  wind tool.
- The inverse, `planetographic_from_planetocentric`, which is a direct evaluation.

**Acceptance.** For `phi_g` = 36.3 on the no-wind reference geoid of Step 5 with Null's set:
`φ_c` = 30.8185 ± 0.001°, ψ = 5.4815 ± 0.001°, converging from the ellipsoid seed 30.8526 in
three or four iterations (handoff §9A.8 lists the iterates). Round trip: the inverse applied to
the result returns 36.3 to 1e-9. For `phi_g` = 36.5 the result is 31.005 (handoff §9A.2), which
demonstrates the sensitivity. With the Step 6 wind included, the value **moves**, by an amount
to be reported; that value, not 30.8185, is the one the reduction will freeze.

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
latitude over the digitized points.

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
   `height_m` from the raw bundle's `table1` group, uncertainty variables present and NaN, and
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
