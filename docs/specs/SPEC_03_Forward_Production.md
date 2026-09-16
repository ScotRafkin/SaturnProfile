# SPEC 03. The forward production and its hydrostatic closure against the source

CASSPIAN Saturn atmosphere reference model. Specification for the coding agent.

Version 0.13, 16 September 2026. Author of record: S. Rafkin. **Status: accepted by the author
at v0.2 (14 September 2026); v0.4 and v0.5 apply the coding agent's pre-execution review of
Step 0 and Step 3, v0.6 the REPORT_03_step0 findings (rulings in §8). Step 0 accepted
(REVIEW_03_step0) and swept at `5cdf07e`; Step 1 accepted (REVIEW_03_step1, `b3efc3d`); Step 2 accepted (REVIEW_03_step2, `994c787`); Step 3
accepted and swept (REVIEW_03_step3, `6ae113e`); Step 4 reviewed and accepted (REVIEW_03_step4;
v0.13 carries the rulings, §8), its acceptance commit and the closure rerun on the clean tree to
follow. This specification closes at that record; SPEC_04 is next.** The amendments in the Appendix have been applied to SPEC_00
(v0.16), SPEC_01 (v0.21) and SPEC_02 (v0.10).
Depends on `SPEC_00_Architecture_and_Data_Files.md` v0.20, the closed
`SPEC_01_Lindal_Tool_Chain.md` v0.28 and the closed `SPEC_02_Refrac_and_Diagnostics.md` v0.11,
which it does not repeat. See §9 for the revision history.

---

## 0. How this specification is to be worked

The protocol of SPEC_01 §0 applies unchanged: one step at a time; the status line at the top
of this document and `STATE.md` read before any step; acceptance script under
`reports/step03_<N>/` (ignored); figures a report cites under `reports/figures/` (committed);
`reports/REPORT_03_step<N>.md` and `reports/REVIEW_03_step<N>.md`; `STATE.md` updated; commit
on acceptance, then the directory sweep of SPEC_00 §8 over `occul_data/` and now also
`forward/`; the full regression (every SPEC_01 and SPEC_02 suite) after any change to `lib`,
`refrac` or `tools`; do not implement ahead. Every angle in `lib` is in radians. **A step that
changes an input file (v0.6, author ruling of 14 September 2026, REPORT_03_step0 decision 1):**
`refrac` and `forward` refuse `-dirty` inputs, so such a step is accepted on a candidate product
built in memory by the acceptance script, with the `-dirty` refusal relaxed inside that script
only and the relaxation named in the output; the registered product is rebuilt on the clean
tree at the sweep, and the suites that run the command-line tools on the on-disk products run
after the sweep and are recorded in the report then. No em dash or
en dash anywhere. No silent choices: a choice this document did not make is a decision in the
report or a question to the author, never a default in the code.

**What this specification covers and what it does not.** The production half of the forward
model, Appendix B5 and B7: the geopotential of a profile from the effective gravity, and
temperature and pressure from refractivity by hydrostatic integration, written as general
library functions that the transfer (SPEC_04) will call at any latitude. Its acceptance is the
**hydrostatic closure**: the production run on the Lindal refractivity at Lindal's own
latitude, on his own tabulated levels, with forward inputs built to match his assumptions,
which must return his tabulated pressure and temperature. That test is a test, not the model;
it is the first thing that can be tested, because nothing else in the forward model is needed
to run it. It also defines the forward run directory and its inputs, the namelist subset the
closure needs, and the forward product kind. It does not specify the transfer (B6), the uniform geopotential
grid (B5.3; SPEC_04), the wind and composition tools for forward runs, the altitude and the
delivered datum (B7.3, B7.4), or the Monte Carlo wrapper. It opens with a Step 0 that corrects
two things the closure exposed in the reduction, because the forward model is not to be
built on a radius grid that is known to be mis-defined or on pressures that carry print
rounding.

**The two findings behind Step 0, stated once.** The reviewing agent ran the hydrostatic
closure independently on the committed `lindal_refractivity.nc` (numpy only, gravity rebuilt from the
embedded G and R) before writing this document. (i) Lindal's altitude h is measured along the
local vertical, the field line of the effective gravity (his Appendix Eqs. 6 and 7,
`dh = dr / cos ψ`). Building the geopotential as the radial component times a radial increment,
which is how the manuscript's B1 and B2 read, leaves a systematic residual of about −5e-3 in
the recovered pressure at every level below the top; building it as the magnitude of the
effective gravity times the field-line increment leaves ±2.6e-3 of scatter and no offset. The
mismatch is `1 / cos ψ − 1 = 4.6e-3` at ψ = 5.49°. (ii) Every row of Table I from 0.20 to
1258.93 mbar lies on the grid `10^(k/100)` mbar within its printed precision (63 of the 65
rows agree to all printed digits; the top decade is printed to two significant figures and
differs from the grid by up to 1.19 percent, at 0.32 mbar). Kind T carries the printed values,
so N in the top decade carries the print rounding; the top-decade closure residual is
±1e-2 on the printed values and ±1e-3 on the grid. The author decided both on 14 September
2026 (§7, decisions 1 and 2).

---

## 1. Step 0: amendments to the reduction, one rebuild, one sweep

**Purpose.** Apply decisions 1 and 2 to the reduction chain so that the anchor product the
closure reads is right, with a single rebuild of `occul_data/lindal/` and a single sweep.
`r0`, `phi_c`, `psi` and everything in the anchoring decision (SPEC_02 decision 9) do not
change; the levels' pressures in the top decade, the number density and refractivity at those
levels, `radius_m` at every level but the anchor, and every hash downstream of the raw bundle
do.

**Deliverable 1: the pressure grid (SPEC_01 v0.21, Steps 2 and 9).**

- `occul_data/lindal/raw/lindal_scalars.toml` gains a table `[pressure_grid]` transcribing
  the inference: `rule = "10^(k/100) mbar"`, `denominator = 100`, `unit = "mbar"`,
  `excluded_printed_values_mbar = [1298.48]` (the bottom row is the end of the data, not a grid
  level, and stays as printed), `basis = "62 of the 65 rows admit exactly one grid value
  within their printed precision, 48 of them printed to four, five or six figures; the three
  rows printed to two figures that admit two grid values (0.20, 0.25, 0.32 mbar) are resolved
  one at a time by the equal spacing of the fixed rows below them; the spacing is 10, 8, 6, 4,
  2 hundredths of a decade in succession"` (v0.5; the v0.3 count of 63 and the v0.4 count of
  49 were both wrong), `value_source = "Table I, inferred by the reviewing agent, 14
  September 2026; not stated by the source. If the source's level list is ever found it
  replaces this inference."`. `raw/notes.md` states the same in prose.
- `lindal_table1.csv` is unchanged: it remains the transcription of what is printed.
- `casspian-lindal-raw` (Step 2) writes `table1/pressure_printed_Pa` (the CSV value converted
  to Pa, `value_source = "Table I"`, unchanged from today's `pressure_Pa`) and
  `table1/pressure_Pa` on the grid, by the rule (v0.4, coding agent's proposal adopted): for
  each row, the set of integers `k` for which `100 · 10^(k/100)` Pa, rounded to the number of
  decimals the CSV prints for that row, equals the printed value. A row with exactly one match
  is fixed. Ambiguous rows are resolved one at a time, starting next to the fixed rows: each
  takes the candidate that continues the spacing in `k` of the two nearest fixed rows on its
  fixed side; if both sides are fixed, both must agree or the row is refused; a row once
  resolved counts as fixed for the next (v0.5). A row still ambiguous after that, or with no
  match and not in `excluded_printed_values_mbar`, is a refusal, not a warning, with the row
  named. A row in `excluded_printed_values_mbar` keeps
  its printed value. No spacing numbers are written into the code; the spacing pattern is what
  the acceptance checks in the report. For Table I: 62 rows have one match; 0.32 mbar takes `k = −50` from
  the spacing of 0.40 and 0.50, then 0.25 takes −60, then 0.20 takes −70; one row is excluded. Per-level flag `pressure_grid_applied` (1 snapped, 0 as
  printed). The rule and the `[pressure_grid]` table go into `history` and into the `scalars`
  group as every other table does. The bundle still interprets nothing it does not carry beside
  the original: both columns are in the file.
- Every downstream tool (`casspian-composition-lindal`, `casspian-lindal-inputs`) reads
  `pressure_Pa`. Kind T's `pressure_Pa` therefore carries the grid values, with a global
  attribute `pressure_grid_rule` quoting the rule and `pressure_printed_Pa` carried as a
  variable with `value_source = "Table I"` (SPEC_00 §6.1 gains this optional variable for
  source profiles). Kind C's levels follow kind T's, as SPEC_02 Step 1 already requires.
- **Note carried in the raw bundle's `history`:** with `k = 0` at 1 mbar the top row is
  `k = −70` (19.9526 Pa). A nearest-integer snap lands 0.32 mbar on `k = −49` (0.3236 mbar)
  instead of `k = −50` (0.3162 mbar), and both round to 0.32; the equal-spacing rule above is
  what settles it, and the acceptance checks the resulting `k` sequence (10 apart to `k = 130`,
  then 8, 6, 4, 2).

**Deliverable 2: Eq. B1 along the local vertical (SPEC_02 v0.10, Step 3; SPEC_00 v0.16 §3.3
and §6.7).**

- `lib.reduction.absolute_radius(h_m, h_ref_m, r0_m, psi_rad)` = `r0 + (h − h_ref) cos ψ`,
  with ψ the frozen anchor tilt. The tabulated altitude is a distance along the local vertical
  from the anchor isobar; its radial projection is `cos ψ` of it. The variation of ψ along the
  profile (5.476° at the bottom to 5.546° at the top, measured) changes `cos ψ` by 8.7e-5 at the
  top, and the projected radius there by about 12 m once the variation is integrated along the path, and is neglected with this number stated in the
  SPEC_00 §6.7 note and in `reduction_record` (`radius_projection_rule = "cos psi at the
  anchor"`, `radius_projection_residual_m`).
- `height_above_anchor_isobar_m` remains `h − h_ref`, the field-line distance; its `long_name`
  says so.
- The latitude drift of the profile along the field line, `dφ_c = tan ψ dr / r`, about
  +0.027° at the top level and −0.010° at the bottom (28 km of horizontal displacement at the
  top, a quarter of the 0.109° label uncertainty in `phi_c`), is neglected; the profile is
  placed at `phi_c` at every level. Recorded in `reduction_record`
  (`latitude_drift_neglected_deg = [top, bottom]`) and in the manuscript notes.
- SPEC_00 §6.7: `radius_m(level)` is defined as the planetocentric radius of the level, from
  the anchor radius and the field-line altitude projected by `cos ψ`; `height_above_anchor_
  isobar_m` as the field-line distance. SPEC_00 §3.3 item 4 reads "B3.2: absolute radius per
  level, Eq. B1 with the tilt: `r = r0 + (h − h_ref) cos ψ`".

**Rebuild and sweep.** `casspian-lindal-raw`, then the composition tool, then
`casspian-lindal-inputs`, then `casspian-refrac` on the manifest, in that order, on a clean
tree after the acceptance commit; every product under `occul_data/lindal/` is rebuilt (T, C, D
and W change hash because the raw bundle hash they record changes, even where their content
does not; the keep-by-content rule of SPEC_02 Step 6 decision 4 already compares the
`raw_bundle` attribute as content, so a new raw bundle rebuilds T and D without any change to
the rule, and `input_hashes` stays excluded from the comparison, v0.6, REPORT_03_step0
finding 3). Products that do not read the raw bundle (G, R) are not rebuilt; their content and
hashes are unchanged (v0.6, REPORT_03_step0 decision 5).

**Expected values (measured by the reviewing agent; to be reproduced, not prescribed).**
`pressure_Pa` at the top level 19.9526 Pa, then 25.1189, 31.6228, 39.8107, 50.1187, 63.0957,
79.4328, 100.0; the 1 bar row and the bottom row unchanged (100,000.0 and 129,848.0 Pa). 65
rows snapped, one as printed. The `k` sequence starts at −70 and is 10 apart to `k = 130`
(20 mbar), then 8 apart to 194, 6 apart to 236, 4 apart to 268, 2 apart to 310, with the
excluded row after it (v0.6: the v0.5 endpoints 230 and 270 were wrong, REPORT_03_step0
finding 1). `n` at the top level 1.041934e22 m⁻³ (was 1.04441e22); at the second
level 1.270497e22; at 1 bar unchanged (5.373124e25). `N` moves at each snapped level by the
same fraction as `p`. `radius_m` at the top level 58,801,571.0 m (was 58,802,888.3, a change of
−1,317.2 m), at the bottom 58,412,566.6 m (+478.3 m), at the anchor level exactly `r0`.
`cos ψ = 0.995406`. F2's effective gravity along the profile changes by under 6e-5 (v0.6:
5.3e-5 at the top level, the Newtonian 2 dr/r = 4.5e-5 divided by `g_eff / g_N` plus the
centrifugal part's own change with `r`; REPORT_03_step0 finding 2).

**Acceptance.** The raw bundle carries both pressure variables, the flag, and the `k` sequence
above; a CSV copy with one printed value edited to something no grid value rounds to is
refused with the row named; a copy with the excluded row moved onto the grid is written with
66 snapped rows and the flag showing it (the exclusion list is data, not code). Kind T reads
back with `pressure_grid_rule` and `pressure_printed_Pa`. `radius_m` at the anchor level equals
`r0` exactly and at the top level equals `r0 + (h_top − h_ref) cos ψ` to 1e-6 m; `reduction_
record` carries the projection rule, its residual and the drift. Every acceptance suite the full regression shows affected (at least the SPEC_01 Step 2 and
Step 9 and the SPEC_02 Step 3 and Step 4 suites; the Step 02_1 suite must also read the
REPORT_03_step0 hash table; the coding agent finds the rest in the regression and lists each
in the report with what changed) is updated for the values above, its other checks unchanged,
and the full regression passes. `r0`, `phi_c`, `psi`, both polar
radii, the asymmetry and every `reduction_record` number of SPEC_02 Step 6 are unchanged to the
last figure REPORT_02_step6 prints. New hashes for every product recorded in REPORT_03_step0 in
the row format the Step 02_1 suite reads. Manuscript notes (§7, decision 3) appended to the
list the handoff carries.

---

## 2. Step 1: `lib.geopotential`, pure

**Purpose.** Eq. B2 restated along the profile's own vertical: the geopotential of each
tabulated level from the effective gravity magnitude and the field-line altitude, with the
gauge on a declared isobar. Handoff §4 for the gravity conventions; SPEC_01 Steps 4 and 6.

**The staggering, agreed before the first integral is coded (SPEC_00 §7.3 rule ii).** Levels
are the tabulated levels of the anchor, index `k = 0` at the top, radius decreasing with `k`,
as kind N stores them. Layers are `k + ½` between consecutive levels.

| Quantity | Lives on | Formed by |
|---|---|---|
| `h_k`, `r_k`, `N_k`, `ℛ̄_k`, `m̄_k` | levels | read from kind N |
| `g_k`, `G_φ,k`, `ψ_k` | levels | `lib.gravity.g_eff_vector` at `(r_k, φ_c, u(φ_c))`, `u` from the anchor's wind at the reference level |
| `|g_eff|_k = sqrt(g_k² + G_φ,k²) = g_k / cos ψ_k` | levels | algebra |
| `ΔΦ_{k+½} = ½ (|g_eff|_k + |g_eff|_{k+1}) (h_{k+1} − h_k)` | layers | trapezoid; `|g_eff|` varies by 1.6 percent over the profile and smoothly, the trapezoid error is below 1e-7 of the layer (bound `(Δh)² g'' / (12 g) = (Δh)² / (2 r²)` with `g'' ≈ 6 g / r²`, 2.5e-8 for a 13 km layer) |
| `Φ_k` | levels | `Φ_a = 0` at the gauge level `a`; `Φ_k = Φ_a + Σ` of the layer increments between `a` and `k`, the sum running upward for `k < a` and downward for `k > a`, so that Φ increases upward |
| `n_k = N_k / ℛ̄_k`, `ρ_k = n_k m̄_k` | levels | Eq. B4 and the ideal gas law; the same ℛ̄ that made N, so `n_k` is the reduction's `number_density_m3` to round-off in the closure |
| `I_{k+½} = ∫_{Φ_{k+1}}^{Φ_k} ρ dΦ` | layers | Step 2 |
| `p_k`, `T_k` | levels | Step 2 |

The gauge level `a` is the level whose tabulated pressure equals the namelist's
`gauge_isobar_Pa` exactly, as floats (v0.4: a label names a level; a tolerance would let a
near miss pass silently); Step 3 refuses a gauge isobar that is not a level of the anchor,
naming the nearest level, (no
interpolation of the gauge in this specification; SPEC_04 places the gauge on the uniform
grid). The wind enters `|g_eff|` through `Ω_abs` exactly as in F2, and for the Lindal anchor
`u(φ_c)` = 2.167 m/s, a 5.6e-5 effect.

**Deliverable: `src/casspian/lib/geopotential.py`**, NumPy in and out, no I/O:

- `effective_gravity_magnitude(u, r, phi_c, Omega, GM, J, degrees, R_norm)` returning
  `(g_mag, g, G_phi, psi)` on `r`, from `lib.gravity.g_eff_vector`.
- `layer_increments(g_mag, h)` returning `ΔΦ_{k+½}` by the trapezoid.
- `geopotential(dPhi, gauge_index)` returning `Φ_k` with `Φ = 0` at the gauge index.
- `geopotential_along_profile(...)`: the three in sequence, for `refrac` style orchestration.

**Expected values (measured; Lindal anchor after Step 0, gauge at the 100 mbar level).**
`|g_eff|` 9.892477 m/s² at the top level, 10.005625 at the anchor, 10.047103 at the bottom;
`ψ` 5.546374° at the top, 5.475835° at the bottom; `Φ` 2.852355e6 m²/s² at the top level and
−1.043743e6 at the bottom; the first layer `ΔΦ` = −1.2666e5 m²/s² (levels are top down, so the
increment from level 0 to level 1 is negative). These are values under the field-line rule; the
radial rule would give 2.839126e6 and −1.038963e6, and the difference is the 4.6e-3 of §0.

**Acceptance.** Uniform gravity (`J = 0`, `Ω = 0`, `u = 0`, `r` large enough that `GM / r²`
varies by under 1e-12 across the column) returns `Φ_k − Φ_a = g (h_k − h_a)` to 1e-12 relative;
`Φ_a` is exactly 0.0; `Φ` is strictly monotonic in `h`; the values above reproduced to the
figures printed; `|g_eff|` at the anchor equals `lib.gravity.g_eff_radial` there divided by
`cos ψ` to round-off; `geopotential` refuses a gauge index off the grid. F5 is not rendered
here (no product yet); it renders in Step 4.

---

## 3. Step 2: `lib.hydrostatic`, pure

**Purpose.** Eq. B5 with `m̄ / ℛ̄` inside the integral, the boundary at the top level, and
Eq. B6 for the temperature. SPEC_00 §7.3 rule (i) applied: density is near exponential in Φ
across a layer, so the layer integral is the exact integral of the log-linear interpolant,
never the trapezoid of endpoint values.

**Deliverable: `src/casspian/lib/hydrostatic.py`**, NumPy in and out, no I/O:

- `layer_mass(rho, Phi)`: `I_{k+½} = (ρ_{k+1} − ρ_k) (Φ_k − Φ_{k+1}) / ln(ρ_{k+1} / ρ_k)`,
  with the limit `ρ_k (Φ_k − Φ_{k+1})` taken when `|ln(ρ_{k+1} / ρ_k)| < 1e-10` (the two
  expressions agree to 1e-10 there); refuses a non-positive density or a layer of zero
  thickness.
- `pressure_from_top(p_b, I)`: `p_0 = p_b`, `p_{k+1} = p_k + I_{k+½}`, a cumulative sum from
  the top, which is the only direction this specification integrates (B7.2).
- `temperature(p, N, R_bar)`: `T = p ℛ̄ / (k_B N)`, Eq. B6, `k_B` from `lib.constants`.
- `density(N, R_bar, m_bar)`: `ρ = (N / ℛ̄) m̄`, Eq. B4 with the mass; `m̄` in kg per molecule
  from `mean_molar_mass_kg_mol` and the Avogadro constant of `lib.constants`.

**Why the rule, in numbers (measured at Table I's geopotential spacing).** The exact-exponential
layer integral recovers an isothermal column to 1e-15 and a column whose temperature falls
linearly with height, from 140 K at the bottom level to 80 K at the top, linear in Φ (v0.8,
REPORT_03_step2 finding 1: the other orientation gives 3.5e-4 and 6.5e-3), to 6e-4; the
trapezoid gives 8.5e-3 on both,
which is larger than the effect Step 0 corrected and would hide it.

**Acceptance.** An isothermal column (uniform composition, `T` = 100 K, `p_b` = 20 Pa) sampled
on the Step 1 `Φ` grid returns `p_k` to 1e-14 relative at every level and `T_k` = 100 K to
1e-14; the linear-`T` column above (140 K at the bottom level, 80 K at the top, linear in Φ)
returns `p` to 6e-4 or better at every level, and the trapezoid alternative, computed in the
acceptance script only, to 8e-3 to 9e-3, both reported;
halving every layer (midpoints inserted with the analytic `ρ`) reduces the linear-`T` error by
a factor of 3.5 to 4.5 (second order); `layer_mass` returns the limit form on a pair of equal
densities and refuses a zero-thickness layer with the layer named; `temperature` applied to the
reduction's own `N`, `ℛ̄` and tabulated `p` returns the tabulated `T` to 1e-12 relative (the
SPEC_02 Step 3 closure, now through this module).

---

## 4. Step 3: the run directory and its inputs, the run namelist, the product kind, and the reader warning

**Purpose.** Everything the closure reads and writes on disk, before the orchestration: the
`forward/<run>/` directory with its own input files, built by the tools under the run prefix
(SPEC_00 §2.3); `lib.control.read_run_namelist` for the SPEC_00 §7.2 subset the closure needs;
the forward product kind; and the `input_hashes` warning of SPEC_00 §8 for every derived kind
(carried from SPEC_02 decision 8).

**Deliverable 0 (v0.9): the season identifier across the chain.** Every file the model reads
or writes carries its season (SPEC_00 v0.17 §5; `docs/CASSPIAN_Seasonal_Design_Note.md` §2 and
§8.1; author decision of 15 September 2026 that all data, input or output, carry a season, and
that everything is valid at an instant). In this step:

- `lib.schema` requires `epoch` on every kind and exactly one of `solar_longitude_deg` (with
  `solar_longitude_source`) and `season_absent_meaning = "uniform"`; `epoch_note` optional.
  `epoch` is parsed: an ISO 8601 date, or on kind `profile` the SPEC_00 §5 string for a season
  declared without a date (v0.11); in this chain every file's date is 1981-08-26 (v0.12). The raw bundle's `role` is always `reduction`: it is the
  transcription of a source and has no build-file section (v0.11, REPORT_03_step3 decision 4).
- The transcription and the reduction chain are amended as SPEC_01 v0.25 and SPEC_02 v0.11 say
  (the `Ls` of the Voyager 2 ingress, 18.2° with the sub-solar latitude 8.06° N, transcribed
  with its source; G and R `uniform`; W with the imaging season; the reduction
  composition `uniform`; T, D and N with the observation's season), **every file of the chain
  carrying `epoch = "1981-08-26"`, the date of the occultation** (author decision of
  15 September 2026, v0.12: one date for the whole chain; a file whose content has no date of
  its own, the harmonic set, the rotation system, the wind curve, the composition, carries the
  date of the observation it serves, and its source's own dating goes in `epoch_note` when it is
  not already stated elsewhere in the file; SPEC_00 v0.19 §5, SPEC_01 v0.28). The G and R
  build-file sections carry `epoch` beside `role` and the tools write it; the wind's data
  properties file `[epoch]` says the same date with its existing note, and the chain is
  rebuilt once by the procedure of §0 for a step that changes an input, with the new hashes
  recorded in this step's report in the row format the Step 02_1 suite reads. Gravity and
  rotation are rebuilt this time, since their content changes.
- The run's inputs carry seasons the same way (deliverable 1): the wind the imaging season,
  gravity and rotation `uniform`, the composition `uniform`; all four dated 1981-08-26 (v0.12).
- The namelist `[run]` gains `solar_longitude_deg`, required, and the optional `date`
  (deliverable 2); in closure mode the run's season must equal the anchor's exactly.
- Kind `profile` carries the run's `solar_longitude_deg`, its `epoch` (the declared date, or
  the SPEC_00 string for a season declared without a date) and `anchor_solar_longitudes_deg`
  (deliverable 3).

No logic acts on the season beyond the closure equality; the identifier is there so nothing
later has to guess.

Two more things ride on the same rebuild (v0.10, §8 rulings 1 and 3):

- **`role` is a required key in every section of every build control file**, `"reduction"` or
  `"forward"`, written into the product's `role` global, and into `composition_role` by the
  composition tool. No default in code: a section without it is refused. `lindal_build.toml`
  gains `role = "reduction"` in every section; the run's build file says `"forward"`. Since the
  chain is rebuilt in this step anyway, the hash changes this causes cost nothing (SPEC_01
  v0.26).
- **Every writer records paths relative to the directory of the file it writes** (SPEC_00
  v0.18 §5): `input_hashes`, `control_file`, `raw_bundle`, `raw_sources`, `master_table`,
  `latitude_conversion_inputs`, `decomposition_geometry`, and any other attribute that names a
  file, with `/` as the separator. The absolute paths some tools wrote were a defect against
  SPEC_00 §5, and SPEC_00 §7 contradicted §5; §7 is corrected. A helper in `lib.io` does the
  conversion once for every writer.

**Deliverable 1: the run directory and its inputs.** `forward/lindal_closure/` holds the
namelist, `lindal_closure_build.toml` (the run's build control file, SPEC_00 §2.2.2), `inputs/`
written by the tools, and `output/` written by the model; nothing else. **Every input kind
the forward model reads is a file of its own under the run prefix, distinct from the file of
the same kind the reduction used, whatever its content.** The reduction's files under
`occul_data/lindal/` and the copies embedded in kind N are what the source assumed and are
never read by `forward`; the files under `inputs/` are what this run assumes. For the closure
the two are meant to agree, and in closure mode the model checks that they do (deliverable 2).
All four carry `role = "forward"` and `profile_or_run = "lindal_closure"`:

- **Kind W, `lindal_closure_wind.nc`:** the same Ingersoll and Pollard curve with the Smith
  points as its uncertainty, on the same latitude and pressure grids, altitude independent,
  the trivial decomposition (`u_cylindrical = u_total`, `u_shear = 0`, `decomposition =
  "trivial_altitude_independent"`), zero at the poles. The shear kernel of SPEC_04 will be
  zero on it, by declaration. Built by `casspian-wind-from-curve`.
- **Kind G, `lindal_closure_gravity.nc`:** Null et al. (1981) at 60,000 km with the Jacobson
  2006 GM, the set the reduction used, from the same `data_static` transcription. Built by
  `casspian-gravity-file`.
- **Kind R, `lindal_closure_rotation.nc`:** System III, from the same `data_static` entry.
  Built by `casspian-rotation-file`.
- **Kind C, `lindal_closure_composition.nc`:** 0.94 H2 and 0.06 He of the remainder after the
  Table I ammonia column, `lindal1985` refractivities, on the anchor's levels, with
  `composition_role = "forward"` (the attribute SPEC_00 §6.2 gives this kind) and the same
  closure declaration. Built by `casspian-composition-lindal`.

No kind D: `forward` refuses a geodesy file (SPEC_00 decision 3); the reference surface comes
from G, R, W and the frozen `r0(phi_c)` in kind N.

```
forward/lindal_closure/
    lindal_closure.toml
    lindal_closure_build.toml
    inputs/
        lindal_closure_composition.nc    kind C, the composition the closure assumes
        lindal_closure_gravity.nc        kind G
        lindal_closure_rotation.nc       kind R
        lindal_closure_wind.nc           kind W
    output/
        lindal_closure_profile.nc        kind profile (deliverable 4)
        figures/
```

They are built by the existing tools, each pointed at the same sources the Lindal reduction
used, with the run's prefix, titles and output paths declared in `lindal_closure_build.toml`
(the same section vocabulary as `lindal_build.toml`, minus `[stage_two]`, which has no
forward counterpart). A thin driver `casspian-run-inputs <run_build.toml>` runs the four sections in
order (gravity, rotation, wind, composition), so a run's inputs are made by one command; it
contains no physics and no value. The composition tool for this run is `casspian-composition-
lindal` with its output redirected, since the generic composition tool does not exist yet
(the forward-inputs specification); its `raw_bundle` pointer is a path outside the run
directory, which the build file may name because a build file is a tool control file, not the
namelist. The files will be content-identical to the reduction's copies apart from prefix,
title, role and writer attributes, and their hashes will differ; that is intended. Disk is
cheap and the distinction between what the source assumed and what a run assumes is not.

`forward/README.md` is rewritten to describe the layout. `.gitignore` already excludes
`forward/*/output/`; the namelist and the build file are committed. Whether `inputs/*.nc` is
committed is the author's call (v0.4, §8 ruling 5): until it is made, the coding agent adds
`forward/*/inputs/*.nc` to `.gitignore` so that the repository's practice (no netCDF product
committed yet, SPEC_00 §2.4 pending) stays uniform, and the four inputs are regenerable from
the build file with their hashes recorded in the Step 3 report.

**Deliverable 2: `lib.control.read_run_namelist(path)` and `load_run_inputs(namelist)`.**
The vocabulary accepted in this specification (SPEC_00 §7.2 v0.16 adds `mode`, `p_b_rule`,
`[diagnostics]`, and `product`, and states the closure rules):

```toml
[run]
name        = "lindal_closure"
description = "Hydrostatic closure of the Lindal reduction at its own latitude under matching forward inputs"
mode        = "closure"          # this specification implements only closure; SPEC_04 adds transfer
solar_longitude_deg = 18.2       # v0.9: the run's season; in closure mode equal to the anchor's
date        = "1981-08-26"       # v0.13: the closure closes on the occultation, so it carries its date

[[anchors]]                      # exactly one entry in closure mode
slug   = "lindal"                # must equal profile_or_run inside the file, or refuse
path   = "../../occul_data/lindal/lindal_refractivity.nc"
weight = 1.0
measurement_uncertainty_scale = 1.0

[inputs]                         # all under inputs/, all prefixed with the run name (SPEC_00 section 7.2)
composition = "inputs/lindal_closure_composition.nc"
gravity     = "inputs/lindal_closure_gravity.nc"
rotation    = "inputs/lindal_closure_rotation.nc"
wind        = "inputs/lindal_closure_wind.nc"

[hydrostatic_boundary]
p_b_rule     = "anchor_profile_top"   # the tabulated pressure at the anchor's top level, read
                                      # from the thermo group embedded in kind N; the only rule
                                      # accepted in closure mode. A general run may instead
                                      # declare p_b_Pa (SPEC_00 section 7.2); exactly one of the two.
p_b_location = "top_of_anchor_profile"

[isobars]
gauge_isobar_Pa = 1.0e4          # where Phi = 0; must be a tabulated level of the anchor in closure mode

[output]
directory = "output"
product   = "lindal_closure_profile.nc"

[diagnostics]                    # as SPEC_00 section 7.1
figures = true
format  = "png"
dpi     = 150
```

Rules. `[inputs]` is required and every rule of SPEC_00 §6.6 and §8 on the run's inputs is
enforced at load: the wind's rotation system and rate equal the run's kind R; the wind covers
the anchor's pressure range and latitude; the components sum to the total; the poles are
zero; each file carries the run prefix and `role = "forward"`; kind C carries
`composition_role = "forward"`. In `closure` mode, additionally: the run's four inputs must be
**content-identical to the anchor's embedded copies** (compared through
`xarray.Dataset.identical`, root and every group, after dropping the writer globals
(`created_by`, `created_at`, `casspian_git_commit`, `casspian_version`, `history`), the naming
and role attributes (`title`, `profile_or_run`, `role`, `composition_role`), and every
provenance attribute that carries a path or a hash of another file (`input_hashes`,
`control_file`, `raw_bundle`, `raw_sources`, `master_table`, `master_table_hash`,
`latitude_conversion_inputs`, `decomposition_geometry`, and any attribute whose value contains
`sha256:`); the report lists the attributes dropped per file (v0.10, §8 ruling 2)), otherwise the run
is refused with the first differing variable or attribute named, because a closure under
inputs that differ from the source's is not a closure; the composition's levels are the
anchor's levels; no `[target]` (the target is the anchor's `phi_c`; a present `[target]` is
refused); `[grid]`, `[numerics]`, `[estimation]` and `datum_isobar_Pa` are refused as **not
implemented in this specification**, a message distinct from the unknown-key refusal, so that
a namelist written for SPEC_04 fails for the right reason. Every refusal of SPEC_00 §7.2 and
§8 that applies (geodesy key, kind mismatch, `-dirty` anchor or input, slug mismatch, a path
outside the run directory other than the anchor) is enforced. `load_run_inputs` reads the
anchor as kind N and the four inputs under their kinds with `lib.io.read` and returns them
with their hashes; the closure comparison is a function of its own, `check_closure_inputs`,
so the report can show what it compared.

**Deliverable 3: the product kind `profile` (SPEC_00 v0.16 §5 and new §6.8).** Written by
`forward`, role `forward`, `casspian_kind = "profile"`, file `<run>_profile.nc`, one latitude,
dimension `level`, `latitude_planetocentric_absent_meaning = "point"` (the SPEC_00 §5 form,
named for the dimension; v0.11) with `latitude_planetocentric_deg` a scalar. Variables, all `provenance = "modeled"` unless noted:

| Variable | Units | Notes |
|---|---|---|
| coordinate `geopotential_m2s2(level)` | m²/s² | `positive = "up"`; strictly monotonic; the gauge stated in `gauge_isobar_Pa` |
| `radius_m(level)` | m | copied from the anchor (`derived`) |
| `height_above_anchor_isobar_m(level)` | m | copied (`measured`) |
| `refractivity(level)` | 1 | copied (`derived`); in closure mode identical to the anchor's |
| `mean_refractivity_m3(level)`, `mean_molar_mass_kg_mol(level)` | m³, kg/mol | from the run's kind C and its species group, on the anchor's levels; in closure mode equal to the anchor's to round-off |
| `number_density_m3(level)` | m⁻³ | Eq. B4 |
| `pressure_Pa(level)` | Pa | Eq. B5; the hydrostatic pressure; attribute `boundary_pressure_Pa` on the variable and as a global |
| `temperature_K(level)` | K | Eq. B6 |
| `pressure_tabulated_Pa(level)`, `temperature_tabulated_K(level)` | Pa, K | closure mode only: the anchor's tabulated values copied for the comparison (`derived`) |
| uncertainty companions of every modeled variable | | present and NaN, `uncertainty_kind = "1sigma"`, `uncertainty_method = "not propagated"`, `uncertainty_terms_unstated` listing every declared term of the anchor (composition, anchor radius, label latitude) with the reason: a first-order companion carrying only the terms that pass through the integral would be read as a total; the Monte Carlo wrapper supplies the product's uncertainty (decision 5) |

Globals (v0.9): `epoch` (the run's declared date, or the SPEC_00 §5 string when only a season
was declared), `solar_longitude_deg` (the run's), `solar_longitude_source = "declared in the
namelist"`, and `anchor_solar_longitudes_deg`, the seasons of the anchors in `[[anchors]]`
order. Scalars: `latitude_planetocentric_deg`, `psi_deg` (copied), `gauge_isobar_Pa` (`index`),
`gauge_level_index` (`index`), `boundary_pressure_Pa` (`index` in closure mode, since it
names the anchor's top level value), `boundary_level_index`. Groups: `anchors/<slug>`, the
kind N file verbatim (variables, attributes and its own groups); `inputs/composition`,
`inputs/gravity`, `inputs/rotation`, `inputs/wind`, the run's four input files verbatim;
`namelist` (attributes `text` and `sha256`, resolved form); `production_record` (attributes: `mode`, the closure input comparison
and its result, the gauge level, `p_b` and its rule, the projection and drift notes copied from the anchor's `reduction_record`, the
layer rule names of Steps 1 and 2, the closure statistics of Step 4, `codata_release`,
`casspian_version`, `casspian_git_commit`). Global `input_hashes` lists the anchor file, the four
inputs and the namelist. `lib.io.read(path, "profile")` returns a `DataTree` and validates; the schema
version is 1.

`tools/plots` dispatches on kind `profile`: F5 (`Φ` against `p` with `h` on a twin axis; the
across-latitude panel of SPEC_02's F5 is dropped until SPEC_04 builds the surface) and F6
(`pressure_Pa / pressure_tabulated_Pa − 1` against `p`, with the boundary level and the gauge
level marked, and the print-rounding envelope of the anchor drawn as a band: half a unit of
the printed pressure column's precision, which after Step 0 is zero at the snapped levels and
3.9e-6 at the off-grid bottom row, plus the height rounding term, half a unit of the printed
height column's precision divided by the local scale height; the precision of a column is the
finest any value of it shows, 0.01 mbar and 0.1 km for Table I, since the last printed figure
of one value cannot be recovered from a float (v0.11, REPORT_03_step3 finding 5). The band is
the leading-order budget of §5, not a bound. The gauge marker on F5 and F6 is labeled with the
declared `gauge_isobar_Pa`, the tabulated value that defines the gauge, not the hydrostatic
pressure found at that level (v0.11). F6 renders
only when the tabulated companions are present; F5 always. The F5 and F6 placeholders written
against `geopotential_m2s2` and `pressure_hydrostatic_Pa` on kind N are removed; kind N
renders F1 to F4 and reports nothing skipped (SPEC_02 decision 10 superseded; §7 decision 4).

**Deliverable 4: the `input_hashes` warning for every derived kind.** `lib.io.read` warns, for
any kind whose file carries `input_hashes`, when a listed path, resolved against the directory
of the file being read, names a file that exists and whose hash differs, and says nothing when
the path does not resolve to a file (a copied file is not a defect). With every writer now
recording paths relative to its product (deliverable 0), this is one rule for every kind, and a
file moved with its inputs keeps its check while a file moved alone reads silently. Kind N keeps
its refusal on a missing recorded hash. This is the SPEC_00 §8 rule as corrected at v0.18,
implemented once for every kind.

**Acceptance.** `casspian-run-inputs forward/lindal_closure/lindal_closure_build.toml` writes
the four inputs, each validating under its kind with the run prefix and `role = "forward"`,
and each content-identical to the anchor's embedded copy under the comparison above (reported
per file with the attributes dropped listed); their hashes differ from the reduction's (also
reported). The Lindal closure namelist loads and resolves; copies with `[inputs]` missing,
`[target]`, `[grid]`, an unknown key, a second `[[anchors]]` entry, a slug that does not
match, `p_b_Pa` alongside `p_b_rule`, a gauge isobar of 9,000 Pa, a `-dirty` anchor, an input
without the run prefix, a namelist that names a geodesy file, and a composition copy with
`x_He` raised and `x_H2` lowered by 0.01 so that the sum rule still holds (and, if the kind C
reader's closure-declaration check trips first, its declared share edited to match, so that
the only refusal left is the closure comparison, with `x_He` named; the report says which
edits were needed) are each refused with the right message (twelve cases, each message
quoted in the report; v0.4), plus (v0.9) a namelist whose `solar_longitude_deg` differs from
the anchor's, refused as not a closure (thirteen cases); every product of the rebuilt chain and
every run input validates with the season attributes, a copy with neither `solar_longitude_deg`
nor `season_absent_meaning` and a copy with both are each refused, kind N's season equals kind
T's, and the transcribed pair satisfies `sin δ_s = sin ε sin Ls` to 0.05° (v0.9); the SPEC_02
Step 5 suite's "F5 and F6 skipped" check is retired in this step, since kind N now reports
nothing skipped; a synthetic
kind `profile` file built in memory validates, reads back as a `DataTree`, and is refused when
`pressure_tabulated_Pa` is present without `temperature_tabulated_K`, when the coordinate is
not monotonic, or when a modeled variable lacks its NaN companion; `casspian-plots` on the
synthetic file renders F5 and F6 with the envelope band and on `lindal_refractivity.nc` renders
F1 to F4 and reports nothing skipped; a kind T copy whose recorded raw bundle hash is edited
warns on read when the raw bundle sits at the recorded relative path beside it (the copy is
placed in a scratch directory with a `raw/` subdirectory holding the bundle), and the same copy
alone in a directory without the bundle reads silently (v0.10); every product of the rebuilt
chain and every run input carries only relative paths, checked by a scan for a drive letter or
a leading `/` in every path-bearing attribute.

---

## 5. Step 4: `forward.production`, the closure run, the product, F5 and F6

**Purpose.** The production (B7) as a function that will serve any latitude: given a
refractivity profile on its levels, the run's composition, gravity, rotation and wind, a gauge
and a boundary pressure, return geopotential, pressure and temperature. Its first use, and this
step's acceptance, is the closure run: the production at the source latitude on the source's
levels, the product written, the figures rendered, and the closure measured against the stated
envelope.

**Deliverable: `src/casspian/forward/production.py`** with `produce(profile, inputs, gauge,
p_b)` (pure orchestration of Steps 1 and 2 on arrays and the run's inputs, no file access,
usable by SPEC_04 at a transferred latitude unchanged) and `run(namelist_path)`, plus the
entry point `casspian-forward <namelist>` (thin: takes the path, dispatches on `mode`, calls
one function, returns; SPEC_04 adds the transfer mode to the same entry point). In closure
mode, in order: read the namelist, the anchor and the run's inputs, and check the inputs
against the anchor's (Step 3); `u(φ_c)` from the run's wind at its reference level; `|g_eff|`,
`ψ` and `Φ` on the anchor's levels with the gauge at the declared isobar (Step 1); `ℛ̄` and
`m̄` on the levels from the run's kind C and its species group; `ρ` from `N`, `ℛ̄`, `m̄`; `p_b`
from the anchor's top level by the rule; `p` by the layer sums from the top and `T` by Eq. B6
(Step 2); the product written with the anchor and the inputs embedded; the figures rendered
when the namelist asks. The closure statistics written to
`production_record`: the residual `pressure_Pa / pressure_tabulated_Pa − 1` at every level as
an attribute array, its maximum magnitude above 2 mbar (the ten top levels), between 2 and 10
mbar, between 10 and 100 mbar, and below 100 mbar excluding the bottom row, the bottom row's
value, and the mean below 10 mbar (v0.13: the bins stated in pressure, a level on an edge
belonging to the deeper bin, so that "the top decade" has one reading and every level falls
under a bound; REPORT_03_step4 finding 1).

**Expected values (measured by the reviewing agent on the committed product with Step 0's
changes applied in memory; to be reproduced, not prescribed).** Residual
`p_hydro / p_tab − 1`: within ±9.8e-4 above 2 mbar (the ten top levels, which is what the
reviewing agent's "top decade" measured), within ±2.2e-3 between 2 and 10 mbar (2.12e-3 at
6.31 mbar, REPORT_03_step4), within ±2.6e-3 between 10 and 100 mbar, within ±1.7e-3 below 100
mbar excluding the bottom row, −3.3e-3 at the bottom row (the off-grid last level, at the
X-band extinction limit), mean below 10 mbar −4.9e-4. Below about 300 mbar the residual is a
trend, negative and growing with depth, not scatter; it sits with the systematic terms below
and is a statement about the source's own integration that the manuscript can make once. `T_hydro /
T_tab − 1` equals the pressure residual at every level to round-off (they are one number:
`N = p_tab ℛ̄ / (k_B T_tab)` and `T_hydro = p_hydro ℛ̄ / (k_B N)`). The residual budget the
envelope rests on, per level and uncorrelated between levels: height rounding ±50 m gives
±1.3e-3 in `ln p` at a 39 km scale height, the dominant term; temperature rounding ±0.05 K
gives ±4e-4 to 6e-4 in the layer density; pressure print rounding is zero at the snapped
levels after Step 0 and 3.9e-6 at the bottom row; discretization a few 1e-4 (Step 2).
Systematic terms: the harmonic set Lindal actually used against Null 1981, bounded near 1e-5 in
`g` by the Campbell and Anderson comparison in manuscript B2; ammonia in `m̄`, up to 4.7e-4 at
the bottom if his mass excluded it; his own integration numerics, unknown. **The negative
control (v0.7, REPORT_03_step1 finding 1):** the same run with `Φ` built as the radial
component times the tabulated altitude increment, `g_k (h_{k+1} − h_k)`, which is the
manuscript's present reading of B1 and B2 and the construction §0's finding exposed, computed
in the acceptance script only with no code path in `forward`; it gives a residual with mean
−5.1e-3 below 10 mbar and −7.8e-3 at the bottom row (one factor of `1 / cos ψ`). The fully
radial reading, `g_k (r_{k+1} − r_k)` with `r` the projected radius, is reported beside it
(two factors; mean about −9.6e-3, bottom about −1.2e-2, measured by the reviewing agent) and
must fail the same bound.

Also in this step (v0.11, REPORT_03_step3 finding 3): the kind C reader enforces the
`source_statement` rule its schema declares for role `reduction` (`conditional_globals`), which
no code read; no product changes and nothing is rebuilt, since the composition tool already
writes the attribute in both roles.

**Acceptance.** `casspian-forward forward/lindal_closure/lindal_closure.toml` writes
`output/lindal_closure_profile.nc`, which reads back as kind `profile` with the anchor and the
four inputs embedded byte-identical (compare through `xarray.Dataset.identical` on every
group); `mean_refractivity_m3` and `mean_molar_mass_kg_mol`, formed from the run's kind C,
equal the anchor's to 1e-14 relative, and `number_density_m3` likewise; `pressure_Pa` at the top level
equals `p_b` exactly and `p_b` equals the anchor's top tabulated pressure; `Φ` is zero at the
gauge level; the residual is within 1.5e-3 in magnitude above 2 mbar (the ten top levels),
within 3e-3 at every level from 2 mbar down excluding the bottom row, within 4e-3 at the
bottom row, and its mean below 10 mbar is within ±1.5e-3 (v0.13: the bins of the statistics
above, every level under one bound); the temperature residual equals the pressure
residual to 1e-12 at every level; the negative control's mean below 10 mbar is below −4e-3, and so is the fully radial
reading's (both must fail the mean bound, which is what shows the bound has teeth; v0.7); F5 and F6 are rendered
by the driver and by `casspian-plots` by hand, byte-identical apart from the footer; the
residual is drawn in F6 against the envelope band, and the levels at which it lies outside the
band are counted and reported with their values (v0.11: the band is the leading-order budget,
±50 m over the scale height, 1.3e-3 at 39 km, while the expected residual between 10 and
100 mbar is ±2.6e-3, so inclusion in the band is a description, not an acceptance; the numeric
bounds above are the acceptance and do not move); the author views
F5 and F6 (attached to the report) and accepts them by eye. A residual outside the ranges
above is a finding with the measured value and the reason, not a loosened check.

---

## 6. What comes after

SPEC_04, the transfer: the uniform geopotential grid of SPEC_00 §7.3 with the re-indexing of
B5.3 and its interpolation error measured against this specification's closure on the
tabulated levels; the reference surface at other latitudes from Eq. B3 with the run's own G, R,
W and the frozen `r0(phi_c)`; the shear kernel, the isobar tracing and the transfer (B6); the
composition-on-pressure fixed point of B7.1 when a chosen composition arrives on its own grid;
altitude and the delivered datum (B7.3, B7.4) and the kind T `model_output` instance, each
step with its own unit acceptance. The end-to-end tests are a specification of their own
after it (SPEC_05): the transfer identity tests (a target equal to the source must return the
input exactly; a transfer to another latitude and back must return it to the discretization
error of the tracing), the uniform composition scaling against Lindal's Fig. 5 caption
(pressure as `m̄ / ℛ̄`, temperature as `m̄`, with `p_b` scaled the same way), and the first runs
under chosen inputs. Builds and tests are kept in separate specifications (author, 14
September 2026). Before or beside SPEC_04, the forward-inputs specification (generic wind tool, composition scenario
generator). Then the Monte Carlo wrapper, which is where the product of this specification
gets an uncertainty.

---

## 7. Decisions made in this document

1. **(Author, 14 September 2026) The model's vertical is the local vertical, and the
   geopotential of a source profile is the line integral of the effective gravity magnitude
   along it.** Lindal's altitude is measured along the field line (his Eqs. 6 and 7); the
   closure confirms it independently of that reading (§0). Eq. B1 gains the projection
   `cos ψ` and is fixed now, in the reduction, with one rebuild (Step 0), because the forward
   model is not to be built on a radius grid known to be mis-defined; `r0` and `phi_c` do not
   move. The latitude drift along the field line is neglected with its size recorded.
2. **(Author, 14 September 2026) The tabulated pressures are the grid `10^(k/100)` mbar.** The
   printed values are rounded representations of it and are kept beside the grid values at
   every stage; the inference is declared in the transcription with its basis and its caveat,
   and the off-grid bottom row stays as printed (Step 0).
3. **Manuscript notes carried from this specification** (appended to the handoff's list):
   the notation list should stop saying z and r are interchangeable; B1 gains `cos ψ` and a
   sentence on the field-line altitude; B2 is restated along the local vertical with the
   magnitude of the effective gravity, or the radial statement is kept with the levels placed
   on the field line, but not the present mixture; the latitude drift along the field line is
   mentioned once with its size; B3.1 or B2 states the pressure grid inference and that Table I
   is printed to two figures in its top decade; the Draft 8 passage in B3.2 arguing for polar
   anchoring is rewritten to decision 9 of SPEC_02.
4. **The forward product kind is `profile`**, with `pressure_Pa` the hydrostatic pressure and
   `pressure_tabulated_Pa` beside it in closure mode; SPEC_02 decision 10's placeholder names
   on kind N are removed and `tools/plots` renders F5 and F6 for kind `profile` only.
5. **Uncertainty companions on the product are NaN with the reason stated.** A companion that
   carried only the composition term through the integral, while the wind, anchor and label
   terms sit outside it, would be read as a total. The Monte Carlo wrapper is where the product
   gets an error bar; `uncertainty_terms_unstated` says so.
6. **The closure is a `forward` run, not a special case, and it is a test, not the model.**
   A run directory, a namelist in the SPEC_00 §7.2 vocabulary with `mode = "closure"`, the same
   entry point and the same `produce` function SPEC_04 will call at a transferred latitude.
   In closure mode the target is the anchor's latitude by definition, so the namelist carries
   no `[target]`.
7. **(Author, 14 September 2026) A forward run's inputs are its own files**, built by the tools
   under the run prefix into `inputs/`, never the reduction's copies or the anchor's embedded
   groups, even when their content is identical, as it is for the closure. What the source
   assumed and what a run assumes are different things and live in different files; the
   closure mode checks that they agree and refuses when they do not.
8. **No regridding in this specification.** The closure runs on the source's tabulated
   levels so that its residual budget contains no interpolation term; the uniform grid and its
   interpolation error belong with the transfer (SPEC_04), measured against this result.
9. **The temperature check and the pressure check are one check**, and the specification says
   why rather than pretending to two.
10. **The radial construction is kept as a negative control in the acceptance**, not in the
   code, so that the closure bound is known to detect the error it was written to detect.

---

## 8. Rulings on the coding agent's pre-execution review of v0.3 (14 September 2026)

1. **Three rows admit two grid values, and the v0.3 text was self-contradictory.** Correct
   on both counts, verified: 0.20 (`k` −70 or −69), 0.25 (−61 or −60), 0.32 (−50 or −49); 62
   rows admit one. The coding agent's rule is adopted as written in Step 0 deliverable 1: a
   row with one match is fixed, a row with several takes the one equally spaced with its fixed
   neighbors, anything still ambiguous is refused, and no spacing numbers go into the code.
   **v0.5:** the v0.4 wording refused Table I, because the three ambiguous rows are adjacent
   and the top row has no fixed neighbor within two rows; the coding agent's second wording
   (resolve one at a time from the fixed side, a resolved row counting as fixed) is adopted.
2. **The "63 of the 65" count was wrong, and so was v0.4's 49.** Verified: of the 62 rows
   with one match, 48 are printed to four, five or six figures (12, 30 and 6), 10 to three and
   4 to two. The basis string carries 48; the transcription takes the v0.5 text.
3. **Twelve refusal cases, eleven listed.** The missing one is the geodesy key, as guessed;
   added.
4. **The `x_He + 0.01` case would fail the sum rule first.** Correct; the case is now a
   compensated edit, and if the closure-declaration check of kind C trips before the closure
   comparison, the declared share is edited too, so that the case tests what it is meant to.
5. **Committing `inputs/*.nc`.** A question for the author, since it would be the first netCDF
   committed to the repository. Until answered, ignored by git and regenerable, hashes in the
   report.
6. **Gauge match.** Exact float equality, refusing with the nearest level named.
7. **Affected suites.** The list was incomplete by design of the coding agent's own finding;
   the acceptance now says every suite the regression shows affected, listed in the report,
   with the Step 02_1 hash-table read and the Step 5 skipped-figures check called out.

**Rulings on REPORT_03_step0 (REVIEW_03_step0, 14 September 2026).** Finding 1, the `k`
run endpoints: correct, verified against the raw bundle (6 apart to 236, 4 apart to 268);
restated here and in SPEC_01 v0.24. Finding 2, the gravity bound: correct, 5.284e-5 reproduced
independently; the bound is 6e-5 with the decomposition stated. Finding 3, keep-by-content:
correct, the `raw_bundle` attribute already does it; the sentence is restated. Finding 4, the
ammonia fill: recorded, no action. Finding 5 and decision 1: the in-memory candidate is the
rule for any step that changes an input, now in §0. Decision 5: gravity and rotation kept.
Decisions 2, 3, 4, 6, 7: accepted as reported.

**Rulings on REPORT_03_step1 (REVIEW_03_step1, 14 September 2026).** Finding 1: correct; the
Step 4 negative control as written at v0.6 was the two-factor construction, while the numbers
quoted were the one-factor construction measured before Step 0. The control is the one-factor
construction, the manuscript's mixture, because that is the error the closure was written to
detect; the two-factor reading is reported beside it and must fail too (Step 4 restated).
Finding 2: noted, no action; the bound holds with a factor of three to spare. Decisions 1 to
6 accepted as reported.

**Rulings on REPORT_03_step2 (REVIEW_03_step2, 15 September 2026).** Finding 1: correct; the
linear-T column is the one falling with height, 140 K at the bottom level and 80 K at the top,
which is what the reviewing agent measured; restated in Step 2. Finding 2: no action.
Decisions 1 to 5 accepted; the `expm1` form (decision 1) is the better evaluation of the same
expression and is recorded as the module's form.

**Rulings on the coding agent's pre-execution review of Step 3 (v0.10, 15 September 2026).**
1. **`role = "forward"` cannot be written by the tools.** Correct. Option (b) is adopted, not
   the recommended (a): `role` is a required key in every build-file section with no default,
   because a default in code is what the standing rule forbids, and because the chain is
   rebuilt in this step for the season identifier anyway, so the hash changes (b) would have
   caused are already being paid. `composition_role` follows the same key.
2. **The closure comparison would refuse identical inputs on `latitude_conversion_inputs`.**
   Correct; the drop list was incomplete. It is now stated as a rule (every attribute that
   carries a path or a hash of another file, plus the writer, naming and role attributes) with
   the known members enumerated, and the report lists what was dropped per file.
3. **The `input_hashes` warning cannot pass its acceptance as specified.** Correct, and the
   cause is a contradiction in SPEC_00: §5 asked for relative paths and §7 said the resolved
   absolute form is recorded, and the tools followed one or the other. Option (b) is adopted,
   not the recommended (a): every writer records paths relative to its product, SPEC_00 §7 is
   corrected, and the chain is rebuilt with the season identifier in this step, so the rebuild
   costs nothing extra. The reasons: an absolute Windows path in a product is a defect for a
   run directory that must be self-contained and for the service deployment of SPEC_00 §11;
   resolving relative paths against "the repository root the file sits in" would have added a
   second convention and a root-finding rule to the reader; and the acceptance then tests the
   rule as written. The writer changes are the cost, and they are the fix of a defect the
   review found, which is what the standing rule about fixes says to do.
4. **The rest, as the coding agent proposed:** the run's wind build section names the raw
   bundle and the run's own gravity and rotation files (the reference geoid for a forward run
   under other harmonics is a forward-inputs specification question, not this step's); the
   composition comparison covers the `species` group; the four run inputs are loaded in the
   acceptance through the §0 relaxation and rebuilt clean at the sweep; `forward/*/inputs/*.nc`
   goes into `.gitignore`.

**Rulings on REPORT_03_step3 (REVIEW_03_step3, 15 September 2026).** Finding 1: correct, the
contradiction was SPEC_01 v0.25's "as already recorded" against SPEC_00 §5's ISO date. The
author's decision (v0.12, superseding the v0.11 wording that asked for the sources' own dates):
every file of the chain carries the occultation date, 1981-08-26, with the source's own dating
in `epoch_note` where it is not already stated; one date, no per-source casework. Applied before
the acceptance commit so the chain is rebuilt once (deliverable 0, SPEC_00 v0.19, SPEC_01
v0.28). Finding 2: as
intended. Finding 3: a reader defect, enforced at Step 4 without a rebuild. Finding 4: accepted;
whether "the reviewing agent" is named in data attributes is a question to the author. Finding
5: correct, the column-wise precision is the envelope's (deliverable 3). Finding 6: deferred to
the Step 4 figures, with the gauge marker labeled by the tabulated gauge pressure. Decisions 1
to 11 accepted; the raw bundle's constant role, the `latitude_planetocentric_absent_meaning`
form and the every-difference refusal are recorded. The reviewing agent's own finding on the
Step 4 band clause is restated in §5. Ruling 5 of this section (committing `inputs/*.nc`)
remains open.

**Rulings on REPORT_03_step4 (REVIEW_03_step4, 16 September 2026).** Finding 1: correct; "the
top decade" of the expected values was the ten top levels, above 2 mbar, and the acceptance's
partition left 2 to 10 mbar under no bound; the bins are now stated in pressure (§5), and check
5 passes on the numbers as measured. Finding 2: F5 is left for the author's figure pass; no
restyling in this step. Finding 3: correct; SPEC_00 v0.20 §7.2 says relative paths and the
resolved namelist lives in the product's `namelist` group, no separate file (decision 6).
Finding 4: adopted; the closure namelist declares `date = "1981-08-26"` (Step 3 deliverable 2),
so the closure product carries the occultation date like every other file of the run; applied
before the acceptance commit and the acceptance rerun. Finding 5: recorded. Decisions 1 to 9
accepted as reported.

**Rulings of 15 September 2026 (v0.9).** The season identifier is added to Step 3 as
deliverable 0 rather than as a step of its own, because Step 3 already touches the schema for
kind `profile` and already builds the run's inputs, and the rebuild of the chain uses the §0
procedure established at Step 0. Nothing else from the seasonal design note enters SPEC_03 or
SPEC_04; the retrieval leg and the seasonal propagator are later specifications.

## 9. Revision history

| Version | Date | Change | Cause |
|---|---|---|---|
| 0.1 | 2026-09-14 | First draft: Step 0 (pressure grid, B1 projection, rebuild), Steps 1 to 4 (geopotential, hydrostatic, namelist and product kind, the closure) with the staggering table, the residual budget and the negative control; decisions 1 to 6, 8 to 10; amendments to SPEC_00, SPEC_01 and SPEC_02 listed for application at acceptance | handoff §8 and §11; the reviewing agent's independent closure of 14 September 2026; author decisions of 14 September 2026 |
| 0.2 | 2026-09-14 | "Closure" renamed the hydrostatic closure throughout, and stated to be a test of the production, not the model; the transfer identity tests placed in SPEC_04 after the transfer exists; `forward/production.py` with a `produce` function the transfer will reuse; the run gets its own `inputs/` built by the tools under the run prefix with `casspian-run-inputs`, and closure mode checks them against the anchor's embedded copies (decision 7); `[inputs]` required in the namelist as SPEC_00 §7.2 always said; §6: end-to-end tests placed in a separate SPEC_05 after the SPEC_04 build | author markup of v0.1 |
| 0.3 | 2026-09-14 | Accepted by the author; Step 0 proceeds; amendments applied to SPEC_00 v0.16, SPEC_01 v0.21 and SPEC_02 v0.10 | author acceptance of v0.2 |
| 0.13 | 2026-09-16 | REVIEW_03_step4: the closure statistics and the acceptance bounds stated in pressure bins (above 2 mbar, 2 to 10, 10 to 100, below 100 excluding the bottom row, bottom row, mean below 10 mbar); expected values gain the 2 to 10 mbar figure and the note on the deep trend; the closure namelist declares `date`; §8 rulings on REPORT_03_step4; status line: Step 4 accepted, SPEC_03 closes at its record; dependency SPEC_00 v0.20 | REVIEW_03_step4 |
| 0.12 | 2026-09-15 | Author decision on REVIEW_03_step3 finding 1: every file of the Lindal chain and of the closure run carries `epoch = "1981-08-26"`, the occultation date, the source's own dating in `epoch_note`; G and R build-file sections carry `epoch`; dependencies SPEC_00 v0.19, SPEC_01 v0.28 | author, 15 September 2026 |
| 0.11 | 2026-09-15 | REVIEW_03_step3: G and R `epoch` as ISO dates with the prose as `epoch_note`, `epoch` parsed by the reader (SPEC_01 v0.27); the raw bundle's role always `reduction`; `latitude_planetocentric_absent_meaning` on kind `profile`; the F6 envelope by column precision and stated to be the budget, not a bound; the gauge marker labeled by the tabulated value; Step 4: the band clause restated as a count, the kind C `source_statement` rule enforced; status line: Step 3 reviewed; §8 rulings on REPORT_03_step3 | REVIEW_03_step3 |
| 0.10 | 2026-09-15 | Coding agent's pre-execution review of Step 3: required `role` in build files (no default), the closure drop list as a rule, paths relative to the product for every writer with the reader resolving against the file's directory, acceptance restated; dependencies SPEC_00 v0.18, SPEC_01 v0.26 | coding agent's pre-execution review |
| 0.9 | 2026-09-15 | Step 3 deliverable 0: the season identifier across the chain (schema, transcription, rebuild, run inputs, namelist `solar_longitude_deg`, kind `profile` globals), with the closure equality and the thirteenth refusal case; dependencies to SPEC_00 v0.17, SPEC_01 v0.25, SPEC_02 v0.11; status line: Step 2 accepted at `994c787` | the seasonal design note of 15 September 2026; author decision |
| 0.8 | 2026-09-15 | Step 2: the linear-T test column's orientation stated (falling with height); status line: Step 2 accepted; §8 rulings on REPORT_03_step2 | REVIEW_03_step2 |
| 0.7 | 2026-09-14 | Step 4: the negative control defined as the one-factor construction `g_k (h_{k+1} − h_k)`, the fully radial reading reported beside it; status line: Step 1 accepted; §8 rulings on REPORT_03_step1 | REVIEW_03_step1 |
| 0.6 | 2026-09-14 | Step 0: `k` run endpoints 236 and 268; gravity bound 6e-5 with its decomposition; keep-by-content through `raw_bundle`; G and R not rebuilt; §0: the in-memory candidate rule for a step that changes an input; §8 rulings on REPORT_03_step0 | REVIEW_03_step0 |
| 0.5 | 2026-09-14 | Step 0: the snap rule restated once more so that adjacent ambiguous rows resolve one at a time from the fixed side; basis count corrected to 48; SPEC_01 v0.23 carries the same | coding agent's second reading of v0.4 |
| 0.4 | 2026-09-14 | Step 0: the grid snap rule restated (equal spacing with fixed neighbors, no spacing numbers in code), the basis count corrected to 62; Step 1: gauge match exact; Step 3: twelfth refusal case (geodesy key), compensated composition edit, `inputs/*.nc` ignored pending the author's decision, the Step 5 skipped-figures check retired; Step 0 acceptance: every affected suite; §8 rulings added, revision history moved to §9; SPEC_01 v0.22 carries the same snap rule | coding agent's pre-execution review |

---

## Appendix. Amendments to earlier specifications, applied at v0.3

Applied to SPEC_00 v0.16, SPEC_01 v0.21 and SPEC_02 v0.10 on 14 September 2026; kept here as
the record of what changed and why. The coding agent commits the author's spec edits in their
own commit before Step 0, as the procedure requires.

**SPEC_00 v0.16.** §3.3 item 4: "B3.2: absolute radius per level, Eq. B1 with the tilt,
`r = r0 + (h − h_ref) cos ψ`; the tabulated altitude is measured along the local vertical".
§5 `casspian_kind` gains `profile`. §6.1: optional `pressure_printed_Pa` and global
`pressure_grid_rule` for a source profile whose pressures are declared on a grid. §6.7:
`radius_m` and `height_above_anchor_isobar_m` defined as in Step 0 deliverable 2; the F5 and
F6 sentence of §6.7 (v0.10) removed. New §6.8, kind `profile`, as Step 3 deliverable 3. §2.3: the run directory gains the build control file and `figures/` under `output/`; every
input kind a run reads is its own file under the run prefix, never a reduction file or an
embedded copy, stated as a rule. §7.2:
`mode`, `p_b_rule` (exactly one of `p_b_rule` and `p_b_Pa`), `[diagnostics]`, `product`;
closure rules as Step 3 deliverable 2. §7.3: the geopotential and hydrostatic integrations
are closed by SPEC_03 Steps 1 and 2 (trapezoid on the smooth gravity magnitude; exact
log-linear layer integral of density); the vertical grid, latitude grid, the other five
integrations and the outer loop remain open for SPEC_04. §8: the `input_hashes` warning
implemented for every derived kind (SPEC_03 Step 3). §10 gains decisions 1, 2, 4 and 5 above
by reference.

**SPEC_01 v0.21.** Step 2: the raw bundle's `table1` group carries `pressure_printed_Pa`,
`pressure_Pa` on the declared grid, and `pressure_grid_applied`, under the `[pressure_grid]`
table of `lindal_scalars.toml`; acceptance values as Step 0 above (top row 19.9526 Pa, 65
snapped rows, the `k` sequence). Step 8: the composition tool reads `pressure_Pa`. Step 9: kind
T carries `pressure_grid_rule` and `pressure_printed_Pa`; the keep-by-content rule compares
`input_hashes` too.

**SPEC_02 v0.10.** Step 3, B1: `absolute_radius` takes `psi` and projects; acceptance: `radius_m`
at the top level `r0 + (h_top − h_ref) cos ψ`; `n` at the top level 1.041934e22 m⁻³ and at the
second level 1.270497e22; the 1 bar and 794.33 mbar checks unchanged. Step 4:
`reduction_record` gains `radius_projection_rule`, `radius_projection_residual_m`,
`latitude_drift_neglected_deg`. Step 5: F5 and F6 move to kind `profile`; kind N reports
nothing skipped. Decision 10 superseded by SPEC_03 decision 4. Decision 9 unchanged: every
number in it is unchanged by Step 0.
