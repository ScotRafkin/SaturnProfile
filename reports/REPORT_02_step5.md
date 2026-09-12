# REPORT 02, Step 5. The standard diagnostics, `casspian.tools.plots`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.7, Step 5,
against `SPEC_00_Architecture_and_Data_Files.md` v0.14, sections 3.5 and 7.1. Status: accepted with two
figure changes (REVIEW_02_step5), applied; see section 7. **Seven of seven acceptance checks pass.** The figures are attached below for
the check the specification reserves for the author.

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
product: the figures are regenerable and ignored by git. The attached copies under
`reports/figures/` are the committed record.

## 0. Before this step, per REVIEW_02_step4

| Commit | What |
|---|---|
| `ed076e9` | SPEC_02 Step 4 with the v0.7 change: the anchor rule spread by a full fixed point rerun under each rule. The rerun gives `r0` 58,537.683 and 58,502.116 km and `phi_c` 30.801965 and 30.807925 deg, all within the v0.7 expected values. The acceptance passes 9 of 9. SPEC_00 v0.14 and SPEC_02 v0.7 were committed with it. |
| `4a9aca2` | The sweep: `lindal_refractivity.nc` rebuilt clean at `ed076e9`, SHA-256 `327f8141...`, recorded in REPORT_02_step4 section 7. |

The status lines of SPEC_02 v0.7 ("Step 5 proceeds") and `STATE.md` were read before starting.

---

## 1. What was built

**`src/casspian/tools/plots/`**, a generic tool tied to no profile.

- **`render.py`** holds `render(path, out_dir=None, format="png", dpi=150)` and the entry point
  **`casspian-plots <file.nc> [--out DIR] [--format png|pdf] [--dpi N]`**.
  - `render` reads `casspian_kind` and dispatches.
    - Kind N gets F1 to F6, those its fields allow.
    - A kind W, T or C file gets its single figure view.
    - Every other kind is refused by name. The CLI exits 2 and writes nothing.
  - It returns a `RenderResult` with each figure written, the combined PDF, each figure skipped
    and why, each footer, and the arrays each figure drew, so an acceptance can check a panel
    against an independent computation.
- **`style.py`** is the one style: fonts, line widths, the color per quantity, the provenance
  markers, layout margins, the anchor isobar and latitude lines, the pressure axis (log,
  decreasing upward, mbar on the label) and the footer. It is applied through
  `matplotlib.rc_context`, so a caller's settings are untouched. There is no seaborn and no
  external stylesheet.
- **`figures_product.py`** draws F1 to F6 from the product and its embedded inputs.
- **`figures_inputs.py`** holds the shared panels (temperature, mole fractions, mean properties,
  wind by latitude) and the single figure views of W, T and C.

**`refrac/product.py`**: `casspian-refrac` renders the standard figures after writing kind N
when the manifest's `[diagnostics] figures = true`. The figures go to `figures/` beside the
product in the declared format and dpi, and it prints what was written and what was skipped.
matplotlib is imported only then.

**`.gitignore`** gains `occul_data/*/figures/`. **`pyproject.toml`** gains `casspian-plots`.

**The figures, as drawn for Lindal:**

- **F1, inputs.**
  - `T(p)`; Lindal states no uncertainty, and the title says so.
  - `u(p)` at `phi_c`, with an inset of `u(phi)` on the reference level and `phi_c` marked.
  - `x_i(p)` on a log axis, NH3 provenance by marker.
  - `m_bar(p)` and `R_bar(p)` on twin axes.
- **F2, gravity along the profile.** Computed by `lib.gravity` on the embedded G, R and W at
  `phi_c` and `radius_m`:
  - `g_N` and `g_eff`;
  - the difference split into the rigid centrifugal part and the wind part, on separate
    scales, since the wind part is 6.1e-4 m/s2 against 1.16;
  - `G_phi` and `psi`.
- **F3, across latitude on the anchor isobar.** Computed by `lib.geoid` and `lib.gravity`:
  - `g_N` and `g_eff` on the wind geoid;
  - the no wind and wind geoids, with `r0` from the file marked and both polar radii annotated;
  - the dynamical height, and `psi(phi)`.
- **F4, the product.**
  - `n(p)`; its uncertainty is not stated, and the title says so.
  - `N(p)` with its band.
  - `N` against `radius_m`.
  - The recovered temperature as a fractional difference. The axis is at least +-2e-12 and
    widens to any larger value, with the +-1e-12 bound shaded.
- **F5, geopotential** and **F6, hydrostatic closure** are skipped for kind N as it stands, and
  named with the field each needs.

## 2. Decisions

1. **The acceptance runs `casspian-refrac` on a copy of the manifest and inputs.** The committed
   `lindal_reduction.toml` has no `[diagnostics]` section. Adding one would change its SHA-256.
   That hash is recorded in REPORT_01_step9, in the product's `input_hashes` and
   `reduction_record`, and checked by the Step 1 acceptance, and the product would have to be
   rebuilt. So check 1 appends `[diagnostics]` to a copy beside copies of the six inputs, as
   the earlier refusal cases do. **For the author:** whether the Lindal manifest should carry
   `[diagnostics] figures = true` is a choice for the next time the product is rebuilt anyway.

2. **How "identical apart from the generation time" is checked.** The generation time is drawn
   alone in a band of 2.5 percent of the figure height at the very bottom (`FOOTER_TIME_BAND`),
   and nothing else is drawn there. The acceptance compares the two renderings as pixel arrays
   with that band masked. The PNG bytes cannot be compared directly, since the time is pixels in
   the image. Everything above the band is identical, and the band itself differs, as it
   should. With the same timestamp passed to both renderings the files are byte identical.

3. **Figure names.**
   - Kind N: `<prefix>_diag_F<k>_<name>.<format>`, with `name` one of `inputs`,
     `gravity_profile`, `gravity_latitude`, `product`, `geopotential`, `hydrostatic`.
   - An input view: `<prefix>_diag_<kind>.<format>`.
   - **The combined PDF is written for kind N only.** The specification states it for the kind
     N set. A single figure view is its own one page file.

4. **The wind view's second panel is `u(p)` at five latitudes** (0, +-30, +-60 deg). A kind W
   file offered alone has no `phi_c`, so the F1 panel "u(p) at phi_c" has no latitude to use.
   The kind T and C views likewise draw no anchor isobar, which only a manifest names.

5. **`--format` and `--dpi` are on the CLI** beside the specified `--out`, defaulting to png and
   150, the manifest example's values. `render` takes both, so the CLI exposes them.

6. **F5 and F6 are implemented for the per level fields the specification names** and were
   exercised on a synthetic copy of the product carrying them. What is not implemented:
   - F5's "Phi on the anchor isobar against latitude";
   - F6's "declared boundary pressure and its level".

   Both need field names SPEC_03 has not defined. F6 marks a boundary if the field carries a
   `boundary_pressure_Pa` attribute; that name is a placeholder until SPEC_03 fixes it.

7. **F2 interpolates the wind to the product's levels linearly in log pressure**, between the
   kind W pressure nodes, at `phi_c` linear in latitude as Step 2 does. For Lindal the wind is
   altitude independent, so any vertical rule gives the same value.

8. **Layout was checked by eye before the acceptance was final.** The first rendering had a
   panel title colliding with an axis label, legends over the NH3 points and the `g_eff`
   minimum, polar annotations below the axes, and the wind part of F2 lying on zero on a shared
   axis. None of this is caught by an acceptance check. All were fixed and the figures below
   are the result.

## 3. Acceptance results

Run by `reports/step02_5/accept_step02_5.py`; full output in `reports/step02_5/output.txt`.

| Check | Measured |
|---|---|
| 1. `casspian-refrac` with `figures = true` | **Pass.** Exit 0. Writes F1 to F4 as PNG and `lindal_diag.pdf` with 4 pages. Prints `skipped F5_geopotential: needs geopotential_m2s2(level)` and `skipped F6_hydrostatic: needs pressure_hydrostatic_Pa(level)`. |
| 2. `casspian-plots` by hand writes identical figures | **Pass.** All four figures are 1650 x 1275 px and identical above the 32 px time band; the band differs. |
| 3. Every figure carries the footer | **Pass.** File name, `casspian_git_commit` and SHA-256 prefix on every figure, the generation time below them; ink present in both footer bands of every PNG. |
| 4. The wind view, and the gravity file refused | **Pass.** `casspian-plots lindal_wind.nc` writes `lindal_diag_wind.png`. `casspian-plots lindal_gravity.nc` exits 2 with `casspian_kind 'gravity' has no diagnostic figures` and writes nothing. |
| 5. F4's recovered temperature below 1e-12 | **Pass.** max `|T_recovered / T - 1|` = 4.4e-16 over 66 levels. The panel's array equals an independent evaluation. |
| 6. F2's `g_eff` at 1 bar against `lib.gravity` directly | **Pass.** Panel 9.994279004116628 m/s2, direct 9.994279004116628 m/s2, difference 0.0. The direct call reads the files on disk, not the embedded copies. |
| 7. Beyond the specification | **Pass.** The default `figures/` beside the product is ignored by git. The T and C views render. A copy of the product carrying synthetic `geopotential_m2s2` and `pressure_hydrostatic_Pa` renders F1 to F6 with nothing skipped. |

## 4. The figures, for the author's review

F1 to F4 were rendered by `casspian-plots` from the committed product
(`occul_data/lindal/lindal_refractivity.nc`, `ed076e9`, clean). The three views were rendered
from the committed input files.

| Figure | File |
|---|---|
| F1, inputs | `reports/figures/step02_5_lindal_diag_F1_inputs.png` |
| F2, gravity along the profile | `reports/figures/step02_5_lindal_diag_F2_gravity_profile.png` |
| F3, gravity and shape across latitude | `reports/figures/step02_5_lindal_diag_F3_gravity_latitude.png` |
| F4, the product | `reports/figures/step02_5_lindal_diag_F4_product.png` |
| wind view | `reports/figures/step02_5_lindal_diag_wind.png` |
| thermo view | `reports/figures/step02_5_lindal_diag_thermo.png` |
| composition view | `reports/figures/step02_5_lindal_diag_composition.png` |

What the figures show that no number in the acceptance does:

- **F2:** at `phi_c` the wind contributes 6.1e-4 m/s2 to `g_N - g_eff`, against 1.16 m/s2 from
  rigid rotation. That is 0.05 percent, as expected for 2 m/s.
- **F3:**
  - The dynamical height is 65 km at the anchor latitude and 127 km at the equator.
  - The wind geoid carries the 28.7 km polar asymmetry, with the south pole higher.
  - `g_eff` is flat between about -5 and +5 deg, where the wind dips to 430 m/s between its two
    490 m/s peaks.
- **F1:**
  - The mean molar mass is constant above the tabulated ammonia and rises by 0.055 percent into
    the deepest levels.
  - The mean refractivity falls by the matching factor, since the `lindal1985` set carries NH3
    at zero.

## 5. Regression

The ten SPEC_01 suites and the SPEC_02 Step 1 to 4 suites were rerun after this step's changes (`reports/step02_5/regression.txt`). All pass: 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 6, 7, 9, 9 for SPEC_02 Steps 1 to 4. The Step 4 suite drives `casspian-refrac`, whose entry point gained the diagnostics hook; its manifest carries no `[diagnostics]`, so it renders nothing, as before.

No em dash or en dash appears in any file written in this step.

## 6. Next step

SPEC_02 is complete with this step once it is accepted. Nothing further is started.

---

## 7. Addendum after review: the two figure changes (SPEC_02 v0.8)

Added 12 September 2026, per `REVIEW_02_step5.md`. Code only in `tools/plots/figures_product.py`.

1. **F4, the refractivity panel.** The fractional uncertainty `refractivity_uncertainty /
   refractivity` is now drawn as a second curve on a linear top axis, in percent, labeled
   "fractional uncertainty of N (percent)", with the +-1 sigma band kept. For Lindal it reads
   2.33 percent at every level; the axis runs from zero to 1.5 times the largest value, so the
   curve sits two thirds of the way across. When a later product carries a visible band, both
   read correctly.
2. **F1, the wind inset.** Its horizontal axis is labeled "planetocentric lat. (deg)".

The Step 5 acceptance was rerun after the change: seven of seven, every measured value as in
section 3. F1 to F4 were rerendered from the committed product and the four copies under
`reports/figures/` replaced.

**The sweep.** After the acceptance commit, `lindal_refractivity.nc` was rebuilt on a clean tree.
The Step 4 regression suite had rewritten it from the working tree, so it carried
`4a9aca2...-dirty`. This hash supersedes `327f8141...` of REPORT_02_step4 section 7.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_refractivity.nc` | refractivity | `a91f75a804ba4bf905be6df3c05049b773cdfbbb532e6aa1f0871b528db3a740` |
