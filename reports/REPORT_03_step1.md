# REPORT 03, Step 1. `lib.geopotential`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 14 September 2026. Specification: `SPEC_03_Forward_Production.md` v0.6, Step 1.
Status: **accepted (REVIEW_03_step1, 14 September 2026).** Finding 1 restated the Step 4 negative
control at SPEC_03 v0.7: the one-factor construction is the control, the two-factor reading is
reported beside it and must fail the same bound. Finding 2 noted, no action.

**Nine of nine acceptance checks pass**: the six the specification asks for (checks 1 to 6) and
three beyond it (4b, 7, 8). Every expected value of the specification is reproduced to the
figures printed.

## 0. Before this step

| Commit | What |
|---|---|
| `73ae948` | REVIEW_03_step0, SPEC_03 v0.6, SPEC_01 v0.24, `STATE.md`. |
| `2149b64` | The Step 0 acceptance commit. |
| `5cdf07e` | The Step 0 sweep record: the chain rebuilt clean at `2149b64`, hashes in REPORT_03_step0 section 7, Steps 02_1, 02_4, 02_5 and 02_6 on the clean products, `STATE.md` accepted. |

The status line of SPEC_03 v0.6 ("Step 1 proceeds after the Step 0 acceptance commit and sweep")
and `STATE.md` were read before starting, and the author said go. The Step 1 acceptance runs on
the clean product of the Step 0 sweep (`lindal_refractivity.nc`, commit `2149b64`, SHA-256
`64c5d01a...`) with the gravity, rotation and wind the reduction used, loaded through
`lib.control`; no relaxation was needed.

---

## 1. What was built

**`src/casspian/lib/geopotential.py`**, NumPy in and out, no I/O, angles in radians, as the
specification's deliverable lists:

- `effective_gravity_magnitude(u, r, phi_c, Omega, GM, J, degrees, R_norm)` returns
  `(g_mag, g, G_phi, psi)` on `r`, from `lib.gravity.g_eff_vector`; `u` and `phi_c` broadcast to
  the shape of `r`.
- `layer_increments(g_mag, h)` returns `dPhi_{k+1/2} = (|g_eff|_k + |g_eff|_{k+1}) / 2 (h_{k+1} −
  h_k)`, the trapezoid.
- `geopotential(dPhi, gauge_index)` returns `Phi_k`, exactly zero at the gauge, with every other
  level the sum of the layer increments between the gauge and itself, taken outward from the gauge
  (decision 1).
- `geopotential_along_profile(u, r, h, phi_c, gauge_index, Omega, GM, J, degrees, R_norm)` runs the
  three in sequence and returns a frozen `GeopotentialProfile` with `Phi`, the layer increments,
  `|g_eff|`, `g`, `G_phi`, `psi` and the gauge index.

The module docstring states the staggering of the specification's table, the field-line rule and
its reason, and the gauge convention. No other file under `src/` changed.

## 2. Acceptance results

Run by `reports/step03_1/accept_step03_1.py`; full output in `reports/step03_1/output.txt`. The
Lindal column: 66 levels, gauge level 29 at 10,000.0 Pa, `phi_c` 30.80556842739218 deg, `u(phi_c)`
2.16708095051181 m/s from the wind file's reference level by `refrac.anchor.wind_of_latitude`,
equal to `u_at_anchor_ms` in the product's record.

| Check | Measured |
|---|---|
| 1. Uniform gravity returns `Phi_k − Phi_a = g (h_k − h_a)` to 1e-12 relative | **Pass.** `J = 0`, `Omega = 0`, `u = 0`, `r = 1e18 m + (h − h_a)`, `GM = 1e37`, so `GM / r^2` varies by 7.8e-13 across the column. Largest relative departure over the 65 levels off the gauge **2.9e-13**; `Phi_a` = 0.0. |
| 2. `Phi_a` is exactly 0.0 | **Pass.** 0.0 on the Lindal column and on the uniform one. |
| 3. `Phi` strictly monotonic in `h` | **Pass.** The sign of every one of the 65 layer increments equals the sign of the height increment; `Phi` decreases with level index, increasing upward. |
| 4. The expected values to the figures printed | **Pass, all eight.** `\|g_eff\|` top 9.892476651 (9.892477), anchor 10.005624644 (10.005625), bottom 10.047103177 (10.047103) m/s²; `psi` top 5.546373691 (5.546374), bottom 5.475834863 (5.475835) deg; `Phi` top 2,852,354.97 (2.852355e6), bottom −1,043,742.65 (−1.043743e6) m²/s²; first layer −126,655.82 (−1.2666e5). |
| 4b. Beyond the specification: which radial reading gives the specification's radial values | **Pass.** The radial component times the altitude increment, `g_k (h_{k+1} − h_k)`, gives top **2,839,126.33** and bottom **−1,038,963.41**, the specification's 2.839126e6 and −1.038963e6; the field-line result exceeds it by 4.659e-3 and 4.600e-3, against `1 / cos psi − 1` = 4.616e-3 at the anchor. The radial component times the radial increment, `g_k (h_{k+1} − h_k) cos psi`, gives 2,826,081.99 and −1,034,189.90, a departure of 9.30e-3 and 9.24e-3, which is `1 / cos² psi − 1` = 9.253e-3 (finding 1). |
| 5. `\|g_eff\|` at the anchor equals `g_eff_radial / cos psi` to round-off | **Pass.** 10.00562464418385 both ways, relative difference **0.0**. `psi` at the anchor level 5.494431541512493 deg, equal to the product's frozen `psi_deg`. |
| 6. `geopotential` refuses a gauge index off the grid | **Pass.** Refused: −1, 66 and 71 ("off the grid of 66 levels (0 to 65); the gauge must be a level"); 2.5, `True` and `np.float64(29.0)` ("not an integer level index"). `np.int64(29)` is accepted and gives the same result as `29`. |
| 7. Beyond the specification: `layer_increments` refuses a malformed column | **Pass.** A zero thickness layer and a reversed layer ("h is not strictly monotonic", with the layer named), a shape mismatch, a single level, a NaN height, each refused with its message. |
| 8. Beyond the specification: the trapezoid error per layer is below 1e-7, as the staggering table states | **Pass.** Against Simpson with `\|g_eff\|` at the midpoint radius, the largest per layer error is **2.85e-8**, on the thickest layer (13.30 km, layer 2), against the table's estimate `(dh)² / (2 r²)` = 2.56e-8 (the measured error is 1.11 times the estimate; finding 2). Effect on `Phi` at the top +0.057 m²/s² (2.0e-8), at the bottom −0.0019 m²/s² (1.9e-9). `\|g_eff\|` varies by 1.56 percent over the profile. |

F5 is not rendered at this step, as the specification says.

## 3. Findings

1. **The two radial constructions in SPEC_03 are not the same construction.** Step 1 quotes the
   radial rule as 2.839126e6 and −1.038963e6, "and the difference is the 4.6e-3 of §0". Those
   numbers are the radial component times the **tabulated altitude increment**, `g_k (h_{k+1} −
   h_k)` (check 4b), the pre-Step 0 reading of B1 and B2 in which `h` was taken as radial. Step 4
   defines its negative control as "the radial component and the radial increment (`g_k`,
   `(h_{k+1} − h_k) cos ψ`)". After Step 0 the radial increment carries `cos ψ` already, so that
   construction departs from the field-line geopotential by `1 / cos² ψ − 1` = 9.25e-3, twice the
   4.6e-3, at every level. The negative control's expected mean of −5.1e-3 below 10 mbar and
   −7.8e-3 at the bottom row look like the one-factor construction measured before Step 0, when
   `r_{k+1} − r_k` equalled `h_{k+1} − h_k`. Nothing in Step 1 depends on it. Raised now so that
   Step 4's negative control is settled before that step: either `g_k (h_{k+1} − h_k)` (the
   manuscript's present mixture, one factor) or `g_k (r_{k+1} − r_k)` (radial throughout, two
   factors). No hydrostatic closure was run at this step, so the residual either construction
   gives, and whether each fails the −4e-3 bound on the control's mean, is not measured here.

2. **The trapezoid error estimate is low by about 11 percent; the stated bound holds.** The
   staggering table estimates the layer error as `(dh)² / (2 r²)`, from `g'' ≈ 6 g / r²`, and states
   it below 1e-7. Measured against Simpson it is 1.11 times the estimate on the thickest layers,
   2.85e-8 at most, so the "below 1e-7" statement holds with a factor of 3.5 to spare. Why the
   estimate is low was not investigated. No change proposed beyond noting it.

## 4. Decisions

1. **The gauge sum runs outward from the gauge.** `Phi_k` for `k > a` is the cumulative sum of
   the increments from `a` down to `k`, and for `k < a` the negative cumulative sum from `a` up to
   `k`, as the staggering table words it, rather than one cumulative sum over the whole column
   with its value at `a` subtracted. `Phi_a` is then exactly zero by construction (not by
   cancellation), and no level near the gauge inherits the round-off of the far end of the
   column.
2. **`layer_increments` refuses what the specification did not mention**: fewer than two levels,
   arrays that are not one dimensional on the same levels, non-finite values, and heights that are
   not strictly monotonic (a zero thickness or reversed layer), each named. Step 2's `layer_mass`
   refuses a zero thickness layer; the geopotential refuses it first.
3. **`geopotential` refuses a gauge index that is not an integer**, including a float that
   happens to be integral (`np.float64(29.0)`) and `True`, as well as a negative index, which
   NumPy would otherwise accept silently by counting from the end. A NumPy integer is accepted.
4. **`effective_gravity_magnitude` broadcasts `u` and `phi_c` to the shape of `r`** and returns
   `(g_mag, g, G_phi, psi)` in the specification's order, although `lib.gravity.g_eff_vector`
   returns `(g, G_phi, magnitude, psi)`.
5. **`geopotential_along_profile` returns a frozen dataclass** carrying the levels and layers
   together, so that Step 4 can write the product and F5 from one object; the specification named
   the function but not its return type.
6. **`u(phi_c)` in the acceptance** is read from the wind file's reference level by
   `refrac.anchor.wind_of_latitude`, the function the reduction uses, and checked against
   `u_at_anchor_ms` in the product's record. How the forward model obtains `u` is Step 4's.

## 5. Regression

A new `lib` module, so the full regression was run (`reports/step03_1/run_regression.sh`,
results in `reports/step03_1/regression.txt`, each suite's output beside it). **All pass:** 6, 14,
8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 6, 7, 9, 9, 7, 7 for SPEC_02 Steps 1 to 6, on the clean
products of the Step 0 sweep. No `-dirty` relaxation was needed: Step 1 changes no input file.

Steps 02_4 to 02_6 rewrite `lindal_refractivity.nc` and its figures from the working tree, and
Step 02_5 recopies the committed Step 5 report figures. The swept product and figures were copied
aside before the suites and restored after them: the product's SHA-256 is `64c5d01a...` before and
after, and `reports/figures/` was restored from git. The Step 03_0 acceptance was not rerun; it
builds a candidate against the Step 0 `before/` copies, which no longer describe the chain, and
Step 1 changes nothing it tests.

No em dash or en dash appears in any file written in this step.

**Sweep** (REVIEW_03_step1 order of work, item 3), after the acceptance commit `b3efc3d`: no input
or product changed in this step. The `-dirty` scan over `occul_data/lindal/` finds every product
clean and unchanged from REPORT_03_step0 section 7: `lindal_raw.nc`, `lindal_thermo.nc`,
`lindal_geodesy.nc`, `lindal_wind.nc`, `lindal_composition.nc` and `lindal_refractivity.nc` carry
`2149b64`; `lindal_gravity.nc` and `lindal_rotation.nc` carry `ece58d2`, kept at Step 0
(REVIEW_03_step0 decision 5). SHA-256 prefixes 07e0727c, 5705f1dc, 51c545a6, 26fe83a3, 6af9dbfb,
64c5d01a, 6db0129c, a5c72017, each equal to section 7 of REPORT_03_step0. Nothing to rebuild.

## 6. Next step

Step 1 is accepted (`b3efc3d`). Step 2, `lib.hydrostatic`, starts when the author says go.
