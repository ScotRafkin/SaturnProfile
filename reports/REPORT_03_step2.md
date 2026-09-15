# REPORT 03, Step 2. `lib.hydrostatic`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 14 September 2026. Specification: `SPEC_03_Forward_Production.md` v0.7, Step 2.
Status: **accepted (REVIEW_03_step2, 15 September 2026).** Finding 1 restated at SPEC_03 v0.8:
the linear-T column falls with height, 140 K at the bottom level and 80 K at the top.

**Eight of eight acceptance checks pass**: the five the specification asks for (checks 1 to 5)
and three beyond it (6 to 8). As first run, check 2 read against the other orientation of the
linear-T column and failed on the trapezoid band; the specification's numbers belong to a
temperature falling with height (finding 1, decision 4).

## 0. Before this step

| Commit | What |
|---|---|
| `dd9c83b` | REVIEW_03_step1, SPEC_03 v0.7 (the Step 4 negative control restated), `STATE.md`. |
| `b3efc3d` | The Step 1 acceptance commit, `lib.geopotential`. |
| `aa1dfea` | The Step 1 record: no product changed, `STATE.md` accepted. |

The status line of SPEC_03 v0.7 ("Step 2 proceeds after the Step 1 acceptance commit") and
`STATE.md` were read before starting, and the author said go. The corrected REVIEW_03_step1 the
author mentioned is not in the working tree: the committed file still reads "still carries
`2149b64`" for every product. Nothing in Step 2 depends on it.

---

## 1. What was built

**`src/casspian/lib/hydrostatic.py`**, NumPy in and out, no I/O, constants from `lib.constants`:

- `density(N, R_bar, m_bar_kg_mol)`: `rho = (N / R_bar) m_bar`, Eq. B4, with `m_bar` in kg/mol as
  kinds N and C carry it, divided by the Avogadro constant.
- `layer_mass(rho, Phi)`: the exact integral of the log-linear density across each layer, in the
  form `rho_k (Phi_k − Phi_{k+1}) expm1(x) / x` with `x = ln(rho_{k+1} / rho_k)`, identical to the
  specification's expression (decision 1), with the limit `rho_k (Phi_k − Phi_{k+1})` where
  `|x| < 1e-10` (`LIMIT_LOG_RATIO`).
- `pressure_from_top(p_b, I)`: `p_0 = p_b`, `p_{k+1} = p_k + I_{k+1/2}`.
- `temperature(p, N, R_bar)`: `T = p R_bar / (k_B N)`, Eq. B6.

The module docstring states the layer rule and its reason and the top down direction. No other file
under `src/` changed.

## 2. Acceptance results

Run by `reports/step03_2/accept_step03_2.py`; full output in `reports/step03_2/output.txt`. The
geopotential grid is the Step 1 grid on the swept product (`2149b64`, 66 levels, `Phi` from
2.852355e6 to −1.043743e6 m²/s², gauge level 29), from the accepted `lib.geopotential` with the
reduction's gravity, rotation and wind. The test columns hold the product's top level composition
uniform (`M` = 2.13508332e-3 kg/mol, `R_bar` = 4.836272e-30 m³), take `p_b` = 20 Pa at the top,
and are built from their analytic pressure, with `N` formed from it and `rho` through the module's
own `density`.

| Check | Measured |
|---|---|
| 1. Isothermal column (`T` = 100 K, `p_b` = 20 Pa): `p` to 1e-14 and `T` = 100 K to 1e-14 relative at every level | **Pass.** Scale geopotential `k_B T / m` = 3.894210e5 m²/s², `p` from 20 to 4.427e5 Pa. Largest `\|p / p_analytic − 1\|` **1.1e-15** (level 56); largest `\|T / 100 − 1\|` **1.1e-15** (level 50). The trapezoid alternative, computed in the script only: **8.475e-3**. |
| 2. Linear-T column: `p` to 6e-4 or better; the trapezoid alternative 8e-3 to 9e-3, both reported | **Pass** on the column falling with height, 140 K at the bottom to 80 K at the top: log-linear **5.987e-4** (level 6), trapezoid **8.785e-3**. The other orientation, reported beside it: log-linear 3.484e-4, trapezoid 6.507e-3, outside the band (finding 1). `T` is linear in `Phi`; `p` from 20 to 2.257e5 Pa. |
| 3. Halving every layer reduces the linear-T error by 3.5 to 4.5 | **Pass.** 5.9870e-4 on the Step 1 grid, 1.4983e-4 with every layer halved (midpoints in `Phi`, the analytic density there): ratio **3.996**; at the bottom level 3.999. The other orientation: 3.995. |
| 4. `layer_mass` returns the limit form on equal densities and refuses a zero-thickness layer with the layer named | **Pass.** Equal densities 2e-3 kg/m³ across 1e5 m²/s²: 200.0, the limit form exactly. A log ratio of 5e-11, inside the limit: −2.5e-11 relative to the exact form. Refused: a zero thickness layer ("does not decrease across layer 6 (levels 6 and 7, ...)"), a column that is not top down ("across layer 0 ..., and 64 more layers"), a zero density ("levels [3]"), a negative density ("levels [4]"). |
| 5. `temperature` on the reduction's own `N`, `R_bar` and tabulated `p` returns the tabulated `T` to 1e-12 | **Pass.** Largest `\|T / T_tab − 1\|` **4.4e-16** over 66 levels. |
| 6. Beyond the specification: `density` against the reduction's number density | **Pass.** `density(N, R_bar, m_bar)` equals `n m_bar / N_A` from `number_density_m3` to 2.2e-16; `rho` from 3.694e-5 to 0.2282 kg/m³. |
| 7. Beyond the specification: the `expm1` form against the literal difference form | **Pass.** Relative to `rho_k dPhi expm1(x) / x`: at log ratio 1e-2 the module 0.0, the literal form −1.0e-14; at 1e-6, −3.3e-16 and −2.1e-11; at 1e-8, +2.2e-16 and +2.1e-9; at 1e-9, 0.0 and **−8.7e-8** (decision 1). |
| 8. Beyond the specification: `pressure_from_top` refusals; `p_0 = p_b` | **Pass.** `p_b` = 0 and NaN refused; a negative layer mass refused with the layer named; `p_0` equals `p_b` exactly. |

## 3. Findings

1. **"Falls linearly from 140 to 80 K" has two readings, and the specification's numbers belong to
   one.** Read top down (140 K at the top, 80 K at the bottom), the log-linear integral gives
   3.48e-4 and the trapezoid 6.51e-3, outside the stated 8e-3 to 9e-3. Read as falling with height
   (140 K at the bottom, 80 K at the top), they are 5.99e-4 and 8.79e-3, the specification's "6e-4"
   and "8.5e-3". The check runs on the second reading; both are reported. Proposed: the
   specification says "falling with height, from 140 K at the bottom level to 80 K at the top,
   linear in `Phi`".
2. **The isothermal column recovers to 1.1e-15**, which the specification's section on the rule
   quotes as 1e-15 and its acceptance bounds at 1e-14. No action.

## 4. Decisions

1. **The layer integral is computed as `rho_k dPhi expm1(x) / x`.** It is the specification's
   `(rho_{k+1} − rho_k) dPhi / ln(rho_{k+1} / rho_k)` exactly, since `rho_{k+1} − rho_k = rho_k
   expm1(x)`, but the literal form divides the difference of two nearly equal densities by a small
   logarithm and loses digits as the log ratio shrinks: −8.7e-8 at a log ratio of 1e-9, far above
   the 1e-10 at which the specification switches to the limit (check 7). The limit branch is kept
   at `|x| < 1e-10` as specified. The log ratios of the Table I layers themselves were not examined.
2. **Refusals the specification did not list.** `layer_mass` refuses a non-finite density or
   geopotential and, beyond a zero thickness layer, any layer across which the geopotential does
   not decrease, so a bottom up column is refused rather than integrated with negative masses;
   the message names the first such layer and counts the rest. `pressure_from_top` refuses a
   boundary pressure that is not positive and finite and any layer mass that is not positive and
   finite, naming the layers.
3. **`density` takes the mean molar mass in kg/mol**, the unit and variable kinds N and C carry
   (`mean_molar_mass_kg_mol`), and converts to kg per molecule inside, as the specification
   describes, so that no caller forms `m_bar` per molecule separately.
4. **The linear-T test column falls with height** (finding 1), linear in `Phi` between 140 K at the
   bottom level and 80 K at the top; its pressure is the closed form
   `p = p_b (T / T_top)^(−m / (k_B beta))` with `beta = dT/dPhi`.
5. **The test columns' composition is the product's top level composition held uniform**, and
   their geopotential grid is the Step 1 grid on the swept product; halving inserts the arithmetic
   midpoint in `Phi` of every layer with the analytic density there, and the error is compared at
   the original levels.

## 5. Regression

A new `lib` module, so the full regression was run (`reports/step03_2/run_regression.sh`,
results in `reports/step03_2/regression.txt`, each suite's output beside it). **All pass:** 6, 14,
8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 6, 7, 9, 9, 7, 7 for SPEC_02 Steps 1 to 6, on the clean
products of the Step 0 sweep.

Steps 02_4 to 02_6 rewrite `lindal_refractivity.nc` and its figures and Step 02_5 recopies the
committed Step 5 report figures, so the swept product and figures were copied aside and restored:
SHA-256 `64c5d01a...` before and after; `reports/figures/` restored from git. The Step 1 acceptance
was not rerun, since Step 2 changes no module it uses; its grid is rebuilt inside the Step 2
acceptance from the same accepted `lib.geopotential`.

No em dash or en dash appears in any file written in this step.

**Sweep** (REVIEW_03_step2 order of work, item 3), after the acceptance commit `994c787`: no input
or product changed in this step. The `-dirty` scan over `occul_data/lindal/` finds every product
clean at the commit it was swept at: `lindal_raw.nc`, `lindal_thermo.nc`, `lindal_geodesy.nc`,
`lindal_wind.nc`, `lindal_composition.nc` and `lindal_refractivity.nc` at `2149b64`,
`lindal_gravity.nc` and `lindal_rotation.nc` at `ece58d2`; SHA-256 prefixes 07e0727c, 5705f1dc,
51c545a6, 26fe83a3, 6af9dbfb, 64c5d01a, 6db0129c, a5c72017, each equal to REPORT_03_step0 section 7.
Nothing to rebuild.

## 6. Next step

Step 2 is accepted (`994c787`). Step 3 proceeds; the author said go.
