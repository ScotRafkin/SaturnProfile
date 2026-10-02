# REPORT 05, Step 2. The column march, faster (decision R, first part)

Coding agent, 1 October 2026. Specification: `docs/specs/SPEC_05_Shear_Experiments.md` v0.6, §2a.
Working tree on `main` at `99931dc`; nothing of this step is committed.

**Result.** No value moved. After fix 1 and again after fix 2, all 518 variables of the 18
registered files, rebuilt as candidates from the working tree, are bit-identical to the committed
files (`tests/step05_2/accept_step05_2.py`, 3 of 3). The full regression, all 28 suites with
`step04_4`, passes at every reference count in 95 minutes on Windows, against about 8 hours before.
A production transfer run takes 12.7 s, against 94 s; `lib.mesh.build_columns` takes 0.05 s of it,
against 78 s. No registered file changes, so none is recommitted.

## 1. What changed

| File | Change |
|---|---|
| `src/casspian/lib/gravity.py` | Fix 1. `validated_degrees`: the degree list checked once and marked, so later calls skip the check. One pass of the Legendre recurrence for `P_l` and `dP_l/dx` together (it ran twice). `HarmonicFactors` / `harmonic_factors`: the latitude-only factors (`sin`, `cos`, `cos^2`, `P_l`, `dP_l/dx`, `(l + 1) J_l`) computed once. `g_eff_from_factors`: the arithmetic of `g_newton`, `G_phi_newton`, `omega_abs`, `g_eff_radial` and `G_phi_eff` term for term and in their order, with the series once for both components and `Omega_abs` once (each ran twice per evaluation). `g_eff_vector` now goes through it. |
| `src/casspian/lib/mesh.py` | Fix 2. `_march`: RK4 in `Phi` with `r` and `z` as arrays over every latitude, one gravity evaluation per stage for every column; the node loop stays, the column loop goes. `column` is its one-latitude case; `build_columns` marches every latitude at once and takes an optional `wind_all`. `stage_geopotentials`: every `Phi` the march evaluates, in its own arithmetic. Signatures kept, so no accepted suite or helper changes. |
| `src/casspian/forward/transfer.py` | Fix 2. `IsobarMap.ln_p_table`, the all-columns form of `ln_p_at`. `_wind_table`: the wind of every column at every stage `Phi`, from one `ln p` table and one `WindField.wind_at` call, passed to `build_columns` as `wind_all`; a point outside the wind grid is refused with the latitude, the geopotential and the pressure named (deliverable 4). The degree list validated once where the gravity file is read. |
| `src/casspian/forward/production.py`, `src/casspian/refrac/anchor.py` | The degree list validated once where the gravity file is read (fix 1). |
| `tests/step05_2/accept_step05_2.py` | The acceptance, below. |
| `tests/run_regression.sh` | The row `step05_2/accept_step05_2 3`. |

**Why the arithmetic is the same.** Fix 1 repeats each formula's operations in their order on
values that are computed once instead of many times. Fix 2 puts the same operations on arrays.
Three things could still have moved a last bit, and each was measured before the march was written:

- The sum over harmonic degrees: a sum over one point's terms (a 1-D reduction) and a sum down the
  first axis of an array (row by row) both associate left to right for 3 and 6 terms in this NumPy
  (2.4.2); 2,000 columns of random terms, every one identical.
- `exp`, `log`, `sin`, `cos`, integer `power`, `hypot` and `arctan2` on arrays and on single values:
  identical on 20,000 values each, in the ranges the model uses, on this machine.
- `np.interp` with an array of points and with one: the same interval and the same slope formula;
  the transfer product's 177 variables confirm it.

On another platform the array routines may round differently from the single-value ones; the 1e-14
rule of REVIEW_05_step0 is the bound there.

## 2. The acceptance

`tests/step05_2/accept_step05_2.py` copies `data_static/`, the Lindal reduction's control files
and raw sources and the two runs' control files into `reports/step05_2/<stage>/tree/`, rebuilds the
18 registered files there in the sweep's order with `_refuse_dirty_commit` relaxed for the script
only, and compares every variable of every group with the committed files. Nothing under
`occul_data/` or `forward/` is written.

| Check | After fix 1 (`output_fix1.txt`) | After fix 2 (`output_fix2.txt`) |
|---|---|---|
| 1. 18 candidates built | 18 of 18, 149.3 s | 18 of 18, 137.6 s |
| 2. every variable bit-identical, or within 1e-14 | 518 of 518 bit-identical | 518 of 518 bit-identical |
| 3. a refusal in the all-columns march names latitude and level | (written with fix 2) | refused: 1,418,412 Pa outside the grid's 1e6 Pa, first met at 9.000000 deg, geopotential -150,000 m2/s2 |

Both records are under `reports/step05_2/`. The fix 1 run had checks 1 and 2 only; check 3 tests
code that fix 2 adds.

## 3. Timings

**One production transfer run** (10 N, 0.05 deg by 5e4 m2/s2, figures off, cProfile, Windows):

| | Before | After fix 1 | After fix 2 |
|---|---|---|---|
| The run | 94 s | 39 s | 12.7 s |
| `lib.mesh.build_columns` | 78 s | 25 s | 0.05 s |
| `lib.geoid.through_anchor` (reference surface) | 8 s | 4.4 s | 4.2 s |
| `lib.control.load_run_inputs` (reading the inputs) | 7 s | 8 s | 7.3 s |

**The regression, every suite, Windows, through the driver.** Before is proof G of SPEC_05 Step 0
(25 suites) and `step04_4`'s own acceptance run; after is this step's full run.

| Suite | Before (s) | After (s) | | Suite | Before (s) | After (s) |
|---|---|---|---|---|---|---|
| `step1/accept_step1` | 3 | 3 | | `step02_6` | 665 | 299 |
| `step1/verify_review_changes` | 2 | 2 | | `step03_1` | 12 | 9 |
| `step2` | 2 | 2 | | `step03_2` | 16 | 14 |
| `step3` | 3 | 4 | | `step03_3` | 181 | 149 |
| `step4` | 0 | 0 | | `step03_4` | 78 | 68 |
| `step5` | 69 | 43 | | `step04_0` | 320 | 286 |
| `step6` | 50 | 29 | | `step04_1` | 36 | 28 |
| `step7` | 147 | 98 | | `step04_2` | 1999 | 530 |
| `step8` | 3 | 3 | | `step04_3` | 4631 | 1532 |
| `step9` | 5 | 6 | | `step04_4` | 15645 | 887 |
| `step02_1` | 11 | 12 | | `step04_5` | 3812 | 917 |
| `step02_2` | 87 | 47 | | `step05_1` | 20 | 23 |
| `step02_3` | 196 | 100 | | `step05_2` | | 205 |
| `step02_4` | 334 | 179 | | **all** | **about 8 h** | **95 min** |
| `step02_5` | 323 | 202 | | | | |

28 of 28 at their reference counts; the registered kind N (`b79c60cb...`), the closure inputs and
product, and the transfer inputs and product restored exactly.

## 4. The regression

A change to `lib` reaches everything (SPEC_05 §2a): the full set, above.

## 5. Findings

1. **What is left of the time, measured.** In `step04_3` (1532 s), `step04_4` (887 s) and
   `step04_5` (917 s) the cost is now mostly in the suites' own loops: `step04_3` integrates columns
   one latitude at a time through `lm.column` for its line integrals, and `tests/step04_5/fields.py`
   and `step04_2` build the cylinder-extended wind one latitude at a time over 3,962 latitudes per
   pass. Deliverable 3 leaves those helpers as they are; moving them to `build_columns` with a
   `wind_all` is a test-side change for a later step if wanted.
2. **The kind N build takes 85 s** (`refrac`, through `lib.geoid` and `refrac.anchor`, not the
   column march). It is in every suite that rebuilds kind N (`step02_x`, `step04_0`, `step05_2`) and
   is the next target for decision R, as is reading the inputs (7 s) and the reference-surface march
   (4 s) in a transfer run.
3. **`build_columns` takes an optional `wind_all`.** Existing calls are unchanged; the transfer is
   the only caller that gives it.
4. **Bit-identity off this platform is not measured.** The array and single-value math routines
   agree on this machine; on Linux they may not, and the 1e-14 rule covers it. No Linux run is made
   (the author's direction, REPORT_05_step0 section 10).

## 6. What is committed on acceptance

The five source files of section 1, `tests/step05_2/accept_step05_2.py`, `tests/run_regression.sh`
and this report. No registered file.
