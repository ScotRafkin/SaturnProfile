# REVIEW 03, Step 1. `lib.geopotential`

Review of `reports/REPORT_03_step1.md` against SPEC_03 v0.6 Step 1. Reviewer: S. Rafkin,
14 September 2026. **Disposition: accepted.** Finding 1 corrects a defect in Step 4's wording,
restated at SPEC_03 v0.7 before that step is reached. Step 2 proceeds after the acceptance
commit.

## What was verified independently

From the swept product (`64c5d01a...`) with numpy and a gravity routine written for the
purpose, not the module under review: `|g_eff|` on the levels from the embedded G, R and W at
`phi_c` with `u` = 2.167 m/s, the trapezoid layer increments on the tabulated altitude, and
the gauge at level 29. Top `Φ` 2,852,354.97, bottom −1,043,742.65, first layer −126,655.82
m²/s²: the report's figures to the last digit. The one-factor radial construction gives
2,839,126.33 and −1,038,963.41, the report's check 4b. The code was read: the outward
cumulative sum from the gauge, the integer gauge index refusal, the monotonic height refusal
and the broadcast of `u` and `phi_c` are as the report describes.

## Rulings on the findings

1. **The Step 4 negative control.** The report is right and the specification was wrong at
   v0.6: the words described the two-factor construction (radial component times the projected
   radial increment) while the numbers, −5.1e-3 and −7.8e-3, belong to the one-factor
   construction (radial component times the tabulated altitude increment), which is what the
   reviewing agent measured before Step 0 and what the manuscript's present B1 and B2 amount
   to. The control is the one-factor construction, because that is the error the closure exists
   to catch. The two-factor reading is reported beside it and must fail the same bound; the
   reviewing agent measures its mean below 10 mbar at about −9.6e-3 and its bottom row at
   about −1.2e-2 on the swept product. Restated in SPEC_03 v0.7 Step 4. Nothing in Step 1
   changes.
2. **The trapezoid error estimate.** Noted; the bound of 1e-7 holds with a factor of three to
   spare and the estimate is a bound on the order, not a prediction. No action.

## On the decisions

1. The outward sum from the gauge: accepted, and the better construction; `Φ_a` exactly zero
   by construction is what the acceptance meant.
2. The refusals in `layer_increments` beyond the specification: accepted.
3. The integer gauge index, `np.float64(29.0)` and `True` refused, `np.int64` accepted:
   accepted.
4. The return order `(g_mag, g, G_phi, psi)` against `g_eff_vector`'s order: accepted; the
   specification named this order and the wrapper honors it.
5. The frozen `GeopotentialProfile`: accepted; Step 4 writes the product from it.
6. `u(phi_c)` by `refrac.anchor.wind_of_latitude` in the acceptance: accepted; how `forward`
   obtains it is Step 4's, and it should be the same function, not a second one.

## Order of work

1. Commit the author's documents (SPEC_03 v0.7, this review) in their own commit.
2. Commit the step (the acceptance commit for Step 1). Push.
3. Sweep: no input or product changed in this step; confirm with the `-dirty` scan that every
   product under `occul_data/` still carries `2149b64` and record that line in the report.
4. Set `STATE.md` to accepted with the commit. Step 2 then proceeds when the author says go.
