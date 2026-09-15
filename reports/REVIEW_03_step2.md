# REVIEW 03, Step 2. `lib.hydrostatic`

Review of `reports/REPORT_03_step2.md` against SPEC_03 v0.7 Step 2. Reviewer: S. Rafkin,
15 September 2026. **Disposition: accepted.** Finding 1 is a wording defect in the
specification, restated at SPEC_03 v0.8. Step 3 proceeds after the acceptance commit.

## What was verified independently

The reviewing agent's own linear-T column, written before the specification, was 140 K at the
bottom level and 80 K at the top, linear in Φ on the Table I geopotential grid; it gave 5.99e-4
for the exact-exponential layer integral and 8.78e-3 for the trapezoid, and 1e-15 for the
isothermal column, the report's figures. So the report read the specification the way it was
meant and the words were at fault. The module was read: the `expm1(x) / x` form is the
specification's expression rearranged (`ρ_{k+1} − ρ_k = ρ_k expm1(x)`), the limit branch is at
1e-10 as specified, the direction check refuses anything but a top-down column, and
`temperature` is Eq. B6 with `k_B` from `lib.constants`. Check 5, the tabulated temperature
recovered to 4.4e-16 from the reduction's own N and ℛ̄, is the SPEC_02 Step 3 closure through
the new module, as the specification asked.

## Rulings on the findings

1. **The orientation of the linear-T column.** Correct, and the report's proposed wording is
   adopted: falling with height, 140 K at the bottom level and 80 K at the top, linear in Φ
   (SPEC_03 v0.8 Step 2). Reporting the other orientation beside it was the right call.
2. **1.1e-15 against "1e-15".** No action; the acceptance bound is 1e-14.

## On the decisions

1. The `expm1` form: accepted, and preferred. It evaluates the same expression without the
   cancellation the literal difference suffers at small log ratios (−8.7e-8 at 1e-9, measured),
   and it makes the 1e-10 limit branch a formality rather than a necessity. The specification's
   expression stands as the definition; the module's form is the evaluation.
2. The refusals beyond the specification (non-finite density or geopotential, a column that is
   not top down, a non-positive boundary pressure or layer mass): accepted. A bottom-up column
   refused rather than integrated with negative masses is the right behavior for a function
   that only integrates downward.
3. `density` taking kg/mol and converting inside: accepted; it matches what kinds N and C
   carry.
4. The closed-form pressure of the linear-T column, `p = p_b (T / T_top)^(−m / (k_B β))`:
   accepted; it is the exact solution and the right reference.
5. The test columns on the Step 1 grid with the product's top-level composition held uniform,
   halving by Φ midpoints with the analytic density: accepted.

## A note on the record

The report says the corrected REVIEW_03_step1 was not in the working tree when Step 2 began.
The first write of the correction did not reach the disk; it was written again on
15 September and read back. The corrected file is in the working tree now and is committed
with this review.

## Order of work

1. Commit the author's documents (SPEC_03 v0.8, the corrected REVIEW_03_step1, this review,
   `STATE.md`) in their own commit.
2. Commit the step (the acceptance commit for Step 2). Push.
3. Sweep: no input or product changed; confirm with the `-dirty` scan that every product
   under `occul_data/` still carries the commit it was swept at (`2149b64`; `ece58d2` for
   gravity and rotation) and record that line in the report.
4. Set `STATE.md` to accepted with the commit. Step 3 then proceeds when the author says go.
