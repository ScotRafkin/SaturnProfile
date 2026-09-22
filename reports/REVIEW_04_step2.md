# REVIEW 04, Step 2. The working mesh, the columns, and the cylinder-extended wind

Review of `reports/REPORT_04_step2.md` against SPEC_04 v0.10 Step 2 and decision P. Reviewer:
S. Rafkin, 22 September 2026. **Disposition: accepted with changes**, the regression's table to be
filled in before the acceptance commit. All three failing checks fail on the specification's
numbers, not on the code: two on a check whose closed form was wrong, one on the sampling of the
construction's inversion curve, which the report itself diagnosed. The rulings are in SPEC_04 §13
(v0.11).

## What was verified

The column at `φ_c` on the anchor's own levels reproduces the reviewing agent's `z_lv` (286,715.2
and −104,098.0 m) and `r − r0` (288,051.3 and −104,576.8 m) to 0.05 m, and the ratio of the two
at the top, 1.004660, is the `1 / cos ψ` projection at the anchor's field-line tilt. The
cylinder-extended wind on the anchor's column, 9.4998, 3.7190, 2.1671 and 1.9254 m/s, lies within
1e-2 of the reviewing agent's 9.490, 3.717, 2.168 and 1.926, and the geopotential shift, −322.6
and +18.0 m²/s², within 0.3 of the reviewing agent's −322.3 and +18.0. The fixed point converged
in four passes to 4e-7 m/s; the reference-level identity holds exactly at all 361 latitudes; the
written file passes the kind W schema, the sum identity and the poles.

## Rulings

1. **Finding 1, the uniform gravity check.** The specification's fault: `J = 0, Ω = 0, u = 0`
   leaves `GM / r²`, and the closed form the check should carry is `r = r0 / (1 − r0 Φ / GM)`,
   which the column meets at 8.9e-16. Restated at v0.11.
2. **Finding 2, the cancellation.** Correct; the bound becomes 1e-7 m absolute (seven ulp of a
   58,516 km radius). Restated.
3. **Finding 3, the inversion curve's sampling.** Accepted and generalized: `s_ref` is sampled
   at the mesh's latitude spacing over the file's whole range, not on the 0.5° grid with a few
   latitudes inserted, so no mapped latitude is more than one mesh spacing from a sample. The
   reviewing agent's values were measured on a 0.05° sampling; the chord error of a 0.5° cell
   is what the report measured as a uniform offset of up to 0.01 m/s.
4. **Finding 4, per hemisphere.** Correct, and decision P now says so: a cylinder cuts the
   reference surface once per hemisphere, the wind is not symmetric, and a symmetric deep wind
   is a different hypothesis for a tool to make later. Nothing in the run's span crosses the
   equator.
5. **Finding 5, `edge_order = 2`.** Correct; the rule is stated.
6. **Finding 6, the 10° N bottom level at 0.054 against 0.05.** Not loosened. The departure has
   the size and sign of finding 3's chord error at the latitude the bottom level maps to; the
   check is rerun under ruling 3 against the same values and bound. If it still misses, the
   report says so with the values and the reviewing agent remeasures under the same sampling
   before any bound moves.

Decisions 1 to 11 are accepted as reported; decision 8 is now in the specification's text, and
decisions 9 and 10 are superseded by rulings 3 and 4 as generalized.

## Order of work

1. Apply rulings 1 to 5 in the acceptance script (checks 1a and 1b restated; `s_ref` sampled at
   the mesh spacing over the file's range); rerun; refresh the report's head, section 3 and
   finding 6 with the remeasured 10° N values.
2. Fill the regression table when the run finishes; every accepted suite at its reference count.
3. Commit the author's documents (SPEC_04 v0.11, this review, `STATE.md`) in their own commit.
4. The acceptance commit for Step 2. Push. No product changes, so no sweep; `STATE.md` to
   accepted with the commit. Step 3 proceeds.
