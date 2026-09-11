# REVIEW 01, Step 9. Stage two of the Lindal tool, kinds T and D and the manifest

Review of `reports/REPORT_01_step9.md` against SPEC_01 v0.17 Step 9. Reviewer: S. Rafkin,
11 September 2026. **Disposition: accepted.** SPEC_00 goes to v0.9 and SPEC_01 to v0.18 for
the rulings below; the author has edited `raw/lindal_scalars.toml`. Commit, then run the
directory sweep, then SPEC_01 is complete.

## Rulings on the findings

1. **Fit residual.** Correct: NaN per surface, `fit_residual_range_m = [3000, 4000]` with the
   quotation. SPEC_00 §6.3 now says so, so the next reader expecting a number is told why
   there is not one.

2. **Gravity citation.** Agreed, and done at the source: `lindal_scalars.toml` gains a
   `[gravity_used_by_source]` table (citations, the Appendix sentence as `value_source`, and
   the pole vector references) in the shape of `[wind_used_by_source]`, whose `how_combined`
   line is also brought up to the v0.14 reading. Kind T's `source_gravity_citation` reads the
   new key. Because the raw bundle hashes `lindal_scalars.toml`, this edit changes the bundle
   and therefore every product downstream of it; that is what the sweep in item 4 is for, and
   it is the right time to do it, since nothing is committed yet.

3. **Longitude.** Fine as built: the scalar is the ingress start, the swath sits beside it,
   the same treatment the latitude has.

4. **The rebuild rule.** Agreed, and it is a better rule: after every acceptance commit, every
   product under `occul_data/` carrying `-dirty` is rebuilt, whichever step wrote it. SPEC_00 §8
   and SPEC_01 §0 now say so. The check that found the two Step 3 products dirty is the
   enforcement and stays in the acceptance script of every later step.

5. **Manifest tolerance in degrees.** Accepted as written. The manifest is a file a person
   edits, `lib` is radians without exception, and `refrac` converts at its boundary; SPEC_02
   will state that conversion explicitly so the inconsistency stays deliberate.

## One addition the report could not have known to make

The manifest's `[geoid]` section predates SPEC_01 v0.16 and does not name the anchor rule.
The reduction has to freeze `r0` the same way Step 7 built the surface, so SPEC_00 v0.9 adds
`anchor_rule = "mean_polar_radius"` to the manifest vocabulary. Add the line to
`lindal_reduction.toml` (the tool writes it) before the rebuild.

## Actions

- Commit Step 9. Run the sweep: rebuild the raw bundle (Step 2 tool, for the scalars edit),
  then every product under `occul_data/lindal/` in dependency order (G, R, W, C, T, D, the
  manifest), so that all carry the clean commit; rerun the Step 9 acceptance to confirm none
  is dirty and every hash in the manifest matches.
- Report the sweep as a short addendum to REPORT_01_step9 with the final commit hash and the
  list of products with their hashes. That addendum closes SPEC_01.
- SPEC_02 (`refrac`) is next and is not to be started before it exists.
