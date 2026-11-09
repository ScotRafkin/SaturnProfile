# REVIEW 20, Step 0. Sections, and the Voyager profiles as anchors

Reviewing agent, 9 November 2026. Report: `reports/REPORT_20_preexecution.md`. Specification:
`docs/specs/SPEC_20_Sections_and_Anchors.md` v0.2, ruled to v0.3.

**Accepted. Proceed with Step 1 as the report's §4 orders it.**

- **The anchors take SPEC_18's line** with one new source-specific module and no change to a shared
  tool or to the model. Lindal's files are read, never written, so acceptance 1 is met by
  `git status`, without a rebuild in a copy.
- **The radius finding is right.** `refrac` uses a profile's heights only above its 100 mbar level, and
  the anchor's absolute radius comes from kind D through the geoid march. The Fig. 9 point would only
  move the radius at which gravity is evaluated, a change under 0.01 K. Below materiality.
- **The section as a byproduct** of two marches per anchor, with the production redone at each
  latitude, agrees with separate runs to 1e-3 K, well inside acceptance 3. Reading the isobar label
  alone would not (0.2 to 0.3 K). Accepted.

**Rulings on §3:**
1. **Voyager 1 heights:** hydrostatic from 100 mbar, gravity at kind D's radius; Fig. 9 a check only.
2. **Composition:** Lindal's C, 94/6 with his ammonia, on each anchor's own levels. One composition
   for every input (SPEC_19). The ammonia's effect on these profiles is far below materiality.
3. **Season:** reuse Lindal's G and R, as SPEC_18 did. Each Voyager 1 file's comment states its own
   date, 13 November 1980; no new key and no solar longitude, since nothing reads it.
4. **sigma_T:** NaN in kind T, as SPEC_18; the figures' bands use Lindal's Fig. 4 curve.
5. **Sections:** commit each control file; ignore the product, which rebuilds in minutes and nothing
   reads.
6. **`tests/step04_0/sweep.py`:** leave it; the new anchors rebuild from their own build files.

**Order of work:** Step 1 as §4 of the report. No regression. Report, then commit after review.
