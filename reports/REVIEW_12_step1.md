# REVIEW 12, Step 1. The Voyager 2 egress profile, and the first comparison

Reviewing agent, 8 October 2026. Report: `reports/REPORT_12_step1.md`. Specification:
`docs/specs/SPEC_12_Lindal_Egress.md` v0.2.

**Accepted once the figure changes below are made and the author has viewed them. SPEC_12 closes
with this step.**

- **The digitizing.** The tool reproduces the CSV. The ingress check gives 1.15 K RMS in temperature
  and 1.74 percent in pressure. Below 50 km the trace runs 1 to 2.4 K cool.
  - There the lapse rate is near 0.8 K/km, so a vertical misplacement of 1 to 2 km is a 1 to 2 K
    error. The altitude ticks themselves depart from a straight line by up to 1.15 km.
  - The error belongs to the figure's axis, so the two curves share it.
- **The comparison.** The transfer does worse than the untransferred anchor everywhere except the
  deep troposphere:

  | Level | Anchor - egress | Delivered - egress |
  |---|---|---|
  | 110 mbar | +0.3 K | +3.9 K |
  | 290 mbar | +1.1 K | +2.4 K |
  | 730 mbar | +4.3 K | -0.9 K |

  - **From 20 to 200 mbar** the delivered profile is up to about 4 K too warm. It starts from the
    radio's northern value and adds the north-to-south change the IRIS gradients imply. At 110 mbar
    IRIS has the south about 5 K warmer than the north (81.4 to 86.6 K); the radio has no change
    (84.0 to 83.8 K).
  - **Above 10 mbar** it follows the anchor, while the egress is up to about 15 K colder near
    0.3 mbar. The wind has no data there, and the relaxed shear leaves no gradient to transfer.
  - **Below 300 mbar** it tracks the egress within about 1 K. This is where the transfer improves on
    the anchor.

**Rulings.**

1. **The IRIS temperatures in `comparison.png`.** `fit`'s value at 36.3 N and at 31.2 S
   planetographic, at 110, 290 and 730 mbar, drawn as markers, one style per latitude. The six
   values also go into acceptance 3's output. The IRIS data have a gap from -27.7 to -38.5 deg, so
   at 31.2 S the fit spans it; the dots on both sides agree.
2. **Labels.** Every latitude in the figures is given in one convention, planetographic as Lindal
   gives his, with planetocentric in parentheses.
3. **A difference panel** in `comparison.png`, against pressure, with the egress's sigma_T as a
   band:
   - delivered minus egress;
   - anchor minus egress;
   - the IRIS value at 31.2 S minus egress, as markers.
4. **The changes of §3 go into acceptance 3's output**, beside the absolute differences.
5. **No correction of the egress by its partner.** The CSV stays as traced, and its offset is
   stated by the ingress check.

Order of work: the figure and output changes; the author's view; then commit and push. No
regression.
- the tool;
- the CSV;
- `tests/step12_1/`;
- the run's control files;
- the report and this review;
- the SPEC_12 row in `STATE.md`.
