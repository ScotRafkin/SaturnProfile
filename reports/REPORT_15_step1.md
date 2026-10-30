# REPORT 15, Step 1. Before acceptance

Coding agent, 30 October 2026. Specification: `docs/specs/SPEC_15_Radio_Targets.md` v0.2. Working
tree on `main` at `3472bf4`, uncommitted. No regression (SPEC_15 §0).

Acceptance: `tests/step15_1/accept_step15_1.py` **5 of 5**; output in
`reports/step15_1/output.txt`, figures and tables beside it.

**In short.**
- **The tracing.** Lindal's Figs. 5 and 6 are traced by `tools/lindal/radio_profiles.py`.
  Their curves loop, curl, cross and touch, so no one scan direction can separate them. Each curve
  is therefore traced along a guide of control points placed by hand on the scan and coded in the
  module, and every point is snapped across the guide to the centre of the ink it lies on (§1).
  Calibration: the ticks lie within 0.0033 and 0.0085 dex in pressure and 0.09 K in temperature of
  their fits.
- **The calibration traces.** Fig. 5's 36.3 N against Lindal's Table I: RMS 0.95 K. Its 31.2 S
  against the Fig. 4 egress of SPEC_12: RMS 1.51 K. Both are 1.5 to 3 K warm from 600 to 1200 mbar.
  That is the figure, not the trace: Table I's points fall a steady 6 to 7 px below the centre of
  the drawn solid line there (§2).
- **Pioneer 11.** Against Kliore's Table 1 as published, Fig. 6's stippled trace is 3.31 K colder on
  average. Against the table rescaled to 94/6 by Lindal's rule it is +0.16 K. Lindal rescaled it,
  and the rescaling of SPEC_15 §1 reproduces his (§3).
- **The targets** in planetocentric latitude: 74.2 S is -71.0680, 3.0 S is -2.3736 and 11.9 S is
  -9.5054. The last is Kliore's 9.5 S, which confirms §2 ruling 1 independently. The conversion
  reproduces the anchor's 30.8056 and SPEC_12's -26.0610 exactly.
- **The comparison** at four targets is §5's table. Neither hypothesis fits the three new
  profiles. At 3.0 S and 11.9 S, above 30 mbar, the observed profiles carry structure of 10 to 30 K
  that the transfer does not make. The untransferred anchor misses them by the same amounts, so
  between 36.3 N and these latitudes the transfer changes little above 30 mbar.

## 1. Tracing (deliverable 1)

**Calibration.** Page 5 is rendered at 600 dpi. The frame and the ticks give a least-squares line of
log10 p in row through the ticks at 1 to 1000 mbar and the frame lines at 0.1 and 10^4 mbar, and one
of T in column through the ticks every 10 K and the frame lines at 80 and 160 K.
- Fig. 5: ticks within 0.0033 dex and 0.087 K of their fits.
- Fig. 6: ticks within 0.0085 dex and 0.058 K.

**Tracing along guides.** A first tracer followed each curve on its own, steering by the ink ahead.
It wandered onto label leaders, circled near the solid's nose, and cut across the dashed curve's
curls. It was replaced by a deterministic one:
- each curve has a guide, a polyline of control points placed by hand on the scan, close together in
  bends and sparse where the curve is smooth;
- most points were placed from row and column scans of the ink, so they sit on the line to a pixel;
- the guide is resampled every 3 px, and each point is snapped across the guide's direction to the
  centre of the ink run it lies on, or of one within 4 px;
- where there is none (a dash's gap) or the run is wider than 1.8 line widths (a crossing, a leader,
  a dot), the guide's point stands.

The points are binned in 0.01 dex of pressure (median), and gaps are filled linearly in log p. Every
hand decision is the guide's points in `CURVES`, and the module's comments say where each curve is
ambiguous. Overlays: `overlay_fig5.png` and `overlay_fig6.png`.

**Where the scan cannot separate two curves**, both traces follow the shared ink:
- 36.3 N and 31.2 S: 4.4 to 9.3 mbar and 72 to 447 mbar;
- 74.2 S and 36.3 N: 0.43 to 0.49, 0.60 to 0.65, 1.86 to 2.29 and 288 to 339 mbar;
- 74.2 S and 31.2 S: 288 to 347 and 537 to 617 mbar;
- 3.0 S and 11.9 S: 1.20 to 1.26 mbar.

On the lower branch of Fig. 5, the dashed curve lies half merged under the solid from 200 to 330
mbar, then joins the dash-dot's band until about 600 mbar. The guide follows that reading of the
scan (comments at `_FIG5_DASHED`).

**Error bars,** half the bar's length from cap to cap, at the row where it is drawn:

| Curve | Top | Bottom |
|---|---|---|
| 74.2 S | 0.141 mbar, ±9.81 K | 1040 mbar, ±6.04 K |
| 3.0 S | 0.402 mbar, ±9.86 K | 1168 mbar, ±4.41 K |
| 36.3 N | 0.189 mbar, ±9.89 K | 1223 mbar, ±4.29 K |
| 31.2 S | 0.211 mbar, ±9.94 K | 1375 mbar, ±4.04 K |

The 74.2 S top bar's centre, 133.1 K, is 0.3 K from its curve's end dot; the other centres were
not checked against their dots. Only the 74.2 S and 3.0 S bars are written to files. The 74.2 S trace's last bin is 1000 mbar, so its bottom bar, drawn at 1040 mbar,
is written on that row, 0.017 dex from where it is drawn.

**Known small faults.**
- The 3.0 S trace jogs about 2 px where its label's arrowhead touches the line near 300 mbar.
- One stippled dash near 1.4 mbar is passed at its end rather than through its middle.

Both are within a line width.

## 2. The calibration traces (check 1)

| Trace | Reference | Levels | RMS | Largest | Levels over 1 K |
|---|---|---|---|---|---|
| Fig. 5, 36.3 N | Lindal Table I | 64 | 0.95 K | -2.13 K at 0.4 mbar | 18 |
| Fig. 5, 31.2 S | `voyager2_egress.csv` (Fig. 4) | 386 | 1.51 K | -4.07 K at 9.26 mbar | 190 |

**Over 1 K, explained.**
- **600 to 1200 mbar, both traces, +1.5 to +3 K.** For each Table I level, its predicted place on
  Fig. 5 and the solid line's centre in that column were measured. Below 600 mbar the drawn line is
  a steady 6 to 7 px (0.018 dex, 4 percent in pressure) above Table I's points. On a line this
  shallow that is about 1.8 K at fixed pressure. Above 500 mbar the offset is within ±3 px. The
  trace is on the line's centre, so the drawn curve departs from Table I at depth. The 31.2 S curve
  shows the same against Fig. 4. Any Fig. 5 profile likely carries the same bias at depth,
  74.2 S included.
- **0.4 mbar on 36.3 N, -2.1 K**, at the solid's nose, where it and the dash-dot touch.
- **4 to 9 mbar on 31.2 S, up to -4.1 K**, where the dashed curve coincides with the solid (§1), so
  the trace carries the solid's temperature there.
- **200 mbar on 36.3 N, -1.05 K**, inside the coincidence of 72 to 447 mbar.

## 3. Pioneer 11 (deliverable 2, check 2)

`occul_data/kliore1980/pioneer11_exit.csv` is the table as published: 128 rows, the uncertainty
ranges on every fifth row, 26 of them.
- **Transcription.** Taken from the PDF's text layer column by column, two OCR slips repaired
  (`88. l`, `3I. l`), and each row checked against the scanned page: radius falling, refractivity
  and pressure rising.
- **The two misprints**, corrected and marked on their rows: 63.7 K at 142.5 mbar is 93.7, and the
  pressure range "-02." at 1.8 mbar is -0.2.
- **Kliore's caption** gives the latitude as 9.5 S (SPEC_15 §0 says about 9.7).

**Rescaling.** `rescale_composition` reads the molar masses and Lindal's refractivities (H2 136, He
35) from `species_master.toml`. From 90/10 to 94/6 by number: temperature ×0.9641 (-3.59 percent),
pressure ×0.9341, mean refractivity ×1.0321.

| Fig. 6 stippled trace minus | Mean | RMS | Largest |
|---|---|---|---|
| Kliore as published | -3.31 K | 3.56 K | -4.76 K at 24.8 mbar |
| Kliore rescaled | +0.16 K | 1.74 K | +7.76 K at 1.68 mbar |

So Lindal rescaled the table, by this rule or one equivalent to it. The largest rescaled difference
is at the table's top row, where Kliore's own range is -11.5 to +21.1 K.

## 4. Data files and runs (deliverables 3 to 6)

**Data files.**
- `occul_data/lindal/voyager1_ingress.csv` (74.2 S): 385 levels, 0.145 to 1000 mbar.
- `occul_data/lindal/voyager1_egress.csv` (3.0 S): 345 levels, 0.417 to 1148 mbar.
- `occul_data/kliore1980/pioneer11_exit.csv`.

Each file's top comment gives the source, the figure or table, the latitude as published and the
epoch. `radio_profiles.read_csv` reads past the comment lines, which `numpy.genfromtxt` alone would
take as the header.

**Targets.** `tests/step15_1/target_latitudes.py` applies the reduction's fixed point, as SPEC_12
did: 74.2 S becomes -71.06804, 3.0 S becomes -2.37362, and 11.9 S becomes -9.50537.

**Runs.**
- Made with `casspian-new-run --from`.
- **The namelists** differ from their parents in `[target]` only.
- **The build files** differ by a comment giving the observation, the epoch, and that the epoch
  difference is not modeled.

Pressure identity, largest over the run's levels:

| Hypothesis | 74.2 S | 3.0 S | 11.9 S |
|---|---|---|---|
| 1 | 0.057 K | 0.053 K | 0.069 K |
| 2 | 0.180 K | 0.147 K | 0.254 K, at 3.99 mbar |

The 31.2 S runs are SPEC_14's.

## 5. For the author

Delivered minus observed (K):

| Target | Profile | 730 | 290 | 110 | 10 | 4.5 | 1 | 0.3 mbar |
|---|---|---|---|---|---|---|---|---|
| 31.2 S (26.06 S) | hypothesis 1 | +0.40 | +2.32 | +5.00 | -1.90 | +0.10 | +5.24 | +15.10 |
| | hypothesis 2 | +0.40 | +2.30 | +5.09 | +0.20 | -0.77 | -2.17 | +2.88 |
| | anchor, untransferred | +4.31 | +1.06 | +0.24 | -3.68 | -0.91 | +4.88 | +14.85 |
| 74.2 S (71.07 S) | hypothesis 1 | +2.52 | +0.12 | +0.40 | -8.10 | -6.30 | +4.21 | +10.76 |
| | hypothesis 2 | +2.52 | +0.11 | +0.50 | -5.88 | -7.26 | -3.70 | -2.29 |
| | anchor, untransferred | +5.56 | -0.23 | -4.78 | -10.12 | -7.54 | +3.66 | +10.33 |
| 3.0 S (2.37 S) | hypothesis 1 | +0.08 | +2.21 | +0.76 | -10.30 | +6.84 | +16.48 | NaN |
| | hypothesis 2 | +0.08 | +2.20 | +0.84 | -8.27 | +6.32 | +10.04 | NaN |
| | anchor, untransferred | +5.29 | +3.26 | -2.27 | -10.97 | +6.99 | +17.29 | NaN |
| 11.9 S (9.51 S) | hypothesis 1 | NaN | NaN | +2.38 | +16.86 | +27.36 | NaN | NaN |
| | hypothesis 2 | NaN | NaN | +2.52 | +20.06 | +26.09 | NaN | NaN |
| | anchor, untransferred | NaN | NaN | -4.15 | +15.10 | +26.80 | NaN | NaN |

IRIS fit minus observed at 110, 290 and 730 mbar:
- 74.2 S: -2.06, +1.10, +0.96 K;
- 3.0 S: -1.72, +3.57, -1.30 K;
- 11.9 S: +0.02 K at 110 mbar (the table stops at 182 mbar).

**What the table shows:**
- **Below 100 mbar** the transfer helps at every target: it removes the anchor's 4 to 6 K warm bias
  at 730 mbar and brings 110 mbar within about 2.5 K of the observation, except at 31.2 S, where
  the radio and IRIS disagree (SPEC_12).
- **At 74.2 S** hypothesis 2 improves 1 and 0.3 mbar, as at 31.2 S. Between 3 and 15 mbar both
  hypotheses are 6 to 8 K cold, and the profile there is the warmest of Fig. 5.
- **At 3.0 S and 11.9 S** the observed profiles have a warm layer of 10 to 30 K above 30 mbar,
  near 10 mbar at 3.0 S and 3 to 10 mbar at 11.9 S. Lindal attributes it to vertically
  propagating waves. Neither hypothesis makes it, and the anchor untransferred misses it by the
  same amounts. Between 30.8 N and these latitudes the transfer changes little above 30 mbar.
- **The epochs differ** for 74.2 S, 3.0 S and 11.9 S (Voyager 1, November 1980; Pioneer 11,
  September 1979), and the runs stay at the anchor's season. Their differences include an unmodeled
  seasonal term.

Figures: `profiles_312S.png`, `profiles_742S.png`, `profiles_30S.png`, `profiles_119S.png` (temperature,
minus observed, density) and `summary.png` (delivered minus observed, four targets by two
hypotheses).

## 6. Files

- `src/casspian/tools/lindal/radio_profiles.py`: the tracer and its guides, the bars, the
  composition rescaling, the CSV reader. Running it rewrites the two Voyager 1 files.
- `occul_data/lindal/voyager1_ingress.csv`, `voyager1_egress.csv`,
  `occul_data/kliore1980/pioneer11_exit.csv`.
- `forward/lindal_iris_v_sl_ii_110_{74s,03s,12s}/` and `forward/lindal_iris_v_sl_iii_110_{74s,03s,12s}/`:
  the control files only.
- `tests/step15_1/accept_step15_1.py` and `tests/step15_1/target_latitudes.py`.
- `docs/Kliore_etal_1980_JGR85_5857.pdf` and `docs/Kliore_etal_1980_Science207_446.pdf` (deliverable 4),
  with SPEC_15, in the acceptance commit.
