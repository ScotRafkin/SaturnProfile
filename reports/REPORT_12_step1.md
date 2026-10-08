# REPORT 12, Step 1. The Voyager 2 egress profile, and the first comparison

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_12_Lindal_Egress.md` v0.2. Working
tree on `main` at `1139e8d`, uncommitted. Acceptance: `tests/step12_1/accept_step12_1.py`, **4 of 4**
(three reported items, and a reproduction check beyond the specification); output in
`reports/step12_1/output.txt`, figures `overlay.png`, `ingress_check.png` and `comparison.png`. No
regression, as §0 directs.

**In short.**
- **The egress is digitized.** Its CSV is reproduced byte for byte by the tool.
- **The ingress check.** Traced and converted the same way, the ingress matches Table I to 1.15 K
  RMS in temperature and 1.74 percent RMS in pressure. Below 50 km it runs 1 to 2.4 K cool, where
  Fig. 4 and Table I themselves differ.
- **The comparison.** The run converges in 5 passes, with its pressure identity within 0.074 K. The
  delivered temperature is above the egress by 3.9 K at 110 mbar and 2.4 K at 290 mbar, and below
  it by 0.9 K at 730 mbar. Measured as changes, to cancel the figure's own offset:
  - Fig. 4's egress is no warmer than its ingress at 110 mbar (+0.03 K);
  - the transfer warms by 3.7 K, the change the IRIS gradients ask for;
  - so the radio occultations and the IRIS retrievals disagree on the temperature change between
    these two latitudes (§4).

## 1. The deliverables

- **`src/casspian/tools/lindal/egress_profile.py`.** It runs as `python -m
  casspian.tools.lindal.egress_profile`, with defaults for the PDF and the output. The module
  docstring states every step and every coded decision. In outline:
  - **Rendering.** Page 4 at 600 dpi by Ghostscript.
  - **Frames.** Found from the fully inked rows and columns.
  - **Calibration.**
    - Temperature is a least-squares line through its 7 ticks, residuals at most 0.06 K.
    - σT is a line through its 10 K tick and the panel's frame lines.
    - Altitude is piecewise linear between its 9 ticks. The ticks depart from one straight line by
      up to 1.15 km (3.5 px), far more than a tick's position is uncertain, so the scan's scale
      varies along the axis.
  - **Masks.** The labels, the 36.3 N leader and its arrowhead, and the 31.2 S leader are blanked
    by coded boxes, a band and a seed pixel.
  - **The tracing**, in four altitude zones, each scanned across the curves: rows above 172 km,
    columns from 172 to 140, rows from 140 to 70, columns below 70.
    - **Zone ends.** The 172 km boundary sits just below where egress and ingress cross, so the two
      curves keep their order within each zone.
    - **Which track is dashed** is read from the figure per zone. As a check, the dashed track should
      be found in more separate pieces; it is in all four zones (pieces 14 against 2, 4 against 1, 5
      against 3, 8 against 3).
    - **Runs.**
      - Runs 3 px or less apart are one line.
      - A run 1.3 line widths wide or more, with both tracks predicted within it, is the two lines
        side by side.
      - A track unseen for 30 lines is held at its last position, and a lone run then goes to the
        other track.
      - Each scan starts 5 km outside its zone.
  - **Binning.** Each curve's samples are the median of each 1 km bin. Empty bins, and spikes more
    than 2 K from the median of their six neighbors, are filled from those neighbors: 4 of 402 on
    the solid curve, 49 of 391 on the dashed and 9 of 397 on σT.
  - **To pressure** (§1 deliverable 2, §2 item 3), with
    `d ln p / dh = -m_bar |g_eff| / (R T)` integrated from 1 bar at 0 km:
    - the profile's own temperature;
    - Lindal's composition, giving 2.1351 to 2.1363 amu;
    - `lib.gravity`'s effective gravity under Lindal's gravity, rotation and wind;
    - the radius `r = r_1bar + h cos(psi)`, the reduction's own placement of an altitude, on
      Lindal's 1 bar ellipse from `lindal_geodesy.nc`. That puts 1 bar at 58,978.7 km at the egress
      latitude and 58,535.8 km at the ingress latitude.

    Each profile's planetocentric latitude is its label converted as the anchor's is: -26.06095 deg
    for the egress (`tests/step11_2/target_latitude.py`) and 30.80557 deg for the ingress, the
    anchor's own.
- **`occul_data/lindal/voyager2_egress.csv`:** 391 levels every 1 km from -16 to 374 km, 1337.7 to
  0.2247 mbar. σT is empty at -16 km, below the σT curve's lowest point.
- **The run** `lindal_iris_ii_110` is rebuilt and rerun at -26.06095353762054 deg. Its target was
  set as §1 deliverable 4 states and committed at `33d8340`. The outer loop takes 5 passes to
  8.3e-9, the mesh is 1142 by 82, and the pressure identity is at most 0.074 K.

## 2. Acceptance

0. **Reproduced** (beyond the specification): the tool writes the committed CSV byte for byte.
1. **The ingress check** (`ingress_check.png`). The solid curve is converted the same way at 36.3 N
   (30.81 N planetocentric) and compared with Table I at its 66 altitudes, from -14.1 to 376.7 km.
   - **Temperature:** RMS 1.15 K, largest -2.39 K at -14.1 km. It is within about 0.5 K from 50 to
     150 km, and 1 to 2.4 K cool below 50 km. Above 150 km it varies between -0.5 and +1.7 K.
   - **Pressure:** RMS 1.74 percent, largest +5.0 percent at 363.9 km (0.25 mbar). It is within
     1 percent below 150 km and grows above, since the integral carries every temperature offset
     below it.
   - **The mean molar mass** at 1 bar is 2.1352 amu, against Lindal's 2.135.
   - **The cool troposphere is Fig. 4 itself.** Drawn over the scan, Table I's points sit about half
     a line above the solid curve between 15 and 45 km, while the trace follows the line's center.
     So the figure and the table differ there by up to about 2 K, and that is the digitizing and
     conversion error the egress carries (SPEC_12 §1 acceptance 1).
2. **The overlay** (`overlay.png`): both curves and σT traced over the scan, for the author.
3. **The comparison** (`comparison.png`), in three panels against pressure (REVIEW_12_step1
   rulings 1 to 4):
   - temperature: the delivered profile from the run, the egress from the CSV with its σT band, the
     anchor, and the IRIS fit (`iris_temperatures.fit` at its defaults) at 36.3 N (30.81 N
     planetocentric) as open circles and at 31.2 S (26.06 S planetocentric) as filled squares;
   - the difference from the egress, with its σT band: delivered minus egress, anchor minus
     egress, and the IRIS fit at 31.2 S minus egress as markers;
   - mass density, the egress's at 2.135 amu.

   | Level | Delivered | Egress (σT) | Delivered - egress | Anchor - egress | IRIS 36.3 N | IRIS 31.2 S | IRIS 31.2 S - egress | Density, delivered / egress - 1 |
   |---|---|---|---|---|---|---|---|---|
   | 110 mbar | 87.69 K | 83.77 K (2.27) | **+3.92 K** | +0.25 K | 81.42 K | 86.63 K | +2.85 K | -4.3 % |
   | 290 mbar | 95.54 K | 93.17 K (2.29) | **+2.37 K** | +1.06 K | 96.00 K | 96.59 K | +3.42 K | -2.4 % |
   | 730 mbar | 116.36 K | 117.30 K (2.73) | **-0.94 K** | +4.32 K | 119.30 K | 114.92 K | -2.38 K | +0.8 % |

   The output also gives §3's changes. Every latitude in the figures and the output is
   planetographic, with planetocentric in parentheses (ruling 2).

## 3. The comparison as changes, to cancel the figure's offset

The egress carries the figure's own offset against Table I, which is up to about 1 K at these
pressures. Fig. 4's two curves share it, so the change from ingress to egress read off the figure
is free of it. That change is set beside the transfer's change from the anchor to the target. The
transfer's change is from the run's product and the anchor's embedded thermo; the figure's from the
two traces, converted the same way.

| Level | Fig. 4: egress - ingress | Transfer: delivered - anchor | Difference | IRIS fit: 31.2 S - 36.3 N | Traced ingress - Table I |
|---|---|---|---|---|---|
| 110 mbar | +0.03 K | +3.67 K | +3.64 K | +5.21 K | -0.28 K |
| 290 mbar | -0.48 K | +1.31 K | +1.79 K | +0.59 K | -0.57 K |
| 730 mbar | -3.46 K | -5.26 K | -1.80 K | -4.38 K | -0.86 K |

These are in acceptance 3's output (REVIEW_12_step1 ruling 4).

## 4. For the author

- **The radio and IRIS disagree on the change between these latitudes.** At 110 mbar Lindal's two
  radio profiles, 36.3 N (30.81 N planetocentric) and 31.2 S (26.06 S planetocentric), are the same temperature to 0.03 K. The IRIS
  fit, from which the wind's shear is built, is about 5 K warmer at the southern latitude.
  - The IRIS fit's change from 36.3 N to 31.2 S is +5.21 K at 110 mbar. SPEC_11 Step 2 measured
    +5.27 K to 31.49 S, its earlier target.
  - The transfer carries that IRIS change faithfully: +3.7 K here. The rest is mostly the
    slope-and-value difference of REVIEW_11_step2's part A.
  - At 730 mbar both instruments see a cooling to the south: the radio -3.5 K, the transfer -5.3 K.
  - Whether the difference at 110 and 290 mbar is the IRIS retrieval, the radio occultation, or the
    wind's construction is the question this comparison raises. It is not settled here.
- **The temperature offset is within the egress's own σT** (2.3 to 2.7 K) at all three levels: the
  run sits 1.7, 1.0 and 0.3 σT from the egress. The change-based difference of 3.6 K at 110 mbar is
  larger than the figure's digitizing error (about 1 K there) and comparable to σT.
- **The ingress check's cool troposphere.** Fig. 4's ingress runs 1 to 2.4 K below Table I between
  the surface and 50 km. Per REVIEW_12_step1 ruling 5 the egress is not corrected by its partner;
  the CSV stays as traced, and the ingress check states its offset.

## 5. Notes

- **Not committed:** the run's inputs and product, as for every unregistered run. The rerun used the
  relaxation of the dirty refusal for its own process, because the tool is uncommitted.
- **`docs/Lindal_et_al_1985_AJ90_1136.pdf`** is in the repository already; the tool reads it.
