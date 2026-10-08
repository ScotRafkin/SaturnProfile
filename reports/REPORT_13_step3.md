# REPORT 13, Step 3. Before acceptance

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_13_SL_Wind.md` v0.5, §2a. Working tree
on `main` at `94a3090`, uncommitted, with Steps 1 and 2's changes. No regression.

Acceptance, rerun: `tests/step13_1/accept_step13_1.py` **3 of 3** and
`tests/step13_2/accept_step13_2.py` **2 of 2**; output in `reports/step13_1/output.txt` and
`reports/step13_2/output.txt`.

**In short.**
- **Item 1, the smooth window.** `fit(window = "smooth")` takes SPEC_10's adaptive width, makes it
  smooth in latitude and never narrower. The derivative's spikes are gone (figure
  `shear_slope_value.png`). The fit's RMS about the points is unchanged to 0.001 K.
- **What the smooth window changes.** Inside the IRIS gap (-27.7 to -38.5) the wider window
  changes the fitted temperatures themselves. At 31.2 S the IRIS fit is 86.26 K at 110 mbar
  instead of 86.63, and 116.21 K at 730 mbar instead of 114.92. The IRIS change from 36.3 N to
  31.2 S becomes +4.84, +0.51 and -3.09 K, against +5.21, +0.59 and -4.38 with the adaptive window.
- **Part (A) is now its linearization alone.** (A1), the plain integral against the change in
  value, is +0.003, +0.013 and -0.029 K at the three levels.
- **Item 2.** `profiles.png` shows the IRIS fit's values at 36.3 N and 31.2 S at the three levels,
  for each window, and the 31.2 S value minus the egress.
- **Item 3** is as reported below. SPEC_00 v0.24 makes its three schema changes official.

## Item 1. A smooth window

### The construction

`iris_temperatures.fit` gains `window`: `"adaptive"` (SPEC_10, the default) or `"smooth"`. The
smooth width (`smooth_width`) is built in three steps:
1. SPEC_10's adaptive width `w_a` is taken on a fine grid of 0.05 deg over the data's span
   (`SMOOTH_STEP_DEG`). The adaptive computation is moved unchanged into `_adaptive_width`.
2. Its running maximum is taken over ±4 deg (`SMOOTH_HALFWIDTH_DEG`).
3. That is averaged with a Hann kernel of the same half-width, `cos^2(pi d / 2H)`, renormalized
   where the data's span cuts it.

**Never narrower.** Every latitude within `H` of `phi0` has a running maximum of at least
`w_a(phi0)`, and the kernel's support is exactly that interval. So the smooth width is at least the
adaptive width at every latitude, by construction; it was not tuned to be. The acceptance measures
the smallest margin, smooth minus adaptive, as 0 at all three levels: the two are equal where the
adaptive width is at its peak or at the plain 4 deg.

**Smooth.** The Hann kernel and its slope vanish at its ends, so the width has a continuous
derivative.

This differs from the suggested Gaussian and smooth maximum, which the specification allows. A
smooth maximum with the adaptive width keeps the adaptive width's corners wherever that width is
the larger, and those corners are what made the spikes. The running maximum puts the smooth width
on or above the adaptive one everywhere, so the corners never return.

The fit at each latitude is otherwise unchanged; only its width differs. The gradient remains the
derivative of the fitted temperatures (`derivative_of_value`), with no filter applied afterward.

### Wiring

- **`lindal_wind.construct`** takes `window` and passes it to `fit`. A wind built with the smooth
  window says so in its `method` attribute.
- **The case `lindal_iris`** takes an optional `window`, `"adaptive"` when absent. So
  `lindal_iris_ii_110` is unchanged; acceptance 0 shows that it reproduces exactly.
- **`lindal_iris_v_ii_110` and `lindal_iris_v_sl_ii_110`** carry `window = "smooth"` in `[shear]`,
  with a comment and a line in their namelists' comments. They were rebuilt and rerun in place.

### Acceptance 1 (`reports/step13_1/shear_slope_value.png`)

The figure has three rows: the two windows' widths, then the gradient and the shear for three
cases. The cases are the slope with the adaptive window (SPEC_11), the derivative with the adaptive
window, and the derivative with the smooth window.

| Level | Fit RMS, adaptive | Fit RMS, smooth | Widest window | dT/dphi largest from the slope, adaptive / smooth | Shear largest from the slope, adaptive / smooth |
|---|---|---|---|---|---|
| 110 mbar | 0.449 K | 0.449 K | 8.53 deg | 1.42 / 1.27 K/deg | +39.6 / +36.5 m/s |
| 290 mbar | 0.496 K | 0.497 K | 8.55 deg | 1.32 / 0.99 K/deg | +21.0 / +15.5 m/s |
| 730 mbar | 0.706 K | 0.707 K | 8.45 deg | 1.41 / 0.96 K/deg | +31.2 / +21.2 m/s |

The fit's RMS is about its 234, 220 and 243 points. K/deg is K per planetocentric degree.

- **The spikes are gone.** The adaptive width rises and falls abruptly, within a degree or two, at
  -13, near -48 and across the IRIS gap. The smooth width spreads each
  rise over about ±4 deg. The smooth-window derivative follows the slope where the adaptive one
  spiked.
- **The largest differences from the slope that remain** are broad, not spikes:
  - at 110 mbar, near -48 and -13 deg, where the slope itself has narrow features;
  - between -10 and +5 deg at 110 and 730 mbar, the band bridged from its edge values, as before.
- **The fit is as good.** Its RMS about the points changes by at most 0.001 K.

### What the smooth window changes in the fitted temperatures

Where the window is wider, the fit averages over more latitude, and in the IRIS gap that moves the
value. The target, 31.2 S (31.11 S by the wind file's conversion), is inside the gap.

| Level | 31.2 S, adaptive | 31.2 S, smooth | Change from 36.3 N, adaptive | Change, smooth |
|---|---|---|---|---|
| 110 mbar | 86.63 K | 86.26 K | +5.21 K | +4.84 K |
| 290 mbar | 96.59 K | 96.51 K | +0.59 K | +0.51 K |
| 730 mbar | 114.92 K | 116.21 K | -4.38 K | -3.09 K |

The values at 36.3 N are the same with either window. The largest change in the fitted value
anywhere on the grid is 1.14 K at -47.8 deg (110 mbar), 1.13 K at -48.0 (290 mbar) and 1.31 K at
-31.3 (730 mbar), where the adaptive width rises abruptly. So the smooth window changes the temperatures the wind stands
for, not only their derivative. With either window the fit spans the gap from the dots on its two
sides (REVIEW_12_step1 ruling 1).

## Item 2. `profiles.png`

`reports/step13_2/profiles.png` draws:
- the IRIS fit's values at 36.3 N (30.81 N planetocentric) as open circles and at 31.2 S (26.06 S
  planetocentric) as filled squares, at 110, 290 and 730 mbar, for each window;
- in the difference panel, the 31.2 S value minus the egress.

Acceptance 2's output gives the six values for each window. The adaptive window's values are
SPEC_12's: 81.42, 96.00 and 119.30 K at 36.3 N, and 86.63, 96.59 and 114.92 K at 31.2 S.

## The rerun of Step 2 (acceptance 2 of §2a)

`lindal_iris_ii_110` is read as SPEC_12 left it. The other two were rebuilt with the smooth window
and rerun:

| Run | Outer loop | Pressure identity |
|---|---|---|
| `lindal_iris_v_ii_110` | 5 passes | at most 0.067 K |
| `lindal_iris_v_sl_ii_110` | 6 passes | at most 0.064 K |

**The round trip.** Each run is compared with its own fit, so the two new runs with the smooth
window's change.

| Level | Run | Delivered | IRIS change | Difference | (A) | (A1) | (A2) | (B) | (C) | (D) |
|---|---|---|---|---|---|---|---|---|---|---|
| 110 mbar | slope, adaptive, IP | +3.67 | +5.18 | **-1.50** | -1.21 | -1.152 | -0.061 | -0.11 | -0.09 | -0.09 |
| | value, smooth, IP | +4.71 | +4.81 | **-0.11** | -0.13 | +0.003 | -0.137 | +0.19 | -0.11 | -0.06 |
| | value, smooth, SL | +4.76 | +4.81 | **-0.06** | -0.13 | +0.003 | -0.137 | +0.26 | -0.11 | -0.07 |
| 290 mbar | slope, adaptive, IP | +1.31 | +0.61 | **+0.70** | +0.31 | +0.308 | +0.000 | +0.44 | -0.21 | +0.16 |
| | value, smooth, IP | +1.19 | +0.53 | **+0.66** | +0.01 | +0.013 | -0.002 | +0.62 | -0.15 | +0.18 |
| | value, smooth, SL | +1.26 | +0.53 | **+0.73** | +0.01 | +0.013 | -0.002 | +0.72 | -0.15 | +0.15 |
| 730 mbar | slope, adaptive, IP | -5.26 | -4.38 | **-0.88** | -0.07 | +0.013 | -0.082 | -0.88 | +0.18 | -0.10 |
| | value, smooth, IP | -4.06 | -3.14 | **-0.91** | -0.07 | -0.029 | -0.043 | -0.85 | +0.14 | -0.14 |
| | value, smooth, SL | -3.91 | -3.14 | **-0.77** | -0.07 | -0.029 | -0.043 | -0.67 | +0.14 | -0.16 |

All values are in K. IP is Ingersoll and Pollard, SL Sanchez-Lavega.

- **(A) is the linearization alone.** (A1) is within 0.03 K at every level; with the adaptive
  window it was 0.15 K at 730 mbar, where the target sits on a bend of the fit.
- **The band, (B), is now the largest part:** +0.6 to +0.7 K at 290 mbar and -0.7 to -0.9 K at
  730 mbar.

**Against the egress** (31.2 S):

| Level | slope, adaptive, IP | value, smooth, IP | value, smooth, SL | Egress sigma_T |
|---|---|---|---|---|
| 110 mbar | +3.91 K | +4.95 K | +5.00 K | 2.27 K |
| 290 mbar | +2.37 K | +2.25 K | +2.32 K | 2.29 K |
| 730 mbar | -0.95 K | +0.26 K | +0.40 K | 2.73 K |

- **110 mbar.** The runs still carry the IRIS warming, now +4.8 K. The radio shows none.
- **730 mbar.** With the smooth window the IRIS cooling to the south is smaller (-3.1 K against
  -4.4), and the runs end 0.3 to 0.4 K above the egress instead of 0.95 K below it.

**The winds** (`winds.png`):
- **Derivative with the smooth window, minus slope:** largest -93.0 m/s near 0.01 mbar at -13.6 deg,
  with an RMS of 8.7 m/s, against -98.1 and 8.8 m/s with the adaptive window. The narrow columns at
  about -48, -32 and -13 deg remain, with similar strength. The gradient there no longer spikes, but
  the derivative still differs from the slope over a few degrees near each of those latitudes (the
  shear panels of `shear_slope_value.png`). Integrated upward through the relaxation, a shear
  difference of 10 to 20 m/s becomes a wind difference of several tens of m/s at the top of the
  grid. So the columns measure how far the derivative and the slope differ there, not the spikes.
- **Sanchez-Lavega minus Ingersoll and Pollard:** the shear's change is largest at -10.5 deg
  (-45.6 m/s). Its RMS is 15.5 m/s within 15 deg of the equator and 1.0 m/s beyond, as before.

## Item 3. Data property files carry only values that matter

**In short.**
- **The files.** `data_static/winds/sanchezlavega2000.toml` and `smith1982_fig4.toml` are written as
  the Appendix gives them, except that `[frame]` and `[source]` are moved to the comments, because no
  code reads them.
- **The demands are dropped.** The wind tool and the schema no longer ask for the keys that went to
  the comments.
- **No value changes.** Both cloud winds, rebuilt, are identical in every value to the builds before
  this change, and both pass the schema. Only descriptive attributes are gone.

### 1. What the code reads

Only `casspian-wind-from-curve` (`tools/wind/build_wind.py`) reads these files. It now reads four
keys:
- `observation_level.value_Pa`, which selects the observed column;
- `latitude.convention`, which becomes `source_latitude_convention`, required of kind W;
- `epoch.value` and `epoch.solar_longitude_deg`, which SPEC_00 requires as values.

None of `rotation_system`, `positive_direction` or `wind_source` is read to compute or select
anything:
- the rotation system comes from the run's rotation file;
- `positive_direction` was never read;
- `wind_source` was only copied into a file attribute that nothing reads.

So `[frame]` and `[source]` are moved to the comments. The frame is already stated in the
Appendix's comments of both files. The Sanchez-Lavega file gains one comment line, "Wind source:
sanchezlavega2000."

### 2. Dropped from the wind tool

- `EPOCH_KEYS` no longer includes `solar_longitude_source`.
- `epoch_note` is no longer read.
- `observation_level.value_source`, `observation_level.justification` and
  `latitude.convention_source` are no longer read.
- For a tabulated source, the required `[source]` table (`TABLE_SOURCE_KEYS`: `citation`,
  `wind_source`, `value_source`, `uncertainty_long_name`, `uncertainty_method`, `method`,
  `caption`) is removed.

**What the wind files no longer carry, as a result:**
- **Global attributes, both winds:** `epoch_note`, `observation_level_justification`,
  `solar_longitude_source`, `source_latitude_convention_source`.
- **Variable attributes, both winds:** `value_source` on `reference_level_pressure_Pa`.
- **Sanchez-Lavega wind only:**
  - the global attributes `method`, `wind_source` and `source_caption_curve`;
  - `value_source` on `u_total_ms`;
  - `uncertainty_method` on `u_total_uncertainty_ms`.

  Its uncertainty's `long_name` is now a fixed text in the tool: "standard deviation of the
  measurements in the latitude bin".
- **The Ingersoll and Pollard wind** keeps the texts the tool itself writes: `method`, `source`,
  `wind_source`, the captions, and its variables' `value_source` and `uncertainty_method`. They
  were never read from its property file.

### 3. Dropped from the schema, for kind W only

- **`method`** and **`observation_level_justification`** are removed from kind W's required global
  attributes.
- **`solar_longitude_source`** is no longer required beside `solar_longitude_deg` on a kind W file.
  The rule is unchanged for every other kind.

SPEC_00 v0.24 makes these three official for kind W.

### 4. One demand kept, on other content

`source` is required of every kind (`AUTHOR_REQUIRED_GLOBALS`). The Sanchez-Lavega wind filled it
from `citation`, which is now a comment. It now takes the table file's own first line and the
file's name: "Saturn Zonal Winds from Sanchez-Lavega et al. (2000) (vasavada_saturn_winds.txt)".
This needed no change to the schema.

### 5. Checks

- **The Ingersoll and Pollard cloud wind**, rebuilt from `lindal_iris_ii_110`'s `[wind]` into the
  scratchpad, has all 12 variables identical in value. It has 34 global attributes, against 38
  before: the four listed in §2.
- **The Sanchez-Lavega cloud wind**, rebuilt from `tests/step13_1/sl_wind_build.toml`, is
  identical in every value to the build before this change and passes the schema.
- **Nothing in a run reads the dropped attributes.**
  - The run's loader compares only `solar_longitude_deg`, the anchor's against the namelist.
  - `lindal_wind` writes its own `method` and `observation_level_justification` on the shear
    wind.
- **No other code reads these files.** `lindal/build_inputs.py` reads a `convention_source` of
  its own, from the raw bundle.

**In the rerun.** Step 1's acceptance check 0 compares the rebuilt Ingersoll and Pollard wind
with `lindal_iris_ii_110`'s input file built before this change. It allows exactly the dropped
attributes and names them: `epoch_note`, `observation_level_justification`,
`reference_level_pressure_Pa.value_source`, `solar_longitude_source` and
`source_latitude_convention_source`. Everything else is identical.

### Files

- **Rewritten:**
  - `data_static/winds/sanchezlavega2000.toml`
  - `data_static/winds/smith1982_fig4.toml`
- **Changed:**
  - `src/casspian/tools/wind/build_wind.py`
  - `src/casspian/tools/wind/curve.py` (`Table.header`)
  - `src/casspian/lib/schema.py`

## Files for items 1 and 2

- **Changed:**
  - `src/casspian/tools/lindal/iris_temperatures.py` (`window`, `smooth_width`, `_adaptive_width`)
  - `src/casspian/tools/lindal/lindal_wind.py`
  - `src/casspian/tools/wind/shear.py`
  - the two new runs' control files
  - `tests/step13_1/accept_step13_1.py`
  - `tests/step13_2/accept_step13_2.py`

## For the SPEC_13 commit

SPEC_11 v0.5, REVIEW_11_step2 and SPEC_00 v0.24 go into the SPEC_13 acceptance commit.
