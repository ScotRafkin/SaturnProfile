# SPEC 13. The gradient from the fitted temperatures, the Sanchez-Lavega wind, and the second run

Version 0.5, 8 October 2026. Author of record: S. Rafkin. **Status: accepted by the author, §3 ruled;
v0.4 adds Step 3 (§2a) on the Step 1 and 2 reports; v0.5 gives the two data property files in full
(Appendix).**
Depends on SPEC_11 and SPEC_12.

## 0. Purpose

Two changes to the wind, each in its own run:

1. The IRIS wind is built from the derivative of the fitted temperatures rather than from the fit's
   local slope (REVIEW_11_step2, carried item 1). The wind then agrees with the temperatures it
   stands for.
2. A second cloud wind from Sanchez-Lavega, Rojas and Sada (2000), Icarus 147, 405
   (`docs/Sanchez-Lavega_etal_cleanPDF.pdf`), in a second run identical to the first except for
   that wind.

The table is `data_static/winds/vasavada_saturn_winds.txt`, which matches the paper's Table II row
for row. Voyager 1 and 2 only, System III. `u_rms` is the standard deviation of the measurements in
each 0.5 deg bin (their Eq. 6), the same kind of number as the uncertainty of the Ingersoll and
Pollard wind.

Protocol as SPEC_04 §0. No regression.

## 1. Step 1. The winds

1. **The gradient.** A parameter `gradient` of the case `lindal_iris`: `"slope"`, `fit`'s
   `gradient_K_per_deg` as SPEC_11 built it and the default, so the existing run is unchanged; or
   `"value"`, the derivative in latitude of `fit`'s `value_K`. Everything else in the construction
   is unchanged, and `fit` itself is unchanged.
2. **The Sanchez-Lavega cloud wind**, through the wind tool's curve path, as the Ingersoll and
   Pollard wind is built:
   - the planetographic column and `u`, with `u_rms` as the uncertainty; the two rows at +/-90
     marked "fake data added here" are dropped;
   - **the ring gap** (-2.1 to -10.7 deg planetographic) bridged by PCHIP between its data edges,
     with no reflection;
   - the small gaps by PCHIP; the poles by the tool's polar rule.

**Acceptance (`tests/step13_1/accept_step13_1.py`, under `reports/step13_1/`).** For the author:

1. the shear at the three levels from the slope (SPEC_11) and from the derivative, against latitude;
2. the two cloud winds against latitude, with the table's points and `u_rms`.

## 2. Step 2. The runs

Two new runs, each made by `casspian-new-run` from the one before, so that each differs from it in
one thing. `lindal_iris_ii_110` is kept as SPEC_11 left it.

1. **`lindal_iris_v_ii_110`**, from `lindal_iris_ii_110`, with `gradient = "value"`.
2. **`lindal_iris_v_sl_ii_110`**, from `lindal_iris_v_ii_110`, with `[wind]` set to the
   Sanchez-Lavega wind.

Every kind file a run is driven by stays in that run's build.

**Acceptance (`tests/step13_2/accept_step13_2.py`, under `reports/step13_2/`).**

1. **The round trip** of each run, as SPEC_11 Step 2's, reported. With the gradient from the fit's
   value, part (A) of REPORT_11_step2 should vanish.
2. **For the author:**
   - the delivered temperature and density against pressure from the three runs, with the egress
     profile and its sigma_T (SPEC_12);
   - `u_total(phi, p)` of the three runs' winds, and the difference of each from the one before.

## 2a. Step 3. Before acceptance (author, 8 October 2026, on REPORT_13_step1 and step2)

1. **A smooth window.** The derivative's spikes come from the fit's window widening abruptly where
   the data thin. They are not physical and are removed at their source.
   - `fit` gains an option for a window whose width is a smooth function of latitude and never
     narrower than SPEC_10's adaptive width. A suggestion: smooth the adaptive width in latitude
     with a Gaussian of a few degrees, then take a smooth maximum with the adaptive width. The
     coding agent may propose another way that meets the requirement.
   - The default stays SPEC_10's adaptive window, so `lindal_iris_ii_110` is unchanged.
   - The case `lindal_iris` takes the option as a parameter. `lindal_iris_v_ii_110` and
     `lindal_iris_v_sl_ii_110` use it and are rebuilt and rerun in place.
   - Not a filter applied to the gradient afterward: the gradient must stay the derivative of the
     fitted temperatures.
2. **`profiles.png`** shows the IRIS fit's values at 36.3 N and 31.2 S at the three levels, as in
   SPEC_12's figure.
3. **Data property files carry only values that matter.** In
   `data_static/winds/sanchezlavega2000.toml` and `smith1982_fig4.toml`, a key stays only if code
   uses its value to compute or select something, or SPEC_00 requires it as a value (the epoch and
   the solar longitude). Everything descriptive moves to the comments at the top: the `*_source`
   keys, `justification`, `enters_reduction`, `epoch_note`, `citation`, `method`, `caption`,
   `uncertainty_long_name`, `uncertainty_method`. Where the tool or the schema now demands one of
   them, the demand is dropped, and the report lists what was dropped. The two files are given in
   full in the Appendix and are written exactly so. If the code reads none of
   `rotation_system`, `positive_direction` or `wind_source`, that key goes to the comments too.

**Acceptance (`reports/step13_1/` and `reports/step13_2/`, rerun).**

1. The gradient and the shear, adaptive and smooth windows, as `shear_slope_value.png`, for the
   author. The fit's RMS about the points with each window.
2. Step 2's acceptance, rerun with the smooth window.

## 3. Ruled (author, 8 October 2026, on REPORT_13_preexecution)

1. **The curve path takes the table** through three new control values, added beside the existing
   ones (the author's leave under "replace, never add"):
   - `curve_format = "table"`: the four columns read as one segment, the +/-90 rows dropped;
   - `gap_rule = "pchip_bridge"`: one PCHIP through every row, with no reflection and no join;
   - `uncertainty_source = "table"`: `u_rms`, interpolated linearly in latitude.

   The Ingersoll and Pollard path and its values are unchanged. Only the tool that writes kind W
   changes; the model does not.
2. **The uncertainty in the gaps.** NaN in the ring gap (-2.1 to -10.7 deg planetographic), where the
   scatter was never measured. Interpolated across the small gaps. Wind values in the ring gap are
   `interpolated`; the rows are `observed`; the polar caps are `extrapolated`.
   - Before the runs, the coding agent reports whether anything in a run reads the wind's
     uncertainty, and what a NaN there does to the delivered profile.
   - How NaN uncertainties are used later (the Monte Carlo, the assimilation) is deferred.
3. **The shear changes with the cloud wind.** `construct` solves for the shear given the wind,
   through the meridional term `(cos phi / r)(du/dphi)_p` and the isobar geometry. So between
   `lindal_iris_v_ii_110` and `lindal_iris_v_sl_ii_110` the reference, the shear and the total all
   change, the shear mostly near the equator. Only the IRIS input to the shear is the same. If the
   implementation would hold the shear fixed, that is a fault to report.
4. **The rest of the pre-execution report stands:**
   - the gradient by central differences of `value_K`, one-sided at the ends;
   - the fit on the cloud wind's latitude grid;
   - one line in `vasavada_saturn_winds.note.md` citing SPEC_13 for the epoch, the frame and
     `u_rms`;
   - the Sanchez-Lavega wind built into the run's own inputs.

## Appendix. The two data property files

`data_static/winds/sanchezlavega2000.toml`:

```toml
# Sanchez-Lavega, A., Rojas, J. F., and Sada, P. V. 2000, Icarus 147, 405, Table II, as received in
# vasavada_saturn_winds.txt (see its note). Mean zonal wind of the loess-smoothed Voyager 1 and 2
# cloud-tracking points in each 0.5 deg bin (their Eq. 5); u_rms is the standard deviation of the
# measurements in each bin (their Eq. 6), linear in latitude between rows, NaN in the ring gap and
# poleward of the rows (SPEC_13 ruling 2).
#
# Latitude planetographic (their abstract). System III, eastward positive (their section 2).
# Observation level: 1 bar, a project assignment as for Ingersoll and Pollard; the source states no
# pressure. A run with the IRIS shear reassigns the cloud wind to 398 mbar (SPEC_11).
# Epoch: the profile combines Voyager 1 (Nov 1980) and Voyager 2 (Aug 1981) imaging and carries the
# date of the Voyager 2 occultation it serves. Ls 18.2 deg computed for that date (15 Sep 2026, JPL
# approximate elements); to be replaced by the ephemeris value.

[observation_level]
value_Pa = 1.0e5

[latitude]
convention = "planetographic"

[frame]
rotation_system = "System III"
positive_direction = "eastward"

[epoch]
value = "1981-08-26"
solar_longitude_deg = 18.2

[source]
wind_source = "sanchezlavega2000"
```

`data_static/winds/smith1982_fig4.toml`:

```toml
# Ingersoll and Pollard (1982) Fig. 5 curve is the wind; the Smith et al. (1982) Fig. 4 points are
# the uncertainty. Both digitized beside this file.
#
# Latitude planetographic (Ingersoll and Pollard Fig. 3 caption; Lindal et al. 1985 Appendix).
# System III, eastward positive (Smith et al. Fig. 4 caption).
# Observation level: 1 bar, a project assignment; neither source states a pressure. The reduction
# wind is altitude independent, so the level sets only which column is flagged observed.
# Epoch: the curve combines Voyager 1 (Nov 1980) and Voyager 2 (late Aug 1981) imaging and carries
# the date of the Voyager 2 ingress occultation it serves. Ls 18.2 deg computed for that date
# (15 Sep 2026, JPL approximate elements); to be replaced by the ephemeris value.

[observation_level]
value_Pa = 1.0e5

[latitude]
convention = "planetographic"

[frame]
rotation_system = "System III"
positive_direction = "eastward"

[epoch]
value = "1981-08-26"
solar_longitude_deg = 18.2
```

## 4. Revision history

- v0.1, 8 October 2026: first draft.
- v0.2, 8 October 2026: the gradient from the derivative of the fitted temperatures (author,
  8 October), as a case parameter; two new runs, each one change from the last.
- v0.3, 8 October 2026: accepted; §3 ruled on the pre-execution report (the table's options, NaN in
  the ring gap, the shear changes with the cloud wind).
- v0.4, 8 October 2026: Step 3: a smooth window for the fit, the IRIS values in `profiles.png`,
  descriptive keys of the data property files moved to comments.
- v0.5, 8 October 2026: the two data property files given in full (Appendix).
