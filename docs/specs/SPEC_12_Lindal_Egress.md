# SPEC 12. The Voyager 2 egress profile, and the first comparison

Version 0.2, 8 October 2026. Author of record: S. Rafkin. **Status: accepted by the author, §2 ruled.**
Depends on SPEC_11 (the run `lindal_iris_ii_110`).

## 0. Purpose

Digitize the Voyager 2 egress temperature profile (Lindal et al. 1985, `docs/Lindal_et_al_1985_AJ90_1136.pdf`)
and compare it with the profile the transfer delivers at the egress latitude. This is the first
comparison against an independent observation: a different instrument, the same day.

Protocol as SPEC_04 §0. No regression: a new tool, and no change to the model.

## 1. Step 1

**Deliverables.**

1. `src/casspian/tools/lindal/egress_profile.py`. Render Lindal's Fig. 4 at 600 dpi, calibrate the
   axes by their ticks, and trace the dashed curve (31.2 S), the solid curve (36.3 N) and the
   sigma_T curve. Any decision by hand is coded in the script.
2. **Altitude to pressure.** Fig. 4's altitude is relative to each profile's 1 bar level. Pressure
   follows by hydrostatic integration from 1 bar with the profile's own temperature and Lindal's
   assumptions, the ones the anchor's reduction uses:
   - composition: `occul_data/lindal/lindal_composition.nc` (2.135 amu);
   - gravity: `lib.gravity` with `lindal_gravity.nc`, `lindal_rotation.nc` and `lindal_wind.nc`,
     at the profile's latitude and altitude;
   - radius: Lindal's 1 bar surface from `lindal_geodesy.nc` (60,268 km equatorial, 54,364 km
     polar) at the profile's planetocentric latitude, with the altitude added above it (§2 item 3).
3. `occul_data/lindal/voyager2_egress.csv`, committed: `altitude_km`, `pressure_mbar`,
   `temperature_K`, `temperature_uncertainty_K` (sigma_T).
4. **The target.** In `lindal_iris_ii_110`, `[target]` is set to Lindal's egress label, 31.2 S
   planetographic, converted to planetocentric as the anchor's label is, and the run is rerun. Only
   the latitude read from the transfer changes.

**Acceptance (`tests/step12_1/accept_step12_1.py`, under `reports/step12_1/`).**

1. **The ingress check.** The solid curve, traced and converted the same way at 36.3 N, against
   Lindal's Table I: temperature and pressure at the table's altitudes. The RMS differences are
   reported. They are the digitizing and conversion error of the egress profile.
2. **The overlay.** The traced curves over the scan, viewed by the author.
3. **The comparison:** `lindal_iris_ii_110`'s delivered temperature and density against pressure,
   with the egress profile and its sigma_T, and the anchor's for reference. The egress density is
   computed from its pressure and temperature with 2.135 amu.

## 2. Ruled (author, 8 October 2026)

1. **The target is Lindal's label, 31.2 S**, as the anchor uses his 36.3 N. The runs targeted
   -26.4 deg planetocentric, which is -31.49 planetographic by the model's rule, the far end of the
   egress swath (31.2 to 31.5 S).
2. **Fig. 4, not Fig. 5.** Fig. 5 plots temperature against pressure, but its three curves crowd
   together from 200 mbar to 1 bar. Fig. 4 has two curves on a linear altitude scale and carries
   sigma_T. The conversion uses Lindal's composition and gravity and is checked on the ingress.
3. **The radius is Lindal's 1 bar surface**, for the egress and the ingress check alike, not a
   run's delivered radius. The egress profile is data and depends on Lindal's inputs only, so every
   run is compared with the same profile. The two radii differ by a few km, a few hundredths of a
   percent in gravity.

## 3. Revision history

- v0.1, 8 October 2026: first draft.
- v0.2, 8 October 2026: accepted; §2 ruled (the target at Lindal's label, Fig. 4); the conversion
  with Lindal's composition, gravity and 1 bar radius.
