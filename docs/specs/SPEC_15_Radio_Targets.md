# SPEC 15. Three more radio profiles, and the two hypotheses at four targets

Version 0.3, 30 October 2026. Author of record: S. Rafkin. **Status: closed, Step 1 accepted
(REVIEW_15_step1).**
Depends on SPEC_12 and SPEC_14.

## 0. Purpose

Hypotheses 1 and 2 have been compared at one target, the Voyager 2 egress at 31.2 S. Three more
radio profiles of the same era become targets for both:

- **Voyager 1 ingress at 74.2 S**, Lindal et al. (1985) Fig. 5
  (`docs/Lindal_et_al_1985_AJ90_1136.pdf`);
- **Voyager 1 egress near 3.0 S**, Lindal et al. Fig. 6;
- **Pioneer 11 exit at 11.9 S**, Kliore et al. (1980, JGR 85, 5857) Table 1, which tabulates the
  profile with its uncertainties. Lindal's Fig. 6 redraws it.

The targets test what 31.2 S could not:
- **3.0 S** lies in the equatorial band, where the shear is the PCHIP bridge and not computed;
- **11.9 S** lies just outside the band;
- **74.2 S** lies in the polar taper. Its order relative to the other profiles changes near
  3 mbar, which a latitude pattern times f(p) cannot reproduce.

Lindal rates the Voyager 1 profiles less reliable than Voyager 2's: scintillation in the
equatorial jet, ring attenuation, and long horizontal paths. Hypothesis 2 was tuned to 31.2 S, so
these three are its first independent test.

Protocol as SPEC_04 §0. No regression: a new tool and new runs, no change to the model.

## 1. Step 1

**Deliverables.**

1. **`src/casspian/tools/lindal/radio_profiles.py`.** Render Lindal's Figs. 5 and 6 at 600 dpi,
   calibrate the axes by their ticks (pressure logarithmic, 0.1 to 10^4 mbar; temperature 80 to
   160 K), and trace:
   - Fig. 5: 36.3 N (solid), 31.2 S (dashed), 74.2 S (dash-dot);
   - Fig. 6: 3.0 S (solid), and 11.9 S (stippled) as a check only.

   Where two curves cannot be told apart in the scan, they coincide within the line's width, and
   both take the shared trace. The report lists those pressure ranges. The horizontal error bars
   are read as temperature uncertainties at the pressures where they are drawn, assigned to the
   curve they sit on. Any decision made by hand is coded in the script. Both figures were computed
   for 94 percent hydrogen and 6 percent helium by number, the anchor's mean molar mass of
   2.135 amu, so nothing is rescaled.
2. **Kliore et al. Table 1, transcribed as published:** equatorial radius, refractivity,
   temperature, pressure, and the temperature and pressure uncertainty ranges where tabulated.
   Two misprints are corrected, each with a comment: 63.7 K at 142.5 mbar is 93.7 K, and the
   pressure range "-02." at 1.8 mbar is -0.2. Each correction is checked against the scanned
   page.
   - **The composition.** The table is computed for 90 percent hydrogen and 10 percent helium.
     Where the comparison reads it, it is rescaled to the anchor's mixture by Lindal's rule (Fig. 5
     caption): temperature in proportion to the mean molecular mass, pressure in proportion to the
     mean molecular mass over the mean molecular refractivity. The rescaling takes the two mixtures'
     values from `data_static/species_master.toml`, in one function, and is about -3.6 percent in
     temperature. The file keeps the published values.
   - **The uncertainties** are Kliore's ranges, which are extremes from the baseline and the
     initial temperature, not one sigma. They are rescaled with the temperatures. Kliore's
     additional +/-3 K for composition is not added, since the rescaling replaces that assumption.
3. **Data files**, committed:
   - `occul_data/lindal/voyager1_ingress.csv` (74.2 S) and `voyager1_egress.csv` (3.0 S), with
     columns `pressure_mbar`, `temperature_K` and `temperature_uncertainty_K`, the last NaN except
     where a bar is drawn;
   - `occul_data/kliore1980/pioneer11_exit.csv`, the table as published.

   Each file's top comment gives the source, the figure or table, the latitude as published and
   the epoch. There is no key for any of them. The calibration traces and the Fig. 6 check trace
   are not committed.
4. **Reference papers.** Both Kliore et al. (1980) papers go into `docs/`, beside Lindal et al.:
   the JGR paper and the Science report (Science 207, 446), which is the preliminary analysis.
5. **The targets,** converted to planetocentric as SPEC_12 converted 31.2 S: 74.2 S, 3.0 S and
   11.9 S planetographic.
6. **Six runs,** made with `casspian-new-run --from`. Each differs from its parent in `[target]`
   only, and the build file's comment gives the observation's epoch:
   - from `lindal_iris_v_sl_ii_110` (hypothesis 1's base): `lindal_iris_v_sl_ii_110_74s`, `_03s`
     and `_12s`;
   - from `lindal_iris_v_sl_iii_110` (hypothesis 2): `lindal_iris_v_sl_iii_110_74s`, `_03s` and
     `_12s`.

**Acceptance (`tests/step15_1/accept_step15_1.py`, under `reports/step15_1/`).**

1. **The calibration traces.** The Fig. 5 trace of 36.3 N against Lindal's Table I, and of 31.2 S
   against `voyager2_egress.csv` (SPEC_12, Fig. 4). RMS and largest difference in kelvins, each
   with the pressure where it occurs. A largest difference over 1 K is explained.
2. **Pioneer 11, two ways.** The Fig. 6 trace of 11.9 S against Kliore's table, published and
   rescaled. This shows whether Lindal rescaled it, and how well the rescaling agrees.
3. **The tracing overlays,** for the author: each figure with every trace drawn over it.
4. **The pressure identity** of each run, in kelvins, level by level.
5. **For the author,** at each of the four targets (31.2 S from SPEC_14's two runs):
   - the delivered temperature and density against pressure, both hypotheses, with the anchor,
     the observed profile and its uncertainties, and the IRIS fit where IRIS has data, as
     `profiles.png` of SPEC_14;
   - delivered minus observed, both hypotheses and the untransferred anchor, at 730, 290, 110, 10,
     4.5, 1 and 0.3 mbar where the profile has data; NaN where it does not;
   - one summary figure, delivered minus observed against pressure, four targets by two
     hypotheses.

## 2. Ruled (author, 30 October 2026)

1. **Latitude conventions.** Lindal labels in planetographic latitude. Kliore et al. give the
   Pioneer 11 exit at about 9.7 S, which is Lindal's 11.9 S in planetocentric, so Kliore's
   latitudes are planetocentric and the label is settled. 74.2 S sits beside 36.3 N and 31.2 S in
   Fig. 5 and is planetographic.
2. **Epochs** go in file and code comments only, with no key or parameter: Voyager 1 on
   13 November 1980; Pioneer 11 on 1 September 1979. The runs stay at the anchor's season
   (26 August 1981), and the difference is not modeled.
3. **Error bars:** the drawn bars, and Kliore's tabulated ranges, with NaN elsewhere.

## 3. Revision history

- v0.1, 30 October 2026: first draft.
- v0.2, 30 October 2026: accepted; Pioneer 11 from Kliore et al. (1980) Table 1, rescaled to the
  anchor's composition by Lindal's rule; §2 ruled.
- v0.3, 30 October 2026: closed.
