# CASSPIAN status, 8 October 2026

Reviewing agent. Supersedes `casspian_status_2026-10-03.md` for the wind work.

## Closed

- **SPEC_08:** kind W redefined. The reference wind and its level are required. The total and the
  shear are optional.
- **SPEC_09:** the IRIS temperatures at 110, 290 and 730 mbar digitized (Conrath and Pirraglia 1983
  Fig. 1).
- **SPEC_10:** `fit`, local linear, 4 deg FWHM.
- **SPEC_11:** the Lindal-only wind, case (ii), and the run `lindal_iris_ii_110`.
- **SPEC_12:** the Voyager 2 egress digitized from Lindal Fig. 4 and converted with Lindal's
  composition, gravity and 1 bar radius. The target is set at Lindal's label, 31.2 S (26.06 S
  planetocentric).
- **SPEC_13:** the gradient from the derivative of the fitted temperatures, a smooth fit window,
  and the Sanchez-Lavega cloud wind. Two runs. Accepted pending the author's look at Step 3's
  figures; committed at `cf33cd6` with SPEC_00 v0.24 (descriptive text in comments) and SPEC_11
  v0.5.

## What the first comparison showed

Delivered minus egress, and the untransferred anchor minus egress:

| Pressure | Delivered - egress | Anchor - egress |
|---|---|---|
| 0.3 mbar | +15.0 K | +14.7 K |
| 1 mbar | +5.1 K | +4.9 K |
| 10 mbar | -2.3 K | -3.7 K |
| 100 mbar | +3.9 K | +0.3 K |
| 300 mbar | +2.1 K | +0.9 K |
| 700 mbar | -0.7 K | +4.2 K |

- **The transfer works.** Below about 300 mbar it beats the anchor.
- **The tropopause region.** From 50 to 300 mbar the delivered profile is too warm. It carries the
  IRIS north-south change (+5.2 K at 110 mbar), which the two radio profiles do not show (+0.03 K).
  The IRIS minus radio offset flips sign between the two latitudes at 110 mbar, while at 730 mbar
  it is the same at both. Composition cannot cause a sign flip.
- **Aloft.** The relaxed shear leaves the delivered profile equal to the anchor's. The radio's
  egress minus ingress changes sign near 4 to 5 mbar: warmer to the south below, colder above.

## What SPEC_13 showed

Delivered minus egress at 31.2 S:

| Level | `lindal_iris_ii_110` | `lindal_iris_v_ii_110` | `lindal_iris_v_sl_ii_110` | Egress sigma_T |
|---|---|---|---|---|
| 110 mbar | +3.9 K | +5.0 K | +5.0 K | 2.3 K |
| 290 mbar | +2.4 K | +2.3 K | +2.3 K | 2.3 K |
| 730 mbar | -1.0 K | +0.3 K | +0.4 K | 2.7 K |

- **The round trip closes:** within 0.11 K at 110 mbar. The model's remaining part is the
  equatorial bridge, 0.6 to 0.9 K at 290 and 730 mbar.
- **The cloud wind barely matters for temperature.** Sanchez-Lavega changes the delivered profile
  by at most 0.14 K at the three levels, though it changes the wind by up to 100 m/s.
- **The IRIS gap.** 31.2 S lies in the IRIS data gap (-27.7 to -38.5 deg). The fit window alone
  moves the IRIS value there by up to 1.3 K (730 mbar).
- **110 mbar** is now the IRIS and radio disagreement in full, about 5 K.

## Decisions (author, 8 October)

- **Each wind is a hypothesis** with its own kind W file and its own run. Earlier hypotheses are kept,
  not overwritten. Hypothesis 1 is case (ii), `lindal_iris_ii_110`.
- **The shear aloft has always been a guess.** The relaxation's start, rate and limit have no data
  behind them. A shape informed by the radio is acceptable, and the egress then checks consistency,
  not independence.
- **The assimilation later pulls toward Lindal's radio egress** (as kind T), starting from
  hypothesis 1, to show the delivered profile moving toward the observation.
- **The regression rule** (judgment of risk), push after each accepted commit, and lean
  documentation are unchanged.

## Next

1. **Fix what the regression found (8 October, on `cf33cd6`, SPEC_13's commit).** Three things, no
   new physics. Detail is in the coding agent's handoff,
   `docs/CASSPIAN_CodingAgent_Handoff_2026-10-22.md` §4.
   - The closure run is refused: its comparison counts five descriptive attributes that SPEC_00
     v0.24 removed. Compare values and ignore descriptive attributes; then run step05_3 checks 4
     to 7, which never ran.
   - SPEC_00 §2.2's listing lacks `iris_temperatures.csv` and `voyager2_egress.csv` (step9 check 6,
     failing since SPEC_09).
   - The suites step09_1 to step13_2 passed by hand but are not in `tests/run_regression.sh`. The
     full set takes about 2 h 56 min, so it runs in chunks.
2. **SPEC_14, after SPEC_13's results:** hypothesis 2, case (iii) "reversed". Above 110 mbar the
   110 mbar shear is multiplied by f(p):
   - f = 1 to p_hold, first guess 15 mbar;
   - linear in ln p through zero at p_cross, first guess 4.5 mbar, continuing at the same slope;
   - held above 0.2 mbar.

   Its own run and kind W. The base is `lindal_iris_v_sl_ii_110`. Watch the equatorial winds
   aloft.
3. **More radio targets:** Lindal Fig. 6 (Voyager 1 egress near 3.0 S, Nov 1980, inside the
   equatorial band; Pioneer 11 near 11.9 S, 1979) and the 74.2 S profile of Fig. 5, digitized in
   one step. Hypotheses 1 and 2 then compared at several targets. The 74.2 S profile changes order
   near 3 mbar, like ingress and egress near 4 to 5 mbar. A pattern times f(p) shear cannot change
   the latitude pattern with height. Lindal rates the Voyager 1 profiles less reliable.
4. **Further wind hypotheses, saved for later:**
   - **Hypothesis 3:** shear from a CIRS cross-section of the same season (Ls 0 to about 16 deg,
     Aug 2009 to Nov 2010, before the storm at 33 to 40 N), nadir and limb up to about 0.01 mbar.
     Interannual change is small at midlatitudes (Fletcher et al. 2016) and large in the tropics:
     about 10 K at the tropical tropopause (Li et al. 2013) and 7.1 K at 2.1 mbar on the equator
     (Sinclair et al. 2014). So CIRS shear is used outside about +/-15 to 20 deg, with larger
     uncertainty or other data across the equator.
   - **A seasonal wind prior** for seasons without observations. This reopens the Seasonal Design
     Note's decision that no shear is derived from a propagated thermal field.
5. **Later:**
   - the IRIS profile catalog;
   - stratospheric IRIS data near 1 mbar;
   - the dry adiabat for normal hydrogen (about 0.85 to 0.9 K/km);
   - the Monte Carlo;
   - convective adjustment that mixes entropy, momentum and composition together;
   - cylinders at depth.
