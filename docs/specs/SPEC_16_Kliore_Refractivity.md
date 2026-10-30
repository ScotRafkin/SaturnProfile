# SPEC 16. Pioneer 11 refractivity against the transferred Lindal anchor

Version 0.3, 30 October 2026. Author of record: S. Rafkin. **Status: closed, Step 1 accepted
(REVIEW_16_step1).**
Depends on SPEC_15.

## 0. Purpose

Refractivity is what an occultation measures. Temperature and pressure follow from it through an
assumed composition. Kliore et al. (1980, JGR 85, 5857) Table 1 tabulates Pioneer 11's refractivity
with radius, from 1.8 to 182 mbar (`occul_data/kliore1980/pioneer11_exit.csv`). This step compares
it with the refractivity the model delivers at 11.9 S (9.51 S planetocentric) from the Lindal
anchor. It also checks that our refractivity arithmetic reproduces Kliore's own.

It prepares Kliore as a second profile (the next specification) and needs no new run: the two
SPEC_15 runs at 11.9 S are used, `lindal_iris_v_sl_ii_110_12s` and `lindal_iris_v_sl_iii_110_12s`.

Protocol as SPEC_04 §0. No regression: acceptance scripts only, no change to the model or tools.

**No new keys, columns or parameters.** Kliore's oblateness, his latitude and the two readings of
the radius mapping are constants in the acceptance script, each with its source in a comment.
Nothing is added to the data file, a control file or a product.

## 1. Step 1

**Deliverables**, all in `tests/step16_1/accept_step16_1.py` and its output.

1. **The self-check.** Kliore's temperature and pressure as published, with his mixture (90 percent
   hydrogen, 10 percent helium by number) and the species refractivities of
   `data_static/species_master.toml`, through the reduction's own refractivity function. The result
   against his tabulated refractivity column, in n-units ((n - 1) x 10^6), level by level.
2. **Kliore's radius, back to the local radius.** His column is the measured radius mapped to the
   equator with a constant oblateness of 0.096 (his text): r_local = r_eq (1 - f) /
   sqrt((1 - f)^2 cos^2 phi + sin^2 phi), with f = 0.096. He gives the latitude as 9.5 S, and his
   latitudes are planetocentric (SPEC_15 §2). Whether he applied the mapping with planetocentric or
   planetographic latitude is not stated. Both readings are carried; at this latitude they differ
   by about 80 km.
3. **Refractivity against radius,** which needs no composition. Delivered ln N against
   `radius_m` from each run, and Kliore's ln N against each reading of r_local.
   - The radius offset that best aligns them, by least squares in ln N over Kliore's levels, and
     the residual after it.
   - The offset itself is a result. It combines the radius mapping, the Pioneer 11 orbit, and the
     model's reference surface at 9.51 S. For scale, Kliore puts the 100 mbar level at about
     60,332 km equatorial radius.
4. **Refractivity against pressure,** which needs the composition. Kliore's pressure rescaled to
   the anchor's mixture by SPEC_15's rule (`radio_profiles.rescale_composition`). His refractivity
   is not rescaled: it is measured. Delivered N against Kliore's N at the rescaled pressures.

**Acceptance (under `reports/step16_1/`).** No bounds: this step measures.

1. **The self-check:** RMS and largest difference, in n-units and percent, with the pressure where
   the largest occurs. A difference beyond the table's rounding (0.1 n-units) is explained. A
   systematic offset means our species refractivities differ from Kliore's, and the report says by
   how much.
2. **For the author:**
   - `refractivity_radius.png`: ln N against radius, both runs and both readings of Kliore's
     radius, before and after the offset; a second panel, the residual in percent;
   - `refractivity_pressure.png`: N against pressure, both runs and Kliore rescaled, with Kliore's
     pressure ranges; a second panel, delivered over observed minus one, in percent;
   - a table, for each run and each reading: the radius offset, the residual RMS in ln N after it,
     and the RMS in pressure space, all in percent; and delivered over observed at 10, 30, 100 and
     150 mbar.

## 2. Ruled (author, 30 October 2026)

1. **Both readings of the radius mapping** are carried, and the fitted offset absorbs the frame.
   The pressure comparison is the one that does not depend on it.
2. **Next:** if the self-check agrees and the comparison is sensible, the next specification
   reduces Kliore's table to a kind N and runs it as a weight-0 validation anchor beside Lindal.
   The comparison logic stays in this suite until the assimilation needs an observation operator.
3. **Information goes in comments,** not keys, unless code requires it as a value.

## 3. Revision history

- v0.1, 30 October 2026: first draft.
- v0.2, 30 October 2026: accepted; §2 ruled; constants in comments, no new keys.
- v0.3, 30 October 2026: closed; the planetocentric reading for what follows; the 17.1 mbar
  misprint corrected to 17.5.
