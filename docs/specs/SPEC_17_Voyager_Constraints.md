# SPEC 17. What the Voyager-epoch observations require of the wind

Version 0.3, 5 November 2026. Author of record: S. Rafkin. **Status: closed, Step 1 accepted
(REVIEW_17_step1).**
Depends on SPEC_15 and SPEC_16.

## 0. Purpose

The wind is two parts: the shear, constrained by temperatures, and the reference wind, which with
gravity and rotation sets the shape of the isobars. This step measures what the Voyager-epoch
observations ask of each, as numbers a fitted wind can use. It fits nothing.

**Near-epoch (author, 30 October 2026).** The Voyager 1 profiles (November 1980) are used with the
Voyager 2 anchor (August 1981) as near-epoch, since nine months of season changes the wind little.
Once time propagation exists they are treated as their own times. Pioneer 11 (September 1979) is
two years away and enters nothing: it is shown for reference only.

Protocol as SPEC_04 §0. No regression: a small data file, one function beside the balance, and
acceptance scripts. No change to the model.

## 1. Step 1

### Part A. The isobar geometry: the 100 mbar surface

1. **Lindal Fig. 9, right panel, digitized.** Five points: the altitude of the 100 mbar surface
   above Lindal's reference geoid, against planetocentric latitude. The reference geoid has a
   polar radius of 54,438 km and the shape set by the external gravity and uniform System III
   rotation (10 h 39 m 22.4 s), with no wind. Rendered at 600 dpi and calibrated by the ticks, as
   in SPEC_15. File: `occul_data/lindal/isobar_100mbar.csv`, columns
   `latitude_planetocentric_deg`, `altitude_km` and `observation` (`voyager1`, `voyager2` or
   `pioneer11`, which the code reads to leave Pioneer out). Source, definition of the reference
   geoid and epochs go in the top comment. Lindal gives no error bars on the points; his geoid fits
   them to 3 to 4 km.
2. **The model's side.** The same reference geoid, from `lib.geoid` with no wind, the run's gravity
   and rotation, and the 54,438 km polar radius. The delivered radius of the 100 mbar surface minus
   it, at each point's latitude:
   - at the anchor, from the anchor's own product;
   - at the three Voyager targets, from the SPEC_15 runs of both hypotheses, at the point's own
     latitude if the product carries the isobar's geometry there; otherwise at the target's
     latitude, with the latitude difference reported.

   100 mbar is the gauge isobar of the Lindal-line runs (`gauge_isobar_Pa = 1.0e4`), so this is
   the isobar whose shape the transfer carries directly.

### Part B. The temperatures: the shear, by latitude segment

3. **The weight.** The model's balance (the relation `construct` uses, `tools/lindal/lindal_wind.py`),
   linearized in the shear with the wind held fixed, gives on each isobar

     Δ(d ln T / dφ) = w(φ, p) Δs(φ, p),

   with w from the run's own gravity, radius, absolute rotation and pressure gradient. A function
   beside `construct` returns w on the run's grid. It is the linearized observation operator the
   fitted wind will use. Each hypothesis uses its own run's w.
4. **Segment constraints.** Every path starts at the anchor, 30.81 N planetocentric. With the
   misfit m_k(p) = ln(T_obs / T_delivered) at each target k, at each anchor level, the observed
   profile taken to the delivered levels linearly in ln p, and NaN where it has no data:
   - segment 1, 30.81 N to 2.37 S: m(3.0 S);
   - segment 2, 2.37 S to 26.06 S: m(31.2 S) minus m(3.0 S);
   - segment 3, 26.06 S to 71.07 S: m(74.2 S) minus m(31.2 S).

   For each segment and level: the required integral of w Δs, the integral of w over the segment,
   and their ratio, the segment-mean Δs in m/s per unit ln p. Segment 1 crosses the equator, where
   w changes sign; its weight integral is reported so the reader can judge the ratio. Segment 3
   carries the Fig. 5 depth bias (REVIEW_15_step1: Lindal's drawn line 1.5 to 3 K warm from 600 to
   1200 mbar), which is marked there and not corrected.
5. **The linearization, tested.** In a copy of the tree, not a registered run: the base wind of
   `lindal_iris_v_sl_ii_110` plus a known Δs, uniform over segment 2's latitudes, tapered smoothly
   to zero over 2 deg at each end, and over 10 to 30 mbar, tapered over half a scale height at each
   end. The wind is integrated from the reference level as the tool does. Run to the three
   Voyager targets at two amplitudes, c and 2c, with c chosen to move the delivered temperature by
   about 2 K. Predicted against delivered change at each target and level.

**Acceptance (`tests/step17_1/accept_step17_1.py`, under `reports/step17_1/`).** No bounds: this
step measures.

1. **The linearization:** at each target, delivered change over predicted in the perturbed band,
   for c and 2c. The departure from linear, and from one, is explained, in particular the share
   from the wind's own change in the balance, which the linearization holds fixed.
2. **For the author:**
   - `isobar_100mbar.png`: Lindal's points and the delivered altitude above the reference geoid
     at each, both hypotheses, with Lindal's best-fit geoid if it can be traced; Pioneer marked
     and not compared;
   - a table of delivered minus observed 100 mbar altitude, in km, at each Voyager point;
   - `weight.png`: w(φ, p) on the base run's grid;
   - `segments.png`: the segment-mean Δs against pressure, three segments by two hypotheses; the
     equatorial profiles above 30 mbar marked as carrying wave structure (SPEC_15);
   - a table of the segment-mean Δs at 730, 290, 110, 30, 10, 4.5 and 1 mbar.

## 2. Notes (accepted with the specification)

1. **The perturbation's shape.** Uniform over a segment and a pressure band is chosen because it
   is easy to predict and to read. It tests the operator, and is not itself a hypothesis.
2. **Next.** The fitted wind (the next specification) solves for Δs on a smooth basis, against
   these segment constraints, the IRIS gradients and the 100 mbar geometry together, weighted by
   their errors. It needs the ruling on the 110 mbar radio and IRIS disagreement first.

## 3. Revision history

- v0.1, 30 October 2026: first draft.
- v0.2, 30 October 2026: accepted; the misfit's interpolation, each run's own weight, the Fig. 5
  depth bias marked.
- v0.3, 5 November 2026: closed.
