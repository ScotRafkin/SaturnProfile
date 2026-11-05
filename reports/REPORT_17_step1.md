# REPORT 17, Step 1. Before acceptance

Coding agent, 30 October 2026. Specification: `docs/specs/SPEC_17_Voyager_Constraints.md` v0.2.
Working tree on `main` at `9e01db4`, uncommitted. No regression (SPEC_17 §0): one function added
beside `construct`, a data file, and scripts. No change to the model; no registered run changed.

Acceptance: `tests/step17_1/accept_step17_1.py` **2 of 2**; output in `reports/step17_1/output.txt`
and `run.log`, figures and `summary.json` beside it. The six SPEC_15 runs at the three Voyager
targets are read as they are. The perturbed runs are copies under `reports/step17_1/perturbed/`
(git-ignored).

**In short.**
- **The 100 mbar surface.**
  - At the anchor, the delivered surface lies 0.8 km above Lindal's point, so the reference geoid
    is reproduced.
  - At the three Voyager targets it lies higher than Lindal's points: +3.3 km at 3.0 S
    (2.37 S planetocentric), +8.4 km at 31.2 S (26.06 S) and +7.1 km at 74.2 S (71.07 S).
  - Lindal's own geoids fit his points to 3 to 4 km. The two southern points are outside that.
  - The two hypotheses differ by 0.1 km at most: at 100 mbar, below the 110 mbar top level, their
    shear is the same.
- **The weight** w is 0.02 to 0.05 s/m per radian. It changes sign at the equator and is largest
  near 60 S, from 30 to 200 mbar.
- **The segment constraints, from 110 to 1300 mbar:**
  - segment 3 (26.06 S to 71.07 S) asks for at most +1.6 m/s per unit ln p;
  - segment 2 (2.37 S to 26.06 S) asks for -7 m/s per unit ln p at 110 mbar;
  - segment 1 (30.81 N to 2.37 S) asks for +1 to +3 m/s per unit ln p.
- **Above 30 mbar,** segments 1 and 2 ask for 10 to 20 m/s per unit ln p, changing sign with
  height. That is the equatorial wave structure of SPEC_15, read as shear.
- **The linearization holds to 2 to 3 percent.**
  - Delivered over predicted is 0.979 at c and 0.969 at 2c, with c = 3.5 m/s per unit ln p
    (+2.0 K at 31.2 S).
  - At zero amplitude, about 0.5 percent is the wind's own change, which the operator holds
    fixed, and 0.6 percent is unexplained.
  - A further 1 percent per 2 K of signal is second order, consistent with w's temperature held
    fixed.
  - 3.0 S, north of the perturbed segment, does not move (1.6e-6 in ln T).

## 1. Fig. 9, digitized (Part A item 1)

`tests/step17_1/digitize_fig9.py` renders journal page 1142 at 600 dpi, as SPEC_15 rendered Figs. 5
and 6.
- **The frame.** The top and bottom rows, the right column and the zero-altitude line are found
  from the ink.
- **Altitude:** a least-squares line through the 15 ticks, 0 to 140 km, RMS 0.24 km.
- **Latitude:** a line through the 19 ticks, 90 N to 90 S, RMS 0.27 deg. The ticks are drafted
  unevenly, by up to 5 px (0.5 deg).
- **The dots** are what survives an opening of the ink by an 11 px disc, at their centroids.
- **The stippled best-fit geoid** is scanned row by row, for the figure only.

`occul_data/lindal/isobar_100mbar.csv`:

| observation | latitude (planetocentric) | altitude (km) | the occultation's latitude (SPEC_15) |
|---|---:|---:|---:|
| voyager2 (anchor) | +30.82 | 61.0 | +30.81 |
| voyager1 | -2.96 | 121.3 | -2.37 |
| pioneer11 | -9.54 | 103.3 | -9.51 (Kliore's 9.5 S) |
| voyager2 | -26.20 | 59.7 | -26.06 |
| voyager1 | -71.84 | 13.3 | -71.07 |

The top comment holds the source, the reference geoid's definition, each point's occultation and
epoch, and Lindal's 3 to 4 km fit. No new keys.

**Three dots sit off their occultations' latitudes,** by -0.59, -0.14 and -0.77 deg. Interpolating
between the neighboring ticks alone gives the same latitudes, so this is where Lindal drew them,
not the calibration. The anchor and Pioneer 11 dots land within 0.03 deg. The drawn latitudes are
recorded as digitized, and the comparison handles the difference (§2).

Lindal's text says the Pioneer 11 point is from radio tracking, not occultation (p. 1141). The
comment says so; it enters nothing.

## 2. The 100 mbar surface (Part A item 2)

**The reference geoid** is `lib.geoid.reference_geoid`:
- through 54,438 km at the pole;
- the run's gravity (Null et al. 1981) and System III rotation (38,362.4 s, Lindal's
  10 h 39 m 22.4 s);
- no wind;
- the equatorial radius is 60,244.0 km.

**The delivered surface** is the run's `reference_surface`: the gauge isobar, 1.0e4 Pa, marched by
Eq. B3 along the path from the anchor. At each target it equals the delivered radius at 100 mbar to
about 20 m. At the anchor, it is the anchor's own `anchor_isobar_radius_m`.

| point | digitized | compared at | observed (km) | hypothesis 1 (km) | hypothesis 2 (km) | delivered minus observed (km) | slope (km/deg) |
|---|---:|---:|---:|---:|---:|---:|---:|
| anchor, 36.3 N | +30.82 | +30.81 | 61.0 | 61.8 | 61.8 | +0.8 | |
| 3.0 S, Voyager 1 egress | -2.96 | -2.37 | 121.3 | 124.6 | 124.6 | +3.3 | +0.85 |
| 31.2 S, Voyager 2 egress | -26.20 | -26.06 | 59.7 | 68.1 | 68.2 | +8.4, +8.5 | +1.63 |
| 74.2 S, Voyager 1 ingress | -71.84 | -71.07 | 13.3 | 20.4 | 20.4 | +7.1 | -0.19 |

**The comparison latitude.** Each digitized latitude lies just beyond its run's path (each path ends
0.05 deg past the target), so all three are compared at the target's latitude, as §1 item 2 says.
The slope is that of the delivered altitude there. Slope times the latitude difference is what the
drawn latitude would add: -0.5, -0.2 and +0.1 km. That leaves +2.8, +8.2 and +7.3 km. The
latitude does not explain the differences.

**What the figure shows** (`isobar_100mbar.png`):
- **Reference geoid reproduced.** The anchor agrees to 0.8 km, within Lindal's 3 to 4 km.
- **North of the equator,** the delivered surface follows Lindal's traced best-fit geoid closely.
- **At the equator,** it is about 3 km higher.
- **From 20 S to 60 S,** it runs 5 to 10 km outside Lindal's curve. The steps near 37 S and 52 S
  appear in both, displaced.
- **Lindal's geoid** comes from an altitude-independent wind (Smith et al. 1982; Ingersoll and
  Pollard 1982). Ours comes from the Sanchez-Lavega cloud wind with the IRIS shear between
  398 mbar and 100 mbar.
- **So the 100 mbar geometry asks for a smaller reference wind effect south of the equator.** It
  constrains the reference wind and the shear below 100 mbar together, which is what the fitted
  wind will weigh. This step does not separate them.
- **Pioneer 11** is drawn hollow and labeled "reference only".

## 3. The weight (Part B item 3)

`lindal_wind.balance_weight(wind, composition, gravity, rotation, pressure_Pa, temperature_K,
equatorial_radius_m=None)` sits beside `construct`.
- **What it computes:** construct's balance linearized in the shear with the wind held fixed:
  `w = -2 Omega_abs r D / g = 2 Omega_abs r m_bar (sin phi - (cos phi / r)(dr/dphi)_p) / (R T)`.
  g cancels.
- **The geometry is construct's,** level by level: the Eq. B3 surface under that level's wind,
  anchored at the gravity file's normalization radius, and `Omega_abs` under the same wind. So `w`
  at a level is exactly `1 / (ds/d(d ln T/dphi))` of construct's solution.
- **The temperature** is the run's own along its path. The isobars' ln N with the target's mean
  refractivity and each isobar's label reproduce the delivered temperature at the target to
  0.06 K.

Each hypothesis's w comes from its own 74.2 S run, whose path crosses all three segments. Values:
- 0.034 at 30 N and -0.050 at 70 S, at 100 mbar;
- 0.020 and -0.031 at 1 mbar;
- through zero at the equator, as `sin phi`.

The two hypotheses' w differ only above 110 mbar, through their temperature: their segment
integrals by up to 4 percent.
`weight.png` shows w over 30.8 N to 71.1 S.

**Cost.** Each level is one march of the wind geoid, about 8 s, so one run's w on 66 levels takes
9 minutes. The whole acceptance run took about 40 minutes.

## 4. The segment constraints (Part B item 4)

**Method.**
- **The misfit:** m_k = ln(T_obs / T_delivered) on each run's levels, the observed profile
  interpolated linearly in ln p, NaN outside its data. The levels are the 66 anchor levels, the
  same index in every run.
- **The segment's requirement** is the misfit differenced as §1 item 4 sets out.
- **The weight integral** is taken from the segment's north end to its south end, with the endpoints
  interpolated on the wind's 0.5 deg grid.

Segment-mean Δs, in m/s per unit ln p, hypothesis 1 / hypothesis 2 (all three rows are in
`output.txt`):

| mbar | 730 | 290 | 110 | 30 | 10 | 4.5 | 1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| segment 1, 30.81 N to 2.37 S | +0.05 / +0.04 | +2.56 / +2.55 | +0.95 / +1.04 | +5.8 / +8.4 | -11.5 / -9.2 | +7.9 / +7.2 | +20.6 / +12.6 |
| segment 2, 2.37 S to 26.06 S | -0.47 / -0.47 | +0.18 / +0.18 | -6.94 / -6.97 | +9.1 / +9.0 | -12.9 / -13.2 | +10.6 / +11.0 | +19.7 / +19.9 |
| segment 3, 26.06 S to 71.07 S | -0.74 / -0.74 | +0.78 / +0.78 | +1.60 / +1.60 | +1.0 / +1.0 | +2.0 / +2.0 | +2.0 / +2.0 | +0.4 / +0.5 |

**Reading the table.**
- **The weight integrals:**
  - segment 1 -0.006 to -0.010 s/m, small, because w passes through zero in it;
  - segment 2 +0.004 to +0.007;
  - segment 3 +0.021 to +0.034.
  Segment 1's ratio is therefore the least robust: a given misfit asks for 3 to 5 times segment 3's
  shear.
- **Below 110 mbar,** the hypotheses agree, since their shears are the same there.
  - Segment 2 asks for -7 m/s per unit ln p at 110 mbar. That is the 31.2 S misfit of
    SPEC_15 (delivered 5 K warm at 110 mbar) less 3.0 S's.
  - Segment 3 asks for +0.8 to +1.6 from 290 to 110 mbar.
  - Segment 3's value at 730 mbar, -0.74, is inside the Fig. 5 depth bias (1.5 to 3 K warm from 600
    to 1200 mbar, REVIEW_15_step1). It is marked and not corrected; a bias of 1.5 to 3 K there is
    -0.4 to -0.9 m/s per unit ln p, so the value is consistent with the bias alone.
- **Above 30 mbar,** segments 1 and 2 swing between -13 and +21 with height. Each swing is the
  equatorial profiles' layering. Segment 2's swings are 3.0 S's with the sign reversed, since it is
  31.2 S minus 3.0 S. Read as shear, this is not a wind a smooth basis can or should fit; both
  panels shade it (`segments.png`).
- **Segment 3 above 30 mbar** has no equatorial profile in it. It asks for a steady +1 to +2.
- **Hypothesis 2 versus hypothesis 1.** Hypothesis 2 asks less of segment 1 at 1 mbar (+12.6
  against +20.6) and the same of segment 3. Its reversal was scaled to 31.2 S, which segment 1's
  path does not reach.

## 5. The linearization, tested (Part B item 5)

**The perturbation** (`tests/step17_1/perturb.py`):
- Δs = c A(phi) B(p);
- A is 1 over segment 2, falling to zero as a raised cosine over the 2 deg inside each end, so
  that nothing reaches segments 1 or 3;
- B is 1 from 10 to 30 mbar, falling as a raised cosine over half a scale height outside each end,
  6.07 and 49.5 mbar, so that Δs is uniform over the band as §1 item 5 states;
- Δu is Δs integrated in ln p from the 398 mbar reference level, as the tool integrates it: zero
  below 49.5 mbar, and `-c A (ln 3 + 0.5)` above 6.07 mbar.

**For the author:** the taper sides, inside the segment in latitude and outside the band in
pressure, are my reading of "tapered at each end"; the specification does not say which side.

**The copies.**
- Each copy is a run directory with:
  - the registered run's namelist, the anchor's path made to reach the repository's data;
  - its composition, gravity and rotation;
  - the base wind, from `lindal_iris_v_sl_ii_110`, with Δu added to `u_total_ms` and `u_shear_ms`.
- A zero-amplitude copy reproduces the registered run exactly: temperature and radius differ by 0.0.
- Each run takes 44 s.
- c = 3.5 m/s per unit ln p moves 31.2 S by +2.0 K in the band; 2c by +4.0 K.

**Predicted:** the integral of `w Delta s` from the anchor to the target, with the base run's w.

| target | amplitude | delivered over predicted, mean [range] | whole balance over w Δs | delivered over whole balance |
|---|---|---|---:|---:|
| 31.2 S | c | 0.9789 [0.9635, 0.9888] | 0.9948 | 0.9841 |
| 31.2 S | 2c | 0.9685 [0.9544, 0.9777] | 0.9944 | 0.9740 |
| 74.2 S | c | 0.9774 [0.9625, 0.9870] | 0.9946 | 0.9827 |
| 74.2 S | 2c | 0.9671 [0.9535, 0.9759] | 0.9943 | 0.9727 |

- **Linearity:** Δ(2c) / (2 Δ(c)) is 0.9894 [0.9887, 0.9906] at both targets.
- **3.0 S:** largest |Δ ln T| 1.6e-6 at c and 3.3e-6 at 2c. Nothing reaches north of the segment.
- **Above the band,** where Δs is zero but Δu is not, the delivered change is +3e-5 in ln T, and
  the whole balance gives the same. The wind change there costs nothing measurable.

**The departure from one, explained.** Extrapolated to zero amplitude, delivered over predicted is
2(0.9789) - 0.9685 = 0.989. That 1.1 percent splits into two parts:
- **The wind's own change: 0.5 percent.**
  - This is the balance evaluated whole, `-(2 Omega_abs r / g)(s D + (cos phi / r)(du/dphi)_p)`,
    under the base and the perturbed winds with the temperature held. Its change, over `w Delta s`,
    is 0.9946.
  - The change is linear in c (0.9948 at c, 0.9944 at 2c).
  - It is the terms the linearization holds fixed: Δu in `(du/dphi)_p` at the tapers, and in
    `Omega_abs`, `g` and the surface.
- **The rest at first order: 0.6 percent, not explained here.** Delivered over the whole balance
  extrapolates to 0.994. The transfer integrates its own kernel, S/g along isobars traced in
  geopotential, with the geometry interpolated on its mesh. The balance here uses the Eq. B3 surface
  at each pressure. The difference between the two is the likely source, but this step does not
  isolate it.

**Second order: 1.0 percent per c.** That is 0.979 against 0.989 at c, and 0.969 at 2c. It is what
holding w's temperature fixed would give.
- Along segment 2 the perturbation warms the band by up to 2.1 percent at c.
- w goes as 1/T, so the operator overstates the change by about half the mean warming along the
  path.

For the fitted wind: the operator is good to about 1 percent at zero amplitude, and to 3 percent for
a 4 K correction. For larger corrections, one relinearization would bring it under 1 percent.

## 6. Files

New, uncommitted:
- `src/casspian/tools/lindal/lindal_wind.py`: `balance_weight` (one function, 40 lines);
- `occul_data/lindal/isobar_100mbar.csv`;
- `tests/step17_1/`: `accept_step17_1.py`, `digitize_fig9.py` and `perturb.py`;
- `reports/step17_1/` (git-ignored):
  - `output.txt`, `run.log` and `summary.json`;
  - `isobar_100mbar.png`, `weight.png`, `segments.png` and `linearization.png`;
  - `perturbed/`, 49 MB, the seven copies;
- this report.

From the author, untracked or modified, not touched:
- `docs/specs/SPEC_17_Voyager_Constraints.md` and `docs/specs/STATE.md`;
- `docs/specs/SPEC_18_Kliore_Kind_N.md` and `docs/casspian_status_2026-10-22.md`.
