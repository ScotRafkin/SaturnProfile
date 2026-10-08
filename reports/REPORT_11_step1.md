# REPORT 11, Step 1. The Lindal-only wind

Coding agent, 7 October 2026. Specification: `docs/specs/SPEC_11_Lindal_Wind.md` v0.3, Step 1 only.
Working tree on `main` at `e5a8d09`, uncommitted. Acceptance: `tests/step11_1/accept_step11_1.py`,
**5 of 5**; output in `reports/step11_1/output.txt`, figures and the four winds in
`reports/step11_1/`. No regression, as §0 directs.

**In short.**
- **Nothing blocks the step.** Five choices §1 leaves open are made and stated in §1 below.
- **The tool** builds the four winds in 5 or 6 iterations, each settling below 0.1 m/s. Each reads as
  kind W and passes the loaders' check.
- **Against the paper.** The 150 mbar shear drawn as Conrath and Pirraglia's Fig. 3 reproduces
  their figure wherever they drew it.
- **For the author** (§3):
  - the cubic bridge overshot inside the equatorial band and is now PCHIP (v0.4 §2 item 4);
  - the 5 deg band leaves the `1/sin(phi)` growth visible out to about 12 deg;
  - the held case reaches 660 m/s at 1 Pa.
- **Step 2.** A proposal for how the winds enter a run's build is in §4.

## 1. The tool, and the choices §1 leaves open

`src/casspian/tools/lindal/lindal_wind.py`: `construct` (no file) and `build` (a `[lindal_wind]`
control section), run as `python -m casspian.tools.lindal.lindal_wind <control>`. The module
docstring gives the method in full.

**The balance.** The relation the transfer uses: along an isobar `d ln N/dphi = K = S/g`. The
composition is uniform in latitude, which the run's is, so with `ln N = const - ln T` at fixed
pressure:

```
(d ln T/dphi)_p = -S/g,  S = 2 Omega_abs r (du/dZ)_R  (Eq. A15)
```

The derivatives are taken to fixed pressure as `lib.kernel` takes them. Solved for
`s = (du/dln p)_phi`:

```
s = [ -g (d ln T/dphi)_p / (2 Omega_abs r) - (cos(phi)/r) (du/dphi)_p ] / D
D = (dln p/dr) [ sin(phi) - (cos(phi)/r) (dr/dphi)_p ]
```

`D` vanishes like `sin(phi)`, the singularity of ruling 2. The second term in the numerator is the
meridional term: a wind that varies only with latitude has a temperature gradient under A15, so
the shear is what the data's gradient needs beyond it.

Choices made, none blocking:

1. **A fifth input, the run's composition.** `dln p/dr = -g m_bar/(R T)` needs the mean molar mass.
   It is read from the run's kind C file at each level and latitude, as the forward code reads it.
   `T` is the IRIS fit's own value at the level.
2. **The isobars' geometry.** At each level, `r(phi)` is the Eq. B3 surface under that level's wind
   (`lib.geoid.wind_geoid`). It is anchored at the gravity file's normalization radius,
   60000 km, because none of §1's inputs gives an isobar's radius. `g` and `Omega_abs` are
   `lib.gravity`'s under that wind. Moving the anchor 300 km either way changes the computed shear
   by at most 0.47, 0.29 and 0.26 m/s per scale height at the three levels (acceptance 4). That is
   about 1 percent of the shear's range.
3. **The equatorial band** is taken in planetocentric latitude, `|phi_c| < 5 deg`, the latitude of
   the model's balance axis.
4. **The output grid.** The cloud wind's 361 latitudes, and 121 pressures from 1 Pa to 1 MPa at 20
   per decade. Each pressure is rounded to six significant digits as the build files write them,
   so 39810.7 Pa is a node exactly.
5. **The flags and the uncertainty.**
   - `value_provenance` is the cloud wind's at 398 mbar and `parameterized` elsewhere.
   - The uncertainty is the cloud wind's column, carried unchanged to every pressure.
   - The globals that would be false are replaced: `vertical_structure`, `method`, `source`,
     `observation_level_Pa`, its justification, and `coverage_pressure_Pa`. The rest are the cloud
     wind's.

The four winds are built by the acceptance, under `reports/step11_1/<case>/`, from the
`lindal_transfer` run's inputs. They are not committed: Step 2 builds them in each run (§4). No
console entry is added (§4).

## 2. Acceptance

0. **Each wind can run the model.**

   | Case | Top level | Above it | Iterations | Last change | `u_total` |
   |---|---|---|---|---|---|
   | i_110 | 110 mbar | held | 6 | 0.066 m/s | -264 to 658 m/s |
   | ii_110 | 110 mbar | relaxed | 5 | 0.081 m/s | -16 to 537 m/s |
   | i_150 | 150 mbar | held | 6 | 0.017 m/s | -268 to 660 m/s |
   | ii_150 | 150 mbar | relaxed | 5 | 0.027 m/s | -15 to 537 m/s |

   For every case:
   - the wind reads as kind W;
   - it passes the loaders' check (`u_total` present, both poles exactly zero);
   - it is array-equal to its construction;
   - its shear at 398 mbar is exactly zero;
   - its three parts sum to round-off.
1. **The shear at the three levels** (`shear_levels.png`).
   - **Range of the computed shear:** -37 to 74 m/s per scale height at the top level, -32 to 35 at
     290 mbar, and -13 to 32 at 730 mbar. Every extreme lies within 12 deg of the equator.
   - **Placing the top level** at 150 rather than 110 mbar changes the top level's shear by at most
     about 1 m/s per scale height.
   - **Held against relaxed:** the level shears differ by 0.007 m/s per scale height at most, within
     the iteration's 0.1 m/s tolerance.
2. **As Fig. 3** (`shear_150mbar_as_fig3.png`). This is `-du/dz` in m/s per scale height, which is
   `du/dln p`, against planetographic latitude, beside the scan of their figure. It reproduces the
   figure:
   - about -37 at 15 N;
   - the +13 to +25 spike near 47 N, and the swings between 40 and 70 N;
   - -25 near -15;
   - the small structure poleward of -45.

   Their balance is geostrophic, with `f = 2 Omega sin(phi)`. The model's adds `Omega_abs`, the
   isobar's tilt and the meridional term, and these make the differences between 5 and 15 deg. They
   drew nothing there. We have +50 to +74 near 5 to 10 N and S, where they have +27 at 10 N.
3. **The winds** (`wind_<case>.png`): `u_total(phi, p)` and `u(p)` at 60, 30, 0, -30 and -60 deg.
4. **Beyond the specification:** the sensitivity to the geometry's anchor, in §1 item 2.

## 3. For the author

- **The bridge, replaced by PCHIP (v0.4 §2 item 4).** The first filing's cubic, matched in value
  and slope, peaked at about 100 m/s per scale height near 2 N at the top level, above any computed
  value (the edge was 75 and climbing), and sagged to 14 between edges of 25 and 34 at 290 mbar.
  PCHIP through the computed values stays between the two edge values:
  - at the top level it falls from about 74 at 5 N to 49 at -5;
  - at 290 mbar it rises from about 19 to 29;
  - at 730 mbar it falls from about 24 to 19.

  The figures and every number in this report are from the rerun with PCHIP.
- **The band's width.** `1/sin(phi)` is 5.7 times its value at 30 deg when at 5 deg, and 2.9 times
  when at 10 deg. The computed shear between 5 and 12 deg is where the largest values are, and where
  Conrath and Pirraglia drew nothing. The band stays at 5 deg (v0.4 §2 item 5): the shear there is 6
  to 10 times its estimated error.
- **The held case aloft.** Held at its 110 mbar value to 1 Pa, about 9 scale heights, the shear
  carries the wind to between -264 and +658 m/s. The relaxed case stays between -16 and +537 m/s.
  The transfer reads only the levels its anchor spans, so this matters to Step 2 only through the
  top of the anchor's profile.
- **Beyond the data's ends.** The taper starts at the last dot of each level: about 74 N
  planetocentric (76.6 to 76.8 N planetographic) and about 89 S. The northern taper over 16 deg is
  long against the structure there. The southern one is half a degree.

## 4. Proposal for Step 2: how the winds enter a run's build

A run's build is five sections in a fixed order, gravity, rotation, wind, shear and composition
(`casspian-run-inputs`). `[wind]` writes the source wind, the cloud wind this tool takes as input,
and `[shear]` writes the wind the namelist reads. The least change is to make this wind a case of
`[shear]`:

1. **A case `lindal_iris` in `casspian-wind-shear`**, whose function is `lindal_wind.construct`. Its
   parameters are `temperatures`, `composition`, `gravity` and `rotation` (paths, so added to the
   tool's `PATH_KEYS`), `top_level_Pa` and `above_top`. `source` is `[wind]`'s output, unchanged.
   The tool's refusals stand: a parameter the case does not use is refused as a typo, as now.
2. **`casspian-run-inputs` builds `[composition]` before `[shear]`**, since this case reads the
   run's composition. The composition reads no wind, so the order changes nothing else. Every
   accepted suite that builds a run would be rerun under that change.
3. **Each run** is made by `casspian-new-run <name> --from forward/lindal_transfer`. Its `[shear]`
   section is set to the case with its two parameters, and its `[target]` to -26.4 deg.
   - The four run names would follow the cases: `lindal_iris_i_110`, `lindal_iris_ii_110`,
     `lindal_iris_i_150` and `lindal_iris_ii_150`.
   - Nothing else in the build file or the namelist changes.

**What stays as it is:**
- `[wind]` and the other tools;
- every existing run and experiment, whose `[shear]` cases are unchanged;
- the model.

`lindal_wind.build` and its `[lindal_wind]` section would then be redundant and could be removed,
leaving `construct` as the one path, and no console entry is needed.

**The alternative** is a sixth build section for this tool. It changes the driver's fixed set of
sections and the check that refuses any other section, so it is more change for the same result.

## 5. Notes

- **A private import.** The acceptance calls `lib.control._admit_wind` directly, to check each wind
  against the loaders' own test without building a run.
- **Run time.** About 3 minutes a case, nearly all of it the geoid marches: three per iteration,
  one per level. The acceptance builds each case twice (through `build`, and through `construct`
  for its figures) plus two sensitivity cases, about 30 minutes.
